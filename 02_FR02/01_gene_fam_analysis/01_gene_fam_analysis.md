# Wide functional scan (HUMAnN3 gene families) — association with incident liver disease in FINRISK 2002 and correlation with gut E. coli

Genome-wide scan of gut functional features (HUMAnN3 UniRef90/EC gene families) for association with incident liver disease, followed by testing which associated features track gut E. coli abundance. This is a DISCOVERY analysis. Validated by the targeted pap-operon KO analysis in SIMPLER.

Exposure:
- Features: HUMAnN3 UniRef90/EC gene families (untargeted, genome-wide scan)
- Measures: detection (presence/absence) AND per-SD scaled abundance

Outcome is liver disease (ICD-10 K70–K77):
- Primary: incident (prospective) 
- Secondary: incident (primary, Cox), prevalent + ever (secondary, logistic)

Data:
- `health_data.RData`: health outcomes + covariates 
- `genefamilies_ec_uniref90_with_names_relab.tsv`: HUMAnN3 gene-family relative abundance (UniRef90 + EC names; RPK)
- `phyloseq.rds`: (Kraken2 / GTDB r214) taxonomic profile; E. coli abundance extracted as compositional log10(relative abundance + 1)

Analytical decisions:
- Covariates: age, sex, BMI, alcohol, recent antibacterial use, any diabetes (PREVAL_DIAB), sequencing depth (ReadDepthMillions) and batch
- Multiple testing: FDR (BH) 

All the data is available on CSC SD Desktop at the folder `/home/camigazo/Workspace/data`

```R
# ==============================================================================
# 0. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse)
library(survival)
library(ggstats)
library(patchwork)
options(width = 200, stringsAsFactors = FALSE)

setwd("/home/camigazo/")
DATA   <- "./data" # data lives in ../data, SD Desktop
OUTDIR <- "01_gene_fam_analysis"; dir.create(OUTDIR, showWarnings = FALSE)

# ==============================================================================
# 1. Load data  (FULL cohort with HUMAnN3 profiles)
# ==============================================================================
load(file.path(DATA, "health_data.RData")) 
ff <- health$followup; bl <- health$baseline; cv <- health$covariates
metadata <- ff[, c("Barcode","PREVAL_LIVERDIS","INCIDENT_LIVERDIS","LIVERDIS_AGEDIFF", "PREVAL_IBD", "PREVAL_DIAB","INCIDENT_DIAB", "BL_USE_RX_J01", "BL_USE_RX_A02BC", "BL_USE_RX_C10AA")]
metadata <- merge(metadata, bl[, c("Barcode","BL_AGE","MEN","BMI","ALKI2_FR02", "CURR_SMOKE", "FIBER_TOTAL", "KOULGR")], by="Barcode", all.x=TRUE)
metadata <- merge(metadata, cv[, c("Barcode","batch_name","ReadDepthMillions")], by="Barcode", all.x=TRUE)
dim(metadata) # 7098 x 19

# ---- min analytic sample: no prevalent liver, complete covariates ----------------
# we will do models later with Tier 1: BL_USE_RX_A02BC, BL_USE_RX_C10AA, CURR_SMOKE and Tier 2: KOULGR, PREVAL_IBD, FIBER_TOTAL
metadata <- metadata %>%
  filter(PREVAL_LIVERDIS == 0) %>%
  mutate(PREVAL_DIAB = factor(PREVAL_DIAB), MEN = factor(MEN), batch_name = factor(batch_name)) %>%
  filter(if_all(c(BL_AGE, MEN, BMI, BL_USE_RX_J01, ALKI2_FR02, batch_name, ReadDepthMillions), ~ !is.na(.)))

cat("analytic sample:", nrow(metadata), "| incident LD cases:", sum(metadata$INCIDENT_LIVERDIS == 1), "\n")
# analytic sample: 6735 | incident LD cases: 133 

# =========================================================================
# 2. HUMAnN3 features — detection (5–95%) + per-SD abundance (>=5%, no cap)
# =========================================================================
gene_fam <- read.csv(file.path(DATA, "genefamilies_ec_uniref90_with_names_relab.tsv"), sep="\t") %>%
  column_to_rownames("X..Gene.Family") %>% t() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  mutate(SampleID = gsub(".R1.trimmed.filtered_Abundance.RPKs", "", SampleID),
         SampleID = gsub("^X", "", SampleID),
         SampleID = gsub("\\.", "-", SampleID)) %>%
  filter(SampleID %in% metadata$Barcode)

# long form, simplify names (keep dictionary)
gl <- reshape2::melt(gene_fam); colnames(gl) <- c("SampleID","Gene.Family","RPK")
gene_dic <- data.frame(Gene.Family = unique(gl$Gene.Family), Gene.Name = gsub("(.*): .*", "\\1", unique(gl$Gene.Family)))
gl <- gl %>% mutate(Gene = gsub("(.*): .*", "\\1", Gene.Family))

# ---- prevalence filters (cohort-wide) -----------------------------------
N    <- n_distinct(gl$SampleID)
prev <- gl %>% filter(RPK > 0) %>% count(Gene, name = "n_det") %>% mutate(prev = n_det / N)

# detection: two-sided 5–95% (binary predictor needs enough 0s AND 1s)
FEAT_DET <- prev %>% filter(prev >= 0.05, prev <= 0.95) %>% pull(Gene)
cat(sprintf("detection features (5-95%%): %d (of %d detected)\n", length(FEAT_DET), nrow(prev)))
#detection features (5-95%): 3551 (of 9591 detected)

# abundance: lower bound only — continuous variation persists in ubiquitous
# genes, so the 95%% cap has no rationale here
FEAT_ABZ <- prev %>% filter(prev >= 0.05) %>% pull(Gene)
cat(sprintf("abundance features (>=5%%):  %d (of %d detected)\n", length(FEAT_ABZ), nrow(prev)))
#abundance features (>=5%):  5255 (of 9591 detected)

# ---- detection matrix (presence/absence) — FEAT_DET ---------------------
det_mat <- gl %>% filter(Gene %in% FEAT_DET) %>%
  mutate(d = as.integer(RPK > 0)) %>%
  select(SampleID, Gene, d) %>% pivot_wider(names_from = Gene, values_from = d, values_fill = 0)

# ---- per-SD scaled abundance: log10(RPK+1) then scale — FEAT_ABZ --------
abz_mat <- gl %>% filter(Gene %in% FEAT_ABZ) %>%
  mutate(a = log10(RPK + 1)) %>%
  select(SampleID, Gene, a) %>% pivot_wider(names_from = Gene, values_from = a, values_fill = 0) %>%
  mutate(across(-SampleID, ~ as.numeric(scale(.x))))

# ---- unscaled log10 abundance (for E. coli correlation / scatter) — FEAT_ABZ ----
ablog_mat <- gl %>% filter(Gene %in% FEAT_ABZ) %>%
  mutate(a = log10(RPK + 1)) %>%
  select(SampleID, Gene, a) %>% pivot_wider(names_from = Gene, values_from = a, values_fill = 0)

# join exposure matrices to metadata (Barcode == SampleID)
M <- metadata %>%
  left_join(det_mat %>% rename_with(~ paste0(.x, "_det"), -SampleID), by = c("Barcode"="SampleID")) %>%
  left_join(abz_mat %>% rename_with(~ paste0(.x, "_abz"), -SampleID), by = c("Barcode"="SampleID")) %>%
  mutate(INCIDENT_LIVERDIS = as.numeric(as.character(INCIDENT_LIVERDIS)))

# how many incident events?
sum(M$INCIDENT_LIVERDIS == 1) # 133

# ==============================================================================
# 3. WIDE SCAN  (Cox per feature, detection and abundance)
#   FDR(BH) within each scan
# ==============================================================================
COV <- "BL_AGE + MEN + BMI + BL_USE_RX_J01 + ALKI2_FR02 + PREVAL_DIAB + batch_name + ReadDepthMillions"

scan_cox <- function(features, suffix, data) {
  map_dfr(features, function(g) {
    col <- paste0(g, suffix)
    m <- tryCatch(coxph(as.formula(paste0("Surv(LIVERDIS_AGEDIFF, INCIDENT_LIVERDIS) ~ ", col, " + ", COV)),
                        data = data), error = function(e) NULL)
    if (is.null(m) || !col %in% rownames(summary(m)$coef)) return(NULL)
    cf <- summary(m)$coef[col, ]
    z  <- tryCatch(cox.zph(m), error = function(e) NULL)                                   # <-- added
    tibble(gene = g,
           HR = exp(cf["coef"]), lo = exp(cf["coef"] - 1.96 * cf["se(coef)"]), hi = exp(cf["coef"] + 1.96 * cf["se(coef)"]),
           p = cf["Pr(>|z|)"], N = m$n, n_events = m$nevent,
           ph_term   = if (!is.null(z) && col %in% rownames(z$table)) z$table[col, "p"] else NA,   # <-- added
           ph_global = if (!is.null(z)) z$table["GLOBAL", "p"] else NA)                            # <-- added
  }) %>% mutate(measure = sub("_", "", suffix), FDR = p.adjust(p, method = "fdr"))
}

res_det <- scan_cox(FEAT_DET, "_det", M)
cat(sprintf("detection scan: %d features, %d at FDR<0.05\n", nrow(res_det), sum(res_det$FDR < 0.05)))
# detection scan: 3551 features, 612 at FDR<0.05

res_ab <- scan_cox(FEAT_ABZ, "_abz", M)
cat(sprintf("abundance scan: %d features, %d at FDR<0.05\n", nrow(res_ab), sum(res_ab$FDR < 0.05)))
# abundance scan: 3551 features, 994 at FDR<0.05

res <- bind_rows(res_det, res_ab) %>%
  left_join(gene_dic, by = c("gene" = "Gene.Name")) %>% # attach full descriptive name
  arrange(measure, p)
write.csv(res, file.path(OUTDIR, "wide_scan_cox_612det_994abz.csv"), row.names = FALSE)

intersect(res_ab  %>% filter(FDR < 0.05) %>% pull (gene), 
          res_det %>% filter(FDR < 0.05) %>% pull (gene)) %>% 
  length() # 403

# significant features = union of DETECTION and ABUNDANCE scans
sig_det <- res_det %>% filter(FDR < 0.05)        # detection-significant
sig_abz <- res_ab  %>% filter(FDR < 0.05)        # abundance-significant
sig     <- bind_rows(sig_det, sig_abz)           # one row per (gene x scan); tagged by `measure`
sig_genes <- unique(sig$gene)                    # deduplicated — correlation/counts are gene-level

cat(sprintf("\nsignificant: %d detection + %d abundance = %d associations across %d unique features\n", nrow(sig_det), nrow(sig_abz), nrow(sig), length(sig_genes)))
#significant: 612 detection + 994 abundance = 1606 associations across 1203 unique features

# ==============================================================================
# 4. E. coli CORRELATION  (do the significant features track gut E. coli?)
#   phyloseq E. coli, compositional log10(rel+1)
# ==============================================================================
library(phyloseq)
ps <- readRDS(file.path(DATA, "phyloseq.rds"))
ps <- subset_taxa(ps, Domain == "Bacteria")
ps <- transform_sample_counts(ps, function(x) x / sum(x))
ecoli <- subset_taxa(ps, Species == "Escherichia coli") %>% psmelt() %>%
  transmute(Barcode = gsub("\\.", "-", Sample), E_coli_log10 = log10(Abundance + 1)) %>%
  filter(Barcode %in% M$Barcode)

cat(sprintf("E. coli extracted: %d samples (analytic N = %d) | non-zero: %d | median(>0) = %.4f\n",
            nrow(ecoli), nrow(M), sum(ecoli$E_coli_log10 > 0),
            median(ecoli$E_coli_log10[ecoli$E_coli_log10 > 0])))
# E. coli extracted: 6735 samples (analytic N = 6735) | non-zero: 2857 | median(>0) = 0.0010

# Spearman: each SIGNIFICANT feature's log10 abundance vs E. coli, BH-FDR
sig_ablog <- ablog_mat %>% select(SampleID, any_of(sig_genes)) %>%
  inner_join(ecoli, by = c("SampleID" = "Barcode"))
corr <- map_dfr(sig_genes, function(g) {
  if (!g %in% names(sig_ablog)) return(NULL)
  ct <- suppressWarnings(cor.test(sig_ablog[[g]], sig_ablog$E_coli_log10, method = "spearman"))
  tibble(gene = g, rho = unname(ct$estimate), p_cor = ct$p.value)
}) %>% mutate(FDR_cor = p.adjust(p_cor, method = "fdr"))

results_final <- sig %>%
  left_join(gene_dic, by = c("gene" = "Gene.Name")) %>%   # adds Gene.Family
  select(gene, gene_name = Gene.Family, measure, HR, lo, hi, p, FDR, ph_global, ph_term) %>%
  left_join(corr, by = "gene")                            # rho/FDR_cor broadcast to both rows of a gene
write.csv(results_final, file.path(OUTDIR, "sig_features_Ecoli_correlation.csv"), row.names = FALSE)

cat(sprintf("of %d unique significant features: %d correlate with E. coli (FDR<0.05)\n", nrow(corr), sum(corr$FDR_cor < 0.05, na.rm = TRUE)))
# of 1203 unique significant features: 1130 correlate with E. coli (FDR<0.05)

# top 20 features by correlation with E. coli (unique features)
results_final %>%
  distinct(gene, .keep_all = TRUE) %>% # gene sig in both scans -> one row
  arrange(desc(rho)) %>% slice_head(n = 20) %>%
  select(gene_name, measure, HR, lo, hi, FDR, rho, ph_global, ph_term)

# ==============================================================================
# 5. FIGURES
# ==============================================================================

# =========================================================================
# 5a. Forest plot — detection AND abundance HRs per feature (Fig 1b styling)
#     blue = detection, red = abundance; filled = FDR<0.05, open = n.s.
# =========================================================================
library(stringr)

# ---- label shortener: trim boilerplate, keep meaning ----
shorten <- function(x) {
  x %>%
    str_replace("Protein of unknown function \\((DUF[0-9]+)\\)", "\\1") %>%
    str_replace("Uncharacterized protein conserved in bacteria \\((DUF[0-9]+)\\)", "\\1") %>%
    str_replace("N-terminal domain", "N-term") %>%
    str_replace("C-terminal domain", "C-term") %>%
    str_replace(" domain$",  "") %>%
    str_replace(" family$",  "") %>%
    str_replace(" protein$", "") %>%
    str_replace("Pili and flagellar-assembly chaperone, ", "") %>%
    str_squish() %>%
    str_trunc(30)
}

# hand-curated overrides for features that must read perfectly
manual <- c(
  "PF13953" = "PapC usher C-term",
  "PF13954" = "PapC usher N-term",
  "PF02753" = "PapD chaperone C-term",
  "PF00345" = "PapD chaperone N-term"
)

pap_ids <- c("PF13953", "PF13954", "PF02753", "PF00345")

# ---- features to display: top 20 unique features by E. coli correlation ----
feat_show <- results_final %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(desc(rho)) %>% slice_head(n = 20) %>% pull(gene)

# ---- pull BOTH scan estimates from `res` (all tested features) ----
#      so every displayed feature shows its detection AND abundance HR,
#      whether or not it was FDR-significant in that scan
both <- res %>%
  filter(gene %in% feat_show) %>%
  mutate(measure    = recode(measure, det = "Detection", abz = "Abundance"),
         measure    = factor(measure, levels = c("Detection", "Abundance")),
         sig        = FDR < 0.05,
         fillvar    = ifelse(sig, as.character(measure), "n.s."),   # white -> open circle
         desc_raw   = sub("^PF[0-9]+: ", "", as.character(Gene.Family)),
         desc_short = ifelse(gene %in% names(manual), manual[gene], shorten(desc_raw)),
         label      = sprintf("%s (%s)", desc_short, gene),
         is_pap     = gene %in% pap_ids)

# ---- order y-axis by HR detection ----
ord  <- both %>% filter(measure == "Detection") %>%
  arrange(HR) %>% pull(label)
both <- both %>% mutate(label = factor(label, levels = ord))

# ---- bold pap rows on the y-axis, in the factor order ----
lab_face <- ifelse(levels(both$label) %in% both$label[both$is_pap], "bold", "plain")

# ---- colours: keep blue (detection) / red (abundance) ----
cols  <- c("Detection" = "#00798C", "Abundance" = "#D1495B")
dodge <- position_dodge(width = 0.6)

p <- ggplot(both, aes(HR, label, colour = measure, group = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = fillvar), shape = 21, size = 2.6, stroke = 0.7, position = dodge) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  scale_colour_manual(values = cols, name = "Measure") +
  scale_fill_manual(values = c(cols, "n.s." = "white"), guide = "none") +
  guides(colour = guide_legend(override.aes = list(shape = 16, linetype = 0, size = 3))) +
  labs(x = "Hazard ratio (95% CI)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
    axis.line.x  = element_line(colour = "grey20", linewidth = 0.4),
    axis.ticks.x = element_line(colour = "grey20", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    axis.text.y  = element_text(size = 9, hjust = 0, colour = "grey15", face = lab_face),
    axis.text.x  = element_text(colour = "grey15"),
    axis.title.x = element_text(colour = "grey15"),
    legend.position    = c(0.88, 0.92),
    legend.background   = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    legend.title = element_text(face = "bold", size = 10),
    legend.key   = element_blank(),
    plot.margin  = margin(6, 14, 6, 6)
  )

ggsave(file.path(OUTDIR, "Forest_top20_det_abz.pdf"), p, width = 6, height = 6)

# =========================================================================
# 5b. E. coli scatter for the four pap 
# =========================================================================
pap_pfams <- c("PF13953","PF13954","PF02753","PF00345")
pap_short <- c(PF13953 = "PapC usher (C-terminal)",
               PF13954 = "PapC usher (N-terminal)",
               PF02753 = "PapD chaperone (C-terminal)",
               PF00345 = "PapD chaperone (N-terminal)")

pap_scat <- map_dfr(pap_pfams, function(g) {
  ablog_mat %>% select(SampleID, all_of(g)) %>%
    inner_join(ecoli, by = c("SampleID" = "Barcode")) %>%
    transmute(gene = pap_short[g], gene_ab = .data[[g]], E_coli_log10)
}) %>% filter(E_coli_log10 > 0, !is.na(gene_ab)) %>%
  mutate(gene = factor(gene, levels = pap_short))   # fix panel order

# rho + P (P underflows -> reported as "< 0.001"); parsed annotation label
pap_stat <- pap_scat %>% group_by(gene) %>%
  summarise(rho   = cor(gene_ab, E_coli_log10, method = "spearman"),
            p_raw = cor.test(gene_ab, E_coli_log10, method = "spearman", exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(p_txt = ifelse(p_raw < 0.0001, "< 0.0001", sprintf("= %.3f", p_raw)),
         label = sprintf("italic(rho)==%.2f*','~italic(P)*'%s'", rho, p_txt))

# accent: red ties this to the abundance story in 5a (this panel IS an
# abundance relationship). swap to "#2c7fb8" (blue) or "grey20" (neutral).
accent <- "#d7301f"

p_pap <- ggplot(pap_scat, aes(E_coli_log10, gene_ab)) +
  geom_point(alpha = 0.18, size = 0.55, colour = "grey45", stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, colour = accent, fill = accent,
              alpha = 0.15, linewidth = 0.7) +
  geom_text(data = pap_stat, parse = TRUE, hjust = 0, vjust = 1,
            aes(x = -Inf, y = Inf, label = label),
            size = 2.8, colour = "grey15") +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +   # right margin so points/CI don't hit edge
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +   # top margin so rho/P clears the data
  labs(x = expression(log[10]("relative abundance of "*italic("E. coli")*" + 1")),
       y = expression(log[10]("relative abundance of "*italic("pap")*" gene + 1"))) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.line   = element_line(colour = "grey20", linewidth = 0.4),
    axis.ticks  = element_line(colour = "grey20", linewidth = 0.4),
    axis.text   = element_text(colour = "grey15", size = 7.5),
    axis.title  = element_text(colour = "grey15"),
    strip.background = element_blank(),                                  # no grey strip boxes
    strip.text       = element_text(face = "bold", size = 8.5, hjust = 0, colour = "grey15"),
    panel.spacing    = unit(10, "pt"),                                   # space between panels
    plot.margin      = margin(6, 10, 6, 6)
  )

ggsave(file.path(OUTDIR, "Scatter_pap_vs_Ecoli_4panel.pdf"), p_pap, width = 5.5, height = 5, device = cairo_pdf)

# =========================================================================
# 5c. HR vs E. coli correlation — separate detection and abundance plots
# =========================================================================
upec      <- c("PF13953","PF13954","PF02753","PF00345")
pap_short <- c(PF13953 = "PapC usher (C-term)",    PF13954 = "PapC usher (N-term)",
               PF02753 = "PapD chaperone (C-term)", PF00345 = "PapD chaperone (N-term)")

# split the union back into the two scans; each row carries that scan's HR
prep <- function(meas) {
  results_final %>% filter(measure == meas) %>%
    mutate(ecoli_sig = FDR_cor < 0.05,
           is_pap    = gene %in% upec,
           pap_label = ifelse(is_pap, pap_short[gene], NA_character_),
           fillcat   = ifelse(ecoli_sig, "corr", "ns"))
}
det_dat <- prep("det")
det_dat$ecoli_sig %>% table()
#FALSE  TRUE 
#8   604

abz_dat <- prep("abz")

# one builder, parameterised by point colour
make_panel <- function(d, point_col) {
  ggplot(d, aes(HR, rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_point(aes(fill = fillcat), shape = 21, colour = point_col,
               size = 1.4, stroke = 0.5) +
    ggrepel::geom_text_repel(data = ~filter(.x, is_pap), aes(label = pap_label),
                             size = 2.3, colour = "grey15", fontface = "bold",
                             min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf,
                             segment.size = 0.2, segment.colour = "grey60") +
    scale_fill_manual(values = c(corr = point_col, ns = "white"), guide = "none") +
    scale_x_log10() +
    labs(x = "Incident liver disease, HR",
         y = expression("Spearman "*italic(rho)*" with "*italic("E. coli"))) +
    theme_classic(base_size = 10) +
    theme(axis.line  = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3))
}

p_det <- make_panel(det_dat, "#00798C")   # teal     = detection
p_abz <- make_panel(abz_dat, "#D1495B")   # red       = abundance

ggsave(file.path(OUTDIR, "Scatter_HR_vs_rho_detection.pdf"), p_det, width = 6, height = 6)
ggsave(file.path(OUTDIR, "Scatter_HR_vs_rho_abundance.pdf"), p_abz, width = 6, height = 6)

# =========================================================================
# pap genes: full detection & abundance models + covariate validation
# =========================================================================
library(broom)

genes <- c("PF13953","PF13954","PF02753","PF00345")
pap_name <- c(PF13953="PapC usher C-term", PF13954="PapC usher N-term", PF02753="PapD chaperone C-term", PF00345="PapD chaperone N-term")

# --- 1. covariate sets -------------------
COV_BASE <- "BL_AGE + MEN + BMI + BL_USE_RX_J01 + ALKI2_FR02 + PREVAL_DIAB + batch_name + ReadDepthMillions"
COV_ADJ  <- paste(COV_BASE, "+ BL_USE_RX_A02BC + BL_USE_RX_C10AA + CURR_SMOKE + factor(KOULGR) + FIBER_TOTAL")

# complete-case sample for the adjusted covariates (fair base-vs-adjusted base)
adj_vars <- c("BL_AGE","MEN","BMI","BL_USE_RX_J01","ALKI2_FR02","PREVAL_DIAB",
              "batch_name","ReadDepthMillions",
              "BL_USE_RX_A02BC","BL_USE_RX_C10AA","CURR_SMOKE","KOULGR","FIBER_TOTAL")
M_cc <- M %>% filter(if_all(all_of(adj_vars), ~ !is.na(.)))

cat(sprintf("complete-case (adjusted) sample: %d\n", nrow(M_cc)))
# complete-case (adjusted) sample: 6267

# --- 3. fit one pap model, return the pap term only ----------------------
fit_pap_term <- function(g, meas, cov, data, model_lab, sample_lab) {
  expo <- paste0(g, "_", meas)
  if (!expo %in% names(data)) { warning("absent: ", expo); return(NULL) }
  f <- as.formula(sprintf("Surv(LIVERDIS_AGEDIFF, INCIDENT_LIVERDIS) ~ `%s` + %s", expo, cov))
  s <- summary(coxph(f, data = data))
  ci <- s$conf.int[1, ]
  tibble(gene = g, gene_name = pap_name[g], measure = meas,
         model = model_lab, sample = sample_lab,
         HR = unname(ci["exp(coef)"]), lo = unname(ci["lower .95"]),
         hi = unname(ci["upper .95"]), p = unname(s$coefficients[1, 5]),
         n = s$n, events = s$nevent)
}

# --- 4. robustness table: base(full, original), adjusted(CC) ------------
robust <- bind_rows(lapply(genes, function(g) bind_rows(
  fit_pap_term(g, "det", COV_BASE, M,     "Base",     "Full"),   # original model, as in main text
  fit_pap_term(g, "det", COV_ADJ,  M_cc,  "Adjusted", "CC"),     # validation
  fit_pap_term(g, "abz", COV_BASE, M,     "Base",     "Full"),
  fit_pap_term(g, "abz", COV_ADJ,  M_cc,  "Adjusted", "CC")
))) %>%
  mutate(measure = recode(measure, det = "Detection", abz = "Abundance"),
         est = sprintf("%.2f (%.2f–%.2f)", HR, lo, hi))

write.csv(robust, file.path(OUTDIR, "pap_models_robustness.csv"), row.names = FALSE)

robust$p %>% max()
robust$HR %>% summary()

# --- 5. FULL coefficient tables for the adjusted models (supplement) -----
fit_pap_full <- function(g, meas, cov, data) {
  expo <- paste0(g, "_", meas)
  f <- as.formula(sprintf("Surv(LIVERDIS_AGEDIFF, INCIDENT_LIVERDIS) ~ `%s` + %s", expo, cov))
  m <- coxph(f, data = data)
  ph <- tryCatch(cox.zph(m)$table["GLOBAL","p"], error = function(e) NA_real_)
  broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(gene = g, gene_name = pap_name[g],
              measure = recode(meas, det = "Detection", abz = "Abundance"),
              term, HR = estimate, lo = conf.low, hi = conf.high, p = p.value,
              ph_global = ph)
}
full_coef <- bind_rows(lapply(genes, function(g) bind_rows(
  fit_pap_full(g, "det", COV_ADJ, M_cc),
  fit_pap_full(g, "abz", COV_ADJ, M_cc)
))) 
write.csv(full_coef, file.path(OUTDIR, "pap_models_full_coefficients.csv"), row.names = FALSE)

# =========================================================================
# Characteristics table for the 6,735-participant functional-analysis sample
# =========================================================================

# re-attach PREVAL_DIAB_T2 (T2D) from health, keyed by Barcode
metadata <- metadata %>% left_join(health$followup %>% select(Barcode, PREVAL_DIAB_T2), by = "Barcode")  # use real name from grep
cat("T2D-diabetes prevalence column joined; non-missing:", sum(!is.na(metadata$PREVAL_DIAB_T2)), "of", nrow(metadata), "\n")
# T2D-diabetes prevalence column joined; non-missing: 6735 of 6735

# bring in follow-up time (time from baseline to death OR administrative censoring)
metadata <- metadata %>% left_join(health$followup %>% select(Barcode, DEATH_AGEDIFF), by = "Barcode")
cat("Follow-up (DEATH_AGEDIFF) non-missing:", sum(!is.na(metadata$DEATH_AGEDIFF)),
    "of", nrow(metadata), "| median:", round(median(metadata$DEATH_AGEDIFF, na.rm = TRUE), 1), "\n")
#Follow-up (DEATH_AGEDIFF) non-missing: 6735 of 6735 | median: 17.8 

# --- formatting helpers --------------------------------------------------
med_iqr <- function(x) {
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf("%g (%g, %g)", round(q[1], 1), round(q[2], 1), round(q[3], 1))
}
n_pct <- function(x, level = 1) {            # x is 0/1 or logical; % of non-missing
  x <- x[!is.na(x)]
  n <- sum(as.integer(as.character(x)) == level)
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * n / length(x))
}

# --- table (now with depth + follow-up rows) -----------------------------
tab <- tribble(
  ~Characteristic, ~Value,
  "Age at baseline, years, median (IQR)",                          med_iqr(metadata$BL_AGE),
  "Men, n (%)",                                                    n_pct(metadata$MEN, level = 1),
  "BMI, kg/m2, median (IQR)",                                      med_iqr(metadata$BMI),
  "Antibacterial (ATC J01) use ≤6 mo before baseline, n (%)",      n_pct(metadata$BL_USE_RX_J01, level = 1),
  "Alcohol consumption past year, g/week, median (IQR)",          med_iqr(metadata$ALKI2_FR02),
  "Sequencing depth (million reads), median (IQR)",
  sprintf("%.2f (%.2f, %.2f)", median(metadata$ReadDepthMillions, na.rm = TRUE),
          quantile(metadata$ReadDepthMillions, .25, na.rm = TRUE),
          quantile(metadata$ReadDepthMillions, .75, na.rm = TRUE)),
  "Follow-up, years, median (IQR)",                               med_iqr(metadata$DEATH_AGEDIFF),
  "Incident liver disease (ICD-10 K70-K77) over follow-up, n (%)", n_pct(metadata$INCIDENT_LIVERDIS, level = 1),
  "Type 2 diabetes prevalence (ICD-10 E11), n (%)",               n_pct(metadata$PREVAL_DIAB_T2, level = 1),
  "Any diabetes prevalence (ICD-10 E10-E14), n (%)",              n_pct(metadata$PREVAL_DIAB, level = 1)
)
print(as.data.frame(tab), right = FALSE)
# Characteristic                                                Value            
# 1  Age at baseline, years, median (IQR)                          50.1 (38.9, 59.1)
# 2  Men, n (%)                                                    3,072 (45.6%)    
# 3  BMI, kg/m2, median (IQR)                                      26.4 (23.7, 29.4)
# 4  Antibacterial (ATC J01) use ≤6 mo before baseline, n (%)      1,178 (17.5%)    
# 5  Alcohol consumption past year, g/week, median (IQR)           36 (9, 102.9)    
# 6  Sequencing depth (million reads), median (IQR)                0.76 (0.55, 1.18)
# 7  Follow-up, years, median (IQR)                                17.8 (17.8, 17.9)
# 8  Incident liver disease (ICD-10 K70-K77) over follow-up, n (%) 133 (2.0%)       
# 9  Type 2 diabetes prevalence (ICD-10 E11), n (%)                195 (2.9%)       
# 10 Any diabetes prevalence (ICD-10 E10-E14), n (%)               233 (3.5%) 

health$followup %>%
  filter(Barcode %in% metadata$Barcode) %>%
  select(INCIDENT_LIVERDIS, INCIDENT_ALC_HEPAT, INCIDENT_FIBROS_CHIRROS, INCIDENT_K11_ALCOLIV, INCIDENT_K11_TOXLIV,
         INCIDENT_K11_DISLIVOTH, INCIDENT_K11_FIBROCHIRLIV, INCIDENT_K11_OTHINFLIV, INCIDENT_NAFLD) %>%
  mutate(across(everything(), ~ as.numeric(. == 1))) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% t()  # <<< count of cases per column
#INCIDENT_LIVERDIS          133
#INCIDENT_ALC_HEPAT          54
#INCIDENT_FIBROS_CHIRROS     19
#INCIDENT_K11_ALCOLIV        54
#INCIDENT_K11_DISLIVOTH      53
#INCIDENT_K11_FIBROCHIRLIV   19
#INCIDENT_K11_OTHINFLIV       8
#INCIDENT_NAFLD              23
```
