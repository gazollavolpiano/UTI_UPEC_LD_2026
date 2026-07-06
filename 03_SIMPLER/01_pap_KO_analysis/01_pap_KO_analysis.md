# Targeted UPEC fimbrial system (P fimbriae) — pap KO association with liver disease in SIMPLER (Clinical Microbiomics KEGG ontology profiles)

Test if the best-characterised UPEC fimbrial system [P fimbriae (pap operon)] is associated with liver disease in SIMPLER with the Clinical Microbiomics KEGG ontology (KO) profiles

Exposure: 
- Genes: pap operon (papA..papK), per-gene AND combined operon
- Measures: detection (presence/absence) AND per-SD scaled abundance

Outcome is liver disease (ICD-10 K70–K77):
- Primary: incident (prospective) 
- Secondary: incident (primary, Cox), prevalent + ever (secondary, logistic)

Data:
- `tidy_clinical_data_6147.RData` — health outcomes + covariates, built in `00_data_organization_SIMPLER.md` 
- `/proj/.../simpler_metagenomics_ko_relative_abundances_v4.0.tsv` — KO NON-downsampled (for ABUNDANCE)
- `/proj/.../simpler_metagenomics_ko_ds_relative_abundances_v4.0.tsv` — KO downsampled (for PREVALENCE / detection / composition)
- `/proj/.../simpler_metagenomics_mgs_relative_abundances_v4.0.tsv` — MGS NON-downsampled; extract E. coli = hMGS.00032 (for abundance)

Analytical decisions:
- Tables*: non-downsampled tables for ABUNDANCE models; downsampled tables for PREVALENCE / detection
- Covariates: age, sex, BMI, alcohol, recent antibacterial use, sequencing depth (ReadDepthMillions), diabetes (PREVAL_E10E14) 
- Multiple testing: FDR (BH) 
(* is based on `simpler_metagenomics_qc_report_v4.0.pdf`)

Modelling strategy: batch (sequencing plate) has 150 levels, so it is handled two ways and results compared:
- Cox model (survival): plate NOT adjusted; PH check reported
- Mixed-effects Cox model (coxme): plate as a random intercept; sensitivity comparison confirming plate adjustment does not change estimates

```R
# ==============================================================================
# 0. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse)
library(survival)
library(coxme)
library(ggstats)
options(width = 200, stringsAsFactors = FALSE)

OUTDIR <- "01_pap_KO_analysis"; dir.create(OUTDIR, showWarnings = FALSE)

# ==============================================================================
# 1. Load data
# ==============================================================================
load("tidy_clinical_data_6147.RData"); health <- df; rm(df)
dim(health)   # 6147 x 40

# E. coli abundance (hMGS.00032), NON-downsampled MGS table
E_coli <- read_tsv("/proj/simp2024014/Omicsdataleverans2/metagenomics_CM/simpler_metagenomics_mgs_relative_abundances_v4.0.tsv",
                   show_col_types = FALSE) %>%
  transmute(SIMPKEY, E_coli_log10 = log10(hMGS.00032 + 1))
health <- health %>% left_join(E_coli, by = "SIMPKEY")
health$E_coli_z <- as.numeric(scale(health$E_coli_log10))   # per-SD, for joint models

# KO tables: abundance = NON-downsampled ; detection = downsampled (per QC report)
ko_abund <- read.delim("/proj/simp2024014/Omicsdataleverans2/metagenomics_CM/simpler_metagenomics_ko_relative_abundances_v4.0.tsv")
ko_prev  <- read.delim("/proj/simp2024014/Omicsdataleverans2/metagenomics_CM/simpler_metagenomics_ko_ds_relative_abundances_v4.0.tsv")
ko_annot <- read.delim("/proj/simp2024014/Omicsdataleverans2/metagenomics_CM/simpler_metagenomics_ko_annotations_v4.0.tsv")

# ==============================================================================
# 2. pap KO exposures  (per-gene detection + per-SD abundance; ALL pap genes)
# ==============================================================================
PAP_SYMBOLS <- c("papA","papB","papC","papD","papE","papF","papG","papH","papI","papJ","papK")
pap_lab <- ko_annot %>%
  mutate(gene = sub(";.*$", "", description)) %>%
  filter(gene %in% PAP_SYMBOLS) %>%
  select(ko_id, gene)
KO_PAP <- pap_lab$ko_id
print(pap_lab)   
# papG (K12522) absent from CM KEGG profile
#    ko_id gene
# 1 K12517 papA
# 2 K12518 papC
# 3 K12519 papD
# 4 K12520 papE
# 5 K12521 papF
# 6 K12523 papK

ko_matrix <- function(tbl, kos, fun, suffix) {
  m <- tbl %>% select(SIMPKEY, any_of(kos)) %>% mutate(across(-SIMPKEY, fun))
  names(m)[-1] <- paste0(names(m)[-1], suffix); m
}
ko_det   <- ko_matrix(ko_prev,  KO_PAP, ~ as.integer(.x > 0),               "_det")    # detection: downsampled
ko_abz   <- ko_matrix(ko_abund, KO_PAP, ~ as.numeric(scale(log10(.x + 1))), "_abz")    # per-SD abundance: non-ds
ko_ablog <- ko_matrix(ko_abund, KO_PAP, ~ log10(.x + 1),                    "_ablog")  # unscaled (scatter/corr)

ko_det %>% select(-SIMPKEY) %>% colSums()
# K12517_det K12518_det K12519_det K12520_det K12521_det K12523_det 
#       1535       3682       1600        227       1331       1186 

H <- health %>%
  reduce(list(ko_det, ko_abz, ko_ablog), left_join, by = "SIMPKEY", .init = .) %>%
  mutate(MEN = factor(MEN), PREVAL_E10E14 = factor(PREVAL_E10E14),
         aliquoting_plate = factor(aliquoting_plate))

gene_of <- function(k) pap_lab$gene[match(k, pap_lab$ko_id)]

# Remove prevalent cases
D_INC <-  H %>% filter(PREVAL_LIVERDIS == 0)

# ==============================================================================
# 3. MODELS  (incident Cox; detection AND abundance)
#   PRIMARY: coxph, plate NOT adjusted, PH reported (cox.zph)
#   SENSITIVITY: coxme, plate as a random intercept
#   BASE vs ADJUSTED covariate sets; FULL coefficient tables (FINRISK-style)
#   FDR(BH) within each measure x engine. Diabetes covariate = PREVAL_E10E14
# ==============================================================================
library(broom)

COV_BASE <- "BL_AGE + MEN + BMI + alcohol_gweek + antibiotic_J01_6mo + ReadDepthMillions + PREVAL_E10E14"
COV_ADJ  <- paste(COV_BASE, "+ statin_6mo + ppi_6mo + factor(education) + fibre + current_smoker")

# covariates the ADJUSTED set adds (for the complete-case subset)
adj <- c("BL_AGE", "MEN", "BMI", "alcohol_gweek", "antibiotic_J01_6mo", "ReadDepthMillions", "PREVAL_E10E14")
adj_extra <- c("statin_6mo", "ppi_6mo", "education", "fibre", "current_smoker")

# complete-case sample 
D_INC_cc <- D_INC %>% filter(if_all(all_of(c(adj, adj_extra)), ~ !is.na(.)))
D_INC <- D_INC %>% filter(if_all(all_of(adj), ~ !is.na(.)))

cat(sprintf("\nincident risk set: base %d | adjusted complete-case %d | cases %d -> %d\n", nrow(D_INC), nrow(D_INC_cc), sum(D_INC$INCIDENT_LIVERDIS == 1), sum(D_INC_cc$INCIDENT_LIVERDIS == 1)))
#incident risk set: base 5936 | adjusted complete-case 5656 | cases 41 -> 37

# ---- 3a. pap-term-only extractor (existing fit_inc, drives main results/FDR) ----
COV <- COV_BASE   # main results use the base model on the full risk set
MODELS <- list()
fit_inc <- function(expo_col, measure, gene, engine, cov = COV, data = D_INC) {
  rhs <- paste(expo_col, "+", cov)
  if (engine == "coxme") rhs <- paste(rhs, "+ (1 | aliquoting_plate)")
  fml <- as.formula(sprintf("Surv(LIVERDIS_AGEDIFF, INCIDENT_LIVERDIS) ~ %s", rhs))
  cc <- complete.cases(data[, all.vars(fml)]); N <- sum(cc)
  n_events <- sum(data$INCIDENT_LIVERDIS[cc] == 1)
  ph_term <- NA; ph_global <- NA
  if (engine == "coxme") {
    m <- tryCatch(coxme(fml, data = data), error = function(e) NULL); if (is.null(m)) return(NULL)
    beta <- fixef(m); if (!expo_col %in% names(beta)) return(NULL)
    se <- sqrt(diag(as.matrix(vcov(m)))); names(se) <- names(beta)
    b <- beta[expo_col]; s <- se[expo_col]
  } else {
    m <- tryCatch(coxph(fml, data = data), error = function(e) NULL); if (is.null(m)) return(NULL)
    co <- summary(m)$coefficients; if (!expo_col %in% rownames(co)) return(NULL)
    b <- co[expo_col, "coef"]; s <- co[expo_col, "se(coef)"]
    z <- tryCatch(cox.zph(m), error = function(e) NULL)
    if (!is.null(z)) { if (expo_col %in% rownames(z$table)) ph_term <- z$table[expo_col, "p"]
                       ph_global <- z$table["GLOBAL", "p"] }
  }
  zval <- b / s; p <- 2 * pnorm(-abs(zval))
  MODELS[[paste(engine, measure, gene, sep = "__")]] <<- m
  data.frame(engine, gene, measure, expo_col,
             HR = exp(b), lo = exp(b - 1.96 * s), hi = exp(b + 1.96 * s),
             p = p, N = N, n_events = n_events, ph_term = ph_term, ph_global = ph_global,
             row.names = NULL)
}

rows <- list()
for (eng in c("coxph", "coxme")) for (k in pap_lab$ko_id) {
  rows[[length(rows)+1]] <- fit_inc(paste0(k, "_det"), "detection", gene_of(k), eng)
  rows[[length(rows)+1]] <- fit_inc(paste0(k, "_abz"), "abundance", gene_of(k), eng)
}

res <- bind_rows(rows) %>%
  group_by(measure, engine) %>% mutate(FDR = p.adjust(p, method = "fdr")) %>% ungroup() %>%
  arrange(measure, engine, p)

write.csv(res, file.path(OUTDIR, "pap_results.csv"), row.names = FALSE)
write.csv(res, file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "pap_results.csv"), row.names = FALSE)

# ---- 3b. FULL coefficient tables (every covariate), coxph, base AND adjusted ----
# base on full risk set; adjusted on complete-case sample (so the two are comparable
# via a base-CC fit too). All three written out FINRISK-style.
fit_full <- function(expo_col, measure, gene, cov, data, model_lab) {
  fml <- as.formula(sprintf("Surv(LIVERDIS_AGEDIFF, INCIDENT_LIVERDIS) ~ %s + %s", expo_col, cov))
  m <- tryCatch(coxph(fml, data = data), error = function(e) NULL); if (is.null(m)) return(NULL)
  ph <- tryCatch(cox.zph(m)$table["GLOBAL", "p"], error = function(e) NA_real_)
  broom::tidy(m, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(gene, measure, model = model_lab, term,
              HR = estimate, lo = conf.low, hi = conf.high, p = p.value,
              ph_global = ph, N = m$n, n_events = m$nevent)
}

full_coef <- bind_rows(lapply(pap_lab$ko_id, function(k) bind_rows(
  fit_full(paste0(k,"_det"), "detection", gene_of(k), COV_BASE, D_INC,    "Base"),
  fit_full(paste0(k,"_abz"), "abundance", gene_of(k), COV_BASE, D_INC,    "Base"),
  fit_full(paste0(k,"_det"), "detection", gene_of(k), COV_BASE, D_INC_cc, "Base_CC"),
  fit_full(paste0(k,"_abz"), "abundance", gene_of(k), COV_BASE, D_INC_cc, "Base_CC"),
  fit_full(paste0(k,"_det"), "detection", gene_of(k), COV_ADJ,  D_INC_cc, "Adjusted"),
  fit_full(paste0(k,"_abz"), "abundance", gene_of(k), COV_ADJ,  D_INC_cc, "Adjusted")
)))
write.csv(full_coef, file.path(OUTDIR, "pap_models_full_coefficients.csv"), row.names = FALSE)
write.csv(full_coef, file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "pap_models_full_coefficients.csv"), row.names = FALSE)

# ---- 3c. pap-term robustness view: base(full) vs adjusted(CC), per gene/measure ----
robust <- full_coef %>%
  filter(grepl("_det$|_abz$", term)) %>%   # keep only the exposure row from each model
  mutate(est = sprintf("%.2f (%.2f-%.2f)", HR, lo, hi),
         P = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))) %>%
  filter(model %in% c("Base", "Adjusted")) %>%
  select(gene, measure, model, est, P, N, n_events) %>%
  pivot_wider(names_from = model, values_from = c(est, P, N, n_events))
print(as.data.frame(robust))

write.csv(robust, file.path(OUTDIR, "pap_robustness_base_vs_adjusted.csv"), row.names = FALSE)
write.csv(robust, file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "pap_robustness_base_vs_adjusted.csv"), row.names = FALSE)

# ==============================================================================
# 4. FIGURES
# ==============================================================================

# ---- FIGURE 1: pap-gene HR forest, DETECTION + ABUNDANCE (coxph, NO plate) --
# FINRISK Fig-1b style: blue = detection, red = abundance.
# Filled with the measure colour = FDR < 0.05; white (open ring) = n.s.
fdf <- res %>% filter(engine == "coxph") %>%
  mutate(expo_col = gsub("(K.*)_.*", "\\1", expo_col)) %>%
  mutate(measure = recode(measure, detection = "Detection", abundance = "Abundance"),
         measure = factor(measure, levels = c("Detection", "Abundance")),
         sig     = FDR < 0.05,
         gene = paste0(gene, " (", expo_col, ")"), 
         fillcat = ifelse(sig, as.character(measure), "ns"))   # colour-fill if sig, else white

# order y-axis by abundance HR (one value per gene; abundance is the primary measure)
ord <- fdf %>% filter(measure == "Abundance") %>% arrange(HR) %>% pull(gene)
ord <- unique(c(ord, setdiff(fdf$gene, ord)))
fdf <- fdf %>% mutate(gene = factor(gene, levels = ord))

cols      <- c("Detection" = "#00798C", "Abundance" = "#D1495B")
fill_cols <- c(cols, "ns" = "white")
dodge     <- position_dodge(width = 0.6)

p_forest <- ggplot(fdf, aes(HR, gene, colour = measure, group = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = fillcat), shape = 21, size = 2.6, stroke = 0.7, position = dodge) +
  scale_x_log10() +
  # coord_cartesian(xlim = c(0.5, 3)) +   # <- uncomment if a thin detection CI blows the axis
  scale_colour_manual(values = cols, name = "Measure") +
  scale_fill_manual(values = fill_cols, guide = "none") +
  guides(colour = guide_legend(override.aes = list(shape = 16, linetype = 0, size = 3))) +
  labs(x = "Hazard ratio (95% CI)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4),
        axis.line.x  = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks.x = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(size = 10, colour = "grey15"),
        axis.text.x  = element_text(colour = "grey15"),
        axis.title.x = element_text(colour = "grey15"),
        legend.position = "none")
ggsave(file.path(OUTDIR, "Forest_pap_det_abz.pdf"), p_forest, width = 6, height = 2.5)
ggsave(file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "Forest_pap_det_abz.pdf"), p_forest, width = 6, height = 2.5)

# ---- FIGURE 2: scatter of pap-gene abundance vs E. coli abundance -----------
# styled to match the FINRISK S5 scatter (theme_minimal, grey gridlines, red accent)
pap_short <- setNames(sapply(pap_lab$ko_id, gene_of), pap_lab$ko_id)
scat <- map_dfr(pap_lab$ko_id, function(k)
    D_INC %>% transmute(gene = unname(pap_short[k]), gene_ab = .data[[paste0(k, "_ablog")]], E_coli_log10)) %>%
  filter(E_coli_log10 > 0, !is.na(gene_ab)) %>%
  mutate(gene = factor(gene, levels = unname(pap_short)))

scat$gene %>% table()
# papA papC papD papE papF papK 
# 5676 5676 5676 5676 5676 5676 

scat %>%
filter(gene_ab != 0) %>%
group_by(gene) %>%
summarise(n=n())
#   gene      n
# 1 papA   1864
# 2 papC   3986
# 3 papD   1942
# 4 papE    273
# 5 papF   1627
# 6 papK   1450

stat_lab <- scat %>% group_by(gene) %>%
  summarise(rho   = cor(gene_ab, E_coli_log10, method = "spearman"),
            p_raw = cor.test(gene_ab, E_coli_log10, method = "spearman", exact = FALSE)$p.value,
            .groups = "drop") %>%
  mutate(p_txt = ifelse(p_raw < 0.0001, "< 0.0001", sprintf("= %.3f", p_raw)),
         label = sprintf("italic(rho)==%.2f*','~italic(P)*'%s'", rho, p_txt))

accent <- "#d7301f"   # red, matching FINRISK S5 (abundance story)

p_scat <- ggplot(scat, aes(E_coli_log10, gene_ab)) +
  geom_point(alpha = 0.18, size = 0.55, colour = "grey45", stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, colour = accent, fill = accent,
              alpha = 0.15, linewidth = 0.7) +
  geom_text(data = stat_lab, parse = TRUE, hjust = 0, vjust = 1,
            aes(x = -Inf, y = Inf, label = label),
            size = 2.8, colour = "grey15") +
  facet_wrap(~ gene, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = expression(log[10]("relative abundance of "*italic("E. coli")*" + 1")),
       y = expression(log[10]("relative abundance of "*italic("pap")*" gene + 1"))) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        axis.line   = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks  = element_line(colour = "grey20", linewidth = 0.4),
        axis.text   = element_text(colour = "grey15", size = 7.5),
        axis.title  = element_text(colour = "grey15"),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 8.5, hjust = 0, colour = "grey15"),
        panel.spacing    = unit(10, "pt"),
        plot.margin      = margin(6, 10, 6, 6))
                 
ggsave(file.path(OUTDIR, "Scatter_pap_vs_Ecoli.pdf"), p_scat, width = 7, height = 5.5, device = cairo_pdf)
ggsave(file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "Scatter_pap_vs_Ecoli.pdf"), p_scat, width = 7, height = 5.5, device = cairo_pdf)
stat_lab %>% select(gene, rho, p_raw) %>% print()

# ---- Within-carrier correlation ----------
# Drops the tied zeros (keeps only pap-carriers; E. coli already > 0 in scat).
scat %>%
  filter(gene_ab > 0) %>%
  group_by(gene) %>%
  summarise(n_within   = n(),
            rho_within = cor(gene_ab, E_coli_log10, method = "spearman"),
            p_within   = cor.test(gene_ab, E_coli_log10, method = "spearman",
                                   exact = FALSE)$p.value,
            .groups = "drop") %>%
  arrange(desc(rho_within))
#   gene  n_within rho_within  p_within
# 1 papF      1627      0.740 1.31e-282
# 2 papK      1450      0.735 7.09e-247
# 3 papC      3986      0.719 0        
# 4 papA      1864      0.686 1.01e-259
# 5 papD      1942      0.681 2.42e-264
# 6 papE       273      0.665 2.77e- 36

# ==============================================================================
# Characteristics table for the 5,936-participant functional-analysis sample
# ==============================================================================
metadata <- D_INC
dim(metadata)
cat("any-diabetes prevalence column joined; non-missing:", sum(!is.na(metadata$PREVAL_E10E14)), "of", nrow(metadata), "\n")

med_iqr <- function(x) {
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf("%g (%g, %g)", round(q[1], 1), round(q[2], 1), round(q[3], 1))
}
n_pct <- function(x, level = 1) {
  x <- x[!is.na(x)]
  n <- sum(as.integer(as.character(x)) == level)
  sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * n / length(x))
}
tab <- tribble(
  ~Characteristic, ~Value,
  "Age at baseline, years, median (IQR)",                              med_iqr(metadata$BL_AGE),
  "Men, n (%)",                                                        n_pct(metadata$MEN, level = 1),
  "BMI, kg/m2, median (IQR)",                                          med_iqr(metadata$BMI),
  "Antibacterial (ATC J01) use \u22646 mo before baseline, n (%)",     n_pct(metadata$antibiotic_J01_6mo, level = 1),
  "Alcohol consumption past year, g/week, median (IQR)",               med_iqr(metadata$alcohol_gweek),
  "Incident liver disease (ICD-10 K70-K77) over follow-up, n (%)",     n_pct(metadata$INCIDENT_LIVERDIS, level = 1),
  "Type 2 diabetes prevalence (ICD-10 E11), n (%)",                    n_pct(metadata$PREVAL_E11, level = 1),
  "Any diabetes prevalence (ICD-10 E10-E14), n (%)",                   n_pct(metadata$PREVAL_E10E14, level = 1),
  "ReadDepthMillions", med_iqr(metadata$ReadDepthMillions)
)
print(as.data.frame(tab), right = FALSE)
#   Characteristic                                                Value            
# 1 Age at baseline, years, median (IQR)                          73.4 (71.1, 76.2)
# 2 Men, n (%)                                                    3,012 (50.7%)    
# 3 BMI, kg/m2, median (IQR)                                      26.1 (23.7, 28.7)
# 4 Antibacterial (ATC J01) use ≤6 mo before baseline, n (%)      711 (12.0%)      
# 5 Alcohol consumption past year, g/week, median (IQR)           36.6 (11.1, 74.7)
# 6 Incident liver disease (ICD-10 K70-K77) over follow-up, n (%) 41 (0.7%)        
# 7 Type 2 diabetes prevalence (ICD-10 E11), n (%)                484 (8.2%)       
# 8 Any diabetes prevalence (ICD-10 E10-E14), n (%)               517 (8.7%)       
# 9 ReadDepthMillions                                             52.2 (34.5, 61.6)

# FOLLOW-UP TIME (years), median (IQR) 
metadata %>%
  summarise(
    fu_median = median(LIVERDIS_AGEDIFF, na.rm = TRUE),
    fu_q1     = quantile(LIVERDIS_AGEDIFF, 0.25, na.rm = TRUE),
    fu_q3     = quantile(LIVERDIS_AGEDIFF, 0.75, na.rm = TRUE)
  )
#   fu_median    fu_q1    fu_q3
# 1  6.592745 4.676249 8.438056

# Check on 
metadata$K70K77_ICD_codes %>% table() %>% sort()
# K703K720     K719     K720 K729K746     K740     K751     K754 K754K746     K759     K766     K703     K743     K768     K750     K760     K769     K729     K746 
#        1        1        1        1        1        1        1        1        1        1        2        2        2        3        3        4        6        9 
```