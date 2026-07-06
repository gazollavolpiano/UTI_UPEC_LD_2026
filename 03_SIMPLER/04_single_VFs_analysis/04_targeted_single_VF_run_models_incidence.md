# Targeted single-VF replication in SIMPLER of the 22 FINRISK-significant VFs (HH/LL subset)

Same Cox machinery / covariates / outcome as 02_02_run_models_incidence.md, but each 22 VF is tested individually, two ways:
- detection: binary presence/absence   -> VFs with 5% <= prevalence < 95%
- abundance: log10(RPKM + 1) (0 for non-carriers) -> VFs with prevalence >= 5%

Scenarios per VF (as in 02_02's table): base | base_cc | adjusted.
BH FDR is computed WITHIN each arm x scenario. Primary = base.

CAVEAT: the HH/LL incident risk set has only 30 events, adjusted model has very low events-per-variable, so adjusted per-VF estimates are unstable.

```R
# ---- 0. Libraries & options --------------------------------------------------
library(tidyverse)
library(survival)
options(width = 220, stringsAsFactors = FALSE)

### SIMPLER: SIMPLER workspace + a REAL wharf path (the FINRISK script left WHARF undefined)
OUTDIR <- "/home/cgazolla/Liver_Coli_2026/04_single_VFs_analysis"
WHARF  <- "/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014"
dir.create(OUTDIR, showWarnings = FALSE)

# ---- 1. Load data ------------------------------------------------------------
# host data
load("/home/cgazolla/Liver_Coli_2026/02_VF_aggregated_analysis/SIMPLER_d0_VF_modelling.RData")  # d0
dim(d0) # 6147 30

# VF-level hits keyed by SampleID -> map to SIMPKEY via the metagenomics key
load("/proj/simp2024014/VF_profiles/vfdb_hits_normalized_VF_level.RData")
vf_hits <- diamond_final; rm(diamond_final)
dim(vf_hits) # 102578 8

dic <- read.delim("/proj/simp2024014/Omicsdataleverans2/simpler_metagenomics_key_fastaq_file_v3.0.tsv") %>%
  mutate(Sample = gsub(".*host_removal/(.*)__[1|2].fq.gz", "\\1", Files)) %>%
  select(Sample, SIMPKEY) %>% distinct()
vf_hits <- left_join(vf_hits, dic, by = c("SampleID" = "Sample"))

# ---- 1b. The pre-specified 22-VF target set (FINRISK base FDR<0.05) -----------
TARGET_VFS <- c(
   "VF0213","VF0394","VF1110","VF0333","VF1138","VF0565","VF0404","VF0228",
   "VF0568","VF0572","VF0227","VF0256","VF0112","VF0571","VF0111","VF0560",
   "VF0221","VF0238","VF0113","VF0239","VF0236","VF0237")   # 22 VFs

# ---- 2. Config ---------------------------------------------------------------
D0_KEY   <- "SIMPKEY" 
VF_KEY   <- "SIMPKEY"
VF_READS <- "NumberOfHighHomologyReads"
VF_RPKM  <- "RPKM"

# SIMPLER: covariates as in 02_02 — NO plate in fixed effects (150 levels)
COV_BASE <- "BL_AGE + MEN + BMI + alcohol_gweek + antibiotic_J01_6mo + ReadDepthMillions + PREVAL_E10E14"
COV_ADJ  <- paste(COV_BASE, "+ statin_6mo + ppi_6mo + factor(education) + fibre + current_smoker")

OUT <- list(time = "LIVERDIS_AGEDIFF", event = "INCIDENT_LIVERDIS", restrict = "PREVAL_LIVERDIS == 0")

PREV_LO <- 0.05    # detection-arm lower prevalence bound
PREV_HI <- 0.95    # detection-arm upper prevalence bound

covars_all <- unique(c(all.vars(as.formula(paste("~", COV_BASE))),
                       all.vars(as.formula(paste("~", COV_ADJ)))))
setdiff(covars_all, names(d0)) # good to go

# ---- 3. Analytic sample + per-VF matrices restricted to the 22 targets --------
base_covars <- all.vars(as.formula(paste("~", COV_BASE)))
dd0 <- d0 %>%
  filter(group2 %in% c("HH", "LL")) %>%
  filter(eval(parse(text = OUT$restrict))) %>%
  drop_na(any_of(base_covars))
ids <- dd0[[D0_KEY]]
cat(sprintf("HH/LL incident risk set: %d participants, %d incident cases\n", length(ids), sum(dd0[[OUT$event]] == 1)))
# HH/LL incident risk set: 4548 participants, 30 incident cases

vf_meta <- vf_hits %>% ungroup() %>% distinct(VFID, VF_Name, VFcategory)

vf_long <- vf_hits %>%
  rename(.key = all_of(VF_KEY), .reads = all_of(VF_READS), .rpkm = all_of(VF_RPKM)) %>%
  filter(.key %in% ids, VFID %in% TARGET_VFS) %>%          ### SIMPLER: keep only the 22 targets
  group_by(.key, VFID) %>%
  summarise(Reads = sum(.reads), RPKM = sum(.rpkm), .groups = "drop") %>%
  filter(Reads > 0)

det_mat <- vf_long %>% transmute(.key, VFID, v = 1L) %>%
  pivot_wider(names_from = VFID, values_from = v, values_fill = 0) %>%
  column_to_rownames(".key") %>% as.matrix()
abz_mat <- vf_long %>% transmute(.key, VFID, v = log10(RPKM + 1)) %>%
  pivot_wider(names_from = VFID, values_from = v, values_fill = 0) %>%
  column_to_rownames(".key") %>% as.matrix()

# rows: zero-fill participants with no target VF detected
align_fill <- function(mat, ids) {
  miss <- setdiff(ids, rownames(mat))
  if (length(miss)) mat <- rbind(mat, matrix(0, length(miss), ncol(mat),
                                             dimnames = list(miss, colnames(mat))))
  mat[ids, , drop = FALSE]
}
# cols: ensure ALL 22 targets are present (a target absent from SIMPLER -> all-zero column)
ensure_cols <- function(mat, cols) {
  add <- setdiff(cols, colnames(mat))
  if (length(add)) mat <- cbind(mat, matrix(0, nrow(mat), length(add),
                                           dimnames = list(rownames(mat), add)))
  mat[, cols, drop = FALSE]
}
det_mat <- ensure_cols(align_fill(det_mat, ids), TARGET_VFS)
abz_mat <- ensure_cols(align_fill(abz_mat, ids), TARGET_VFS)
stopifnot(identical(rownames(det_mat), ids), identical(colnames(det_mat), TARGET_VFS))

# ---- 4. Arms within the targeted set -----------------------------------------
prevalence <- colMeans(det_mat)                                  # over the SIMPLER HH/LL subset
FEAT_ABZ <- TARGET_VFS                                           # all 22 tested in abundance (pre-specified)
FEAT_DET <- names(prevalence)[prevalence >= PREV_LO & prevalence < PREV_HI]   # detection only where contrast exists
cat(sprintf("Targets by arm -> abundance: %d (all) | detection: %d (5-95%% prevalent in SIMPLER)\n",
            length(FEAT_ABZ), length(FEAT_DET)))
#Targets by arm -> abundance: 22 (all) | detection: 22 (5-95% prevalent in SIMPLER)

det_df <- as.data.frame(det_mat); det_df[[D0_KEY]] <- rownames(det_mat)
abz_df <- as.data.frame(abz_mat); abz_df[[D0_KEY]] <- rownames(abz_mat)
names(det_df)[names(det_df) != D0_KEY] <- paste0("DET__", names(det_df)[names(det_df) != D0_KEY])
names(abz_df)[names(abz_df) != D0_KEY] <- paste0("ABZ__", names(abz_df)[names(abz_df) != D0_KEY])
dd0 <- dd0 %>% left_join(det_df, by = D0_KEY) %>% left_join(abz_df, by = D0_KEY)

# ---- 5. Fitting (identical to FINRISK) ---------------------------------------
fit_one_vf <- function(dd, expo, cov_fit, cov_cc, want_sd = FALSE) {
  dd$expo  <- expo
  form_cc  <- as.formula(sprintf("Surv(%s, %s) ~ expo + %s", OUT$time, OUT$event, cov_cc))
  form_fit <- as.formula(sprintf("Surv(%s, %s) ~ expo + %s", OUT$time, OUT$event, cov_fit))
  ddc <- dd[complete.cases(dd[, all.vars(form_cc)]), ]
  n <- nrow(ddc); ev <- sum(ddc[[OUT$event]] == 1)
  na_out <- list(est=NA,lo=NA,hi=NA,p=NA,est_sd=NA,lo_sd=NA,hi_sd=NA,ph=NA,phg=NA,n=n,cases=ev)
  if (n == 0 || length(unique(ddc$expo)) < 2) return(na_out)
  fit <- tryCatch(coxph(form_fit, data = ddc), error = function(e) NULL)
  if (is.null(fit)) return(na_out)
  s <- summary(fit)$coefficients
  if (!"expo" %in% rownames(s)) return(na_out)
  b <- s["expo","coef"]; se <- s["expo","se(coef)"]; p <- s["expo","Pr(>|z|)"]
  z   <- tryCatch(cox.zph(fit), error = function(e) NULL)
  ph  <- if (!is.null(z) && "expo" %in% rownames(z$table)) z$table["expo","p"] else NA
  phg <- if (!is.null(z)) z$table["GLOBAL","p"] else NA
  o <- list(est = exp(b), lo = exp(b-1.96*se), hi = exp(b+1.96*se), p = p,
            est_sd = NA, lo_sd = NA, hi_sd = NA, ph = ph, phg = phg, n = n, cases = ev)
  if (want_sd) {
    sdv <- sd(ddc$expo)
    o$est_sd <- exp(b*sdv); o$lo_sd <- exp((b-1.96*se)*sdv); o$hi_sd <- exp((b+1.96*se)*sdv)
  }
  o
}

SCEN <- list(base    = list(fit = COV_BASE, cc = COV_BASE),
             base_cc = list(fit = COV_BASE, cc = COV_ADJ),
             adjusted= list(fit = COV_ADJ,  cc = COV_ADJ))

rows <- list()
for (v in TARGET_VFS) {
  for (arm in c("detection", "abundance")) {
    if (arm == "detection" && !v %in% FEAT_DET) next
    if (arm == "abundance" && !v %in% FEAT_ABZ) next
    expo    <- dd0[[paste0(if (arm == "detection") "DET__" else "ABZ__", v)]]
    want_sd <- arm == "abundance"
    for (sc in names(SCEN)) {
      r <- fit_one_vf(dd0, expo, SCEN[[sc]]$fit, SCEN[[sc]]$cc, want_sd = want_sd)
      rows[[length(rows)+1]] <- data.frame(
        VFID = v, arm = arm, scenario = sc,
        cases = r$cases, n = r$n, prevalence = prevalence[[v]],
        HR = r$est, lo = r$lo, hi = r$hi, p = r$p,
        HR_perSD = r$est_sd, lo_perSD = r$lo_sd, hi_perSD = r$hi_sd,
        ph_term = r$ph, ph_global = r$phg, row.names = NULL)
    }
  }
}
res <- do.call(rbind, rows) %>%
  left_join(vf_meta, by = "VFID") %>%
  relocate(VF_Name, VFcategory, .after = VFID)

# ---- 6. BH FDR within arm x scenario (over the 22 targets); write ------------
res$FDR <- NA_real_
for (a in unique(res$arm)) for (s in unique(res$scenario)) {
  idx <- which(res$arm == a & res$scenario == s & !is.na(res$p))
  res$FDR[idx] <- p.adjust(res$p[idx], "BH")
}
res <- res %>% arrange(arm, scenario, p)

res <- res %>% mutate(VF_Name = recode(VF_Name, "Allantion utilization" = "Allantoin utilization"))

write.csv(res, file.path(OUTDIR, "simpler_singleVF_replication_full.csv"), row.names = FALSE)
write.csv(res, file.path(WHARF, "simpler_singleVF_replication_full.csv"), row.names = FALSE)

res %>% arrange(FDR) %>% head() # nothing is significant

# ---- 7. Forest — ALL 22 targets, base scenario (filled = replicates in SIMPLER) ----
# show every target so non-replicating VFs are visible as open points
fdf <- res %>%
  filter(scenario == "base", VFID %in% TARGET_VFS) %>%
  mutate(
    HRx = ifelse(arm == "abundance", HR_perSD, HR),     # abundance -> per-SD ; detection -> binary
    lox = ifelse(arm == "abundance", lo_perSD, lo),
    hix = ifelse(arm == "abundance", hi_perSD, hi),
    measure = recode(arm, detection = "Detection", abundance = "Abundance"),
    measure = factor(measure, levels = c("Detection", "Abundance")),
    sig     = !is.na(FDR) & FDR < 0.05,
    gene    = paste0(VF_Name, " (", VFID, ")"),
    fillcat = ifelse(sig, as.character(measure), "ns")
  ) %>%
  filter(!is.na(HRx))

# order by ABUNDANCE HR (primary in SIMPLER); detection-only VFs appended
ord_abz <- fdf %>% filter(measure == "Abundance") %>% arrange(HRx) %>% pull(gene)
ord_det <- fdf %>% filter(measure == "Detection", !gene %in% ord_abz) %>% arrange(HRx) %>% pull(gene)
fdf <- fdf %>% mutate(gene = factor(gene, levels = unique(c(ord_det, ord_abz))))

cols      <- c("Detection" = "#00798C", "Abundance" = "#D1495B")
fill_cols <- c(cols, "ns" = "white")
dodge     <- position_dodge(width = 0.6)

p_forest <- ggplot(fdf, aes(HRx, gene, colour = measure, group = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lox, xmax = hix), height = 0.25, linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = fillcat), shape = 21, size = 2.6, stroke = 0.7, position = dodge) +
  scale_x_log10() +
  # coord_cartesian(xlim = c(0.5, 3)) +   # <- uncomment if a wide CI blows the axis
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
        axis.text.y  = element_text(size = 9, colour = "grey15"),
        axis.text.x  = element_text(colour = "grey15"),
        axis.title.x = element_text(colour = "grey15"),
        plot.caption = element_text(size = 7, colour = "grey40", hjust = 0),
        legend.position = "top")

n_gene <- nlevels(fdf$gene)
H <- max(3, 0.32 * n_gene + 1)
ggsave(file.path(OUTDIR, "Forest_singleVF_replication_SIMPLER.pdf"), p_forest, width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Forest_singleVF_replication_SIMPLER.svg"), p_forest, width = 6.5, height = 5)

ggsave(file.path(WHARF, "Forest_singleVF_replication_SIMPLER.pdf"), p_forest, width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(WHARF, "Forest_singleVF_replication_SIMPLER.svg"), p_forest, width = 6.5, height = 5)

# ==============================================================================
# 8. NEW: AGGREGATE 22-VF EXPOSURES — two single tests:
#      (i)  SUMMED ABUNDANCE  : log10(total RPKM summed over the 22 + 1)  -> dose
#      (ii) ANY DETECTION     : carries >=1 of the 22 vs carries none     -> binary
# ==============================================================================

# ---- analytic sample --------
dd <- d0 %>%
  filter(eval(parse(text = OUT$restrict))) %>%
  drop_na(any_of(base_covars))
ids <- dd[[D0_KEY]]
cat(sprintf("Sample: %d participants, %d incident cases\n", length(ids), sum(dd[[OUT$event]] == 1)))
# Sample: 5936 participants, 41 incident cases

# ---- build the two aggregate exposures over the 22 targets -------------------
agg22 <- vf_hits %>%
  rename(.key = all_of(VF_KEY), .reads = all_of(VF_READS), .rpkm = all_of(VF_RPKM)) %>%
  filter(.key %in% ids, VFID %in% TARGET_VFS) %>%
  group_by(.key, VFID) %>%
  summarise(RPKM = sum(.rpkm), Reads = sum(.reads), .groups = "drop") %>%
  filter(Reads > 0) %>% 
  group_by(.key) %>%
  summarise(sum_rpkm = sum(RPKM), n_det = n_distinct(VFID), .groups = "drop")

scores <- tibble(.key = ids) %>%
  left_join(agg22, by = ".key") %>%
  mutate(
    sum_rpkm = replace_na(sum_rpkm, 0),
    n_det    = replace_na(n_det, 0),
    sumab22  = log10(sum_rpkm + 1), # (i) summed abundance of the 22 (log10 RPKM)
    any22    = as.integer(n_det > 0) # (ii) detection of ANY of the 22 (>=1 vs 0)
  ) %>%
  rename(!!D0_KEY := ".key")

dd <- dd %>% left_join(scores %>% select(all_of(D0_KEY), sumab22, any22), by = D0_KEY)

cat(sprintf("\nCarriers of >=1 of the 22: %d (%.1f%%) | non-carriers: %d\n", sum(dd$any22 == 1), 100*mean(dd$any22 == 1), sum(dd$any22 == 0)))
# Carriers of >=1 of the 22: 4022 (67.8%) | non-carriers: 1914

cat("summed-abundance (log10) distribution:\n"); print(summary(dd$sumab22[dd$any22 == 1]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009822 0.346322 1.145842 1.264715 1.942757 3.920230 

# ---- fit each exposure x 3 scenarios -----------------------------------------
EXPO <- list(
  list(lab = "Summed abundance of 22 VFs (log10 RPKM)", var = "sumab22", sd = TRUE),   # report per-SD
  list(lab = "Detection of ANY of 22 VFs (>=1 vs 0)",   var = "any22",   sd = FALSE)   # binary HR
)

rows <- list()
for (e in EXPO) for (sc in names(SCEN)) {
  r <- fit_one_vf(dd, dd[[e$var]], SCEN[[sc]]$fit, SCEN[[sc]]$cc, want_sd = e$sd)
  rows[[length(rows)+1]] <- data.frame(
    exposure = e$lab, scenario = sc, cases = r$cases, n = r$n,
    HR = r$est, lo = r$lo, hi = r$hi,                    # binary: present-vs-absent ; abundance: per log10
    HR_perSD = r$est_sd, lo_perSD = r$lo_sd, hi_perSD = r$hi_sd,  # abundance only
    p = r$p, ph_term = r$ph, ph_global = r$phg, row.names = NULL)
}
agg22_res <- do.call(rbind, rows)
print(agg22_res, digits = 3) # read the BASE rows; abundance -> HR_perSD, detection -> HR
#                                  exposure scenario cases    n   HR    lo   hi HR_perSD lo_perSD hi_perSD     p ph_term ph_global
# 1 Summed abundance of 22 VFs (log10 RPKM)     base    41 5936 1.22 0.917 1.62     1.22    0.917     1.62 0.173   0.265     0.356
# 2 Summed abundance of 22 VFs (log10 RPKM)  base_cc    37 5656 1.18 0.869 1.59     1.18    0.869     1.59 0.293   0.220     0.405
# 3 Summed abundance of 22 VFs (log10 RPKM) adjusted    37 5656 1.18 0.872 1.60     1.18    0.872     1.59 0.284   0.227     0.282
# 4   Detection of ANY of 22 VFs (>=1 vs 0)     base    41 5936 1.54 0.726 3.25       NA       NA       NA 0.261   0.355     0.393
# 5   Detection of ANY of 22 VFs (>=1 vs 0)  base_cc    37 5656 1.34 0.626 2.88       NA       NA       NA 0.449   0.392     0.475
# 6   Detection of ANY of 22 VFs (>=1 vs 0) adjusted    37 5656 1.35 0.630 2.90       NA       NA       NA 0.439   0.383     0.328

# ==============================================================================
# Full model output: Summed abundance of 22 VFs (log10 RPKM, exposure rescaled to per-SD), within-carrier dose gradient (carriers only, sumab22 > 0)
# ==============================================================================
base_vars <- all.vars(as.formula(sprintf("Surv(%s, %s) ~ sumab22 + %s", OUT$time, OUT$event, COV_BASE)))

dd_fit <- dd %>%
  filter(sumab22 > 0) %>%                                   # within-carrier dose gradient
  filter(complete.cases(across(all_of(base_vars)))) %>%
  mutate(sumab22_sd = sumab22 / sd(sumab22))                # exposure -> per-SD (carriers-only SD)
form_base <- as.formula(sprintf("Surv(%s, %s) ~ sumab22_sd + %s", OUT$time, OUT$event, COV_BASE))
fit_base  <- coxph(form_base, data = dd_fit)

# --- full coefficient table; sumab22_sd row is now PER SD ---------------------
coef_tbl <- broom::tidy(fit_base, conf.int = TRUE, exponentiate = TRUE) %>%
  transmute(
    term,
    HR       = estimate,          # exposure: per 1 SD ; covariates: per unit / per level
    lo       = conf.low,
    hi       = conf.high,
    se_logHR = std.error,         # SE stays on the log-HR scale (no SE for a ratio)
    z        = statistic,
    p        = p.value
  ) %>%
  mutate(term = recode(term, sumab22_sd = "Summed abundance 22 VFs (per SD)"))

# --- PH diagnostics: exposure term + global (per-variable test) ----------------
zph <- cox.zph(fit_base)
out_tbl <- coef_tbl %>%
  mutate(model      = "summed_abundance_22VF_base_perSD_carriers",
         n          = fit_base$n,
         nevent     = fit_base$nevent,
         ph_exposure= zph$table["sumab22_sd", "p"],
         ph_global  = zph$table["GLOBAL", "p"], .before = 1)

print(out_tbl, digits = 3)
#   model                                         n nevent ph_exposure ph_global term                                HR    lo    hi se_logHR      z      p
# 1 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 Summed abundance 22 VFs (per SD) 1.18  0.836  1.67  0.176    0.942 0.346 
# 2 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 BL_AGE                           0.955 0.890  1.03  0.0359  -1.27  0.204 
# 3 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 MENMen                           0.798 0.369  1.72  0.393   -0.574 0.566 
# 4 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 BMI                              1.03  0.951  1.11  0.0395   0.676 0.499 
# 5 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 alcohol_gweek                    1.00  0.995  1.01  0.00343  0.411 0.681 
# 6 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 antibiotic_J01_6mo1              1.83  0.785  4.27  0.432    1.40  0.161 
# 7 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 ReadDepthMillions                0.992 0.974  1.01  0.00951 -0.820 0.412 
# 8 summed_abundance_22VF_base_perSD_carriers  4022     32      0.0195     0.347 PREVAL_E10E141                   2.59  1.13   5.92  0.423    2.25  0.0246

# ==============================================================================
# Boxplot: summed abundance of the 22 VFs (SD log10 RPKM) by liver-disease, within-carrier dose gradient
# ==============================================================================

# ---- per-participant table ---------------------------------------------------
pt <- dd %>%
  transmute(sumab22, case = .data[["INCIDENT_LIVERDIS" ]]) %>%
  mutate(grp = factor(ifelse(case == 1, "Liver\ndisease", "No liver\ndisease"), levels = c("No liver\ndisease", "Liver\ndisease")))

pt <- pt %>% filter(sumab22 > 0) # carriers only 

# SIMPLER adjusted Cox HR (base model)
hr_label <- "HR 1.18 (0.84\u20131.67), P = 0.3"

# n per group beneath each box
nlab <- pt %>% count(grp) %>% mutate(lab = paste0("n=", n))

case_col <- "#D1495B"; non_col <- "grey55"

p_box <- ggplot(pt, aes(grp, sumab22)) +
#  geom_violin(fill = "grey92", colour = NA, scale = "width", width = 0.9) +
  geom_jitter(data = ~filter(.x, case == 0), aes(colour = "No liver disease"), width = 0.16, size = 0.4, alpha = 0.10) +
  geom_jitter(data = ~filter(.x, case == 1), aes(colour = "Liver disease"), width = 0.08, size = 1.5, alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", colour = "grey25", linewidth = 0.5, alpha = 0.5) +
  geom_text(data = nlab, aes(x = grp, y = -Inf, label = lab), vjust = -0.6, size = 2.9, colour = "grey40", inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = Inf, label = hr_label, vjust = 1.6, size = 3.1, colour = "grey15") +
  scale_colour_manual(values = c("Liver disease" = case_col, "No liver disease" = non_col), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  labs(x = NULL, y = "Summed abundance of 22 VFs\nlog10(summed RPKM + 1)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor    = element_blank(),
        panel.grid.major.y  = element_line(colour = "grey88", linewidth = 0.4),
        axis.line.y  = element_line(colour = "grey20", linewidth = 0.4),
        axis.line.x  = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks.y = element_line(colour = "grey20", linewidth = 0.4),
        axis.text    = element_text(colour = "grey15"),
        axis.title.y = element_text(colour = "grey15"),
        legend.position = "none",
        plot.margin = margin(10, 14, 8, 10))
 
ggsave(file.path(OUTDIR, "Box_VF22_sumabund_by_case_SIMPLER.svg"), p_box, width = 2.5, height = 4.2)
ggsave(file.path(OUTDIR, "Box_VF22_sumabund_by_case_SIMPLER.pdf"), p_box, width = 2.5, height = 4.2, device = cairo_pdf)

ggsave(file.path(WHARF, "Box_VF22_sumabund_by_case_SIMPLER.svg"), p_box, width = 2.5, height = 4.2)
ggsave(file.path(WHARF, "Box_VF22_sumabund_by_case_SIMPLER.pdf"), p_box, width = 2.5, height = 4.2, device = cairo_pdf)
```