# Single-VF association with incident liver disease (FINRISK, HH/LL subset)

Same Cox machinery / covariates / outcome as 02_02_run_models_incidence.md, but each VF is tested individually, two ways:
- detection: binary presence/absence   -> VFs with 5% <= prevalence < 95%
- abundance: log10(RPKM + 1) (0 for non-carriers) -> VFs with prevalence >= 5%

Scenarios per VF (as in 02_02's table): base | base_cc | adjusted.
BH FDR is computed WITHIN each arm x scenario. Primary = base.

CAVEAT: the HH/LL incident risk set has only 47 events, adjusted model has very low events-per-variable, so adjusted per-VF estimates are unstable.

```R
# ---- 0. Libraries & options --------------------------------------------------
library(tidyverse)
library(survival)
options(width = 220, stringsAsFactors = FALSE)

OUTDIR <- "/home/camigazo/04_single_VFs_analysis"
dir.create(OUTDIR, showWarnings = FALSE)

# ---- 1. Load data ------------------------------------------------------------
load("/home/camigazo/02_VF_aggregated_analysis/FINRISK_d0_VF_modelling.RData")  # d0
dim(d0) # 7226 x 29

# gene-level VF hits, keyed by Barcode
load("/home/camigazo/data/vfdb_hits_normalized.RData")
vf_hits <- diamond_final; rm(diamond_final)
dim(vf_hits) # 1007163 x 11

# ---- 2. Config ---------------------------------------------------------------
D0_KEY   <- "Barcode"                   
VF_KEY   <- "SampleID"                    
VF_READS <- "NumberOfHighQualityReads"
VF_RPKM  <- "RPKM_HighQuality"

### FINRISK: covariate sets use `batch` (4 levels, already a factor in d0) 
COV_BASE <- "BL_AGE + MEN + BMI + ALKI2_FR02 + BL_USE_RX_J01 + ReadDepthMillions + PREVAL_DIAB + batch_name"
COV_ADJ  <- paste(COV_BASE, "+ BL_USE_RX_A02BC + BL_USE_RX_C10AA + KOULGR + FIBER_TOTAL + CURR_SMOKE")

OUT <- list(time = "LIVERDIS_AGEDIFF", event = "INCIDENT_LIVERDIS", restrict = "PREVAL_LIVERDIS == 0")

PREV_LO <- 0.05    # lower prevalence bound (both arms)
PREV_HI <- 0.95    # upper prevalence bound (detection arm only)

covars_all <- unique(c(all.vars(as.formula(paste("~", COV_BASE))),
                       all.vars(as.formula(paste("~", COV_ADJ)))))
setdiff(covars_all, names(d0))#ok

# ---- 3. Analytic sample + per-VF detection / abundance matrices --------------
# HH/LL incident risk set (FINRISK universe ~1,764), complete on BASE covariates
base_covars <- all.vars(as.formula(paste("~", COV_BASE)))
dd0 <- d0 %>%
  filter(group2 %in% c("HH", "LL")) %>%
  filter(eval(parse(text = OUT$restrict))) %>%
  drop_na(any_of(base_covars))                       
ids <- dd0[[D0_KEY]]
cat(sprintf("HH/LL incident risk set: %d participants, %d incident cases\n", length(ids), sum(dd0[[OUT$event]] == 1)))
# expect 1837 participants, 47 incident cases (FINRISK HH/LL)

# VF metadata (name/category per VFID) before aggregating
vf_meta <- vf_hits %>% ungroup() %>% distinct(VFID, VF_Name, VFcategory)

# Aggregate gene-level hits to VF level per sample (summed reads / RPKM), restrict to HH/LL
vf_long <- vf_hits %>%
  rename(.key = all_of(VF_KEY), .reads = all_of(VF_READS), .rpkm = all_of(VF_RPKM)) %>%
  filter(.key %in% ids) %>%
  group_by(.key, VFID) %>%
  summarise(Reads = sum(.reads), RPKM = sum(.rpkm), .groups = "drop") %>%
  filter(Reads > 0)

# wide detection (0/1) and abundance (log10(RPKM+1)) matrices
det_mat <- vf_long %>% transmute(.key, VFID, v = 1L) %>%
  pivot_wider(names_from = VFID, values_from = v, values_fill = 0) %>%
  column_to_rownames(".key") %>% as.matrix()
abz_mat <- vf_long %>% transmute(.key, VFID, v = log10(RPKM + 1)) %>%
  pivot_wider(names_from = VFID, values_from = v, values_fill = 0) %>%
  column_to_rownames(".key") %>% as.matrix()

# align to the full analytic roster; participants with no VF row -> all-zero
align_fill <- function(mat, ids) {
  miss <- setdiff(ids, rownames(mat))
  if (length(miss)) mat <- rbind(mat, matrix(0, length(miss), ncol(mat),
                                             dimnames = list(miss, colnames(mat))))
  mat[ids, , drop = FALSE]
}
det_mat <- align_fill(det_mat, ids)
abz_mat <- align_fill(abz_mat, ids)
stopifnot(identical(rownames(det_mat), ids), identical(colnames(det_mat), colnames(abz_mat)))

# ---- 4. Prevalence filters -> feature sets -----------------------------------
prevalence <- colMeans(det_mat)                                  # over the FINRISK HH/LL subset
FEAT_ABZ <- names(prevalence)[prevalence >= PREV_LO]             # abundance arm
FEAT_DET <- names(prevalence)[prevalence >= PREV_LO & prevalence < PREV_HI]  # detection arm (primary in FINRISK)
cat(sprintf("VFs tested -> detection: %d (%.0f-%.0f%%) | abundance: %d (>=%.0f%%)\n", length(FEAT_DET), 100*PREV_LO, 100*PREV_HI, length(FEAT_ABZ), 100*PREV_LO))
#VFs tested -> detection: 33 (5-95%) | abundance: 33 (>=5%)

# attach exposure columns to dd0 (prefixed; CC masks below only ever touch `expo` + covars)
det_df <- as.data.frame(det_mat); det_df[[D0_KEY]] <- rownames(det_mat)
abz_df <- as.data.frame(abz_mat); abz_df[[D0_KEY]] <- rownames(abz_mat)
names(det_df)[names(det_df) != D0_KEY] <- paste0("DET__", names(det_df)[names(det_df) != D0_KEY])
names(abz_df)[names(abz_df) != D0_KEY] <- paste0("ABZ__", names(abz_df)[names(abz_df) != D0_KEY])
dd0 <- dd0 %>% left_join(det_df, by = D0_KEY) %>% left_join(abz_df, by = D0_KEY)

# ---- 5. Fitting --------------------------------------------------------------
# cov_cc defines the complete-case sample; cov_fit defines the fitted model.
#   base     : cov_cc = cov_fit = COV_BASE                 (base covars, full CC)
#   base_cc  : cov_cc = COV_ADJ, cov_fit = COV_BASE        (base covars, adjusted CC sample)
#   adjusted : cov_cc = cov_fit = COV_ADJ                  (adjusted covars, adjusted CC sample)
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
  if (want_sd) {                                # per-SD-of-this-VF rescaling (abundance arm)
    sdv <- sd(ddc$expo)
    o$est_sd <- exp(b*sdv); o$lo_sd <- exp((b-1.96*se)*sdv); o$hi_sd <- exp((b+1.96*se)*sdv)
  }
  o
}

SCEN <- list(base    = list(fit = COV_BASE, cc = COV_BASE),
             base_cc = list(fit = COV_BASE, cc = COV_ADJ),
             adjusted= list(fit = COV_ADJ,  cc = COV_ADJ))

rows <- list()
all_vfs <- union(FEAT_DET, FEAT_ABZ)
for (i in seq_along(all_vfs)) {
  v <- all_vfs[i]
  if (i %% 25 == 0) cat(sprintf("  ... %d / %d VFs\n", i, length(all_vfs)))
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
        HR = r$est, lo = r$lo, hi = r$hi, p = r$p,             # binary HR (det) or per-log10 (abz)
        HR_perSD = r$est_sd, lo_perSD = r$lo_sd, hi_perSD = r$hi_sd,
        ph_term = r$ph, ph_global = r$phg, row.names = NULL)
    }
  }
}
res <- do.call(rbind, rows) %>%
  left_join(vf_meta, by = "VFID") %>%
  relocate(VF_Name, VFcategory, .after = VFID)

res %>% filter(p < 0.05) %>% head(2)
res %>% filter(p < 0.05) %>% pull(VF_Name) %>% unique()

# ---- 6. BH FDR within arm x scenario; write ----------------------------------
res$FDR <- NA_real_
for (a in unique(res$arm)) for (s in unique(res$scenario)) {
  idx <- which(res$arm == a & res$scenario == s & !is.na(res$p))
  res$FDR[idx] <- p.adjust(res$p[idx], "BH")
}
res <- res %>% arrange(arm, scenario, p)

res %>% filter(FDR < 0.05) %>% filter(scenario=="base") %>% pull(VF_Name) %>% unique()
res <- res %>% mutate(VF_Name = recode(VF_Name, "Allantion utilization" = "Allantoin utilization"))

write.csv(res, file.path(OUTDIR, "finrisk_singleVF_incident_full.csv"), row.names = FALSE)

# list any significant VF to be tested 
vf_list_any <- res %>% filter(FDR < 0.05) %>% pull(VFID) %>% unique() 
write.csv(vf_list_any, file.path(OUTDIR, "finrisk_25singleVF_significant.csv"), row.names = FALSE)
vf_list_any
#[1] "VF0394" "VF0333" "VF0568" "VF0227" "VF0404" "VF1138" "VF0228" "VF0560" 
# "VF0571" "VF0112" "VF0565" "VF0239" "VF0213" "VF0111" "VF0572" "VF0256"
# "VF0221" "VF0229" "VF0220" "VF0236" "VF0237" "VF0238" "VF1110" "VF0113"
#[25] "VF0144"

base <- res %>% filter(scenario == "base")

# any-arm significant (union) and per-arm
sig_any  <- base %>% filter(FDR < 0.05) %>% pull(VFID) %>% unique()
sig_det  <- base %>% filter(arm == "detection", FDR < 0.05) %>% pull(VFID) %>% unique()
sig_abz  <- base %>% filter(arm == "abundance", FDR < 0.05) %>% pull(VFID) %>% unique()

length(sig_any) # 22 total distinct VFs 
length(sig_abz); length(sig_det)      # 22 in abundance / 17 in detection
length(intersect(sig_det, sig_abz))   # 17 significant in BOTH
length(setdiff(sig_abz, sig_det))     # 5 abundance-only
length(setdiff(sig_det, sig_abz))     # 0 detection-only

# VFs that become significant ONLY in the adjusted model (added covariates)
adj <- res %>% filter(scenario == "adjusted")
sig_adj_any <- adj %>% filter(FDR < 0.05) %>% pull(VFID) %>% unique()
adj_only <- setdiff(sig_adj_any, sig_any) # new under adjustment
res %>% filter(VFID %in% adj_only) %>% distinct(VFID, VF_Name) # names for the text
#VFID    VF_Name
#1 VF0229 Aerobactin
#2 VF0220 P fimbriae
#3 VF0144    Capsule

# ---- 7. single-VF HR forest — DETECTION + ABUNDANCE (base scenario only) ----------------------------------

# ---- rows to show -----------------------------------------------------------
# Default: VFs reaching FDR < 0.05 in EITHER base arm (the liver-disease-associated set)
sig_vfs <- res %>%
  filter(scenario == "base", FDR < 0.05) %>%
  pull(VFID) %>% unique()

# ---- assemble plot frame ----------------------------------------------------
# Each VF contributes up to two points: detection (binary HR) and abundance (per-SD HR)
fdf <- res %>%
  filter(scenario == "base", VFID %in% sig_vfs) %>%
  mutate(
    HRx = ifelse(arm == "abundance", HR_perSD, HR),     # abundance -> per-SD ; detection -> binary
    lox = ifelse(arm == "abundance", lo_perSD, lo),
    hix = ifelse(arm == "abundance", hi_perSD, hi),
    measure = recode(arm, detection = "Detection", abundance = "Abundance"),
    measure = factor(measure, levels = c("Detection", "Abundance")),
    sig     = !is.na(FDR) & FDR < 0.05,
    gene    = paste0(VF_Name, " (", VFID, ")"),
    fillcat = ifelse(sig, as.character(measure), "ns")   # fill colour if sig, else white
  ) %>%
  filter(!is.na(HRx)) # drop arms that didn't converge

# ---- y-axis order: by DETECTION HR (primary); abundance-only VFs appended ----
ord_det <- fdf %>% filter(measure == "Detection") %>% arrange(HRx) %>% pull(gene)
ord_abz <- fdf %>% filter(measure == "Abundance", !gene %in% ord_det) %>% arrange(HRx) %>% pull(gene)
fdf <- fdf %>% mutate(gene = factor(gene, levels = unique(c(ord_abz, ord_det))))  # high HR at top

# ---- draw -------------------------------------------------------------------
cols      <- c("Detection" = "#00798C", "Abundance" = "#D1495B")
fill_cols <- c(cols, "ns" = "white")
dodge     <- position_dodge(width = 0.6)

p_forest <- ggplot(fdf, aes(HRx, gene, colour = measure, group = measure)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lox, xmax = hix), height = 0.25, linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = fillcat), shape = 21, size = 2.6, stroke = 0.7, position = dodge) +
  scale_x_log10() +
  # coord_cartesian(xlim = c(0.3, 5)) +   # <- uncomment if a wide detection CI blows the axis
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

# height scales with the number of VFs
n_gene <- nlevels(fdf$gene)
H <- max(3, 0.32 * n_gene + 1)
ggsave(file.path(OUTDIR, "Forest_singleVF_det_abz_FINRISK.pdf"), p_forest, width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Forest_singleVF_det_abz_FINRISK.svg"), p_forest, width = 6.5, height = 5)    # editable text for Inkscape (svglite)

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
# Sample: 6735 participants, 133 incident cases

# ---- build the two aggregate exposures over the 22 targets -------------------
TARGET_VFS<-vf_list_any
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
# Carriers of >=1 of the 22: 1754 (26.0%) | non-carriers: 4981

cat("summed-abundance (log10) distribution:\n"); print(summary(dd$sumab22[dd$any22 == 1]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.3079  1.0074  1.3479  1.6548  2.2452  4.1509 

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
# exposure scenario cases    n   HR   lo   hi HR_perSD lo_perSD hi_perSD        p ph_term ph_global
# 1 Summed abundance of 22 VFs (log10 RPKM)     base   133 6735 1.41 1.20 1.66     1.33     1.16     1.53 4.37e-05   0.869     0.244
# 2 Summed abundance of 22 VFs (log10 RPKM)  base_cc   125 6267 1.41 1.19 1.67     1.33     1.16     1.54 6.09e-05   0.735     0.548
# 3 Summed abundance of 22 VFs (log10 RPKM) adjusted   125 6267 1.39 1.18 1.65     1.32     1.15     1.52 1.14e-04   0.705     0.888
# 4   Detection of ANY of 22 VFs (>=1 vs 0)     base   133 6735 1.74 1.20 2.52       NA       NA       NA 3.37e-03   0.716     0.240
# 5   Detection of ANY of 22 VFs (>=1 vs 0)  base_cc   125 6267 1.77 1.21 2.60       NA       NA       NA 3.15e-03   0.647     0.543
# 6   Detection of ANY of 22 VFs (>=1 vs 0) adjusted   125 6267 1.75 1.19 2.56       NA       NA       NA 4.08e-03   0.646     0.888

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
  mutate(#model      = "summed_abundance_22VF_base_perSD_carriers",
         n          = fit_base$n,
         nevent     = fit_base$nevent,
         ph_exposure= zph$table["sumab22_sd", "p"],
         ph_global  = zph$table["GLOBAL", "p"], .before = 1)

print(out_tbl, digits = 3)
# n nevent ph_exposure ph_global term                                 HR     lo    hi se_logHR       z            p
# 1  1754     47       0.589     0.507 Summed abundance 22 VFs (per SD)  1.47   1.13   1.92 0.136     2.84    0.00455    
# 2  1754     47       0.589     0.507 BL_AGE                            0.994  0.971  1.02 0.0117   -0.509   0.610      
# 3  1754     47       0.589     0.507 MENMen                            0.840  0.445  1.59 0.324    -0.537   0.591      
# 4  1754     47       0.589     0.507 BMI                               1.07   1.01   1.13 0.0285    2.39    0.0171     
# 5  1754     47       0.589     0.507 ALKI2_FR02                        1.00   1.00   1.00 0.000667  5.31    0.000000110
# 6  1754     47       0.589     0.507 BL_USE_RX_J011                    0.683  0.298  1.57 0.424    -0.898   0.369      
# 7  1754     47       0.589     0.507 ReadDepthMillions                 1.00   0.815  1.23 0.106     0.0323  0.974      
# 8  1754     47       0.589     0.507 PREVAL_DIAB1                      1.08   0.322  3.64 0.619     0.128   0.898      
# 9  1754     47       0.589     0.507 batch_nameRKL0004_FinRisk49-80    1.23   0.603  2.50 0.362     0.565   0.572      
# 10  1754     47       0.589     0.507 batch_nameRKL001_FinRisk_Manary   1.81   0.891  3.68 0.362     1.64    0.101 

# ==============================================================================
# Boxplot: summed abundance of the 22 VFs (log10 RPKM) by liver-disease, within-carrier dose gradient
# ==============================================================================

# ---- per-participant table ---------------------------------------------------
pt <- dd %>%
  transmute(sumab22, case = .data[["INCIDENT_LIVERDIS" ]]) %>%
  mutate(grp = factor(ifelse(case == 1, "Liver\ndisease", "No liver\ndisease"), levels = c("No liver\ndisease", "Liver\ndisease")))

pt <- pt %>% filter(sumab22 > 0) # carriers only -> within-carrier dose gradient 
dim(pt)

# adjusted Cox HR (base model)
hr_label <- "HR 1.47 (1.13\u20131.92), P < 0.004"

# n per group beneath each box
nlab <- pt %>% count(grp) %>% mutate(lab = paste0("n=", n))

case_col <- "#D1495B"; non_col <- "grey55"

p_box <- ggplot(pt, aes(grp, sumab22)) +
 # geom_violin(fill = "grey92", colour = NA, scale = "width", width = 0.9) +
  geom_jitter(data = ~filter(.x, case == 0), aes(colour = "No liver disease"), width = 0.16, size = 0.4, alpha = 0.10) +
  geom_jitter(data = ~filter(.x, case == 1), aes(colour = "Liver disease"), width = 0.08, size = 1.5, alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", colour = "grey25", linewidth = 0.5, alpha=0.5) +
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

ggsave(file.path(OUTDIR, "Box_VF22_sumabund_by_case_FINRISK.svg"), p_box, width = 2.5, height = 4.2)
ggsave(file.path(OUTDIR, "Box_VF22_sumabund_by_case_FINRISK.pdf"), p_box, width = 2.5, height = 4.2, device = cairo_pdf)

```