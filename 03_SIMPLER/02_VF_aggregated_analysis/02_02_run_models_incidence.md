# Association of microbiome-derived VF predictors with incident liver disease (SIMPLER)

Single-state Cox models testing whether microbiome-derived VF predictors are associated with incident liver disease (ICD-10 K70–K77), using the modelling dataset `SIMPLER_d0_VF_modelling.RData` (prepared in `02_02_run_models_incidence.md`).

Exposures (each continoues VF predictor is tested on both its log10 scale and per-SD of that log10 value):
- VF burden: log10(summed RPKM + 1)
- VF richness: log10(summed richness + 1)
- E. coli abundance: log10(relative abundance + 1)
- VF groups: HH vs LL (high richness/high burden vs low/low)
- Joint "profile vs load" models

Modelling strategy: batch (sequencing plate) has 150 levels, so it is handled two ways and results compared:
- Cox model (survival): plate NOT adjusted; PH check reported
- Mixed-effects Cox model (coxme): plate as a random intercept; sensitivity comparison confirming plate adjustment does not change estimates, uses COV_BASE + sequencing_plate

Two models for Cox model (survival), differing in covariates:
- COV_BASE: "BL_AGE + MEN + BMI + alcohol_gweek + antibiotic_J01_6mo + ReadDepthMillions + PREVAL_E10E14"
- COV_ADJ: COV_BASE "+ statin_6mo + ppi_6mo + factor(education) + fibre + current_smoker"

```R
# ==============================================================================
# 0. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse)
library(survival)
library(coxme)
options(width = 220, stringsAsFactors = FALSE)

OUTDIR <- "/home/cgazolla/Liver_Coli_2026/02_VF_aggregated_analysis"
setwd(OUTDIR)

# ==============================================================================
# 1. Load data
# ==============================================================================
load("SIMPLER_d0_VF_modelling.RData")
dim(d0)# 6147 x 30

# ==============================================================================
# 2. Configuration: covariate sets, exposures, outcome                          # <<< CHANGED (whole section)
# ==============================================================================
PLATE <- "aliquoting_plate"  

# Scenarios differ by COVARIATE SET
COV_BASE <- "BL_AGE + MEN + BMI + alcohol_gweek + antibiotic_J01_6mo + ReadDepthMillions + PREVAL_E10E14"
COV_ADJ  <- paste(COV_BASE, "+ statin_6mo + ppi_6mo + factor(education) + fibre + current_smoker")

EXP <- list(
  list(lab="E. coli abundance (log10)",       var="E_coli_log10",    coef="E_coli_log10",    add=NA),
  list(lab="E. coli abundance (per SD)",       var="E_coli_z",        coef="E_coli_z",        add=NA),
  list(lab="VF richness (log10)",              var="richness_log10",  coef="richness_log10",  add=NA),
  list(lab="VF richness (per SD)",             var="richness_log10_z",coef="richness_log10_z",add=NA),
  list(lab="VF burden (log10)",                var="burden_log10",    coef="burden_log10",    add=NA),
  list(lab="VF burden (per SD)",               var="burden_log10_z",  coef="burden_log10_z",  add=NA),
  list(lab="VF group (HH vs LL)",              var="group2",          coef="group2HH",        add=NA),
  list(lab="VF richness (per SD) + E. coli",   var="richness_log10_z",coef="richness_log10_z",add="E_coli_z"),
  list(lab="VF burden (per SD) + E. coli",     var="burden_log10_z",  coef="burden_log10_z",  add="E_coli_z"),
  list(lab="VF group + E. coli",               var="group2",          coef="group2HH",        add="E_coli_z")
)

OUT <- list(time="LIVERDIS_AGEDIFF", event="INCIDENT_LIVERDIS", restrict="PREVAL_LIVERDIS==0")

# ==============================================================================
# 3. Fitting
#   coxph: base (full sample) | base_cc (base covars on adjusted CC rows) | adjusted (adj covars)
#   coxme: plate random intercept on base covars (plate sensitivity only)
# ==============================================================================
build_form <- function(var, add, cov_str, rand = FALSE) {
  rhs <- if (is.na(add)) var else paste(var, "+", add)
  rhs <- paste(rhs, "+", cov_str)
  if (rand) rhs <- paste(rhs, sprintf("+ (1 | %s)", PLATE))
  as.formula(sprintf("Surv(%s, %s) ~ %s", OUT$time, OUT$event, rhs))
}
cc_mask <- function(dd, form) complete.cases(dd[, all.vars(form)])

stats_coxph <- function(fit, pick, phvar) {
  s <- summary(fit)$coefficients
  if (!pick %in% rownames(s)) return(c(est=NA,lo=NA,hi=NA,p=NA,ph_term=NA,ph_global=NA))
  b <- s[pick, ]
  est <- exp(b["coef"]); lo <- exp(b["coef"]-1.96*b["se(coef)"]); hi <- exp(b["coef"]+1.96*b["se(coef)"])
  ph_term <- NA; ph_global <- NA
  z <- tryCatch(cox.zph(fit), error = function(e) NULL)
  if (!is.null(z)) {
    if (phvar %in% rownames(z$table)) ph_term <- z$table[phvar, "p"]
    ph_global <- z$table["GLOBAL", "p"]
  }
  c(est=unname(est), lo=unname(lo), hi=unname(hi), p=unname(b["Pr(>|z|)"]),
    ph_term=unname(ph_term), ph_global=unname(ph_global))
}
stats_coxme <- function(fit, pick) {
  b_all <- fixef(fit); se_all <- sqrt(diag(as.matrix(vcov(fit)))); names(se_all) <- names(b_all)
  if (!pick %in% names(b_all)) return(c(est=NA,lo=NA,hi=NA,p=NA,ph_term=NA,ph_global=NA))
  b <- b_all[pick]; se <- se_all[pick]
  c(est=unname(exp(b)), lo=unname(exp(b-1.96*se)), hi=unname(exp(b+1.96*se)),
    p=unname(2*pnorm(-abs(b/se))), ph_term=NA, ph_global=NA)
}
mk <- function(lab, term, scenario, engine, cases, n, s) data.frame(
  exposure=lab, term=term, scenario=scenario, engine=engine,
  est=s["est"], lo=s["lo"], hi=s["hi"], p=s["p"], cases=cases, n=n,
  ph_term=s["ph_term"], ph_global=s["ph_global"], row.names=NULL, stringsAsFactors=FALSE)

# terms reported per exposure: the VF/E.coli term, plus the E.coli term for joint models
terms_to_report <- function(e) {
  out <- list(list(term="VF/E.coli", pick=e$coef, phvar=e$var, lab=e$lab))
  if (!is.na(e$add)) out[[2]] <- list(term="E.coli", pick=e$add, phvar=e$add,
                                      lab=paste0(e$lab, " [E.coli term]"))
  out
}

dd0  <- d0[eval(parse(text = OUT$restrict), d0), ]   # incident risk set (PREVAL_LIVERDIS == 0)
dd0$PREVAL_LIVERDIS %>% table()
dim(dd0) # 6147 to 6097
rows <- list()

for (e in EXP) {
  base_form <- build_form(e$var, e$add, COV_BASE)
  adj_form  <- build_form(e$var, e$add, COV_ADJ)
  cme_form  <- build_form(e$var, e$add, COV_BASE, rand = TRUE)

  # analytic samples: base CC (full) and adjusted CC (the shared base_cc / adjusted sample)
  cc_base <- cc_mask(dd0, base_form); dd_base <- dd0[cc_base, ]
  cc_adj  <- cc_mask(dd0, adj_form);  dd_adj  <- dd0[cc_adj, ]
  n_base  <- nrow(dd_base); ev_base <- sum(dd_base[[OUT$event]] == 1)
  n_adj   <- nrow(dd_adj);  ev_adj  <- sum(dd_adj[[OUT$event]]  == 1)

  f_base <- tryCatch(coxph(base_form, data = dd_base), error = function(x) NULL)  # Base: base covars, full sample
  f_bcc  <- tryCatch(coxph(base_form, data = dd_adj),  error = function(x) NULL)  # Base_CC: base covars, adjusted sample
  f_adj  <- tryCatch(coxph(adj_form,  data = dd_adj),  error = function(x) NULL)  # adjusted: adj covars, adjusted sample
  f_cme  <- tryCatch(coxme(cme_form,  data = dd_base), error = function(x) NULL)  # coxme: plate sensitivity (base covars)

  for (tt in terms_to_report(e)) {
    if (!is.null(f_base)) rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "base",    "coxph", ev_base, n_base, stats_coxph(f_base, tt$pick, tt$phvar))
    if (!is.null(f_bcc))  rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "base_cc", "coxph", ev_adj,  n_adj,  stats_coxph(f_bcc,  tt$pick, tt$phvar))
    if (!is.null(f_adj))  rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "adjusted","coxph", ev_adj,  n_adj,  stats_coxph(f_adj,  tt$pick, tt$phvar))
    if (!is.null(f_cme))  rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "coxme",   "coxme", ev_base, n_base, stats_coxme(f_cme, tt$pick))
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUTDIR, "simpler_liver_incident_full.csv"), row.names = FALSE)
write.csv(res, file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "simpler_liver_incident_full.csv"), row.names = FALSE)

# ==============================================================================
# 4. Summary table: base -> base_cc -> adjusted (+ coxme plate sensitivity)
# ==============================================================================
fmt_est <- function(e,l,h) ifelse(is.na(e), "—", sprintf("%.2f (%.2f-%.2f)", e, l, h))
fmt_p   <- function(p) ifelse(is.na(p), "—", formatC(p, format = "g", digits = 3))
fmt_cn  <- function(c,n) ifelse(is.na(n), "—", sprintf("%d / %d", c, n))
fmt_ph  <- function(p) ifelse(is.na(p), "n/a", sprintf("%.3f", p))

main <- res[res$term == "VF/E.coli", ]; expo_order <- sapply(EXP, function(x) x$lab)
getsc <- function(sc) { m <- main[main$scenario == sc, ]; m[match(expo_order, m$exposure), ] }
b <- getsc("base"); bc <- getsc("base_cc"); aj <- getsc("adjusted"); cm <- getsc("coxme")

tab <- data.frame(
  Exposure   = expo_order,
  N_full     = fmt_cn(b$cases, b$n),                       # Base: full analytic sample
  HR_base    = fmt_est(b$est, b$lo, b$hi),  P_base    = fmt_p(b$p),
  N_cc       = fmt_cn(aj$cases, aj$n),                     # complete-case sample (base_cc == adjusted)
  HR_baseCC  = fmt_est(bc$est, bc$lo, bc$hi), P_baseCC = fmt_p(bc$p),
  HR_adj     = fmt_est(aj$est, aj$lo, aj$hi), P_adj    = fmt_p(aj$p),
 # HR_coxme   = fmt_est(cm$est, cm$lo, cm$hi), P_coxme  = fmt_p(cm$p),  # plate random-intercept sensitivity
  stringsAsFactors = FALSE
)

write.csv(tab, file.path(OUTDIR, "simpler_liver_incident_table.csv"), row.names = FALSE)
write.csv(tab, file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "simpler_liver_incident_table.csv"), row.names = FALSE)

# PH check (coxph base & adjusted): exposure-term p / global p
for (i in seq_along(expo_order)) {
  cat(sprintf("  %-32s base %s / %s   adj %s / %s\n", expo_order[i],
              fmt_ph(b$ph_term[i]),  fmt_ph(b$ph_global[i]),
              fmt_ph(aj$ph_term[i]), fmt_ph(aj$ph_global[i])))
}
  # E. coli abundance (log10)        base 0.401 / 0.423   adj 0.251 / 0.304
  # E. coli abundance (per SD)       base 0.401 / 0.423   adj 0.251 / 0.304
  # VF richness (log10)              base 0.962 / 0.485   adj 0.753 / 0.368
  # VF richness (per SD)             base 0.962 / 0.485   adj 0.753 / 0.368
  # VF burden (log10)                base 0.200 / 0.315   adj 0.140 / 0.237
  # VF burden (per SD)               base 0.200 / 0.315   adj 0.140 / 0.237
  # VF group (HH vs LL)              base 0.351 / 0.081   adj 0.290 / 0.239
  # VF richness (per SD) + E. coli   base 0.955 / 0.457   adj 0.752 / 0.287
  # VF burden (per SD) + E. coli     base 0.216 / 0.413   adj 0.157 / 0.313
  # VF group + E. coli               base 0.363 / 0.136   adj 0.302 / 0.317
```