# Association of microbiome-derived VF predictors with incident liver disease (FINRISK)

Single-state Cox models testing whether microbiome-derived VF predictors are associated with incident liver disease (ICD-10 K70–K77), using the modelling dataset `FINRISK_d0_VF_modelling.RData` (prepared in `02_02_run_models_incidence.md`).

Exposures (each continoues VF predictor is tested on both its log10 scale and per-SD of that log10 value):
- VF burden: log10(summed RPKM + 1)
- VF richness: log10(summed richness + 1)
- E. coli abundance: log10(relative abundance + 1)
- VF groups: HH vs LL (high richness/high burden vs low/low)
- Joint "profile vs load" models

Modelling strategy: batch (sequencing plate) has 4 levels, so no need for Mixed-effects Cox model (coxme).
- Cox model (survival): plate adjusted; PH check reported

Two models for Cox model (survival), differing in covariates:
- COV_BASE: "BL_AGE + MEN + BMI + BL_USE_RX_J01 + ALKI2_FR02 + PREVAL_DIAB + batch_name + ReadDepthMillions
- COV_ADJ: COV_BASE "+ BL_USE_RX_A02BC + BL_USE_RX_C10AA + CURR_SMOKE + KOULGR + FIBER_TOTAL"

```R
# ==============================================================================
# 0. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse)
library(survival)
options(width = 220, stringsAsFactors = FALSE)

OUTDIR <- "/home/camigazo/02_VF_aggregated_analysis"
setwd(OUTDIR)

# ==============================================================================
# 1. Load data
# ==============================================================================
load("FINRISK_d0_VF_modelling.RData")
dim(d0)# 7226 x 29

# ==============================================================================
# 2. Configuration: covariate sets, exposures, outcome
# ==============================================================================
# Scenarios differ by COVARIATE SET
COV_BASE <- paste("BL_AGE + MEN + BMI + BL_USE_RX_J01 + ALKI2_FR02 + PREVAL_DIAB + batch_name + ReadDepthMillions")
COV_ADJ  <- paste(COV_BASE, "+ BL_USE_RX_A02BC + BL_USE_RX_C10AA + CURR_SMOKE + KOULGR + FIBER_TOTAL")

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
# Each exposure is now fit as THREE coxph models:
###   base     = base covars on the full analytic sample
###   base_cc  = base covars on the adjusted complete-case sample
###   adjusted = adjusted covars on the same complete-case sample
# ==============================================================================
build_form <- function(var, add, cov_str) {
  rhs <- if (is.na(add)) var else paste(var, "+", add)
  rhs <- paste(rhs, "+", cov_str)
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
dim(dd0) # 7182 
rows <- list()

for (e in EXP) {
  base_form <- build_form(e$var, e$add, COV_BASE)
  adj_form  <- build_form(e$var, e$add, COV_ADJ)
  
  # analytic samples: base CC (full) and adjusted CC (the shared base_cc / adjusted sample)
  cc_base <- cc_mask(dd0, base_form); dd_base <- dd0[cc_base, ]
  cc_adj  <- cc_mask(dd0, adj_form);  dd_adj  <- dd0[cc_adj, ]
  n_base  <- nrow(dd_base); ev_base <- sum(dd_base[[OUT$event]] == 1)
  n_adj   <- nrow(dd_adj);  ev_adj  <- sum(dd_adj[[OUT$event]]  == 1)
  
  f_base <- tryCatch(coxph(base_form, data = dd_base), error = function(x) NULL)  # Base: base covars, full sample
  f_bcc  <- tryCatch(coxph(base_form, data = dd_adj),  error = function(x) NULL)  # Base_CC: base covars, adjusted sample
  f_adj  <- tryCatch(coxph(adj_form,  data = dd_adj),  error = function(x) NULL)  # adjusted: adj covars, adjusted sample
  
  for (tt in terms_to_report(e)) {
    if (!is.null(f_base)) rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "base",    "coxph", ev_base, n_base, stats_coxph(f_base, tt$pick, tt$phvar))
    if (!is.null(f_bcc))  rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "base_cc", "coxph", ev_adj,  n_adj,  stats_coxph(f_bcc,  tt$pick, tt$phvar))
    if (!is.null(f_adj))  rows[[length(rows)+1]] <- mk(tt$lab, tt$term, "adjusted","coxph", ev_adj,  n_adj,  stats_coxph(f_adj,  tt$pick, tt$phvar))
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUTDIR, "finrisk_liver_incident_full.csv"), row.names = FALSE)

# ==============================================================================
# 4. Summary table: base -> base_cc -> adjusted
# ==============================================================================
fmt_est <- function(e,l,h) ifelse(is.na(e), "—", sprintf("%.2f (%.2f-%.2f)", e, l, h))
fmt_p   <- function(p) ifelse(is.na(p), "—", formatC(p, format = "g", digits = 3))
fmt_cn  <- function(c,n) ifelse(is.na(n), "—", sprintf("%d / %d", c, n))
fmt_ph  <- function(p) ifelse(is.na(p), "n/a", sprintf("%.3f", p))

main <- res[res$term == "VF/E.coli", ]; expo_order <- sapply(EXP, function(x) x$lab)
getsc <- function(sc) { m <- main[main$scenario == sc, ]; m[match(expo_order, m$exposure), ] }
b <- getsc("base"); bc <- getsc("base_cc"); aj <- getsc("adjusted")

tab <- data.frame(
  Exposure   = expo_order,
  N_full     = fmt_cn(b$cases, b$n),                       # Base: full analytic sample
  HR_base    = fmt_est(b$est, b$lo, b$hi),  P_base    = fmt_p(b$p),
  N_cc       = fmt_cn(aj$cases, aj$n),                     # complete-case sample (base_cc == adjusted)
  HR_baseCC  = fmt_est(bc$est, bc$lo, bc$hi), P_baseCC = fmt_p(bc$p),
  HR_adj     = fmt_est(aj$est, aj$lo, aj$hi), P_adj    = fmt_p(aj$p),
  stringsAsFactors = FALSE
)

write.csv(tab, file.path(OUTDIR, "finrisk_liver_incident_table.csv"), row.names = FALSE)

# PH check (coxph base & adjusted): exposure-term p / global p
for (i in seq_along(expo_order)) {
  cat(sprintf("  %-32s base %s / %s   adj %s / %s\n", expo_order[i],
              fmt_ph(b$ph_term[i]),  fmt_ph(b$ph_global[i]),
              fmt_ph(aj$ph_term[i]), fmt_ph(aj$ph_global[i])))
}
# E. coli abundance (log10)        base 0.544 / 0.215   adj 0.376 / 0.850
# E. coli abundance (per SD)       base 0.544 / 0.215   adj 0.376 / 0.850
# VF richness (log10)              base 0.402 / 0.278   adj 0.565 / 0.324
# VF richness (per SD)             base 0.402 / 0.278   adj 0.565 / 0.324
# VF burden (log10)                base 0.336 / 0.287   adj 0.514 / 0.354
# VF burden (per SD)               base 0.336 / 0.287   adj 0.514 / 0.354
# VF group (HH vs LL)              base 0.473 / 0.413   adj 0.785 / 0.586
# VF richness (per SD) + E. coli   base 0.401 / 0.238   adj 0.565 / 0.253
# VF burden (per SD) + E. coli     base 0.343 / 0.121   adj 0.526 / 0.137
# VF group + E. coli               base 0.473 / 0.384   adj 0.784 / 0.568

```