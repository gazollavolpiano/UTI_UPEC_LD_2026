# Virulence factor predictor construction and health outcome dataset assembly (FINRISK)

Assemble a single analysis dataframe combining microbiome-derived VF predictors, liver disease outcomes, and harmonized covariates

Predictors:
- `vfdb_hits_normalized.RData`: for VF detection (from `01_virulence_factor_detection_metagenomes.md`)
- `phyloseq.rds`: Phyloseq object with the microbiome data

Outcomes + covariates `health_data.RData`, build with health outcomes contained in:
- Follow-up: `Metagenomics_FR0207_endpoints_FU19_2023-02-21_v2.txt.gz`
- Baseline: `Metagenomics_FR02_phenotypes_2023-02-23.txt.gz`
- PCA from host genotype: `10PC_calculated_by_Owen.csv`
- Batch information: `microbiome_batch_information_KnightLab.tsv`
- Read depth (filtered host-depleted reads): `multiqc_fastqc.txt`
- File with NMR data: `FR02_Microbiome_NMR_2022-06-14.txt` 

Output dataframe:

Predictors: (log transformed and log transformed per SD)
- VF burden: log10(summed RPKM + 1)
- VF richness: log10(summed richness + 1)
- E. coli abundance: log10(relative abundance + 1)
- VF groups: HH vs LL

Outcome related variables: PREVAL_LIVERDIS, INCIDENT_LIVERDIS, LIVERDIS_AGEDIFF

Covariates:
- COV_BASE <- "BL_AGE + MEN + BMI + BL_USE_RX_J01 + ALKI2_FR02 + PREVAL_DIAB + batch_name + ReadDepthMillions
- COV_ADJ: COV_BASE "+ BL_USE_RX_A02BC + BL_USE_RX_C10AA + CURR_SMOKE + KOULGR + FIBER_TOTAL"

All the data is available on CSC SD Desktop at the folder `/home/camigazo/Workspace/data`

```R
# ==============================================================================
# 0. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse)
library(ggExtra)
library(patchwork)
library(phyloseq)
options(width = 200, stringsAsFactors = FALSE)

setwd("/home/camigazo")
DATA   <- "./data"
OUTDIR <- "02_VF_aggregated_analysis"; dir.create(OUTDIR, showWarnings = FALSE)

# ==============================================================================
# 1. Load data
# ==============================================================================

# Health outcomes + covariates
load(file.path(DATA, "health_data.RData"))
ff <- health$followup; bl <- health$baseline; cv <- health$covariates
health <- ff %>%
  select(Barcode, PREVAL_LIVERDIS, INCIDENT_LIVERDIS, LIVERDIS_AGEDIFF, PREVAL_DIAB, BL_USE_RX_J01, BL_USE_RX_A02BC, BL_USE_RX_C10AA) %>%
  left_join(bl %>% select(Barcode, BL_AGE, MEN, BMI, ALKI2_FR02, CURR_SMOKE, KOULGR, FIBER_TOTAL), by = "Barcode") %>%
  left_join(cv %>% select(Barcode, batch_name, ReadDepthMillions), by = "Barcode")
dim(health) # 7098 x 17

# VF detection (VFDB hits)
load(file.path(DATA, "vfdb_hits_normalized.RData"))
vfdb_hits <- diamond_final; rm(diamond_final)
dim(vfdb_hits) # 1007163 x 11

# Microbiome data from phyloseq, get rel. abundance
ps <- readRDS(file.path(DATA, "phyloseq.rds"))
ps <- subset_taxa(ps, Domain == "Bacteria")
ps <- transform_sample_counts(ps, function(x) x / sum(x))
abrel <- psmelt(ps) 
dim(abrel) # 59571144 x 10
abrel$Abundance %>% max() # needs *100, should have added before
abrel$Abundance <- abrel$Abundance * 100

# ==============================================================================
# 2. Build VF predictors (richness, burden, HH/LL groups)
# ==============================================================================
# Aggregate gene-level hits to VF level
vfdb_hits <- vfdb_hits %>%
  group_by(SampleID, VFID, VF_Name, VFcategory) %>%
  summarise(Reads = sum(NumberOfHighQualityReads), RPKM = sum(RPKM_HighQuality)) %>%
  filter(Reads !=0)

# Drop cohort-wide singletons (single detection = noise) and zero-read rows
to_remove <- vfdb_hits %>% group_by(VFID) %>% summarise(prevalence = n()) %>% filter(prevalence == 1) %>% pull(VFID)
length(to_remove) # 26 singletons
vfdb_hits <- vfdb_hits %>% filter(!VFID %in% to_remove, Reads > 0)
n_distinct(vfdb_hits$VFID) # 93 non-singleton VFs
n_distinct(vfdb_hits$SampleID) # 2156 / x samples with VF detection

# ---- Richness: number of distinct VFs per sample ----
richness_vfs <- vfdb_hits %>% ungroup() %>% distinct(SampleID, VFID) %>% count(SampleID, name = "richness") %>%
  mutate(richness = as.numeric(richness), richness_log10 = log10(richness + 1))
median(richness_vfs$richness) # 3
set.seed(123)
km <- kmeans(richness_vfs$richness, centers = 2)
richness_vfs$richness_group <- ifelse(km$cluster == 1, 1, 0)
table(richness_vfs$richness_group) # 1:555 / 0:1601

# ---- Burden: summed RPKM per sample ----
burden_vfs <- vfdb_hits %>% group_by(SampleID) %>% summarise(burden = sum(RPKM)) %>%
  mutate(burden_log10 = log10(burden + 1))
median(burden_vfs$burden) #20.08
set.seed(123)
km <- kmeans(burden_vfs$burden_log10, centers = 2)
burden_vfs$burden_group <- ifelse(km$cluster == 1, 1, 0)
table(burden_vfs$burden_group) # 1:683 / 0:1473

# ---- Combine + cross-classify into High/Low groups (oriented to match FR02) ----
virulence_summary <- left_join(richness_vfs, burden_vfs, by = "SampleID") %>% rename(Barcode = SampleID)
virulence_summary <- left_join(data.frame(Barcode=sample_names(ps)), virulence_summary)
dim(virulence_summary) # 7226 x 7

rich_hi <- virulence_summary$richness_group == 1   # cluster 1 == High 
burd_hi <- virulence_summary$burden_group   == 1

virulence_summary$richness_burden_group <- ifelse(is.na(rich_hi) | is.na(burd_hi), "No VF detected",
                                                  ifelse( rich_hi &  burd_hi, "Richness: High, Burden: High",
                                                          ifelse( rich_hi & !burd_hi, "Richness: High, Burden: Low",
                                                                  ifelse(!rich_hi &  burd_hi, "Richness: Low, Burden: High",
                                                                         "Richness: Low, Burden: Low"))))
print(table(virulence_summary$richness_burden_group))   # HH 535 / HL 20 / LH 148 / LL 1453 / none 5070

virulence_summary$richness_group <- factor(ifelse(rich_hi, "High", "Low"), levels = c("Low","High"))
virulence_summary$burden_group   <- factor(ifelse(burd_hi, "High", "Low"), levels = c("Low","High"))

# Primary contrast: HH vs LL
HH <- "Richness: High, Burden: High"; LL <- "Richness: Low, Burden: Low"
virulence_summary$group2 <- factor(ifelse(virulence_summary$richness_burden_group %in% c(LL, HH),
                                          virulence_summary$richness_burden_group, NA),
                                   levels = c(LL, HH), labels = c("LL","HH"))

virulence_summary <- left_join(cv %>% select(Barcode, ReadDepthMillions), virulence_summary)

# ==============================================================================
# 2b. QC — read depth across VF groups 
# ==============================================================================

# Short labels
lab <- c("Richness: High, Burden: High" = "HH",
         "Richness: High, Burden: Low"  = "HL",
         "Richness: Low, Burden: High"  = "LH",
         "Richness: Low, Burden: Low"   = "LL",
         "No VF detected"               = "No VF detected")

dat <- virulence_summary %>%
  mutate(grp = factor(unname(lab[as.character(richness_burden_group)]),
                      levels = c("No VF detected", "HH", "HL", "LH", "LL")))

# Omnibus Kruskal-Wallis
kruskal.test(ReadDepthMillions ~ grp, data = dat)
#Kruskal-Wallis chi-squared = 724.9, df = 4, p-value < 2.2e-16

# Pairwise Wilcoxon, BH-adjusted
pw <- pairwise.wilcox.test(dat$ReadDepthMillions, dat$grp, p.adjust.method = "BH")

fmt_p <- function(p) ifelse(is.na(p), "",
                            ifelse(p < 2.2e-16, "< 2e-16",
                                   ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                                          formatC(p, format = "f", digits = 5))))

as.data.frame(apply(pw$p.value, c(1, 2), fmt_p), stringsAsFactors = FALSE)
#No VF detected       HH       HL      LH
#HH        < 2e-16                          
#HL       9.55e-14 1.01e-10                 
#LH        0.01614  < 2e-16 1.69e-12        
#LL        < 2e-16  0.59544 1.53e-10 < 2e-16

# Read depth by 4-way richness-burden group (VF-detected only)
dat %>% 
  group_by(grp) %>%
  summarise(median = median(ReadDepthMillions), sd = sd(ReadDepthMillions), n = n())
#grp            median    sd     n
#1 No VF detected  0.686 0.576  5070
#2 HH              1.06  1.57    535
#3 HL              3.70  5.01     20
#4 LH              0.642 0.437   148
#5 LL              1.08  1.50   1453

# ==============================================================================
# 2c. Exposure distributions
# ==============================================================================
library(patchwork)

# Palette keyed to the full group names (ordered HH -> HL -> LH -> LL)
grp_levels <- c("Richness: High, Burden: High", "Richness: High, Burden: Low",
                "Richness: Low, Burden: High",  "Richness: Low, Burden: Low")
vf_cols4 <- c("Richness: High, Burden: High" = "#FA7F72",   # coral
              "Richness: High, Burden: Low"  = "#E8A33D",   # amber
              "Richness: Low, Burden: High"  = "#4FA3A5",   # teal
              "Richness: Low, Burden: Low"   = "#675BC6")   # purple

vf_plot_df <- virulence_summary %>%
  filter(!is.na(richness)) %>% # drop "No VF detected"
  mutate(grp = factor(richness_burden_group, levels = grp_levels))

# k-means cut points (richness on raw scale, burden on log10 scale)
rich_cut <- with(vf_plot_df, mean(c(max(richness[richness_group == "Low"]),
                                    min(richness[richness_group == "High"]))))
burd_cut <- with(vf_plot_df, mean(c(max(burden_log10[burden_group == "Low"]),
                                    min(burden_log10[burden_group == "High"]))))

# Base scatter: larger points, stronger alpha, bigger legend keys
p_joint <- vf_plot_df %>%
  ggplot(aes(richness, burden_log10, colour = grp)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = rich_cut, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = burd_cut, linetype = 2, colour = "grey40") +
  scale_colour_manual(name = "VF group", values = vf_cols4) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(x = "Richness (distinct VFs)",
       y = "Burden, log10(summed RPKM + 1)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")          # legend handled separately below

# Marginal densities (no legend on this object)
p_joint_m <- ggMarginal(p_joint, groupColour = TRUE, groupFill = TRUE,
                        type = "density", alpha = 0.4)

# Pull a clean legend from a twin plot, then lay it beside the marginal plot
legend <- cowplot::get_legend(
  p_joint +
    theme(legend.position = "right",
          legend.title = element_text(size = 13),
          legend.text  = element_text(size = 11),
          legend.key.size = unit(0.9, "cm"))
)

combined <- wrap_elements(p_joint_m) + wrap_elements(legend) + plot_layout(widths = c(4, 1.3))

svg(file.path(OUTDIR, "distribution_richness_burden.svg"), width = 8, height = 6)
print(combined)
dev.off()

# ==============================================================================
# 3. Assemble modelling dataframe (VF + E. coli + health)
# ==============================================================================
d0 <- virulence_summary %>% select(Barcode, 
                                   richness, richness_log10, richness_group,
                                   burden, burden_log10, burden_group,
                                   richness_burden_group, group2, ReadDepthMillions)

# E. coli log10(rel. abundance + 1)
ecoli_df <- abrel %>% filter(Species == "Escherichia coli") %>%  mutate(E_coli_log10 = log10(Abundance + 1), Barcode=Sample) %>% select(Barcode, E_coli_log10)
d0 <- merge(d0, ecoli_df)

# Health outcomes + covariates
length(intersect(d0$Barcode, health$Barcode))  # 7098
d0 <- left_join(d0, health)
nrow(d0) # 7226

# ---- Outcome counts (QC) ----
print(table(d0$PREVAL_LIVERDIS)); print(table(d0$INCIDENT_LIVERDIS))

# ==============================================================================
# 4. Standardised exposures (per-SD of the log10 values)
# ==============================================================================
d0$E_coli_z         <- as.numeric(scale(d0$E_coli_log10))
d0$richness_log10_z <- as.numeric(scale(d0$richness_log10))
d0$burden_log10_z   <- as.numeric(scale(d0$burden_log10))

print(sapply(d0[c("E_coli_z","richness_log10_z","burden_log10_z")],
             function(x) round(c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), nonNA = sum(!is.na(x))), 3)))
#       E_coli_z richness_log10_z burden_log10_z
# mean         0                0              0
# sd           1                1              1
# nonNA     7226             2156           2156

# ==============================================================================
# 5. Select, coerce types, save
# ==============================================================================

str(head(d0,1))
# $ Barcode              : chr "...."
# $ ReadDepthMillions    : num 0.463
# $ richness             : num NA
# $ richness_log10       : num NA
# $ richness_group       : Factor w/ 2 levels "Low","High": NA
# $ burden               : num NA
# $ burden_log10         : num NA
# $ burden_group         : Factor w/ 2 levels "Low","High": NA
# $ richness_burden_group: chr "No VF detected"
# $ group2               : Factor w/ 2 levels "LL","HH": NA
# $ E_coli_log10         : num 0
# $ PREVAL_LIVERDIS      : int 0
# $ INCIDENT_LIVERDIS    : int 0
# $ LIVERDIS_AGEDIFF     : num 10.7
# $ PREVAL_DIAB          : int 1
# $ BL_USE_RX_J01        : int 0
# $ BL_USE_RX_A02BC      : int 0
# $ BL_USE_RX_C10AA      : int 1
# $ BL_AGE               : num 57.9
# $ MEN                  : int 0
# $ BMI                  : num 49.6
# $ ALKI2_FR02           : num 0
# $ CURR_SMOKE           : int 0
# $ KOULGR               : int 2
# $ FIBER_TOTAL          : int NA
# $ batch_name           : Factor w/ 4 levels "Knight_FinRisk_Repool2",..: 2
# $ E_coli_z             : num -0.225
# $ richness_log10_z     : num NA
# $ burden_log10_z       : num NA

# ---- coerce to modelling types ----
d0$richness_burden_group <- factor(d0$richness_burden_group)
d0$MEN   <- factor(d0$MEN, levels = c(0, 1), labels = c("Women","Men"))
d0$PREVAL_DIAB <- factor(d0$PREVAL_DIAB)
d0$BL_USE_RX_J01 <- factor(d0$BL_USE_RX_J01)
d0$BL_USE_RX_C10AA <- factor(d0$BL_USE_RX_C10AA)
d0$BL_USE_RX_A02BC <- factor(d0$BL_USE_RX_A02BC)
d0$CURR_SMOKE <- factor(d0$CURR_SMOKE)
d0$KOULGR <- factor(d0$KOULGR)

save(d0, file = file.path(OUTDIR, "FINRISK_d0_VF_modelling.RData")) 

# ==============================================================================
# 6. Table 4 — VF carriage & E. coli abundance 
# ==============================================================================
# Helpers
fmt_n_pct <- function(n, denom) sprintf("%s (%.1f%%)", format(n, big.mark = ",", trim = TRUE), 100 * n / denom)
med_iqr   <- function(x, d = 1) {
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE)
  f <- function(v) format(round(v, d), big.mark = ",", trim = TRUE, nsmall = d)
  sprintf("%s (%s, %s)", f(q[1]), f(q[2]), f(q[3]))
}

# ---- SIMPLER values (computed from d0) ----
N_sp <- nrow(d0) # 7226
N_detect <- d0 %>% filter(richness !=0) %>% nrow()
det  <- d0$richness_burden_group != "No VF detected"   # VF-detected participants
g    <- d0$richness_burden_group

# E. coli raw rel. abundance + detection, back-transformed from the retained log10 column
# (d0 keeps only E_coli_log10 = log10(relab + 1); these recover the raw scale exactly)
ecoli_relab <- 10^d0$E_coli_log10 - 1
ecoli_det   <- ecoli_relab > 0

sp <- c(
  profiled    = format(N_sp, big.mark = ","),
  seqdepth    = med_iqr(d0$ReadDepthMillions, 2),
  vf_det      = fmt_n_pct(sum(det),  N_sp),
  vf_nodet    = fmt_n_pct(sum(!det), N_sp),
  richness    = med_iqr(d0$richness[det], 0), # distinct VFs, among detected
  burden      = med_iqr(d0$burden[det],   1), # summed RPKM, among detected
  HH          = fmt_n_pct(sum(g == "Richness: High, Burden: High"), N_detect),
  LL          = fmt_n_pct(sum(g == "Richness: Low, Burden: Low"),   N_detect),
  HL          = fmt_n_pct(sum(g == "Richness: High, Burden: Low"),  N_detect),
  LH          = fmt_n_pct(sum(g == "Richness: Low, Burden: High"),  N_detect),
  ecoli_det   = fmt_n_pct(sum(ecoli_det), N_sp),
  ecoli_relab = med_iqr(ecoli_relab[ecoli_det], 3) # among detected
)

sp %>% as.data.frame()
# profiled                   7,226
# seqdepth          0.76 (0.55, 1.17)
# vf_det             2,156 (29.8%)
# vf_nodet           5,070 (70.2%)
# richness               3 (1, 13)
# burden         20.1 (8.7, 184.6)
# HH                   535 (24.8%)
# LL                 1,453 (67.4%)
# HL                     20 (0.9%)
# LH                    148 (6.9%)
# ecoli_det          3,082 (42.7%)
# ecoli_relab 0.238 (0.095, 0.873)

```