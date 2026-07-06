# Virulence factor predictor construction and health outcome dataset assembly (SIMPLER)

Assemble a single analysis dataframe combining microbiome-derived VF predictors, liver disease outcomes, and harmonized covariates

Predictors:
- `vfdb_hits_normalized_VF_level.RData`: for VF detection (from `01_virulence_factor_detection_metagenomes.md`)
- `simpler_metagenomics_mgs_relative_abundances_v4.0.tsv`: MGS relative abundances, for E. coli abundances

Outcomes + covariates:
- `tidy_clinical_data_6147.RData`: health outcomes (from `00_data_organization_SIMPLER.md`; source data in `/proj/simp2024014/Dataleverans`)

Output dataframe:

Predictors: (log transformed and log transformed per SD)
- VF burden: log10(summed RPKM + 1)
- VF richness: log10(summed richness + 1)
- E. coli abundance: log10(relative abundance + 1)
- VF groups: HH vs LL

Outcome related variables: PREVAL_LIVERDIS, INCIDENT_LIVERDIS, LIVERDIS_AGEDIFF

Covariates:
- COV_BASE: "BL_AGE + MEN + BMI + alcohol_gweek + antibiotic_J01_6mo + ReadDepthMillions + PREVAL_E10E14"
- COV_ADJ: COV_BASE "+ statin_6mo + ppi_6mo + factor(education) + fibre + current_smoker"

```R
# ==============================================================================
# 0. GLOBAL CONFIGURATION & LIBRARIES
# ==============================================================================
library(tidyverse)
library(ggExtra)
library(patchwork)
options(width = 200, stringsAsFactors = FALSE)

OUTDIR <- "02_VF_aggregated_analysis"; dir.create(OUTDIR, showWarnings = FALSE)

# ==============================================================================
# 1. Load data
# ==============================================================================
# Health outcomes + covariates
load("tidy_clinical_data_6147.RData")
health <- df; rm(df)
dim(health) # 6147 x 40

# VF detection (gene-level VFDB hits)
load("/proj/simp2024014/VF_profiles/vfdb_hits_normalized_VF_level.RData")
vfdb_hits <- diamond_final; rm(diamond_final)
dim(vfdb_hits) # 102578 x 8

# MGS relative abundances (for E. coli)
abrel <- read.csv("/proj/simp2024014/Omicsdataleverans2/metagenomics_CM/simpler_metagenomics_mgs_relative_abundances_v4.0.tsv", sep = "\t")

# Barcode <-> SIMPKEY key
dic <- read.delim("/proj/simp2024014/Omicsdataleverans2/simpler_metagenomics_key_fastaq_file_v3.0.tsv") %>%
  mutate(Sample = gsub(".*host_removal/(.*)__[1|2].fq.gz", "\\1", Files)) %>%
  select(Sample, SIMPKEY) %>%
  distinct()

# ==============================================================================
# 2. Build VF predictors (richness, burden, HH/LL groups)
# ==============================================================================
# Aggregate gene-level hits to VF level
vfdb_hits <- vfdb_hits %>%
  group_by(SampleID, VFID, VF_Name, VFcategory) %>%
  summarise(Reads = sum(NumberOfHighHomologyReads), RPKM = sum(RPKM)) %>%
  filter(Reads !=0)

# Drop cohort-wide singletons (single detection = noise) and zero-read rows
to_remove <- vfdb_hits %>% group_by(VFID) %>% summarise(prevalence = n()) %>% filter(prevalence == 1) %>% pull(VFID)
length(to_remove) # 38 singletons
vfdb_hits <- vfdb_hits %>% filter(!VFID %in% to_remove, Reads > 0)
n_distinct(vfdb_hits$VFID) # 180 non-singleton VFs
n_distinct(vfdb_hits$SampleID) # 6141 / 6147 samples with VF detection

# ---- Richness: number of distinct VFs per sample ----
richness_vfs <- vfdb_hits %>% ungroup() %>% distinct(SampleID, VFID) %>% count(SampleID, name = "richness") %>%
  mutate(richness = as.numeric(richness), richness_log10 = log10(richness + 1))
median(richness_vfs$richness) # 13
set.seed(123)
km <- kmeans(richness_vfs$richness, centers = 2)
richness_vfs$richness_group <- ifelse(km$cluster == 1, 1, 0)
table(richness_vfs$richness_group) # 1:3459 / 0:2682

# ---- Burden: summed RPKM per sample ----
burden_vfs <- vfdb_hits %>% group_by(SampleID) %>% summarise(burden = sum(RPKM)) %>%
  mutate(burden_log10 = log10(burden + 1))
median(burden_vfs$burden) # 24.46
set.seed(123)
km <- kmeans(burden_vfs$burden_log10, centers = 2)
burden_vfs$burden_group <- ifelse(km$cluster == 1, 1, 0)
table(burden_vfs$burden_group) # 1:4857 / 0:1284

# ---- Combine + cross-classify into High/Low groups (oriented to match FR02) ----
virulence_summary <- left_join(richness_vfs, burden_vfs, by = "SampleID") %>% rename(Barcode = SampleID)
virulence_summary <- left_join(dic, virulence_summary, by = c("Sample" = "Barcode"))
dim(virulence_summary) # 6147 x 8

rich_hi <- virulence_summary$richness_group == 0   # cluster 0 == High (opposite of FR02 labelling)
burd_hi <- virulence_summary$burden_group   == 0

virulence_summary$richness_burden_group <- ifelse(is.na(rich_hi) | is.na(burd_hi), "No VF detected",
                                            ifelse( rich_hi &  burd_hi, "Richness: High, Burden: High",
                                            ifelse( rich_hi & !burd_hi, "Richness: High, Burden: Low",
                                            ifelse(!rich_hi &  burd_hi, "Richness: Low, Burden: High",
                                                  "Richness: Low, Burden: Low"))))
print(table(virulence_summary$richness_burden_group))   # HH 1268 / HL 1414 / LH 16 / LL 3443 / none 6

virulence_summary$richness_group <- factor(ifelse(rich_hi, "High", "Low"), levels = c("Low","High"))
virulence_summary$burden_group   <- factor(ifelse(burd_hi, "High", "Low"), levels = c("Low","High"))

# Primary contrast: HH vs LL
HH <- "Richness: High, Burden: High"; LL <- "Richness: Low, Burden: Low"
virulence_summary$group2 <- factor(ifelse(virulence_summary$richness_burden_group %in% c(LL, HH),
                                          virulence_summary$richness_burden_group, NA),
                                   levels = c(LL, HH), labels = c("LL","HH"))

virulence_summary <- left_join(virulence_summary, health %>% select(Sample, ReadDepthMillions), by = "Sample")

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
#Kruskal-Wallis chi-squared = 117.55, df = 4, p-value < 2.2e-16

# Pairwise Wilcoxon, BH-adjusted
pw <- pairwise.wilcox.test(dat$ReadDepthMillions, dat$grp, p.adjust.method = "BH")

fmt_p <- function(p) ifelse(is.na(p), "",
                     ifelse(p < 2.2e-16, "< 2e-16",
                     ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                            formatC(p, format = "f", digits = 5))))

as.data.frame(apply(pw$p.value, c(1, 2), fmt_p), stringsAsFactors = FALSE)
#    No VF detected       HH      HL      LH
# HH        0.00104                         
# HL        0.00052 4.85e-08                
# LH        0.04101  0.11795 0.02649        
# LL        0.00145  0.00344 < 2e-16 0.20651

# Read depth by 4-way richness-burden group (VF-detected only)
dat %>% 
  group_by(grp) %>%
  summarise(median = median(ReadDepthMillions), sd = sd(ReadDepthMillions), n = n())
#   grp            median    sd     n
# 1 No VF detected   23.8  4.16     6
# 2 HH               52.4 20.0   1268
# 3 HL               55.8 23.2   1414
# 4 LH               42.2 19.0     16
# 5 LL               50.7 20.4   3443

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

svg(file.path("/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014", "distribution_richness_burden.svg"), width = 8, height = 6)
print(combined) 
dev.off()

# ==============================================================================
# 3. Assemble modelling dataframe (VF + E. coli + health)
# ==============================================================================

virulence_summary$Barcode <- virulence_summary$SIMPKEY   # ID column = SIMPKEY
d0 <- virulence_summary %>% select(Barcode, SIMPKEY,
                                   richness, richness_log10, richness_group,
                                   burden, burden_log10, burden_group,
                                   richness_burden_group, group2, ReadDepthMillions)

# E. coli (hMGS.00032): log10(rel. abundance + 1)
ecoli_df <- abrel %>% mutate(E_coli_log10 = log10(hMGS.00032 + 1)) %>% select(SIMPKEY, E_coli_log10)
d0 <- merge(d0, ecoli_df, by = "SIMPKEY")

# Health outcomes + covariates
metadata <- health[, c("SIMPKEY",
                       "PREVAL_LIVERDIS","INCIDENT_LIVERDIS","LIVERDIS_AGEDIFF",
                       "BL_AGE", "MEN", "BMI", "alcohol_gweek", "antibiotic_J01_6mo", "PREVAL_E10E14", "aliquoting_plate",
                       "statin_6mo", "ppi_6mo", "education", "fibre", "current_smoker")]
length(intersect(d0$Barcode, metadata$SIMPKEY))  # 6147
d0 <- merge(d0, metadata, by = "SIMPKEY")
nrow(d0) # 6147 (no duplication)

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
# nonNA     6147             6141           6141

# ==============================================================================
# 5. Select, coerce types, save
# ==============================================================================

str(head(d0,1))
#  $ SIMPKEY              : chr "SIMPxxxx"
#  $ Barcode              : chr "SIMPxxxxxx"
#  $ richness             : num 20
#  $ richness_log10       : num 1.32
#  $ richness_group       : Factor w/ 2 levels "Low","High": 2
#  $ burden               : num 13.7
#  $ burden_log10         : num 1.17
#  $ burden_group         : Factor w/ 2 levels "Low","High": 1
#  $ richness_burden_group: chr "Richness: High, Burden: Low"
#  $ group2               : Factor w/ 2 levels "LL","HH": NA
#  $ ReadDepthMillions    : num 58.4
#  $ E_coli_log10         : num 0.00995
#  $ PREVAL_LIVERDIS      : int 0
#  $ INCIDENT_LIVERDIS    : int 0
#  $ LIVERDIS_AGEDIFF     : num 10.7
#  $ BL_AGE               : num 97.8
#  $ MEN                  : int 0
#  $ BMI                  : num 24.4
#  $ alcohol_gweek        : num NA
#  $ antibiotic_J01_6mo   : int 0
#  $ PREVAL_E10E14        : int 0
#  $ aliquoting_plate     : chr "FH00245445"
#  $ statin_6mo           : int 0
#  $ ppi_6mo              : int 0
#  $ education            : int 1
#  $ fibre                : num NA
#  $ current_smoker       : int NA
#  $ E_coli_z             : num -0.379
#  $ richness_log10_z     : num 0.588
#  $ burden_log10_z       : num -0.703

# ---- coerce to modelling types ----
d0$richness_burden_group <- factor(d0$richness_burden_group)
d0$MEN   <- factor(d0$MEN, levels = c(0, 1), labels = c("Women","Men"))
d0$antibiotic_J01_6mo <- factor(d0$antibiotic_J01_6mo)
d0$PREVAL_E10E14 <- factor(d0$PREVAL_E10E14)
d0$aliquoting_plate <- factor(d0$aliquoting_plate)
d0$ppi_6mo <- factor(d0$ppi_6mo)
d0$statin_6mo <- factor(d0$statin_6mo)
d0$education <- factor(d0$education)
d0$current_smoker <- factor(d0$current_smoker)

save(d0, file = file.path(OUTDIR, "SIMPLER_d0_VF_modelling.RData")) 

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
N_sp <- nrow(d0) # 6147
N_detect <- d0 %>% filter(richness !=0) %>% nrow()
det  <- d0$richness_burden_group != "No VF detected"   # VF-detected participants
g    <- d0$richness_burden_group

# E. coli raw rel. abundance + detection, back-transformed from the retained log10 column
# (d0 keeps only E_coli_log10 = log10(relab + 1); these recover the raw scale exactly)
ecoli_relab <- 10^d0$E_coli_log10 - 1
ecoli_det   <- ecoli_relab > 0

sp <- c(
  profiled    = format(N_sp, big.mark = ","),
  seqdepth    = med_iqr(d0$ReadDepthMillions, 1),
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
# profiled                      6,147
# seqdepth          52.2 (34.5, 61.6)
# vf_det                6,141 (99.9%)
# vf_nodet                   6 (0.1%)
# richness                 13 (4, 29)
# burden            24.5 (14.7, 72.4)
# HH                    1,268 (20.6%)
# LL                    3,443 (56.1%)
# HL                    1,414 (23.0%)
# LH                        16 (0.3%)
# ecoli_det             5,883 (95.7%)
# ecoli_relab    0.018 (0.003, 0.153)

abrel %>% filter(hMGS.00032 !=0) %>% pull(hMGS.00032) %>% summary()
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  0.00010  0.00294  0.01840  0.89440  0.15266 89.01175 
```