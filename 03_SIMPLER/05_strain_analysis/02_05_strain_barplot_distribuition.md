# SIMPLER: E. coli strain (StrainScan) typing & ST prevalence

Code to finish the E. coli strain-level typing analysis the SIMPLER cohort, figuring out which E. coli sequence types (STs) are present in each metagenome and how common each ST is across the cohort. 

- `metadata_ecoli_db_with_STs.RData`: ST, CC, cluster ID, for 1431 genomes of the database, created from `00_Escherichia_coli_typing.md`
- `ezclermont_results.txt` (added to the folder): Clermont phylotypes, created from `00_Escherichia_coli_typing.md`
```R
# ---- 0. Config & libraries ---------------------------------------------------
library(tidyverse)
options(width = 200, stringsAsFactors = FALSE)

OUTDIR <- "05_strain_analysis"; dir.create(OUTDIR, showWarnings = FALSE)
WHARF  <- "/proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014"

STRAINSCAN_DIR <- "/proj/simp2024014/Ecoli_StrainScan/00_default_config/strainscan_output"
PATH_CLUSTER_ST <- "/proj/simp2024014/Ecoli_StrainScan/00_default_config/metadata_ecoli_db_with_STs.RData"
PATH_CLER <- paste0(OUTDIR,"/ezclermont_results.txt")

# ---- 1. Read StrainScan reports ----------------------------------------------
sample_dirs <- list.files(STRAINSCAN_DIR, full.names = TRUE)
N_SAMPLES <- 4548

read_report <- function(d) {
  f <- file.path(d, "final_report.txt")
  if (!file.exists(f)) return(NULL)
  x <- tryCatch(read.delim(f, header = TRUE), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  ab_col <- intersect(c("Relative_Abundance", "Relative_Abundance_Inside_Cluster"), names(x))[1]
  if (is.na(ab_col)) return(NULL)
  tibble(Sample             = basename(d),
         Species            = "Escherichia coli",
         Cluster_ID         = x$Cluster_ID,
         Relative_Abundance = x[[ab_col]])
}

load("/home/cgazolla/Liver_Coli_2026/02_VF_aggregated_analysis/SIMPLER_d0_VF_modelling.RData")
d0 %>%
filter(PREVAL_LIVERDIS==0) %>%
drop_na(group2, BL_AGE, MEN, BMI, alcohol_gweek, antibiotic_J01_6mo, ReadDepthMillions, PREVAL_E10E14) %>%
pull(SIMPKEY) -> samples_to_filt

read.delim("/proj/simp2024014/Omicsdataleverans2/simpler_metagenomics_key_fastaq_file_v3.0.tsv") %>%
  mutate(Sample = gsub(".*host_removal/(.*)__[1|2].fq.gz", "\\1", Files)) %>%
  select(Sample, SIMPKEY) %>%
  distinct() %>%
  filter(SIMPKEY %in% samples_to_filt) %>%
  pull(Sample) -> samples_to_filt

length(samples_to_filt) # 4548

strain_hits <- map(sample_dirs, read_report) %>% bind_rows() %>% filter(Sample %in% samples_to_filt)
n_distinct(strain_hits$Sample) # 1479 samples with >=1 E. coli strain detected

# ---- 2. Map clusters to STs ---------------------------------------
load(PATH_CLUSTER_ST)
cluster_st <- metadata %>% distinct(Cluster_ID = cluster_ID_tmp, ST)

# cluster detected in samples but missing from the DB annotation -> label Unknown
setdiff(unique(strain_hits$Cluster_ID), cluster_st$Cluster_ID)   # C725
cluster_st <- bind_rows(cluster_st, tibble(Cluster_ID = "C725", ST = "Unknown"))

# clusters mapping to >1 ST: collapse to a single comma-separated label
multi_st <- cluster_st %>% count(Cluster_ID, name = "n_st") %>% filter(n_st > 1) %>% pull(Cluster_ID)
length(multi_st)                                                 # 13 clusters
strain_hits %>% filter(Cluster_ID %in% multi_st) %>%
  summarise(n_samples_affected = n_distinct(Sample))             # 71 samples
cluster_st <- cluster_st %>%
  group_by(Cluster_ID) %>%
  summarise(ST = paste(unique(ST), collapse = ","), .groups = "drop")

strain_hits <- strain_hits %>%
  left_join(cluster_st, by = "Cluster_ID") %>%
  distinct(Sample, Species, Cluster_ID, ST, Relative_Abundance)

# ---- 3. ST prevalence across the cohort --------------------------------------
st_prevalence <- strain_hits %>%
  distinct(ST, Sample) %>%                     # a sample may carry several STs
  count(ST, name = "n_samples") %>%
  mutate(prevalence_pct = 100 * n_samples / N_SAMPLES) %>%
  arrange(desc(n_samples))

write.csv(st_prevalence, file.path(OUTDIR, "ST_prevalence_SIMPLER.csv"), row.names = FALSE)
write.csv(st_prevalence, file.path(WHARF, "ST_prevalence_SIMPLER.csv"), row.names = FALSE)

# ---- 4. List of distinct detected STs (for the reference-genome panel) -------
detected_STs <- strain_hits$ST %>%
  strsplit(",") %>% unlist() %>% unique()
detected_STs <- detected_STs[!grepl("Unknown", detected_STs)]
length(detected_STs)      # 187 distinct STs
write.table(detected_STs, file.path(OUTDIR, "STs_detected_SIMPLER.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(detected_STs, file.path(WHARF, "STs_detected_SIMPLER.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---- 5. Summary counts -------------------------------------------------------
n_distinct(strain_hits$Sample)        # 1479 samples with strain detection
n_distinct(strain_hits$Cluster_ID)    # 303 distinct clusters

# ----top-20 STs, phylogroup lineages highlighted ----------
# We want to highlight if the ST belong to B2, F, or D, use metadata
clermont <- read.csv(PATH_CLER, sep="\t", header=FALSE) %>% 
mutate(accessions=gsub("GC._(.........).*","\\1", V1)) %>%
select(accessions, V2) %>%
rename(Clermont=V2)

metadata <- left_join(metadata, clermont)

phylogroup_B2_F_D_STS <- metadata %>% filter(Clermont %in% c("B2", "F", "D")) %>% pull(ST) %>% unique() %>% as.character()

top20 <- st_prevalence %>%
  filter(ST != "Unknown") %>%
  slice_max(prevalence_pct, n = 20, with_ties = FALSE) %>%
  mutate(
    is_B2_F_D = map_lgl(strsplit(ST, ","), ~ any(.x %in% phylogroup_B2_F_D_STS)),
    label   = paste0("ST", gsub(",", "/", ST)),
    label   = fct_reorder(label, prevalence_pct)
  )

p_st <- ggplot(top20, aes(prevalence_pct, label, fill = is_B2_F_D)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#B30000", `FALSE` = "grey75"),
                    labels = c(`TRUE` = "B2, D or F-associated lineage", `FALSE` = "Other"),
                    name   = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Prevalence (% of SIMPLER samples)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor    = element_blank(),
        panel.grid.major.x  = element_line(colour = "grey88", linewidth = 0.4),
        axis.line.x  = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks.x = element_line(colour = "grey20", linewidth = 0.4),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(colour = "grey15"),
        axis.text.x  = element_text(colour = "grey15"),
        axis.title.x = element_text(colour = "grey15"),
        legend.position = "top")

ggsave(file.path(OUTDIR, "Barplot_top20_STs_SIMPLER.pdf"), p_st, width = 5.5, height = 5, device = cairo_pdf)
ggsave(file.path(OUTDIR, "Barplot_top20_STs_SIMPLER.svg"), p_st, width = 5.5, height = 5)

ggsave(file.path(WHARF, "Barplot_top20_STs_SIMPLER.pdf"), p_st, width = 5.5, height = 5, device = cairo_pdf)
ggsave(file.path(WHARF, "Barplot_top20_STs_SIMPLER.svg"), p_st, width = 5.5, height = 5)
```