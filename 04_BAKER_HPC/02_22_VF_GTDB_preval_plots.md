# Prepare prevalence of the 22 signal VFs in reference genomes (GTDB)

Input (Baker HPC): `gtdb_r124_vfdb_setA_diamond_70ident_90scov_filtered_hits.Rdata`: one row per genome x VF DIAMOND hit (≥70% identity, ≥90% subject coverage)

We want to output a plot showing:
- VF prevalence across all reference genomes belonging to the species detected in the FINRISK 2002 metagenomes (list is `finrisk_species_prev.csv`)
- as above, but restricted to the family Enterobacteriaceae

```R
# ==============================================================================
# 0. Config & libraries
# ==============================================================================
library(tidyverse)
options(width = 200, stringsAsFactors = FALSE)

BASE        <- "/labs/workspace/users/camilagv/active_2023/09_september_4_GTDB_VF_analysis"
PATH_HITS   <- file.path(BASE, "master_table_2025/gtdb_r124_vfdb_setA_diamond_70ident_90scov_filtered_hits.Rdata")
PATH_GENOME <- file.path(BASE, "master_table_2025/all_files.txt")
PATH_FINPREV<- file.path("finrisk_species_prev.csv")

# detection threshold = 1% of the FINRISK
MIN_DETECTION_SAMPLES <- 1

# ==============================================================================
# 1. Prepare data
# ==============================================================================

RANKS <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# the 22 FINRISK liver-disease-associated VFs (base model, FDR < 0.05)
TARGET_VFS <- c(
  "VF0213","VF0394","VF1110","VF0333","VF1138","VF0565","VF0404","VF0228",
  "VF0568","VF0572","VF0227","VF0256","VF0112","VF0571","VF0111","VF0560",
  "VF0221","VF0238","VF0113","VF0239","VF0236","VF0237")

# ---- 1. Load reference-genome VF hits ----------------------------------------
load(PATH_HITS) # loads `all_filtered_hits`
ref_hits <- all_filtered_hits; rm(all_filtered_hits)
dim(ref_hits) # 24580200 x 15

# ---- 2. Split GTDB taxonomy (parse each lineage ONCE, then join) -------------
lineages <- ref_hits %>%
  distinct(gtdb_taxonomy) %>%
  separate(gtdb_taxonomy, into = RANKS, sep = ";", remove = FALSE) %>%
  mutate(across(all_of(RANKS), ~ sub("^[a-z]__", "", .)))   # strip d__/p__/.../s__ prefixes

ref_hits <- ref_hits %>% left_join(lineages, by = "gtdb_taxonomy")

# ---- 3. Restrict to species detected in FINRISK (>=x% prevalence) ------------
finrisk_species <- read.csv(PATH_FINPREV, header = TRUE) %>%
  filter(prevalence_pct >= MIN_DETECTION_SAMPLES) %>%
  pull(Species)
length(finrisk_species) # 2259

ref_hits <- ref_hits %>% filter(Species %in% finrisk_species)
n_distinct(ref_hits$Species) # 1809

# ---- 4. Total reference genomes per species (prevalence denominator) ---------
analysed_accessions <- read_lines(PATH_GENOME) %>%
  sub(".*(GCA_[0-9]+\\.[0-9]+).*", "\\1", .) %>%
  unique()

length(analysed_accessions)# 394913

genome_counts <- tibble(ncbi_assembly_accession = analysed_accessions) %>%
  right_join(ref_hits %>% distinct(ncbi_assembly_accession, gtdb_taxonomy),
             by = "ncbi_assembly_accession") %>%
  count(gtdb_taxonomy, name = "genome_number") %>%
  arrange(desc(genome_number))

head(genome_counts, 2)
#   gtdb_taxonomy                                   genome_number
#   .......;s__Escherichia coli                              33849
#   ......;s__Klebsiella pneumoniae                          14975

ref_hits <- ref_hits %>% left_join(genome_counts, by = "gtdb_taxonomy")

# ---- 5. Per-species VF prevalence (distinct genomes carrying each VF) --------
vf_prevalence <- ref_hits %>%
  distinct(vfid, vf_name, vf_category, gtdb_taxonomy, genome_number,
           ncbi_assembly_accession) %>%           # one row per genome x VF
  count(vfid, vf_name, vf_category, gtdb_taxonomy, genome_number,
        name = "n_genomes_with_vf") %>%
  mutate(prevalence_pct = 100 * n_genomes_with_vf / genome_number) %>%
  left_join(lineages, by = "gtdb_taxonomy")        # add Domain..Species back

# ---- 6. Annotations: genome-count group + target-VF flag ---------------------
vf_prevalence <- vf_prevalence %>%
  mutate(
    genome_count_group = cut(genome_number,
                             breaks = c(0, 100, 500, 5000, Inf), right = FALSE,
                             labels = c("<100 genomes", "100-500 genomes","500-5K genomes", ">5K genomes")),
    is_target = vfid %in% TARGET_VFS
  )

# ---- 7. the 22 target VFs only ---------------------------
vf_prevalence_target <- vf_prevalence %>% filter(is_target)

# Correct "Allantion"
vf_prevalence_target$vf_name <- gsub("Allantion","Allantoin", vf_prevalence_target$vf_name)

# for the manuscript: species (>=x% in FINRISK) carrying >=1 of the 22 VFs
n_distinct(vf_prevalence_target$Species) # 104

# ==============================================================================
# 3. Plot VF prevalence across reference genomes
# ==============================================================================
library(ComplexHeatmap)
library(circlize)

# --- combined VF label (name + id) used for the heatmap columns ---------------
vf_prevalence_target <- vf_prevalence_target %>%
  mutate(vf_label = paste0(vf_name, " (", vfid, ")"))

cap_lab <- vf_prevalence_target %>% filter(vfid == "VF0560") %>% distinct(vf_label) %>% pull()

# --- ROW order: families by signal (# distinct target VFs carried), then ------
#     species within family by their own signal (so E. coli leads Enterobacteriaceae)
fam_signal <- vf_prevalence_target %>%
  group_by(Family) %>%
  summarise(fam_nvf = n_distinct(vfid), fam_prev = sum(prevalence_pct), .groups = "drop") %>%
  arrange(desc(fam_nvf), desc(fam_prev))

row_order <- vf_prevalence_target %>%
  group_by(Family, Species) %>%
  summarise(sp_nvf = n_distinct(vfid), sp_prev = sum(prevalence_pct), .groups = "drop") %>%
  mutate(Family = factor(Family, levels = fam_signal$Family)) %>%
  arrange(Family, desc(sp_nvf), desc(sp_prev))

# --- COLUMN order: by carriage within Enterobacteriaceae, Capsule first --------
ent_order <- vf_prevalence_target %>%
  filter(Family == "Enterobacteriaceae") %>%
  group_by(vf_label) %>%
  summarise(nsp = n_distinct(Species), .groups = "drop") %>%
  arrange(desc(nsp))
all_labels <- sort(unique(vf_prevalence_target$vf_label))
col_order  <- unique(c(cap_lab, ent_order$vf_label, all_labels))   # capsule, then by signal, then leftovers

# --- continuous prevalence matrix; absent pairs stay NA (= not detected) -------
prev_mat <- vf_prevalence_target %>%
  reshape2::dcast(Species ~ vf_label, value.var = "prevalence_pct") %>%
  column_to_rownames("Species") %>%
  as.matrix()
prev_mat <- prev_mat[row_order$Species, col_order, drop = FALSE]

# --- colours ------------------------------------------------------------------
col_prev <- colorRamp2(c(0, 25, 50, 75, 100), c("#DEEBF7", "#ABD9E9", "#FDAE61", "#E34A33", "#B30000"))

# genome count -> greyscale so it does not compete with the blue prevalence ramp
gc_cols <- c("<100 genomes" = "#F0F0F0", "100-500 genomes" = "#BDBDBD", "500-5K genomes" = "#737373", ">5K genomes" = "#252525")

# --- right annotation: genome count ---------
gc_grp <- vf_prevalence_target %>% distinct(Species, genome_count_group)
gc_grp <- gc_grp[match(rownames(prev_mat), gc_grp$Species), ]
row_anno <- rowAnnotation(
  "Genome count" = gc_grp$genome_count_group,
  col = list("Genome count" = gc_cols),
  annotation_name_gp = gpar(fontsize = 0)
)

# --- heatmap ------------------------------------------------------------------
ht <- Heatmap(
  prev_mat,
  name              = "Prevalence\n(% of genomes)",
  col               = col_prev,
  na_col            = "white",                          # not detected
  heatmap_legend_param = list(at = c(0, 25, 50, 75, 100)),
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_split         = factor(row_order$Family, levels = fam_signal$Family),
  row_title_rot     = 0,
  row_title_side    = "left",
  row_title_gp      = gpar(fontsize = 10),
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = 8, fontface = "italic"),
  column_names_side = "top",
  column_names_rot  = 90,
  column_names_gp   = gpar(fontsize = 10),
  border            = TRUE,
  rect_gp           = gpar(col = "grey85", lwd = 0.5),
  right_annotation  = row_anno
)

# extra legend explaining the white (NA) cells
na_lgd <- Legend(labels = "Not detected", legend_gp = gpar(fill = "white"), border = "grey60", title = "")

svg("VFDB_GTDB_22VFs.svg", width = 10, height = 14)
draw(ht, annotation_legend_list = list(na_lgd), padding = unit(c(5, 10, 5, 2), "mm"))
dev.off()

```