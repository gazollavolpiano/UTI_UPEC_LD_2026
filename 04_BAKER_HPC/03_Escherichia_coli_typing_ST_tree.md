# E. coli phylogeny — select GTDB genomes for the 197 FINRISK STs (<=10/ST highest CheckM quality) for SKA2 alignment + NJ tree

The list of STs is available at `list_of_252_STs_found_in_FR02_SIMPLER.csv`

```R
library(tidyverse)

# ---- 1. Lineage annotation: ST, clonal complex, Clermont phylogroup ----------
strip_version <- function(x) sub("\\.\\d.*", "", x)     # GCA_000001.1_... -> GCA_000001
st <- read.csv("/labs/workspace/users/camilagv/2024/03_FR02_liver/03_phylogeny_ecoli/01_representative_phylogeny/lineages_and_genome_selection/corrected_output_ecoli_achtman_4_mlst_Seeman.csv", header = FALSE) %>%
  transmute(ACCESSION = strip_version(sub("/labs/sysgen/raw/internal/mericlab/databases/GTDB/bacteria_v214/", "", V1, fixed = TRUE)), ST = as.character(V3))

cc <- read_tsv("/labs/workspace/users/camilagv/2024/03_FR02_liver/03_phylogeny_ecoli/01_representative_phylogeny/lineages_and_genome_selection/sep_23_2024_achtman_escherichia_cc_profiles.csv", show_col_types = FALSE) %>%
  transmute(ST = as.character(ST), clonal_complex) %>%
  distinct() # guard: one CC row per ST

ez <- read_tsv("/labs/workspace/users/camilagv/2024/03_FR02_liver/03_phylogeny_ecoli/01_representative_phylogeny/lineages_and_genome_selection/ezclermont_results.txt", col_names = c("ACCESSION", "Clermont"), show_col_types = FALSE) %>%  mutate(ACCESSION = strip_version(ACCESSION))

lineages <- st %>%
  left_join(cc, by = "ST") %>%
  left_join(ez, by = "ACCESSION")

setdiff(st$ACCESSION, ez$ACCESSION) # empty (Clermont typed for all)
dim(lineages) # 33849 x 4

# ---- 2. Restrict to the 252 STs --------------------------------------
finrisk_STs   <- read_lines("list_of_252_STs_found_in_FR02_SIMPLER.csv")
lineages_fr02 <- lineages %>% filter(ST %in% finrisk_STs)
dim(lineages_fr02) # 17386 x 4

# ---- 3. Add GTDB CheckM quality metadata -------------------------------------
ecoli_meta <- read_tsv(paste0("/labs/mericlab/databases/GTDB/bac120_metadata_r214.tsv"), show_col_types = FALSE) %>%
  filter(grepl("s__Escherichia coli$", gtdb_taxonomy)) %>%
  transmute(ACCESSION = strip_version(ncbi_genbank_assembly_accession),
            ncbi_genbank_assembly_accession, ncbi_isolation_source, checkm_completeness, checkm_contamination, contig_count)

lineages_fr02 <- lineages_fr02 %>% left_join(ecoli_meta, by = "ACCESSION")
nrow(lineages_fr02) # 17386 genomes

# ---- 4. Keep up to 5 highest-quality genomes per ST -------------------------
lineages_fr02$quality_score <- (lineages_fr02$checkm_completeness - 5*lineages_fr02$checkm_contamination)

panel <- lineages_fr02 %>%
  filter(ST != "undefined") %>%
  group_by(ST) %>%
  slice_max(quality_score, n = 5, with_ties = FALSE) %>%   # deterministic top-5
  ungroup()

nrow(panel) # 1044
setdiff(finrisk_STs, panel$ST) # none missing
panel %>% count(ST, name = "n_genomes") %>% count(n_genomes)   # STs by genomes kept
#   n_genomes     n
# 1         1    32
# 2         2    19
# 3         3     9
# 4         4    13
# 5         5   179

# ---- 5. Add VF data --------------------------------------------------------
# the 22 FINRISK liver-disease-associated VFs (base model, FDR < 0.05)
TARGET_VFS <- c(
  "VF0213","VF0394","VF1110","VF0333","VF1138","VF0565","VF0404","VF0228",
  "VF0568","VF0572","VF0227","VF0256","VF0112","VF0571","VF0111","VF0560",
  "VF0221","VF0238","VF0113","VF0239","VF0236","VF0237")
  
load("/labs/workspace/users/camilagv/active_2023/09_september_4_GTDB_VF_analysis/master_table_2025/gtdb_r124_vfdb_setA_diamond_70ident_90scov_filtered_hits.Rdata")
hits <- all_filtered_hits; rm(all_filtered_hits)
dim(hits)  # 24580200 x 15

# Tally the 22 target VFs per genome 
vf_presence <- hits %>%
  filter(vfid %in% TARGET_VFS) %>%
  mutate(acc = sub("\\.[0-9]+$", "", ncbi_assembly_accession)) %>%
  distinct(acc, vfid) %>% # one row per (genome, VF)
  count(acc, name = "n_target_vfs") # how many of the 22 present

# Join onto panel, flag complete carriers
panel <- panel %>%
  left_join(vf_presence, by = c("ACCESSION" = "acc")) %>%
  mutate(
    n_target_vfs = coalesce(n_target_vfs, 0L),  # absent / none of the 22 -> 0
    has_all_22   = n_target_vfs == length(TARGET_VFS)
  )

table(panel$has_all_22)
# FALSE  TRUE 
#  1032    12 

table(panel$n_target_vfs)

panel %>% count(Clermont, has_all_22)
#    Clermont  has_all_22     n
#  1 A         FALSE        289
#  2 B1        FALSE        371
#  3 B2        FALSE        153
#  4 C         FALSE         22
#  5 D         FALSE         86
#  6 E         FALSE         51
#  7 F         FALSE         19
#  8 F         TRUE          12
#  9 G         FALSE         27
# 10 U         FALSE          5
# 11 U/cryptic FALSE          1
# 12 cryptic   FALSE          8

summary(panel$n_target_vfs)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 13.00   17.00   17.00   17.66   19.00   22.00 

# ---- 6. Write outputs --------------------------------------------------------
panel %>% pull(ncbi_genbank_assembly_accession) %>% strip_version() %>% write_lines("Ecoli_accessions_1044_genomes_selected.txt")
write.csv(panel, "metadata_lineages_1044_genomes_selected_252_STs_E_coli.csv", row.names = FALSE)
```

# Whole-genome alignment (SKA2 0.3.7, k=31) + NJ tree (RapidNJ)

Use 1044 E. coli genomes + 8 E. fergusonii outgroup 

```bash
# I added the E. fergusonii genomes to a file Fergusonii_accessions.txt
# gtdb_genomes_paths.txt has the paths to the genomes stored in the HPC

# Resolve genome paths
cat Ecoli_accessions_1044_genomes_selected.txt Fergusonii_accessions.txt > phylogeny_accessions.txt

# -F: match accessions as fixed strings (avoids one accession being a substring
#     of another and pulling extra paths)
grep -F -f phylogeny_accessions.txt gtdb_genomes_paths.txt > phylogeny_paths.txt

# Build the SKA input list (name <tab> path)
while read -r f; do
  base="${f##*/}"
  printf '%s\t%s\n' "${base%.fna}" "$f" >> phylogeny.list
done < phylogeny_paths.txt

# SKA build + align
micromamba activate gubbins

sbatch -p standard --mem=500GB -c 8 --time=04:00:00 --job-name=ska_1 \
  --wrap='ska build -o al_seqs -k 31 -f phylogeny.list --threads 8' # Submitted batch job 2334066

sbatch -p standard --mem=500GB -c 8 --time=10:00:00 --job-name=ska_2 \
  --wrap='ska align -o al_seqs.aln al_seqs.skf --threads 8 --verbose' # Submitted batch job 2334070

# Neighbour-joining tree
rapidnj al_seqs.aln -i fa > NJ_coliseqs.nwk
```

We need to tidy the tree tip labels to match the metadata

```R
library(ape)

# --- 1. Read tree -------------------------------------------------------------
tr <- read.tree("NJ_coliseqs.nwk")
length(tr$tip.label) # 1052 tips
head(tr$tip.label)

# --- 2. Strip tip labels to bare accession ------------------------------------
tr$tip.label <- sub("'", "", tr$tip.label)
tr$tip.label <- sub("^(GC[AF]_[0-9]+).*", "\\1", tr$tip.label)
head(tr$tip.label) 

# --- 3. Write simplified tree -------------------------------------------------
write.tree(tr, "NJ_coliseqs_simplified.nwk")
```

I will use Microreact to visualize the tree and metadata together