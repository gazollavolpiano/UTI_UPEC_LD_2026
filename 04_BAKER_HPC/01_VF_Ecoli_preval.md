# Prepare VF prevalence in E. coli reference genomes (GTDB)

Builds `VF_Ecoli_preval.csv`: for each VFDB virulence factor, the percentage of GTDB E. coli genomes that carry it. We will use it for annotation for the VF heatmap (joined onto the VF panel).

Input (BAKER HPC)
- `gtdb_r124_vfdb_setA_diamond_70ident_90scov_filtered_hits.Rdata`: one row per genome × VF DIAMOND hit (≥70% identity, ≥90% subject coverage)

Output: `VF_Ecoli_preval.csv` with `vfid`, `vf_name`, `vf_category`, `n_genomes` (E. coli genomes carrying the VF), `ecoli_prevalence_pct`.

```R
# ==============================================================================
# 0. Global configuration & libraries
# ==============================================================================
library(tidyverse)
options(width = 200, stringsAsFactors = FALSE)

# ==============================================================================
# 1. Load reference-genome VF hits
#    GTDB genomes x VFDB setA, DIAMOND blastx, retained at >=70% id / >=90% scov
# ==============================================================================
load("/labs/workspace/users/camilagv/active_2023/09_september_4_GTDB_VF_analysis/master_table_2025/gtdb_r124_vfdb_setA_diamond_70ident_90scov_filtered_hits.Rdata")
hits <- all_filtered_hits; rm(all_filtered_hits)
dim(hits)  # 24580200 x 15

# ==============================================================================
# 2. Subset to E. coli genomes
# ==============================================================================
ecoli <- hits %>% filter(str_detect(gtdb_taxonomy, ";s__Escherichia coli$"))
ecoli$gtdb_taxonomy %>% unique() # [1] "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"

n_ecoli_hit_genomes <- n_distinct(ecoli$ncbi_assembly_accession)  # E. coli genomes with >=1 hit
n_ecoli_hit_genomes # 33849, correct

# ==============================================================================
# 3. Per-VF prevalence in E. coli genomes
# ==============================================================================
ecoli_preval <- ecoli %>%
  distinct(vfid, vf_name, vf_category, ncbi_assembly_accession) %>%  # one row per genome x VF
  count(vfid, vf_name, vf_category, name = "n_genomes") %>%
  mutate(ecoli_prevalence_pct = 100 * n_genomes / n_ecoli_hit_genomes)

# ==============================================================================
# 4. Comprehensive table: every detected VF, 0% if absent from E. coli
# ==============================================================================
all_vfs <- hits %>% distinct(vfid, vf_name, vf_category)

VF_Ecoli_preval <- all_vfs %>%
  left_join(ecoli_preval %>% select(vfid, n_genomes, ecoli_prevalence_pct), by = "vfid") %>%
  mutate(
    n_genomes            = replace_na(n_genomes, 0L),
    ecoli_prevalence_pct = replace_na(ecoli_prevalence_pct, 0)
  ) %>%
  arrange(desc(ecoli_prevalence_pct))

dim(VF_Ecoli_preval) # 623 x 5
head(VF_Ecoli_preval, 10)

# Write
write.csv(VF_Ecoli_preval, "VF_Ecoli_preval.csv", row.names = FALSE)
```