# List taxa present in FINRISK metagenomes

```R
library(tidyverse)
library(phyloseq)

OUTDIR <- "/home/camigazo/05_VF_in_GTDB"
dir.create(OUTDIR, showWarnings = FALSE)

# FINRISK species-level metagenomes (FULL set — every profiled sample)
ps <- readRDS("/home/camigazo/data/phyloseq.rds")
ps <- subset_taxa(ps, Domain == "Bacteria") 
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# count matrix, taxa as rows
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)

# collapse to SPECIES (sum counts per GTDB species), drop unclassified
sp <- as.data.frame(tax_table(ps))$Species
ok <- !is.na(sp) & sp != ""
otu_sp <- rowsum(otu[ok, , drop = FALSE], group = sp[ok])   # species x samples

# prevalence = fraction of samples with reads > 0, across the FULL set
n_samp <- ncol(otu_sp)
prev   <- rowMeans(otu_sp > 0)

species_prev <- tibble(
  Species         = rownames(otu_sp),
  n_detected      = rowSums(otu_sp > 0),
  prevalence      = prev,
  prevalence_pct  = 100 * prev
) %>%
  arrange(desc(prevalence))

write_csv(species_prev, file.path(OUTDIR, "finrisk_species_prev.csv"))
```