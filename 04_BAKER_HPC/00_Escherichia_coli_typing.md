# Escherichia coli typing 

We want lineages defined with for a tree recosntruction with E. fergusonii as an outgroup.: 
- MLST (Achtman), for STs and CCs
- Clermont 2013 PCR phylotypes

First, we need a list of genomes from GTDB metadata.

```R
# library
library(tidyverse)

# read GTDB metadata
metadata <- read.csv("bac120_metadata_r214.tsv", sep="\t", header = TRUE)

# filter for s__Escherichia coli
ecoli <- metadata[grepl("s__Escherichia coli$", metadata$gtdb_taxonomy),] 
dim(ecoli) # 33849   110

# get the outgroup
fergusonii <- metadata[grepl("s__Escherichia fergusonii$", metadata$gtdb_taxonomy),]
dim(fergusonii) # 3 110
sort(fergusonii$contig_count) # 8 genomes with 1 contig, I am gonna add just these ones
fergusonii <- fergusonii[fergusonii$contig_count == 1,]

# write file with accessions for E. coli and E. fergusonii
ecoli$ncbi_genbank_assembly_accession %>% write.table("Ecoli_accessions.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
fergusonii$ncbi_genbank_assembly_accession %>% write.table("Fergusonii_accessions.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

```bash
# grep gtdb_genomes_paths.txt (paths for each genome) with Ecoli_accessions.txt to filter it
grep -f Ecoli_accessions.txt gtdb_genomes_paths.txt > Ecoli_accessions_for_mlst.txt
wc Ecoli_accessions_for_mlst.txt -l # 33849 ok
```

Add the following script to the folder:

```bash
# !bin/bash
# mlst_Seeman.sh
cat Ecoli_accessions_for_mlst.txt | while read line
do
  echo $line
  /micromamba/envs/MLST_Seeman/bin/mlst $line --csv --scheme ecoli_achtman_4 --threads 1 >> ecoli_achtman_4_mlst_Seeman.csv
done
```

Now run with:

```bash
# run with the selected scheme (ecoli_achtman_4):
bash mlst_Seeman.sh

# need to modify files to open them with R
cut -d ',' -f 1-3 ecoli_achtman_4_mlst_Seeman.csv > corrected_output_ecoli_achtman_4_mlst_Seeman.csv

# download data on CC (clonal complexes) for Escherichia coli
wget https://rest.pubmlst.org/db/pubmlst_ecoli_achtman_seqdef/schemes/4/profiles_csv
mv profiles_csv sep_23_2024_achtman_escherichia_cc_profiles.csv
```

Also run ClermonTyping for these genomes.

```bash
# Use ezclermont (0.7.0)
cat Ecoli_accessions_for_mlst.txt | parallel 'ezclermont {} 1>> ezclermont_results.txt  2>> ezclermont_results.log'
```

Now we have the files:

- `corrected_output_ecoli_achtman_4_mlst_Seeman.csv`: with the STs for each genome
- `sep_23_2024_achtman_escherichia_cc_profiles.csv`: with the CCs for each ST
- `ezclermont_results.txt`: with the Clermont phylotypes

# StrainScan: distribution of STs and CCs in the pre-built database of E. coli

Steps:
- list the genomes in the database
- see if we have the genomes in /labs/sysgen/raw/internal/mericlab/databases/GTDB/bac120_metadata_r214.tsv (make sure they are not in list_of_15_bacteria_genomes_missing.txt)
- see if we have ST for the genomes in /labs/sysgen/workspace/users/camilagv/2024/08_LIVERDIS/02_Ecoli_quick_tree/mlst_Seeman/corrected_output_ecoli_achtman_4_mlst_Seeman.csv
- see if we have VFs for the genomes in /labs/sysgen/workspace/users/camilagv/active_2023/09_september_4_GTDB_VF_analysis/vfdb_hits_gtdb.RData
- download missing genomes? --> no, I will ignore them for now (only 2 genomes are missing)
- profile missing STs 
- see the distribution of STs and CCs in the clusters

```R
# library
library(tidyverse)

# read the cluster file (cluster ID | cluster size | reference genomes in the cluster)
hclsMap_95_recls <- read.delim("DB_Ecoli/Tree_database/hclsMap_95_recls.txt", header = FALSE, sep = "\t") %>% rename(cluster_ID = V1, cluster_size = V2, reference_genomes = V3)
head(hclsMap_95_recls)
## how many genomes
hclsMap_95_recls %>% pull(reference_genomes) %>% strsplit(",") %>% unlist() %>% unique() %>% length() # 1433
## how many clusters
hclsMap_95_recls %>% nrow() # 823

# read the file with genomes we have on the HPC cluster
bac120_metadata_r214 <- read.delim("/labs/sysgen/raw/internal/mericlab/databases/GTDB/bac120_metadata_r214.tsv", header = TRUE, sep = "\t")
## how many genomes from the database are in the HPC cluster?
### we need to extract important information from the accessions
accessions_bac120_r214 <- bac120_metadata_r214 %>% pull(ncbi_genbank_assembly_accession)
accessions_bac120_r214 <- gsub("\\..*", "", accessions_bac120_r214) # remove version number
accessions_bac120_r214 <- gsub("GCA_", "", accessions_bac120_r214) # remove GCA_
accessions_bac120_r214[1:10] # looks right
accessions_hclsMap <- hclsMap_95_recls %>% pull(reference_genomes) %>% strsplit(",") %>% unlist() %>% unique() 
accessions_hclsMap <- gsub("GCF_", "", accessions_hclsMap) # remove GCF_
setdiff(accessions_hclsMap, accessions_bac120_r214) %>% length() # 2 genomes are missing from the HPC cluster... ok!

# read the file with STs
mlst_Seeman <- read.delim("/labs/sysgen/workspace/users/camilagv/2024/08_LIVERDIS/02_Ecoli_quick_tree/mlst_Seeman/corrected_output_ecoli_achtman_4_mlst_Seeman.csv", header = FALSE, sep = ",")
accessions_mlst <- mlst_Seeman %>% pull(V1) %>% unique()
accessions_mlst <- gsub(".*GCA_", "", accessions_mlst) # remove GCF_
accessions_mlst <- gsub("\\..*", "", accessions_mlst) # remove version number
accessions_mlst[1:10] # looks right
setdiff(accessions_hclsMap, accessions_mlst) %>% length() # 1165 are not profiled... ok, we need to profile them

# read the file with VFs
# load data with VF results
load("/labs/sysgen/workspace/users/camilagv/active_2023/09_september_4_GTDB_VF_analysis/vfdb_hits_gtdb.RData") 
vfdb_hits_gtdb <- vfdb_hits_gtdb %>% select(assembly_accession, VF_Name, VFID)
dim(vfdb_hits_gtdb) # 5763655       6
accessions_vfdb <- gsub("GC._", "", vfdb_hits_gtdb$assembly_accession) %>% unique()
setdiff(accessions_hclsMap, accessions_vfdb) %>% length() # 2 genomes are missing... ok!

# list files that will be used for ST profiling
accessions_to_profile <- setdiff(accessions_hclsMap, accessions_mlst)
writeLines(accessions_to_profile, "accessions_to_profile.txt")
```

Now do the ST annotation for 1165(-2) genomes:

```bash
# activate enviroment
micromamba activate MLST_Seeman
# use ecoli_achtman_4 and ecoli schemes

# add paths to genomes in accessions_to_profile.txt
find /labs/sysgen/raw/internal/mericlab/databases/GTDB/bacteria_v214/ -type f -name "*.fna" > genomes_for_mlst_list.txt
wc genomes_for_mlst_list.txt -l # 394917
grep -f accessions_to_profile.txt genomes_for_mlst_list.txt > genomes_for_mlst_list_ST.txt ## grep genomes to filter it
wc genomes_for_mlst_list_ST.txt -l # 1163 ... true, 2 genomes are missing!
rm genomes_for_mlst_list.txt  accessions_to_profile.txt
```

Add the following script to the folder:

```bash
# !bin/bash
# mlst_Seeman.sh
cat genomes_for_mlst_list_ST.txt | while read line
do
  echo $line
  /home/cgazollavolpiano/micromamba/envs/MLST_Seeman/bin/mlst $line --csv --scheme ecoli --threads 1 >> ecoli_mlst_Seeman.csv
  /home/cgazollavolpiano/micromamba/envs/MLST_Seeman/bin/mlst $line --csv --scheme ecoli_achtman_4 --threads 1 >> ecoli_achtman_4_mlst_Seeman.csv
done
```

Now run with:

```bash
# run with the schemes:
bash mlst_Seeman.sh

# see number of files done:
wc ecoli_mlst_Seeman.csv -l  # 1163
wc ecoli_achtman_4_mlst_Seeman.csv -l # 1163
tail slurm-*.out

#  The files have a problem to be open on R, lets modify them a bit...
cut -d ',' -f 1-3 ecoli_mlst_Seeman.csv > corrected_output_ecoli_mlst_Seeman.csv
cut -d ',' -f 1-3 ecoli_achtman_4_mlst_Seeman.csv > corrected_output_ecoli_achtman_4_mlst_Seeman.csv
```

Evaluate the distribuition of STs whithin each StrainScan cluster:

```R
# library
library(tidyverse)

# read the cluster file (cluster ID | cluster size | reference genomes in the cluster)
hclsMap_95_recls <- read.delim("/labs/sysgen/workspace/users/camilagv/2024/03_FR02_liver/02_strain_detection/DB_Ecoli/Tree_database/hclsMap_95_recls.txt", header = FALSE, sep = "\t") %>% rename(cluster_ID = V1, cluster_size = V2, reference_genomes = V3)
## format the data to long
hclsMap_95_recls <- hclsMap_95_recls %>% separate_rows(reference_genomes, sep = ",")
## add simplified accessions
hclsMap_95_recls$accessions <- gsub("GCF_", "", hclsMap_95_recls$reference_genomes) # remove GCF_
## plot number of genomes per cluster
svg("n_genomes_per_cluster.svg", width = 8, height = 3)
hclsMap_95_recls %>% group_by(cluster_ID) %>% summarise("n_genomes"=n()) %>% ggplot(aes(x = n_genomes)) + geom_histogram(binwidth = 1) + labs(title = "Number of Genomes per Cluster", x = "Number of Genomes", y = "Frequency")
dev.off()
## who are the ones with a lot of genomes?
hclsMap_95_recls %>% group_by(cluster_ID) %>% summarise("n_genomes"=n()) %>% arrange(desc(n_genomes)) %>% head(10)

# read the files with STs
mlst1 <- read.delim("/labs/sysgen/workspace/users/camilagv/2024/08_LIVERDIS/02_Ecoli_quick_tree/mlst_Seeman/corrected_output_ecoli_achtman_4_mlst_Seeman.csv", header = FALSE, sep = ",")
mlst2 <- read.delim("/labs/sysgen/workspace/users/camilagv/2024/03_FR02_liver/02_strain_detection/mlst_Seeman_1165_genomes/corrected_output_ecoli_achtman_4_mlst_Seeman.csv", header = FALSE, sep = ",")
mlst <- bind_rows(mlst1, mlst2) %>% rename(ST = V3)
## filter only accessions in the database and tidy
mlst$accessions <- gsub(".*GCA_", "", mlst$V1) # remove GCF_
mlst$accessions <- gsub("\\..*", "", mlst$accessions) # remove version number
mlst <- mlst %>% filter(mlst$accessions %in% hclsMap_95_recls$accessions) %>% select(ST, accessions) %>% distinct()
## add CC annotation to ST results and tidy
cc <- read.csv("/labs/sysgen/workspace/users/camilagv/2024/08_LIVERDIS/02_Ecoli_quick_tree/mlst_Seeman/may_7_2024_achtman_escherichia_cc_profiles.csv", header = TRUE, sep = "\t") %>% select(ST, clonal_complex)
cc$ST <- as.character(cc$ST)
mlst <- mlst %>% left_join(cc, by = "ST")

# merge cluster and mlst data
metadata <- left_join(mlst, hclsMap_95_recls) %>% arrange(cluster_ID)
## add "Unknown" to ST and clonal_complex if blank or "-"
metadata$ST[metadata$ST == ""] <- "Unknown"
metadata$ST[metadata$ST == "-"] <- "Unknown"
metadata$ST[is.na(metadata$ST)] <- "Unknown"
metadata$clonal_complex[metadata$clonal_complex == ""] <- "Unknown"
metadata$clonal_complex[metadata$clonal_complex == "-"] <- "Unknown"
metadata$clonal_complex[is.na(metadata$clonal_complex)] <- "Unknown"
dim(metadata) # 1431    6
## convert metadata columns to factors
metadata$ST <- as.factor(metadata$ST)
metadata$clonal_complex <- as.factor(metadata$clonal_complex)
metadata$cluster_ID <- as.factor(metadata$cluster_ID)

# create a table with  STs | cluster (separated by ",")
metadata$cluster_ID_tmp <- paste0("C",metadata$cluster_ID)
cluster_table <- metadata %>% select(cluster_ID_tmp, ST) %>% distinct() %>% group_by(ST) %>% summarise(cluster_ID_tmp = paste(cluster_ID_tmp, collapse = ","))
write.csv(cluster_table, "cluster_table.csv", row.names = FALSE)
## can one cluster have more than one ST? NO
metadata %>% group_by(cluster_ID) %>% summarise("n_distinctST"=n_distinct(ST)) %>% filter(n_distinctST > 1) 
metadata %>% filter(cluster_ID == 791) %>% pull(ST) %>% unique() # check example

# quick check: who are the two genomes missing?
setdiff(hclsMap_95_recls$accessions, metadata$accessions) # 000147855, 012935215
hclsMap_95_recls %>% filter(reference_genomes %in% c("GCF_000147855", "GCF_012935215")) # 000147855 is in cluster 791

# save the metadata
save(metadata, file = "metadata_ecoli_db_with_STs.RData")
```