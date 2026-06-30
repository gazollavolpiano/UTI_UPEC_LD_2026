# Detection of virulence factors 

We are going to run the detection of virulence factors based on diamond and the VFDB database. The VFDB database has 4193 protein sequences curated from bacterial pathogens. The diamond software is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.

Add the loop scrip vfdb_diamond_array.sh to the folder. The scripts uses the following parameters:
* --max-target-seqs 1 --max-hsps 1: limits the search to the single best aligning sequence (--max-target-seqs 1) and the single best high-scoring segment pair (--max-hsps 1) per query sequence.
* --sensitive: uses the sensitive mode of DIAMOND, which is slower but more accurate.
* IMPORTANT: do not use --subject-cover, --id, and --query-cover parameters, they simply do not work as expected. Instead, we are going to filter the results in the R script.

Notes:
* the file `seqkit_stats_results_6147_participants.csv` has the read depth for each sample, calculated with seqkit stats. 
* before starting, I transferred the file `VFDB_setA_pro.dmnd` to the analysis folder, it is the diamond database created with the VFDB protein sequences.

This script accepts a single sample name (via -s) and runs the diamond blastx command on its corresponding file. It checks that the file exists.

```bash
#!/bin/bash
# vfdb_diamond_singlesample.sh
# This script processes a single sample using diamond blastx.
#
# Usage: bash vfdb_diamond_singlesample.sh -s sample_name

# Parse command-line arguments
while getopts "s:" arg; do
  case $arg in
    s) sample="$OPTARG" ;;
    *) echo "Usage: $0 -s sample_name"; exit 1 ;;
  esac
done

if [ -z "$sample" ]; then
  echo "Usage: $0 -s sample_name"
  exit 1
fi

# Construct the sample file path
sample_file=$(echo /metagnomics/*/*/host_removal/${sample}__1.fq.gz)

if [ ! -f "$sample_file" ] || [ ! -r "$sample_file" ]; then
  echo "$(date +"%T") - Error: Sample file $sample_file not found or is not readable."
  exit 1
fi

# Run diamond blastx on the sample
echo "$(date +"%T") - Processing sample $sample with diamond blastx..."
diamond blastx --db /proj/simp2024014/VF_profiles/VFDB_setA_pro.dmnd \
  --query "$sample_file" \
  --out "vfdb_diamond_output/${sample}_vfdb_diamond.txt" \
  --outfmt 6 qtitle stitle pident length sstart send slen qlen \
  --max-target-seqs 1 --max-hsps 1 --sensitive --threads 10

if [ $? -eq 0 ]; then
  echo "$(date +"%T") - Finished processing sample $sample."
else
  echo "$(date +"%T") - Error processing sample $sample."
fi
```

# Run the processing of the results: create vfdb_hits.RData

This step is going to process the results from the diamond analysis. The script is going to create two groups:
* Homologous: subject coverage ≥80%, identity ≥ 50%, query coverage ≥ 50%
* High-Quality: subject coverage ≥80%, identity ≥ 90%, query coverage ≥ 50% --> SELECTED FOR THE ANALYSIS

Add the following R script (virulence_factor_process.R) to the folder:

```R
# !/usr/bin/env Rscript
# virulence_factor_process.R

# call libs
library(tidyverse)
library(dplyr)
library("IRanges",warn.conflicts=FALSE)

# list the files in the folder
files <- list.files(pattern = "vfdb_diamond.txt", full.name=TRUE)

# read each file with result and process statistics
diamond_final <- c() 
for (file in files){
    # get sample name
    sample <- gsub("./(.*)_.*_diamond.txt","\\1",file)
    cat("processing", sample, "\n")

    # read file and add column names
    gene_hits <- tryCatch(read.delim(file, header=F), error=function(e) NULL)
    if(is.null(gene_hits)){next} # if the file is empty, go to the next file
    colnames(gene_hits) <- c("qtitle", "stitle", "pident", "length", "sstart", "send", "slen", "qlen")

    # filter for subject coverage ≥80% 
    gene_hits$qlen_aa <- round(gene_hits$qlen/3) # slen is the length of the read in nucleotides, I need to "convert it to amino acids"
    gene_hits$min_lenght <- round(gene_hits$qlen_aa*0.8) # calculate the minimum length for the alignment (80% of the read length)
    gene_hits$remove <- ifelse(gene_hits$length<gene_hits$min_lenght, 1, 0) # mark the alignments with less than 80% of the read length
    gene_hits <- gene_hits %>% filter(remove==0) # remove alignments with less than 80% of the read length
    if(dim(gene_hits)[1]==0){next} # if the file is empty now, go to the next file

    # add function to calculate the coverage of the gene using IRanges package
    calc_coverage <- function(gene_hits_dataframe){  
      coverage <- data.frame()
      for(gene in unique(gene_hits_dataframe$stitle)){
        tmp <- gene_hits_dataframe %>% filter(stitle==gene) # get the alignments for the gene
        reads <- dim(tmp)[1] # count the number of reads for the gene
        tmp$starts <- apply(tmp[,c("sstart", "send")],1,min) # get the start position for each alignment
        tmp$stops <- apply(tmp[,c("sstart", "send")],1,max) # get the stop position for each alignment
        myrange <- IRanges(start=tmp$starts,end=tmp$stops) # create a range for the gene
        cov_posit <- as.numeric(coverage(myrange)) # get a vector with the coverage for each position (1 if covered, 0 if not covered), it goes until the last position covered (not the last on the gene)
        cov_posit <- sum(cov_posit > 0) # count the number of positions with coverage
        gene_len <- unique(tmp$slen) # get the length of the gene
        gene_cov <- (cov_posit/gene_len)*100 # calc gene coverage
      coverage <- as.data.frame(rbind(coverage, cbind(gene, gene_len, gene_cov, reads)))
    }
    coverage$gene_cov <- as.numeric(coverage$gene_cov)
    coverage$reads <- as.numeric(coverage$reads)
    coverage$gene_len <- as.numeric(coverage$gene_len)
    return(coverage)
    }

    # COMMENTED BECAUSE WE ARE NOT USING THE "less_homology" GROUP ANYMORE
    # # filter data for "less_homology" group
    # ## remove hits with less than 50% of identity (pident)
    # gene_hits_less_homology <- gene_hits %>% filter(pident>=50) 
    # if(dim(gene_hits_less_homology)[1]==0){next} # if the file is empty now, go to the next file (high_homology group will be empty too)
    # ## calculate the gene coverage for each gene 
    # gene_hits_less_homology <- calc_coverage(gene_hits_less_homology)
    # ## filter hits with less than 50% of coverage
    # gene_hits_less_homology <- gene_hits_less_homology[gene_hits_less_homology$gene_cov>=50,]
    # if(dim(gene_hits_less_homology)[1]==0){next} # if the file is empty now, go to the next file (high_homology group will be empty too)
    # ## rename columns and filter columns
    # gene_hits_less_homology <- gene_hits_less_homology %>% select(gene_header=1, GeneLengthAminoAcids=2, NumberOfLessHomologyReads=4)

    # filter data for "high_homology" group
    ## remove hits with less than 90% of identity (pident)
    gene_hits_high_homology <- gene_hits %>% filter(pident>=90)
    ## calculate the gene coverage for each gene
    gene_hits_high_homology <- calc_coverage(gene_hits_high_homology)
    ## filter hits with less than 50% of coverage
    gene_hits_high_homology <- gene_hits_high_homology[gene_hits_high_homology$gene_cov>=50,]

    if(dim(gene_hits_high_homology)[1]==0){
      # COMMENTED BECAUSE WE ARE NOT USING THE "less_homology" GROUP ANYMORE
      # set gene_hits with only the homologous hits (high_homology group will be empty)
      # gene_hits <- gene_hits_less_homology
      ## add NumberOfHighHomologyReads column with 0
      # gene_hits$NumberOfHighHomologyReads <- 0

      next # if the file is empty now, go to the next file

      } else {
        ## rename columns and filter columns
        gene_hits_high_homology <- gene_hits_high_homology %>% select(gene_header=1, GeneLengthAminoAcids=2, NumberOfHighHomologyReads=4)

        # tidy final results
        ## combine data
        # gene_hits <- full_join(gene_hits_less_homology, gene_hits_high_homology, by = c("gene_header", "GeneLengthAminoAcids"))
        gene_hits <- gene_hits_high_homology

        ## add 0 instead of NA
        gene_hits$NumberOfHighHomologyReads[is.na(gene_hits$NumberOfHighHomologyReads)] <- 0
        #gene_hits$NumberOfLessHomologyReads[is.na(gene_hits$NumberOfLessHomologyReads)] <- 0
      }

    # add the VFID column
    gene_hits$VFID <- gsub(".*_\\((VF\\d+)\\).*", "\\1", gene_hits$gene_header)

    # add the sample name
    gene_hits$SampleID <- sample

    # organize columns order
    gene_hits <- gene_hits %>% select(SampleID, VFID, GeneLengthAminoAcids, #NumberOfLessHomologyReads, 
                                      NumberOfHighHomologyReads, gene_header)

    # append data
    diamond_final <- rbind(diamond_final, gene_hits)
}

# save the data
save(diamond_final, file="vfdb_hits.RData")

# print that the script finished
cat("The script finished!")
```
Now we can check the results in R. 

We need to compute the RPKM (Reads Per Kilobase of transcript, per Million mapped reads) for each VF and **genes**.

```R
# Libraries
library(tidyverse)

# Load the data
load("vfdb_diamond_output/vfdb_hits.RData")

# Check the data
dim(diamond_final) # 348843x5
diamond_final$SampleID %>% unique() %>% length() # 6141... 

# What samples are missing?
sample_names <- read_lines("sample_names.txt")
setdiff(sample_names, diamond_final$SampleID) 
# I checked these 6 samples, no hit mets the filtering criteria

# Read file with seqkit stats results
read_depth <- read.csv("seqkit_stats_results_6147_participants.csv")
## filter for forward reads (we are using only the forward reads for the analysis)
read_depth <- read_depth[read_depth$Direction==1,] %>% select(Sample, num_seqs)
## calculate the read depth in millions
read_depth$ReadDepthMillions <- as.numeric(read_depth$num_seqs) / 1e6
read_depth <- read_depth %>% select(Sample, ReadDepthMillions)
summary(read_depth$ReadDepthMillions) 
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 18.76   34.47   52.19   51.09   61.62  358.83 
## make sure all samples are present
setdiff(diamond_final$SampleID, read_depth$Sample) # 0, ok

# read the metadata for the VFs
metadata <- read.csv("VF_metadata.csv", header = T) %>% select(VFID, VF_Name, VFcategory) 
## test if all are present 
setdiff(diamond_final$VFID, metadata$VFID) # 0, ok, good to go

# merge the data (add metadata and read depth to diamond_final)
diamond_final <- left_join(diamond_final, read_depth, by=c("SampleID"="Sample")) # add read depth to the data
diamond_final <- left_join(diamond_final, metadata, by=c("VFID"="VFID")) # add metadata to the data

###############################################################################################
#    RPKM for individual genes (e.g. fimB, fimE, fimA, fimI, fimC, fimD, fimF, fimG, and fimH)
###############################################################################################

# each VFID has multiple DIFFERENT genes (e.g. VF0394 has fliA, motA, cheB etc.) composing it
# in addition, an VFID can have genes with multiple alleles/sequences (e.g. LPS (VF0367) has two acpXL alleles)... we need to concatenate the results for these genes from the same VF (sum their abundances)
## add "Gene" column using the gene_header column 
diamond_final$Gene <- str_extract(diamond_final$gene_header, "(?<=_\\()[^)]+(?=\\))")
## add "Accession" column using the gene_header column
diamond_final$Accession <- str_extract(diamond_final$gene_header, "\\([a-z]{2}\\|[^)]+\\)")
diamond_final$Accession <- gsub("\\(|\\)", "", diamond_final$Accession)
## here you can see an exaple of the VF with genes with multiple alleles
diamond_final %>% ungroup() %>% filter(Gene=="acpXL") %>% select(Gene, Accession, VFID) %>% distinct() 
#    Gene       Accession   VFID
# 1 acpXL gb|WP_002963616 VF0367
# 2 acpXL gb|WP_002963985 VF0367
## you can also have different alleles being used in different VFIDs...
diamond_final %>% ungroup() %>% filter(Gene=="cheA") %>% select(Gene, Accession, VFID) %>% distinct() 
#   Gene       Accession   VFID
# 1 cheA gb|WP_032902679 VF0394
# 2 cheA gb|YP_002343725 VF0114

# calculate RPK (reads per kilobase) - DIVIDE the read counts by the length of each gene in kilobases
diamond_final$GeneLengthNucleotides <- diamond_final$GeneLengthAminoAcids * 3 # get the lenght of the gene in nucleotides
diamond_final$RPK <- diamond_final$NumberOfHighHomologyReads / (diamond_final$GeneLengthNucleotides / 1e3) # calculate RPK for high quality

# calculate the RPKM (RPK/total read count to millions)
diamond_final$RPKM <- (diamond_final$RPK / diamond_final$ReadDepthMillions) 

# calculate TPM (RPK / sum(RPK) * 1e6)
diamond_final <- diamond_final %>% group_by(SampleID) %>% mutate(TPM = (RPK / sum(RPK)) * 1e6) 

# Sum the abundances of the same gene inside the same VFID
# this means correct for VF with genes with multiple alleles
diamond_final <- diamond_final %>%
group_by(SampleID, Gene, VFID,  VF_Name, VFcategory, ReadDepthMillions) %>%
summarise(NumberOfHighHomologyReads = sum(NumberOfHighHomologyReads),
          RPKM = sum(RPKM),
          TPM = sum(TPM))

# save data
save(diamond_final, file="vfdb_hits_normalized_gene_level.RData")

# quick inspection
## sample size with detected VFs
diamond_final %>% ungroup() %>% filter(NumberOfHighHomologyReads>0) %>% select(SampleID) %>% unique() %>% nrow() # 6141 / 6147 = 99.9%
## number of detected VFIDs
diamond_final %>% ungroup() %>% filter(NumberOfHighHomologyReads>0) %>% select(VFID) %>% unique() %>% nrow() # 218
## number of detected genes
diamond_final %>% ungroup() %>% filter(NumberOfHighHomologyReads>0) %>% select(Gene) %>% unique() %>% nrow() # 984

###############################################################################################
#                    RPKM at the VF level (e.g., Type 1 fimbriae, VF0221) 
###############################################################################################

# Sum the abundances at the VFID level
diamond_final <- diamond_final %>%
group_by(SampleID, VFID,  VF_Name, VFcategory, ReadDepthMillions) %>%
summarise(NumberOfHighHomologyReads = sum(NumberOfHighHomologyReads),
          RPKM = sum(RPKM),
          TPM = sum(TPM))

# save data
save(diamond_final, file="vfdb_hits_normalized_VF_level.RData")
```