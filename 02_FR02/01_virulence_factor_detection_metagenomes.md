# Detection of virulence factors 

We are going to run the detection of virulence factors based on diamond and the VFDB database. The VFDB database has 4193 protein sequences curated from bacterial pathogens. The diamond software is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data.

Add the loop scrip vfdb_diamond_array.sh to the folder. The scripts uses the following parameters:
* --max-target-seqs 1 --max-hsps 1: limits the search to the single best aligning sequence (--max-target-seqs 1) and the single best high-scoring segment pair (--max-hsps 1) per query sequence.
* --sensitive: uses the sensitive mode of DIAMOND, which is slower but more accurate.
* IMPORTANT: do not use --subject-cover, --id, and --query-cover parameters, they simply do not work as expected. Instead, we are going to filter the results in the R script.

```bash
#!/bin/bash
#vfdb_diamond_batch.sh

# pull one single line/filename from the batch file with sample names
while getopts "s:" arg; do
  case $arg in
    s) sample_list=$OPTARG;;
  esac
done

cat $sample_list | while read line 
do
  # run diamond with vfdb
  /homes/cvolpian/micromamba/envs/diamond/bin/diamond blastx --db VFDB_setA_pro.dmnd --query /csc/fr_metagenome/microbiome/microbiome2018-02-26/qc/${line}/filtered/${line}.R1.trimmed.filtered.fastq.gz \
  --out "vfdb_diamond_output/${line}_vfdb_diamond.txt" \
  --outfmt 6 qtitle stitle pident length sstart send slen qlen \
  --max-target-seqs 1 --max-hsps 1 --sensitive --threads 4
  # add sample to a log file
  echo "$(date +"%T"), finished $line..." >> vfdb_diamond_log.log
done
```

Now, run the batch script:

```bash
# create a directory to store the results
mkdir vfdb_diamond_output

# create batch files (9 batches)
s=$(wc -l < sample_names.txt); b=8; n=$((s / b)); split --verbose -l $n sample_names.txt --additional-suffix _batch.txt --numeric-suffixes b

# create scripts for each batch (medmem.q)
for file in b0*_batch.txt; do
batch="$(basename -- $file | sed 's/_batch.txt//')"
header='#!/bin/bash \n# run_'${batch}'.sh \n#$ -q medmem.q \n#$ -V \n#$ -N diamond_'${batch}' \n#$ -l mem_total=40g \n#$ -pe make 4 \n#$ -cwd \n#$ -j y \n#$ -o '${batch}'.log \n \n'
exec='sh vfdb_diamond_batch.sh -s '$file' '
echo -e $header $exec > run_${batch}.sh
done

# execute the batches
qsub run_b00.sh # medmem.q, 7892929
qsub run_b01.sh # medmem.q, 7892930
qsub run_b02.sh # medmem.q, 7892931
qsub run_b03.sh # medmem.q, 7892932
qsub run_b04.sh # medmem.q, 7892933
qsub run_b05.sh # medmem.q, 7892934
qsub run_b06.sh # medmem.q, 7892935
qsub run_b07.sh # medmem.q, 7892936
qsub run_b08.sh # medmem.q, 7892937
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
library("IRanges",warn.conflicts=FALSE)

# list the files in the folder
files <- list.files(path= "vfdb_diamond_output", pattern = "_vfdb_diamond.txt", full.name=TRUE)

# read each file with result and process statistics
diamond_final <- c() 
for (file in files){
    # get sample name
    sample <- gsub(".*_diamond_output/(.*)_.*_diamond.txt","\\1",file)
    cat("processing", sample, "\n")

    # read file and add column names
    gene_hits <- tryCatch(read.delim(file, header=F), error=function(e) NULL)
    if(is.null(gene_hits)){next} # if the file is empty, go to the next file
    colnames(gene_hits) <- c("qtitle", "stitle", "pident", "length", "sstart", "send", "slen", "qlen")

    # filter for subject coverage ≥80% (common for both Homologous and High-Quality groups)
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

    # filter data for "Homologous" group
    ## remove hits with less than 50% of identity (pident)
    gene_hits_homologous <- gene_hits %>% filter(pident>=50) 
    if(dim(gene_hits_homologous)[1]==0){next} # if the file is empty now, go to the next file (High-Quality group will be empty too)
    ## calculate the gene coverage for each gene 
    gene_hits_homologous <- calc_coverage(gene_hits_homologous)
    ## filter hits with less than 50% of coverage
    gene_hits_homologous <- gene_hits_homologous[gene_hits_homologous$gene_cov>=50,]
    if(dim(gene_hits_homologous)[1]==0){next} # if the file is empty now, go to the next file (High-Quality group will be empty too)
    ## rename columns and filter columns
    gene_hits_homologous <- gene_hits_homologous %>% select(gene_header=1, GeneLengthAminoAcids=2, NumberOfHomologousReads=4)

    # filter data for "High-Quality" group
    ## remove hits with less than 90% of identity (pident)
    gene_hits_high_quality <- gene_hits %>% filter(pident>=90)
    ## calculate the gene coverage for each gene
    gene_hits_high_quality <- calc_coverage(gene_hits_high_quality)
    ## filter hits with less than 50% of coverage
    gene_hits_high_quality <- gene_hits_high_quality[gene_hits_high_quality$gene_cov>=50,]

    if(dim(gene_hits_high_quality)[1]==0){
      # set gene_hits with only the homologous hits (High-Quality group will be empty)
      gene_hits <- gene_hits_homologous
      ## add NumberOfHighQualityReads column with 0
      gene_hits$NumberOfHighQualityReads <- 0
      } else {
        ## rename columns and filter columns
        gene_hits_high_quality <- gene_hits_high_quality %>% select(gene_header=1, GeneLengthAminoAcids=2, NumberOfHighQualityReads=4)

        # tidy final results
        ## combine data
        gene_hits <- full_join(gene_hits_homologous, gene_hits_high_quality, by = c("gene_header", "GeneLengthAminoAcids"))
        ## add 0 instead of NA
        gene_hits$NumberOfHighQualityReads[is.na(gene_hits$NumberOfHighQualityReads)] <- 0
        gene_hits$NumberOfHomologousReads[is.na(gene_hits$NumberOfHomologousReads)] <- 0
      }

    # add the VFID column
    gene_hits$VFID <- gsub(".*_\\((VF\\d+)\\).*", "\\1", gene_hits$gene_header)

    # add the sample name
    gene_hits$SampleID <- sample

    # organize columns order
    gene_hits <- gene_hits %>% select(SampleID, VFID, GeneLengthAminoAcids, NumberOfHomologousReads, NumberOfHighQualityReads, gene_header)

    # append data
    diamond_final <- rbind(diamond_final, gene_hits)
}

# save the data
save(diamond_final, file="vfdb_hits.RData")

# print that the script finished
cat("The script finished!")
```

Run it with:

```bash
# run the R script
header='#!/bin/bash \n# virulence_factor_process.sh \n#$ -q medmem.q \n#$ -V \n#$ -N virulence_factor_process \n#$ -l mem_total=40g \n#$ -pe make 4 \n#$ -cwd \n#$ -j y \n#$ -o virulence_factor_process.log \n \n'
exec='Rscript virulence_factor_process.R'
echo -e $header $exec > virulence_factor_process.sh
micromamba activate R
qsub virulence_factor_process.sh # medmem.q, 7894344
```

# Calculate RPKM and add metadata: create vfdb_hits_normalized.RData

We need to calculate the RPKM for each gene and also add additional metadata for the genes.

For this, we will need to provide: 
* file with the number of reads for each sample
* file with metadata for the genes (VF_metadata.csv)

This will create the file vfdb_hits_normalized.RData, which is the final file with the results.

```R
# library
library(tidyverse)

# read the data with the hits ("diamond_final")
load("vfdb_hits.RData")

# read the read depth file
read_depth <- read.delim("/csc/fr_metagenome/microbiome/microbiome2018-02-26/qc/multiQC_per_sample/multiqc_data/multiqc_fastqc.txt", header = T, sep = "\t") %>% select(Sample, Total.Sequences)
## filter for forward reads (we are using only the forward reads for the analysis)
read_depth$direction <- gsub(".*.R1.trimmed.filtered.fastq.gz", "forward", read_depth$Sample)
read_depth <- read_depth[read_depth$direction=="forward",] %>% select(Sample, Total.Sequences)
## correct the sample names
read_depth$SampleID <- gsub(".R1.trimmed.filtered.fastq.gz", "", read_depth$Sample) 
## calculate the read depth in millions
read_depth$ReadDepthMillions <- as.numeric(read_depth$Total.Sequences) / 1e6
read_depth <- read_depth %>% select(SampleID, ReadDepthMillions)
median(read_depth$ReadDepthMillions) # 0.750678, as expected
## make sure all samples are present
setdiff(diamond_final$SampleID, read_depth$SampleID) # 0, ok

# read the metadata for the genes
metadata <- read.csv("misc/VF_metadata.csv", header = T) %>% select(VFID, VF_Name, VFcategory) 
## test if all genes are present 
setdiff(diamond_final$VFID, metadata$VFID) # 0, ok, good to go

# merge the data (add metadata and read depth to diamond_final)
diamond_final <- left_join(diamond_final, read_depth, by="SampleID") # add read depth to the data
diamond_final <- left_join(diamond_final, metadata, by=c("VFID"="VFID")) # add metadata to the data

# calculate RPK (reads per kilobase) - DIVIDE the read counts by the length of each gene in kilobases
diamond_final$GeneLengthNucleotides <- diamond_final$GeneLengthAminoAcids * 3 # get the lenght of the gene in nucleotides
diamond_final$RPK_Homologous <- diamond_final$NumberOfHomologousReads / (diamond_final$GeneLengthNucleotides / 1e3) # calculate RPK for homologous
diamond_final$NumberOfHighQualityReads <- as.numeric(diamond_final$NumberOfHighQualityReads)
diamond_final$RPK_HighQuality <- diamond_final$NumberOfHighQualityReads / (diamond_final$GeneLengthNucleotides / 1e3) # calculate RPK for high quality

# calculate the RPKM (RPK/total read count to millions)
diamond_final$RPKM_Homologous <- (diamond_final$RPK_Homologous / diamond_final$ReadDepthMillions) # calculate RPKM for homologous
diamond_final$RPKM_HighQuality <- (diamond_final$RPK_HighQuality / diamond_final$ReadDepthMillions) # calculate RPKM for high quality

# each VFID has multiple DIFFERENT genes (e.g. VF0394 has fliA, motA, cheB etc.) composing it
# in addition, an VFID can have genes with multiple alleles/sequences (e.g. LPS (VF0367) has two acpXL alleles)... we need to concatenate the results for these genes from the same VF (sum their abundances)
## add "Gene" column using the gene_header column 
diamond_final$Gene <- str_extract(diamond_final$gene_header, "(?<=_\\()[^)]+(?=\\))")
## add "Accession" column using the gene_header column
diamond_final$Accession <- str_extract(diamond_final$gene_header, "\\([a-z]{2}\\|[^)]+\\)")
diamond_final$Accession <- gsub("\\(|\\)", "", diamond_final$Accession)
## save a file with the genes, accessions and gene_header just for inspection (or if someone has a question)
diamond_final %>% ungroup() %>% select(Gene, Accession, gene_header) %>% distinct() %>% write.csv("misc/genes_accessions_and_headers.csv", row.names = F)
## here you can see an exaple of the VF with genes with multiple alleles
diamond_final %>% ungroup() %>% filter(Gene=="acpXL") %>% select(Gene, Accession, VFID) %>% distinct() 
## you can also have different alleles being used in different VFIDs...
diamond_final %>% ungroup() %>% filter(Gene=="cheA") %>% select(Gene, Accession, VFID) %>% distinct() 
## sum the abundances of the same gene inside the same VFID
diamond_final <- diamond_final %>%
group_by(SampleID, Gene, VFID,  VF_Name, VFcategory) %>%
summarise(NumberOfHomologousReads = sum(NumberOfHomologousReads),
NumberOfHighQualityReads = sum(NumberOfHighQualityReads),
RPKM_Homologous = sum(RPKM_Homologous),
RPKM_HighQuality = sum(RPKM_HighQuality))

# save data
save(diamond_final, file="vfdb_hits_normalized.RData")
```