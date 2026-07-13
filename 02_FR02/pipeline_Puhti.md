# Using Puhit for metagenomic analysis of FR02 samples

The cluster manual is at https://csc-training.github.io/csc-env-eff/

To log in to Puhti, go to https://user-auth.csc.fi/idp/profile/oidc/authorize?execution=e1s3 and click in open terminal.

In the terminal on browser, do:

```bash
# go to your folder in project_2007341
cd /projappl/project_2007341/camila_mericlab

# install kraken2 (2.1.3) and bracken (2.9) using 
# to use conda we need to create a containerized environment using the Tykky wrapper (https://docs.csc.fi/computing/containers/tykky/)
# create a new directory for the installation in the /projappl/<your_project>/... folder
mkdir metagenome_taxonomy
```

Now we need to add the following file "metagenome_taxonomy_env.yml" in the folder

```yml
channels:
  - conda-forge
  - bioconda
dependencies:
  - kraken2=2.1.3
  - bracken=2.9
  - biopython=1.81
```

Now, let's create the container:

```bash
# to request an interactive session you need to specify the CSC project number
# ask one core for 1 h: 
sinteractive --account project_2007341 --time 01:00:00 --mem 20000 # 20 GB

# start containerize
module purge
module load tykky # load the Tykky module
conda-containerize new --mamba --prefix metagenome_taxonomy metagenome_taxonomy_env.yml

# add the <install_dir>/bin directory to $PATH:
 export PATH="/projappl/project_2007341/camila_mericlab/metagenome_taxonomy/bin:$PATH" 
 ```

Now we need to download the Kraken2 indexes:

```bash
# first, install the HPRC indices 
cd /projappl/project_2007341/camila_mericlab/
mkdir databases
cd databases
mkdir k2_HPRC_20230810 && cd k2_HPRC_20230810
wget https://zenodo.org/record/8339732/files/k2_HPRC_20230810.tar.gz
tar -xvzf k2_HPRC_20230810.tar.gz
```

I need to transfer the GTDB r214 indexes from BAKER HPC to Puhti.

I need to add this is to the scratch folder (/scratch/project_2007341/camila_mericlab/databases) because it can hold the size of the database (128 GB).

```bash
screen -S transfering
srun --pty -p long --mem=20G -t 03-00:00:00 /bin/bash
cd /labs/sysgen/workspace/users/camilagv/kraken_snakemake/databases/gtdb_r214/128gb
scp hash.k2d.gz camigazo@puhti.csc.fi:/scratch/project_2007341/camila_mericlab/databases 
```

Now, on Puhti:

```bash
cd /scratch/project_2007341/camila_mericlab/databases
sbatch --account project_2007341 --time=00-03:00:00 --mem=50GB -c 1 --job-name= --wrap="gunzip hash.k2d.gz" 
squeue -u camigazo # check if it is done
# now send the other files...
```

On BAKER HPC:

```bash
screen -S transfering
srun --pty -p long --mem=20G -t 03-00:00:00 /bin/bash
cd /labs/sysgen/workspace/users/camilagv/kraken_snakemake/databases/gtdb_r214/128gb
scp opts.k2d camigazo@puhti.csc.fi:/scratch/project_2007341/camila_mericlab/databases 
scp taxo.k2d camigazo@puhti.csc.fi:/scratch/project_2007341/camila_mericlab/databases 
scp database150mers* camigazo@puhti.csc.fi:/scratch/project_2007341/camila_mericlab/databases 
cd /labs/sysgen/workspace/users/camilagv/kraken_snakemake/utils
scp extract_kraken_reads.py camigazo@puhti.csc.fi:/projappl/project_2007341/camila_mericlab
```

# Classification of metagenomic sequences with Kraken 2 / Bracken and GTDB r204 

Classify the reads with the human pangenome indexes (HPRC) and extract the non-human reads:

```bash
# go to folder 
cd /projappl/project_2007341/camila_mericlab

# create output directory
mkdir temporary_files

# list samples
ls -1 /scratch/project_2007341/DATA/filtered/*.R1.trimmed.filtered.fastq.gz| xargs -n1 basename | sed 's/.R1.trimmed.filtered.fastq.gz//' > sample_names.txt
wc sample_names.txt -l # 8089, ok!

# create 6 batches ony with files starting with 82...
grep "^82" sample_names.txt > b_82_batch.txt
wc b_82_batch.txt -l # 7235 samples, ok!
split --verbose -l1447 b_82_batch.txt --additional-suffix _batch.txt --numeric-suffixes b
```

Add the following script (step1.sh) to the folder. This script will run kraken2 with the HPRC database:

```bash
#!/bin/bash
#step1.sh

# Function to process a single sample with Kraken2
process_sample() {
    local sample=$1

    # run human read classification 
    kraken2 --db /scratch/project_2007341/camila_mericlab/taxonomic_classification/databases/k2_HPRC_20230810/db/ \
            --threads 4 --paired --gzip-compressed \
            /scratch/project_2007341/DATA/filtered/${sample}.R1.trimmed.filtered.fastq.gz \
            /scratch/project_2007341/DATA/filtered/${sample}.R2.trimmed.filtered.fastq.gz \
            --output temporary_files/${sample}.HPRC.kraken

    # run human read extraction
    python /scratch/project_2007341/camila_mericlab/taxonomic_classification/extract_kraken_reads.py \
    -k temporary_files/${sample}.HPRC.kraken -s1 /scratch/project_2007341/DATA/filtered/${sample}.R1.trimmed.filtered.fastq.gz \
    -s2 /scratch/project_2007341/DATA/filtered/${sample}.R2.trimmed.filtered.fastq.gz -t 0 \
    -o temporary_files/nonhuman_${sample}.1.fastq -o2 temporary_files/nonhuman_${sample}.2.fastq --fastq-output

    # count read number from non-human reads
    wc -l temporary_files/nonhuman_${sample}.1.fastq | awk '{print $1/4}' > temporary_files/nonhuman_read_number_${sample}.txt
    
    # run bacteria classification
    kraken2 --db /scratch/project_2007341/camila_mericlab/databases --threads 4 --output temporary_files/${sample}.128gb_kraken_gtdb_output.txt \
    --report temporary_files/${sample}.128gb_kraken_gtdb_report.txt --paired temporary_files/nonhuman_${sample}.1.fastq temporary_files/nonhuman_${sample}.2.fastq \
    --confidence 0.1 --use-names 

    # remove non-human reads and other files not needed
    rm temporary_files/nonhuman_${sample}.1.fastq temporary_files/nonhuman_${sample}.2.fastq
    rm temporary_files/${sample}.HPRC.kraken
    rm temporary_files/${sample}.128gb_kraken_gtdb_output.txt    

    # run bracken
    bracken -d /scratch/project_2007341/camila_mericlab/databases -i temporary_files/${sample}.128gb_kraken_gtdb_report.txt \
    -o temporary_files/${sample}.128gb_kraken_gtdb_bracken_output.txt -r 150 -l S -t 25

    # remove bracken not needed files
    rm temporary_files/${sample}.128gb_kraken_gtdb_report.txt
    rm temporary_files/${sample}.128gb_kraken_gtdb_report_bracken_species.txt
}

# Parse command line arguments
while getopts "s:" arg; do
  case $arg in
    s) sample_list=$OPTARG;;
    *) echo "Usage: $0 -s sample_list.txt"; exit 1;;
  esac
done

# Process each sample
cat $sample_list | while read sample; do
    process_sample "$sample"
done
``` 

Now run:

```bash
# send slurm array script to cluster
module --force purge
module load tykky
export PATH="/scratch/project_2007341/camila_mericlab/taxonomic_classification/metagenome_taxonomy/bin:$PATH" 

# make directory for results
mkdir temporary_files

# run test batch
sbatch --account project_2007341 --time=00-01:00:00 --mem=132GB -c 4 --job-name=test --wrap="bash step1.sh -s b_test_batch.txt" 
rm b_test_batch.txt b_82_batch.txt

# problem with size... I am gonna move it to the scratch folder (I changed the script above step1.sh)
cd /projappl/project_2007341
mv camila_mericlab taxonomic_classification
mv taxonomic_classification /scratch/project_2007341/camila_mericlab
cd /scratch/project_2007341/camila_mericlab/taxonomic_classification

# send batches
for file in b*_batch.txt; do
batch="$(basename -- $file | sed 's/_batch.txt//')"
sbatch --account project_2007341 --time=03-00:00:00 --mem=132GB -c 4 --job-name=${batch} --wrap="bash step1.sh -s $file" 
done

# organize files after run
mkdir misc
mv slur* misc/
```

The run stopped suddenly because of the time limit. 

I need to know the last sample that was processed to continue from there. 

Add the following script (check_batches.sh) to the folder. This script will check the last processed sample in each batch:

```bash
#!/bin/bash
#check_batches.sh

# Directory where the batch files are located
BATCH_DIR="/scratch/project_2007341/camila_mericlab/taxonomic_classification"

# Directory of temporary files
TEMP_DIR="/scratch/project_2007341/camila_mericlab/taxonomic_classification/temporary_files"

# Loop through each batch file
for batch_file in ${BATCH_DIR}/b*_batch.txt; do
    echo "Checking $batch_file"
    last_processed_sample=""

    # Read each sample in the batch file
    while read sample; do
        # Check if the specific file for this sample exists
        if [ -f "${TEMP_DIR}/nonhuman_read_number_${sample}.txt" ]; then
            last_processed_sample=$sample
        else
            # If the file for the current sample doesn't exist, break the loop
            break
        fi
    done < $batch_file

    if [ -n "$last_processed_sample" ]; then
        echo "Last processed sample in $batch_file: $last_processed_sample"
    else
        echo "No samples processed in $batch_file"
    fi
done
```
Lets remove the processed samples from the batch files:

```bash
#!/bin/bash
#clean_batches.sh

# Processed samples from each batch file
declare -A last_processed_samples
last_processed_samples["/scratch/project_2007341/camila_mericlab/taxonomic_classification/b00_batch.txt"]="xxxxxxxx-7"
last_processed_samples["/scratch/project_2007341/camila_mericlab/taxonomic_classification/b01_batch.txt"]="xxxxxxxx-5"
last_processed_samples["/scratch/project_2007341/camila_mericlab/taxonomic_classification/b02_batch.txt"]="xxxxxxxx-3"
last_processed_samples["/scratch/project_2007341/camila_mericlab/taxonomic_classification/b03_batch.txt"]="xxxxxxxx-7"
last_processed_samples["/scratch/project_2007341/camila_mericlab/taxonomic_classification/b04_batch.txt"]="xxxxxxxx-9"

for batch_file in "${!last_processed_samples[@]}"; do
    last_sample=${last_processed_samples[$batch_file]}
    echo "Processing $batch_file, including last processed sample $last_sample"

    # Temporary file to store the new batch file
    temp_file="${batch_file}.tmp"

    # Flag to start copying lines to new file
    start_copy=false

    while read sample; do
        if [ "$start_copy" = true ]; then
            # Append sample to temporary file
            echo "$sample" >> "$temp_file"
        elif [ "$sample" = "$last_sample" ]; then
            # Include the last processed sample and start copying from here
            echo "$sample" >> "$temp_file"
            start_copy=true
        fi
    done < "$batch_file"

    # Replace the original batch file with the new one
    mv "$temp_file" "$batch_file"
    echo "Updated $batch_file"
done
```

Now let me check again:

```bash
# check if the new files are ok
cat b0*_batch.txt | wc -l # 4259... ok! 7235 - 4259 = 2976
ls temporary_files/nonhuman_read_number_* | wc -l # 2981... 2981-5(batches) = 2976, ok!!!
```

I am also scared of finishing our remaining Billing Units (BU)... So I guess I will run bacteria/archea read classification separately to avoid requesting a lot of memory for the whole script...

```bash
#!/bin/bash
#step1a.sh

# Function to process a single sample with Kraken2
process_sample() {
    local sample=$1

    # run human read classification 
    kraken2 --db /scratch/project_2007341/camila_mericlab/taxonomic_classification/databases/k2_HPRC_20230810/db/ \
            --threads 4 --paired --gzip-compressed \
            /scratch/project_2007341/DATA/filtered/${sample}.R1.trimmed.filtered.fastq.gz \
            /scratch/project_2007341/DATA/filtered/${sample}.R2.trimmed.filtered.fastq.gz \
            --output temporary_files/${sample}.HPRC.kraken

    # run human read extraction
    python /scratch/project_2007341/camila_mericlab/taxonomic_classification/extract_kraken_reads.py \
    -k temporary_files/${sample}.HPRC.kraken -s1 /scratch/project_2007341/DATA/filtered/${sample}.R1.trimmed.filtered.fastq.gz \
    -s2 /scratch/project_2007341/DATA/filtered/${sample}.R2.trimmed.filtered.fastq.gz -t 0 \
    -o temporary_files/nonhuman_${sample}.1.fastq -o2 temporary_files/nonhuman_${sample}.2.fastq --fastq-output

    # count read number from non-human reads
    wc -l temporary_files/nonhuman_${sample}.1.fastq | awk '{print $1/4}' > temporary_files/nonhuman_read_number_${sample}.txt
}

# Parse command line arguments
while getopts "s:" arg; do
  case $arg in
    s) sample_list=$OPTARG;;
    *) echo "Usage: $0 -s sample_list.txt"; exit 1;;
  esac
done

# Process each sample
cat $sample_list | while read sample; do
    process_sample "$sample"
done
``` 

Now run:

```bash
# send slurm array script to cluster
module --force purge
module load tykky
export PATH="/scratch/project_2007341/camila_mericlab/taxonomic_classification/metagenome_taxonomy/bin:$PATH" 

# run test batch (xxxxxxxx-8)
ls temporary_files/*xxxxxxxx-8*
sbatch --account project_2007341 --time=00-00:45:00 --mem=40GB -c 4 --job-name=test --wrap="bash step1a.sh -s b_test_batch.txt" 
ls temporary_files/*xxxxxxxx-8*
# finished in less than 1 min, ok!

# send batches
for file in b*_batch.txt; do
batch="$(basename -- $file | sed 's/_batch.txt//')"
sbatch --account project_2007341 --time=00-35:30:00 --mem=40GB -c 4 --job-name=${batch} --wrap="bash step1a.sh -s $file" 
done

# check if all files were processed
ls temporary_files/*2.fastq | wc -l # 4259, ok!
```

Now send the last part that needs more memory:

```bash
#!/bin/bash
#step1b.sh

# Function to process a single sample with Kraken2
process_sample() {
    local sample=$1

    # run bacteria/archea read classification
    kraken2 --db /scratch/project_2007341/camila_mericlab/databases --threads 4 --output temporary_files/${sample}.128gb_kraken_gtdb_output.txt \
    --report temporary_files/${sample}.128gb_kraken_gtdb_report.txt --paired temporary_files/nonhuman_${sample}.1.fastq temporary_files/nonhuman_${sample}.2.fastq \
    --confidence 0.1 --use-names 

    # run bracken
    bracken -d /scratch/project_2007341/camila_mericlab/databases -i temporary_files/${sample}.128gb_kraken_gtdb_report.txt \
    -o temporary_files/${sample}.128gb_kraken_gtdb_bracken_output.txt -r 150 -l S -t 25

}

# Parse command line arguments
while getopts "s:" arg; do
  case $arg in
    s) sample_list=$OPTARG;;
    *) echo "Usage: $0 -s sample_list.txt"; exit 1;;
  esac
done

# Process each sample
cat $sample_list | while read sample; do
    process_sample "$sample"
done
``` 

Send it with:

```bash
# list files to process
ls temporary_files/*2.fastq | xargs -n1 basename | sed 's/.1.fastq//' | sed 's/.2.fastq//' | sed 's/nonhuman_//' > files_to_process.txt
wc files_to_process.txt -l # 4259, ok!

# create 5 batches 
split --verbose -l852 files_to_process.txt --additional-suffix _batch.txt --numeric-suffixes b

# test with 4 samples
head -n4 b00_batch.txt > test.txt
sbatch --account project_2007341 --time=00-00:30:00 --mem=132GB -c 4 --job-name=test --wrap="bash step1b.sh -s test.txt" 
#it took 19 min, or 4.75 min per sample, so... 852 * 4.75 = 4047 min = 67.45 h / 24 = 2.8 days, ok!

# send batches
module --force purge
module load tykky
export PATH="/scratch/project_2007341/camila_mericlab/taxonomic_classification/metagenome_taxonomy/bin:$PATH" 

for file in b*_batch.txt; do
batch="$(basename -- $file | sed 's/_batch.txt//')"
sbatch --account project_2007341 --time=03-00:00:00 --mem=132GB -c 4 --job-name=${batch} --wrap="bash step1b.sh -s $file" 
done

# check if all files were processed
ls *.128gb_kraken_gtdb_bracken_output.txt | wc -l # 7226-7235=9 are missing... but some simply didnt have any reads assigned to bacteria after human removal
```

Now, I will create the phyloseq object ...

```bash
# to request an interactive session you need to specify the CSC project number
# ask one core for 1 h: 
sinteractive --account project_2007341 --time 01:00:00 --mem 20000 # 20 GB

# start R
module purge
module load r-env
R
```

```R
# library
library(phyloseq)
library(tidyverse)

# list files
files <- list.files(pattern = ".128gb_kraken_gtdb_bracken_output.txt") 

# for each file read and merge
merged_bracken_output <- data.frame()
for(f in files){
    print(f)
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
            select(-matches("kraken_assigned_reads|added_reads|frac|taxonomy_id|taxonomy_lvl")) %>% 
            rename(Species = name) 
    sample<- gsub(".128gb_kraken_gtdb_bracken_output.txt", "", f)
    df$Sample <- sample
    merged_bracken_output <- rbind(merged_bracken_output, df)
}

merged_bracken_output <- merged_bracken_output %>%
                pivot_wider(names_from = Sample, values_from = new_est_reads)

# read the taxonomy auxiliary file
taxonomy <- read.table("gtdb_v214_taxonomy.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# add the full taxonomy to the merged bracken output 
merged_bracken_output <- merge(taxonomy, merged_bracken_output, all.y = TRUE)

# replace NA with 0
# merged_bracken_output[is.na(merged_bracken_output)] <- 0

# write file to taxonomy_tabulated_output.tsv
write.table(merged_bracken_output, file="taxonomy_tabulated_output.tsv", sep = "\t", row.names = FALSE)

# create otu table and tax table
otu <- merged_bracken_output %>% 
select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, Species) %>%
column_to_rownames("Species")

tax <- merged_bracken_output %>% 
select(Domain, Phylum, Class, Order, Family, Genus, Species) 
rownames(tax) <- tax$Species
 
# create phyloseq object
physeq <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE), 
                   tax_table(as.matrix(tax)))

# write phyloseq object
saveRDS(physeq, "phyloseq.rds")

# organize read numbers
path <- "/scratch/project_2007341/camila_mericlab/taxonomic_classification/temporary_files/read_number"
read_depth_nonhuman_files <- list.files(path=path, pattern = "nonhuman_read_number_*")
read_depth_nonhuman <- data.frame()
for (i in 1:length(read_depth_nonhuman_files)) {
  sample <- gsub("nonhuman_read_number_(.*).txt", "\\1", read_depth_nonhuman_files[i])
  cat(sample, "\n")
  read_number <- read.table(paste0(path,"/",read_depth_nonhuman_files[i]))
  read_number <- read_number$V1
  read_depth_nonhuman <- rbind(read_depth_nonhuman, cbind("Sample"=sample, "Non_Human_Read_Number"=read_number))
}

reads_assigned <- as.data.frame(sample_sums(physeq))
colnames(reads_assigned) <- c("Reads_Classified")
reads_assigned$Sample <- rownames(reads_assigned)

read_track <- merge(read_depth_nonhuman, reads_assigned, by = "Sample")
read_track$Percent_Reads_Classified <- read_track$Reads_Classified / as.numeric(read_track$Non_Human_Read_Number) * 100

write.table(read_track, file="sanity_check_read_track.txt", sep = "\t", row.names = FALSE)
```

Now that I finish, I am going to transfer the files to ATLAS to continue with the analysis...

```bash
# go to foder on puhti
cd /scratch/project_2007341/camila_mericlab/taxonomic_classification

# send files to ATLAS 
scp phyloseq.rds cvolpian@atlas.fimm.fi:/csc/fr_metagenome/workspaces/camilagv/2024
scp sanity_check_read_track.txt cvolpian@atlas.fimm.fi:/csc/fr_metagenome/workspaces/camilagv/2024
scp taxonomy_tabulated_output.tsv cvolpian@atlas.fimm.fi:/csc/fr_metagenome/workspaces/camilagv/2024
```

