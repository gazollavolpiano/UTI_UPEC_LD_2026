# Profile E. coli strains with StrainScan

One caveat is that the StrainScan tool does not run with multiple threads to profile the strains (it can only be used to build the database).

```bash
# load conda module
module load conda
source conda_init.sh

# install tool (strainscan 1.0.14)
# conda create -n strainscan -c bioconda strainscan
conda activate strainscan
strainscan -h
conda list # 1.0.14 is installed

# I added DB_Ecoli.zip (1.4G) downloaded from Google drive to the wharf folder on Bianca
# now, let's transfer it to the analysis folder
mv /proj/nobackup/simp2024014/wharf/cgazolla/cgazolla-simp2024014/DB_Ecoli.zip . 
unzip DB_Ecoli.zip
rm DB_Ecoli.zip

# list samples
ls -1 /proj/simp2024014/Omicsdataleverans/metagnomics/*/*/host_removal/*__1.fq.gz| xargs -n1 basename | sed 's/__1.fq.gz//' > sample_names.txt
wc -l sample_names.txt # 6147
```

The problem with BIANCA cluster is that sometimes we get permission denied errors when processing the samples. We gonna create a script that reads the sample names from the list and loops indefinitely. For each sample, if the corresponding file (using your directory and naming pattern) exists and is readable, it submits an sbatch job with the specified resources. Once submitted, that sample is not checked again.

```bash
#!/bin/bash
# check_samples.sh
# This script continuously monitors sample availability.
# When both forward and reverse sample files become accessible,
# it submits a job to process that sample.
#
# Usage: bash check_samples.sh -s sample_names.txt

# Parse command-line arguments
while getopts "s:" arg; do
  case $arg in
    s) sample_list="$OPTARG" ;;
    *) echo "Usage: $0 -s sample_list"; exit 1 ;;
  esac
done

if [ -z "$sample_list" ]; then
  echo "Usage: $0 -s sample_list"
  exit 1
fi

# Read sample names into an array
mapfile -t samples < "$sample_list"

# Declare an associative array to track submitted samples
declare -A submitted

# Continuous monitoring loop
while true; do
  for sample in "${samples[@]}"; do
    # Skip if this sample has already been submitted
    if [[ ${submitted["$sample"]} ]]; then
      continue
    fi

    # Construct the expected sample file paths for forward and reverse reads
    forward_file=$(echo /proj/simp2024014/Omicsdataleverans/metagnomics/*/*/host_removal/${sample}__1.fq.gz)
    reverse_file=$(echo /proj/simp2024014/Omicsdataleverans/metagnomics/*/*/host_removal/${sample}__2.fq.gz)

    # Check if both the forward and reverse sample files exist and are readable
    if [ -f "$forward_file" ] && [ -r "$forward_file" ] && [ -f "$reverse_file" ] && [ -r "$reverse_file" ]; then
      echo "$(date +"%T") - Sample $sample is available (both forward and reverse reads). Submitting job..."
      # Submit a job to process this single sample using strainscan_defaultparam_singlesample.sh
      sbatch -A simp2024014 --partition core --time=00-03:00:00 --mem=78GB --job-name=strain_scan_${sample} --wrap="bash strainscan_defaultparam_singlesample.sh -s $sample"
      # Mark sample as submitted so it is not re-submitted
      submitted["$sample"]=1
    else
      echo "$(date +"%T") - Sample $sample not available yet (missing forward or reverse read)."
    fi
  done
  # Sleep for 5 minutes before checking again 
  sleep 300
done
```

This script accepts a single sample name (via -s) and runs the strainscan command on its corresponding file. It checks that the file exists.

```bash
#!/bin/bash
# strainscan_defaultparam_singlesample.sh
# This script processes a single sample using strainscan.
#
# Usage: bash strainscan_defaultparam_singlesample.sh -s sample_name

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
sample_file=$(echo /proj/simp2024014/Omicsdataleverans/metagnomics/*/*/host_removal/${sample}__1.fq.gz)
sample_file_rev=$(echo /proj/simp2024014/Omicsdataleverans/metagnomics/*/*/host_removal/${sample}__2.fq.gz)

if [ ! -f "$sample_file" ] || [ ! -r "$sample_file" ]; then
  echo "$(date +"%T") - Error: Sample file $sample_file not found or is not readable."
  exit 1
fi

# Run strainscan on the sample
echo "$(date +"%T") - Processing sample $sample with strainscan..."
strainscan -i "$sample_file" -j "$sample_file_rev" -d /proj/simp2024014/Ecoli_StrainScan/DB_Ecoli --output_dir "strainscan_output/${sample}" 

if [ $? -eq 0 ]; then
  echo "$(date +"%T") - Finished processing sample $sample."
else
  echo "$(date +"%T") - Error processing sample $sample."
fi
```

```bash
# create the output directory
mkdir strainscan_output

# activate conda environment
module load conda
source conda_init.sh
conda activate strainscan

# monitoring script so that it continuously checks sample availability and submits jobs
sbatch -A simp2024014 --partition core --time=10-00:00:00 --mem=5GB -c 1 --job-name=check_samples --wrap="bash check_samples.sh -s sample_names.txt"
```