# StrainScan with E. coli pre-built database 

Use the script `strainscan_batch.sh`. It uses some parameters to increase the sensitivity for FINRISK shallow sequeincing data:
- "--strain_prob 1": if this parameter is set to 1, then the algorithm will output the probabolity of detecting a strain (or cluster) in low-depth (e.g. <1x) samples
- "--low_dep 2": this parameter can be set to "1" if the sequencing depth of input data is very low (e.g. < 5x). For super low depth ( < 1x ), you can use "-l 2" (default: -l 0)

```bash
#!/bin/bash
# strainscan_batch.sh

# parse command-line arguments
while getopts "s:" arg; do
  case $arg in
    s) sample_list=$OPTARG;;
    *) echo "Usage: $0 -s sample_list"; exit 1;;
  esac
done

# create output dir if it doesn't exist
mkdir -p strainscan_output

# process each sample in the list
while IFS= read -r line; do
  # run StrainScan
  python StrainScan.py -i "/microbiome2018-02-26/qc/${line}/filtered/${line}.R1.trimmed.filtered.fastq.gz" \
  -j "/microbiome2018-02-26/qc/${line}/filtered/${line}.R2.trimmed.filtered.fastq.gz" \
  -d DB_Ecoli --output_dir "strainscan_output/${line}" --strain_prob 1 --low_dep 2 
  # check if command was successful
  if [ $? -eq 0 ]; then
    echo "$(date +"%T"), finished $line..." >> strainscan_batch.log
  else
    echo "$(date +"%T"), error processing $line" >> strainscan_batch.log
  fi
done < "$sample_list"
```

Execute it with:

```bash
# Use the script with list of sample names 
bash strainscan_batch.sh -s list_samples_batch.txt 
```

Use R to concatenate the results and get the distribution of strains in the samples, including the annotation of the STs and the VF content of the STs. 