
########Step 1: First FastQC (Raw reads quality control)
#!/usr/bin/env bash
# Enable strict error handling and debugging
set -euxo pipefail
set -o errtrace

# Input and output directories
INPUT_DIR="#Directory where input data (*_1.fq.gz / *_2.fq.gz or *_fwd.fq.gz / *_rev.fq.gz is)"
OUTPUT_DIR=""

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FastQC on files starting with "RUB"
for file in $INPUT_DIR/MAR*.fq.gz; do
fastqc -o "$OUTPUT_DIR" "$file"
done
echo "first qc for raw data finished"


#######step2: Running SortMeRNA (rRNA removal)
#!/usr/bin/env bash

set -euxo pipefail
set -o errtrace

######downloading Sortmerna database 
wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.4/database.tar.gz

###indexing sortmerna database so the mapping will be easier and faster
#!/usr/bin/env bash

set -euxo pipefail
set -o errtrace

./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/rfam-5.8s-database-id98.fasta --index 1 
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/rfam-5s-database-id98.fasta --index 1
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-arc-16s-id95.fasta --index 1
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-arc-23s-id98.fasta --index 1
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-bac-16s-id90.fasta --index 1
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-bac-23s-id98.fasta --index 1
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-euk-18s-id95.fasta --index 1
./sortmerna --ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-euk-28s-id98.fasta --index 1



sortmerna_path="/home/mozhgan/miniconda3/envs/bumble/bin/sortmerna"

# Directory containing sample files
input_directory="/media/mozhgan/bgi/halictid_data_andreia/other_halictids/concat_rawreads_240731"
sample_output="/media/mozhgan/bgi//raw/007_marginatum_andreia/02_sortmerna"

# Summary file path
summary_file="${sample_output}/summary.txt"

# Initialize summary file
if [[ ! -f "$summary_file" ]]; then
echo "Sample\tTotal Reads Failing E-value Threshold (%)" > "$summary_file"
fi

# Find samples starting with 'RUB' and handle both forward and reverse files
samples=($(ls ${input_directory} | grep "^MAR" | sed 's/_concatenated_[12].fq.gz//g' | uniq))

# Iterate over samples
for sample in "${samples[@]}"; do
# Construct forward and reverse file names
input_file1="${input_directory}/${sample}_concatenated_1.fq.gz"
input_file2="${input_directory}/${sample}_concatenated_2.fq.gz"

# Check if both files exist
if [[ ! -f "$input_file1" || ! -f "$input_file2" ]]; then
echo "Missing paired files for sample: $sample. Skipping."
continue
fi

# Create rrna_trimmed folder if it doesn't exist
mkdir -p "${sample_output}"

# Delete existing key-value DB directory
rm -rf "/home/mozhgan/sortmerna/run/kvdb"

# Run SortMeRNA with --aligned option to save the aligned log file
"$sortmerna_path" \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/rfam-5.8s-database-id98.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/rfam-5s-database-id98.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-arc-16s-id95.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-arc-23s-id98.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-bac-16s-id90.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-bac-23s-id98.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-euk-18s-id95.fasta \
--ref /home/mozhgan/miniconda3/envs/bumble/bin/sortmerna_ref/silva-euk-28s-id98.fasta \
--reads "$input_file1" \
--reads "$input_file2" \
--out2 \
--fastx \
--idx-dir /home/mozhgan/sortmerna/run/idx \
--other "${sample_output}/${sample}" \
--paired_in \
--threads 111 \
--aligned "${sample_output}/${sample}.aligned.log" \
|& tee "${sample_output}/${sample}.log"

# Extract percentage of reads failing the E-value threshold from the log file
failing_percentage=$(grep "Total reads failing E-value threshold" "${sample_output}/${sample}.aligned.log.log" | awk '{print $NF}' | sed 's/[()]//g')

# Append the sample name and failing percentage to the summary file
echo -e "${sample}\t${failing_percentage}" >> "$summary_file"
done

echo "sortmerna finished"

#step3: FastQC (Post rRNA filtering QC)
#!/usr/bin/env bash

# Enable strict error handling and debugging
set -euxo pipefail
set -o errtrace
fastqc="/usr/local/bin/fastqc"
# Input and output directories
INPUT_DIR="/media/mozhgan/bgi//raw/007_marginatum_andreia/02_sortmerna"
OUTPUT_DIR="/media/mozhgan/bgi//raw/007_marginatum_andreia/03_qc_after_sortmerna"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/MAR*.fq.gz; do
if [[ ! "$file" =~ log ]]; then
fastqc -o "$OUTPUT_DIR" "$file"
fi
done


echo "qc for sortmerna finished."

###step4: Trimmomatic (Adapter and quality trimming)
#!/bin/bash

# Enable error checking and tracing
set -euxo pipefail
set -o errtrace
trimmomatic="/usr/local/bin/trimmomatic"
# Input directory containing all *_fwd.fq.gz and *_rev.fq.gz files
samples_directory="/media/mozhgan/bgi//raw/007_marginatum_andreia/02_sortmerna"

# Output directory for trimmed files and logs
output_directory="/media/mozhgan/bgi//raw/007_marginatum_andreia/04_trimmed"

# Report file to save trimming summary
report_file="$output_directory/trimming_summary.txt"

# Create output directory if it does not exist
mkdir -p "$output_directory"

# Initialize the report file with a header (only if it doesn't already exist)
if [[ ! -f "$report_file" ]]; then
echo -e "Sample\tReads_Before\tReads_After\tPercentage_Trimmed\tPercentage_Remaining" > "$report_file"
fi

# Find all *_fwd.fq.gz files and iterate over them
for file_1 in "$samples_directory"/*_fwd.fq.gz; do
# Derive the base name without _fwd.fq.gz
base_sample=$(basename "$file_1" _fwd.fq.gz)

# Construct the reverse file name
file_2="${samples_directory}/${base_sample}_rev.fq.gz"

# Check if both forward and reverse files exist
if [[ -f "$file_1" && -f "$file_2" ]]; then
# Run Trimmomatic on the pair of fastq.gz files
trimmomatic PE \
-threads 25 \
"$file_1" \
"$file_2" \
"$output_directory/${base_sample}_trimmomatic_1.fq.gz" \
"$output_directory/${base_sample}_unpaired_1.fq.gz" \
"$output_directory/${base_sample}_trimmomatic_2.fq.gz" \
"$output_directory/${base_sample}_unpaired_2.fq.gz" \
ILLUMINACLIP:/home/mozhgan/miniconda3/envs/bumble/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:2:TRUE \
SLIDINGWINDOW:5:5 \
LEADING:25 \
TRAILING:25 \
MINLEN:50 |& tee "$output_directory/${base_sample}.log"

# Calculate reads before and after trimming
pre_trim_count=$(zcat "$file_1" | wc -l | awk '{print $1/4}')
post_trim_count=$(zcat "$output_directory/${base_sample}_trimmomatic_1.fq.gz" | wc -l | awk '{print $1/4}')

# Calculate percentages
percentage_trimmed=$(awk "BEGIN {print (($pre_trim_count - $post_trim_count) / $pre_trim_count) * 100}")
percentage_remaining=$(awk "BEGIN {print ($post_trim_count / $pre_trim_count) * 100}")

# Append the results to the trimming summary
echo -e "${base_sample}\t${pre_trim_count}\t${post_trim_count}\t${percentage_trimmed}\t${percentage_remaining}" >> "$report_file"
else
  echo "Reverse file for sample ${base_sample} does not exist. Skipping."
fi
done


echo "trimmomatic finished"


#step5: FastQC (Post trimming QC)
#!/usr/bin/env bash

# Enable strict error handling and debugging
set -euxo pipefail
set -o errtrace
fastqc="/usr/local/bin/fastqc"
# Input and output directories
INPUT_DIR="/media/mozhgan/bgi//raw/007_marginatum_andreia/04_trimmed"
OUTPUT_DIR="/media/mozhgan/bgi//raw/007_marginatum_andreia/05_qc_after_trimming"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/MAR*_trimmomatic_*.fq.gz; do
if [[ ! "$file" =~ log ]]; then
fastqc -o "$OUTPUT_DIR" "$file"
fi
done

echo "qc for trimmomatic finished"


####step6: HISAT2 (Alignment to reference genome)

#!/usr/bin/env bash

set -euxo pipefail
set -o errtrace

###downloading the genome of L. malachurum
wget https://beenomes.princeton.edu/wp-content/uploads/2021/04/LMAL.zip

###indexing the genome to make the mapping process easier
#you have to be in the same directory as hisat2 is installed

./hisat2-build /media/mozhgan/bgi//raw/genome/LMAR_genome_v2.1.1.fasta mgt_index

# Directories
samtools="/usr/local/bin/samtools"
samples_directory="/media/mozhgan/bgi//raw/007_marginatum_andreia/04_trimmed"
log_directory="/media/mozhgan/bgi//raw/007_marginatum_andreia/06_sam/sam_log"
sam_directory="/media/mozhgan/bgi//raw/007_marginatum_andreia/06_sam"
hisat2_executable="/media/mozhgan/bgi//raw/hisat2-2.2.1/hisat2"
genome_index="/media/mozhgan/bgi//raw/hisat2-2.2.1/mgt_index"

# Create directories if they don't exist
mkdir -p "$log_directory"
mkdir -p "$sam_directory"

# Get the list of samples
#samples=($(find "$samples_directory" -type f -name "*_trimmomatic_1.fq.gz" | sed 's/_trimmomatic_1.fq.gz//' | xargs -n 1 basename))
samples=($(find "$samples_directory" -name "*_trimmomatic_1.fq.gz" ! -name "*aligned*" | sed 's/.*\/\(.*\)_trimmomatic_1\.fq\.gz/\1/' | sort | uniq))


# Run HISAT2 for each sample
for sample in "${samples[@]}"; do
(
  "$hisat2_executable" \
  -p 40 -q -x "$genome_index" \
  -1 "${samples_directory}/${sample}_trimmomatic_1.fq.gz" \
  -2 "${samples_directory}/${sample}_trimmomatic_2.fq.gz" \
  -S "${sam_directory}/${sample}.sam"
) > "${log_directory}/${sample}.log" 2>&1
done

# Generate summary.txt in the log directory
summary_file="${log_directory}/summary.txt"
echo -e "Sample\tMapped Reads\tTotal Reads\tMapping Percentage (%)" > "$summary_file"

for sample in "${samples[@]}"; do
# Extract relevant numbers from the log file
total_reads=$(grep "reads; of these:" "${log_directory}/${sample}.log" | awk '{print $1}')
mapped_reads=$(grep "overall alignment rate" "${log_directory}/${sample}.log" | awk '{print $1}' | sed 's/%//')
mapped_count=$(echo "scale=0; ($total_reads * $mapped_reads) / 100" | bc)

# Calculate percentage, exclude the aligned reads
aligned_reads=$(grep -E "aligned 0 times" "${log_directory}/${sample}.log" | wc -l)
mapped_reads_adjusted=$(($mapped_count - $aligned_reads))

# Append to summary file
echo -e "${sample}\t${mapped_reads_adjusted}\t${total_reads}\t$(echo "scale=2; ($mapped_reads_adjusted / $total_reads) * 100" | bc)" >> "$summary_file"
done

echo "Mapping summary saved in: $summary_file"
echo "sam finished"



#step7: Samtools (SAM to BAM conversion)
#!/usr/bin/env bash

set -euxo pipefail
set -o errtrace

### CODE 1: Convert SAM to BAM ###
samtools="/usr/local/bin/samtools"
# Directory containing input SAM files
input_dir="/media/mozhgan/bgi//raw/007_marginatum_andreia/06_sam"

# Directory for output BAM files
output_dir="/media/mozhgan/bgi//raw/007_marginatum_andreia/07_bam"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over all SAM files in the input directory
for sam_file in "$input_dir"/*.sam; do
# Extract the base name (without extension) for each file
sample=$(basename "$sam_file" .sam)

# Convert SAM to BAM using samtools
samtools view -@ 111 -F 4 -b "$sam_file" > "${output_dir}/${sample}.bam"
done

echo "SAM to BAM conversion completed."

#####step 8: HTSeq-count (Gene expression quantification)

# Specify the paths to BAM files, GFF file, output directory, and log files
bam_dir="/media/mozhgan/bgi//raw/007_marginatum_andreia/07_bam"
gff_file="/media/mozhgan/bgi//raw/genome2/LMAR_OGS_v2.1.1.gff3"
htseq_count="/usr/local/bin/htseq-count"
output_dir="/media/mozhgan/bgi//raw/007_marginatum_andreia/08_htseq_count"
log_dir="$output_dir/logs"
summary_file="$output_dir/mapping_rates.txt"

# Create the output and log directories if they don't exist
mkdir -p "$output_dir"
mkdir -p "$log_dir"

# Initialize the summary file
echo "Sample Name,Total Reads,Mapped Reads,Mapping Rate (%)" > "$summary_file"

# Loop through all BAM files in the BAM directory
for bam_file in "$bam_dir"/*.bam; do
sample=$(basename "$bam_file" .bam)
output_file="${output_dir}/${sample}.count"
log_file="$log_dir/${sample}_log.txt"

# Check if the BAM file exists
if [[ ! -f "$bam_file" ]]; then
echo "BAM file for sample $sample does not exist: $bam_file" | tee -a "$log_file"
continue
fi

# Run htseq-count and capture standard output and error
echo "Processing $bam_file" | tee -a "$log_file"
{ time htseq-count -f bam -s no -t gene -i ID "$bam_file" "$gff_file" > "$output_file"; } 2>&1 | tee -a "$log_file"

# Calculate and log mapping statistics using samtools
total_reads=$(samtools view -c "$bam_file")
mapped_reads=$(samtools view -c -F 4 "$bam_file")
mapping_rate=$(echo "scale=2; ($mapped_reads/$total_reads)*100" | bc)

echo "Total reads: $total_reads" | tee -a "$log_file"
echo "Mapped reads: $mapped_reads" | tee -a "$log_file"
echo "Mapping rate: $mapping_rate%" | tee -a "$log_file"

# Append mapping statistics to the summary file
echo "$sample,$total_reads,$mapped_reads,$mapping_rate" >> "$summary_file"
echo "----------------------------------" | tee -a "$log_file"
done

echo "HTSeq-count run completed."

