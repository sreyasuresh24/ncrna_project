#!/bin/bash -l

###### Run this script on MARS to map single cell data to the S. mansoni genome (v10) using Cellranger
###### Run the script like this (changing resources as required, depending on data):
# sbatch --cpus-per-task=1 --mem=34G /mnt/data/project0014/Sreya/ncrna

############# SLURM SETTINGS #############
#SBATCH --account=project0014   # Account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=cr_mapping   # Descriptive job name of your choice
#SBATCH --output=%x-%j.out      # Output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # Error file name will contain job name + job ID
#SBATCH --partition=nodes       # Which partition to use, default on MARS is “nodes"
#SBATCH --time=0-20:00:00       # Time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=34G               # Memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # Number of nodes to allocate, default is 1
#SBATCH --cpus-per-task=12      # Number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --mail-user=2931382S@student.gla.ac.uk # Email address for notifications
#SBATCH --mail-type=END         # Mail me when my jobs ends
#SBATCH --mail-type=FAIL        # Mail me if my jobs fails

############# LOADING MODULES #############
module load apps/cellranger/7.2.0 

############# CODE ############

# Define reference transcriptome
REFERENCE="/mnt/data/project0014/Reference_data/Smansoni_genome"

# Define output directory
OUTPUT_DIR="/mnt/data/project0014/Sreya/ncrna/cellranger/cellranger_count_Sm"

# Define dir containing fastq files
FASTQ_DIR="/mnt/data/project0014/Sreya/ncrna/output_results/datasets/fastq_data"

#Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

#Loop through each sample directory in the FASTQ directory
for sample_dir in "$FASTQ_DIR"/*; do
    # Ensure the item is a directory
    if [ -d "$sample_dir" ]; then
        sample_id=$(basename "$sample_dir")
        echo "Processing sample: $sample_id"
        echo "$sample_dir"
        #Run cellranger count S. mansoni ref
		cellranger count --id="$sample_id" --transcriptome="$REFERENCE" --sample="$sample_id" --localcores 12 --localmem 34 --include-introns false --fastqs="$sample_dir";

#Check if cellranger count was successful
        if [ $? -eq 0 ]; then
            echo "Sample $sample_id processed successfully."
            # Move the output to the desired directory
            mv "$sample_id" "$OUTPUT_DIR/$sample_id"
        else
            echo "Error processing sample $sample_id" >&2
      fi      
	fi
done


