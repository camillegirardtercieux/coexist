#!/bin/bash
#SBATCH --job-name=Chap_2_22_06_17
#SBATCH --account=agap
#SBATCH --partition=agap_long
#SBATCH --output=/lustre/girardtercieuxc/Chap_2_22_06_17_log.%A_%a.txt
#SBATCH --error=/lustre/girardtercieuxc/Chap_2_22_06_17_err.%A_%a.txt
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-800%25

# Loading singularity module
module load R/packages/4.1.0 
module load cv-standard
module load gcc/8.5.0
module load gdal/3.4.2
module load geos/3.7.2

# echo $SLURM_ARRAY_TASK_ID

# Paths
script="/home/girardtercieuxc/Chap_2/_R/make.R"

# run python script inside container
Rscript $script
# Then use Sys.getenv('SLURM_ARRAY_TASK_ID') in script.R

# End
