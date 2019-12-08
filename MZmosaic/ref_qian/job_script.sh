#!/bin/bash
# Created by the Melbourne Bioinformatics job script generator for SLURM
# Tue Jun 12 2018 09:53:46 GMT+1000 (AEST)

# Partition for the job:
#SBATCH -p sysgen

# The name of the job:
#SBATCH --job-name="qian_ape"

# The project ID which this job should run under:
#SBATCH --account="SG0011"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=120000

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=fengq2@student.unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=1-0:0:00
#SBATCH -e "slurm-%A_%a.out"

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The job command(s):
/vlsci/SG0011/qian-feng/MZmosaic/mosaic -ma -seq combined_ape_454_renamedForMosaic.fasta -aa -tag combined_ape_454_renamedForMosaic20180612 -psum -group 2 db target -target target -del 0.0166765773591 -eps 0.273305573191 -rec 0.014 
