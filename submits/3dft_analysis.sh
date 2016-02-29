#!/bin/sh

#SBATCH --job-name=3dFT_analysis
#SBATCH --partition=pre
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

#SBATCH --time=1-00:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

echo "Date:"
date '+%s'
echo "Github has:"
git rev-parse --verify HEAD
echo "Using ACI / HCP / Slurm cluster."
echo "JobID = $SLURM_JOB_ID"
echo "Using $SLURM_NNODES nodes"
echo "Using $SLURM_NODELIST nodes."
echo "Number of cores per node: $SLURM_TASKS_PER_NODE"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo $@
echo ""

# Executable
python /home/maldonis/3dft/3dft_analysis.py $@

echo "Finished on:"
date '+%s'
