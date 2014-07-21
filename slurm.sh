#!/bin/sh

#SBATCH --job-name=3dIFT
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
echo "Using ACI / HCP / Slurm cluster."
echo "JobID = $SLURM_JOB_ID"
echo "Using $SLURM_NNODES nodes"
echo "Using $SLURM_NODELIST nodes."
echo "Number of cores per node: $SLURM_TASKS_PER_NODE"
echo "Submit directory: $SLURM_SUBMIT_DIR"
echo ""

# Executable
mpiexec inverse3dft $SLURM_JOB_ID paramfile_t1

echo "Finished on:"
date '+%s'
