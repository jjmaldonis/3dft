#!/bin/bash

#$ -N IFT_3D

#$ -q all.q
#$ -pe orte 1

#$ -l h_rt=INFINITY
#$ -cwd
#$ -o $JOB_NAME.$JOB_ID.out
#$ -e $JOB_NAME.$JOB_ID.err
#$ -V

# Write the date and time to outfile.
echo "Submitted on:"
date '+%s'

MPI_HOME=/share/apps/openmpi_intel_20130712/bin
###$MPI_HOME/mpiexec -n $NSLOTS inverse3dft $@
$MPI_HOME/mpiexec -n $NSLOTS 3dft $@

echo "Finished on:"
date '+%s'
