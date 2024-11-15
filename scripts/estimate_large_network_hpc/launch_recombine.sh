#!/bin/bash
#SBATCH --job-name=recombine
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --partition=shared-cpu,shared-bigmem,public-cpu,public-bigmem
#SBATCH --mail-user=lionel.voirol@unige.ch
#SBATCH --mail-type=NONE
#SBATCH --output estimate_station_ngl/outfile/outfile_recombine.out
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.0.0
INFILE=estimate_station_ngl/recombine.R
OUTFILE=estimate_station_ngl/report/recombine.Rout
srun R CMD BATCH $INFILE $OUTFILE
