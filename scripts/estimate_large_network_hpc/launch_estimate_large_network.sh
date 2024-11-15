#!/bin/bash
#SBATCH --partition=shared-cpu,shared-bigmem,public-cpu,public-bigmem
#SBATCH --time=00-00:20:00
#SBATCH --mem-per-cpu=25000 # in MB
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mail-user=lionel.voirol@unige.ch
#SBATCH --job-name=estimate_large_station
#SBATCH --mail-type=NONE
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
module load GCC/9.3.0 OpenMPI/4.0.3 R/4.0.0

INFILE=estimate_station_ngl/estimate_large_network.R
OUTFILE=estimate_station_ngl/report/report_${n}_${SLURM_ARRAY_TASK_ID}.Rout
OUTLOG=estimate_station_ngl/outfile/outfile_${n}_${SLURM_ARRAY_TASK_ID}.out

exec > $OUTLOG 2>&1


srun R CMD BATCH $INFILE $OUTFILE

