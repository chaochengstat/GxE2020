#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH -t 5-10:00:00
#SBATCH -J Ex4
#SBATCH -o log.%A_%a
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --array=1-8 

echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'

module load R

srun Rscript linear_Pskew.R