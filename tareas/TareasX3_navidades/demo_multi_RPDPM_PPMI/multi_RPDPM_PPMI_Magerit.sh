#!/bin/bash
##----------------------- Start job description ----------------------- 
#SBATCH --partition=standard
#SBATCH --job-name=multi_RPDPM_PPMI
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=4096
#SBATCH --time=160:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=carlos.platero@upm.es
#SBATCH --output=out-%j.log
#SBATCH --error=err-%j.log
#SBATCH --chdir=/home/r503/r503739/PPMI/Analysis/multi_RPDPM_PD/
##------------------------ End job description ------------------------ 


# Carga de los módulos necesarios (MATLAB y R)
module load MATLAB/2022b



# Definición del directorio SCRATCH
export SCRATCH="/home/r503/r503739/.matlab/local_cluster_jobs/R2022b"
mkdir -p $SCRATCH/$SLURM_JOB_ID


# Lanzamiento del script
matlab -nodisplay  < /home/r503/r503739/PPMI/Analysis/multi_RPDPM_PD/multi_RPDPM_PPMI_parfor.m

# Limpieza del directorio SCRATCH
rm -rf $SCRATCH/$SLURM_JOB_ID
