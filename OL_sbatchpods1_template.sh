#!/bin/bash
#SBATCH --account=s1189 
#SBATCH --constraint=hasw
##SBATCH --qos=long 
#SBATCH --time=11:59:59 
#SBATCH --output="slurm_output_%x_%j.out"
#SBATCH --error="slurm_error_%x_%j.out"
#SBATCH --ntasks=nttasks --ntasks-per-node=nTasksPerNode
#SBATCH --job-name=OL

. /usr/share/modules/init/bash
#module purge
module load comp/intel/20.0.0.166
module load mpi/impi/20.0.0.166
ulimit -s unlimited
cd /discover/nobackup/syatheen/Grey/THP/noah-da-gp
/usr/local/other/PoDS/PoDS/pods.py -x /discover/nobackup/syatheen/Grey/THP/noah-da-gp/Execfile_OL -n nTasksPerNode
exit 0

 
