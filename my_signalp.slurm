#!/bin/sh
#SBATCH --array=0-9
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=30
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --job-name=python_signalp
#SBATCH --error=SLURM_report/error_report/python_signalp.%J.err
#SBATCH --output=SLURM_report/output/python_signalp.%J.out


module load python/3.8
module load signalp/5.0b

top_folder="MGYG-HGUT-00$SLURM_ARRAY_TASK_ID"
for sec_level_folder in /work/juancui/yinchaohe/geno_fastas/$top_folder/*; do

python HCC_smallprotein_signalp.py $sec_level_folder
done
wait



