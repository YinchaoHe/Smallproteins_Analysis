#!/bin/sh
#SBATCH --array=0-39
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=25
#SBATCH --mem=100gb
#SBATCH --time=168:00:00
#SBATCH --job-name=python_tmhmm
#SBATCH --error=SLURM_report/error_report/python_tmhmm.%J.err
#SBATCH --output=SLURM_report/output/python_tmhmm.%J.out

module load python/3.7
module load emboss/6.6
module load tmhmm/2.0c


top_folder="MGYG-HGUT-00$SLURM_ARRAY_TASK_ID"

for sec_level_folder in /work/juancui/yinchaohe/geno_fastas/$top_folder/*; do
    {
        python HCC_smallprotein_tmhmm_v2.py $sec_level_folder
    }&
done
wait
