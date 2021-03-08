#!/bin/sh
#SBATCH --partition=gpu
#SBATCH --gres=gpu  
#SBATCH --mem=5gb
#SBATCH --time=00:20:00
#SBATCH --job-name=python_tmhmm
#SBATCH --error=SLURM_report/error_report/python_tmhmm.%J.err
#SBATCH --output=SLURM_report/output/python_tmhmm.%J.out

module load python/3.8
pip install Bio
module load emboss/6.6
module load tmhmm/2.0c
python HCC_smallprotein_tmhmm.py --nargs "MGYG-HGUT-001"
echo "I am task '%J' on node `hostname`"


