#!/bin/sh
#SBATCH --partition=gpu
#SBATCH --gres=gpu
#SBATCH --mem=5gb
#SBATCH --time=00:20:00
#SBATCH --job-name=python_signalp
#SBATCH --error=SLURM_report/error_report/python_signalp.%J.err
#SBATCH --output=SLURM_report/output/python_signalp.%J.out

module load python/3.8
module load signalp/4.1
python smallprotein_signalp.py -p HCC_test
echo "I am task '%J' on node `hostname`"


