#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 48:00:00
#SBATCH --mem=50G
#SBATCH -J combo_4
#SBATCH -o combo_4_%j.out
#SBATCH -e combo_4_%j.err



echo "START"

module load python/3.10.1


python3 Run_CombinationAnalysis_4.py


echo "END"




