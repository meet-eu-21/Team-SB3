#!/usr/bin/bash
#SBATCH --mem=200GB
#SBATCH --ntasks=10
#SBATCH --partition=long
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=adrienpauron@gmail.com
module load hmmlearn
module load python
python multi_pipeline.py

