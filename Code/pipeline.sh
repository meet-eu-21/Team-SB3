#!/usr/bin/bash

#SBATCH --mem=200GB
#SBATCH --cpus-per-task=4
module load python
python multi_pipeline.py

