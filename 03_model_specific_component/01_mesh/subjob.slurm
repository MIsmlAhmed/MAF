#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --account=rpp-kshook
#SBATCH --time=3:00:00
#SBATCH --mem=90000M
#SBATCH --job-name=meshflow


##################
#1. load the needed modules
    
module restore scimods

#2. load python environments (that contains all needed packages)
    
source $HOME/virtual-envs/scienv/bin/activate

# run python
python write_mesh.py