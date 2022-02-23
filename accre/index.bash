#!/bin/bash
#SBATCH --mail-user=deborah.j.park@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=90

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon index -t Mus_musculus.GRCm38.cdna.all.fa -i mt_index_b -k 31
