#!/bin/bash
#SBATCH --mail-user=deborah.j.park@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=300


/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-1_2.fastq 
/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s1/NS-1_1_val_1.fq.gz 
/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s1/NS-1_2_val_2.fq.gz
