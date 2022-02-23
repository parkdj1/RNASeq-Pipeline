#!/bin/bash
#SBATCH --mail-user=deborah.j.park@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=800



ml Intel/2017.4.196 IntelMPI Python/3.6.3
source /home/parkdj1/software/trimgalore/bin/activate


gzip -d /scratch/pua_lab/NS-3_1.fastq.gz
gzip -d /scratch/pua_lab/NS-3_2.fastq.gz

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-3_2.fastq /scratch/pua_lab/NS-3_2.fastq --outdir=/scratch/pua_lab/s3/fastqc_3

/home/parkdj1/software/TrimGalore/trim_galore --paired  --fastqc --gzip -o /scratch/pua_lab/s3 /scratch/pua_lab/NS-3_1.fastq /scratch/pua_lab/NS-3_2.fastq

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s3/NS-3_1_val_1.fq.gz /scratch/pua_lab/s3/NS-3_2_val_2.fq.gz --outdir=/scratch/pua_lab/s3/fastqc_3

rm NS-3_1.fastq
rm NS-3_2.fastq

gzip -d /scratch/pua_lab/NS-4_1.fastq.gz
gzip -d /scratch/pua_lab/NS-4_2.fastq.gz

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-4_2.fastq /scratch/pua_lab/NS-4_2.fastq --outdir=/scratch/pua_lab/s4/fastqc_4

/home/parkdj1/software/TrimGalore/trim_galore --paired  --fastqc --gzip -o /scratch/pua_lab/s4 /scratch/pua_lab/NS-4_1.fastq /scratch/pua_lab/NS-4_2.fastq

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s4/NS-4_1_val_1.fq.gz /scratch/pua_lab/s4/NS-4_2_val_2.fq.gz --outdir=/scratch/pua_lab/s4/fastqc_4

rm NS-4_1.fastq
rm NS-4_2.fastq


gzip -d /scratch/pua_lab/NS-5_1.fastq.gz
gzip -d /scratch/pua_lab/NS-5_2.fastq.gz

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-5_2.fastq /scratch/pua_lab/NS-5_2.fastq --outdir=/scratch/pua_lab/s5/fastqc_5

/home/parkdj1/software/TrimGalore/trim_galore --paired  --fastqc --gzip -o /scratch/pua_lab/s5 /scratch/pua_lab/NS-5_1.fastq /scratch/pua_lab/NS-5_2.fastq

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s5/NS-5_1_val_1.fq.gz /scratch/pua_lab/s5/NS-5_2_val_2.fq.gz --outdir=/scratch/pua_lab/s5/fastqc_5


gzip -d /scratch/pua_lab/NS-6_1.fastq.gz
gzip -d /scratch/pua_lab/NS-6_2.fastq.gz

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-6_2.fastq /scratch/pua_lab/NS-6_2.fastq --outdir=/scratch/pua_lab/s6/fastqc_6

/home/parkdj1/software/TrimGalore/trim_galore --paired  --fastqc --gzip -o /scratch/pua_lab/s6 /scratch/pua_lab/NS-6_1.fastq /scratch/pua_lab/NS-6_2.fastq

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s6/NS-6_1_val_1.fq.gz /scratch/pua_lab/s6/NS-6_2_val_2.fq.gz --outdir=/scratch/pua_lab/s6/fastqc_6

gzip -d /scratch/pua_lab/NS-7_1.fastq.gz
gzip -d /scratch/pua_lab/NS-7_2.fastq.gz

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-7_2.fastq /scratch/pua_lab/NS-7_2.fastq --outdir=/scratch/pua_lab/s7/fastqc_7

/home/parkdj1/software/TrimGalore/trim_galore --paired  --fastqc --gzip -o /scratch/pua_lab/s7 /scratch/pua_lab/NS-7_1.fastq /scratch/pua_lab/NS-7_2.fastq

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s7/NS-7_1_val_1.fq.gz /scratch/pua_lab/s7/NS-7_2_val_2.fq.gz --outdir=/scratch/pua_lab/s7/fastqc_7

gzip -d /scratch/pua_lab/NS-8_1.fastq.gz
gzip -d /scratch/pua_lab/NS-8_2.fastq.gz

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/NS-8_2.fastq /scratch/pua_lab/NS-8_2.fastq --outdir=/scratch/pua_lab/s8/fastqc_8

/home/parkdj1/software/TrimGalore/trim_galore --paired  --fastqc --gzip -o /scratch/pua_lab/s8 /scratch/pua_lab/NS-8_1.fastq /scratch/pua_lab/NS-8_2.fastq

/home/parkdj1/software/FastQC/fastqc /scratch/pua_lab/s8/NS-8_1_val_1.fq.gz /scratch/pua_lab/s8/NS-8_2_val_2.fq.gz --outdir=/scratch/pua_lab/s8/fastqc_8
