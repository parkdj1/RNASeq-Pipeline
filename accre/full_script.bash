#!/bin/bash
#SBATCH --mail-user=deborah.j.park@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1000

filePath=/data/pua_lab/Serezani_2020/
start=5072-AB-
ext=.fq
startNum=10
numSamp=18

fqcPath=/home/parkdj1/software/
TGPath=/home/parkdj1/software/
salmonPath=/home/parkdj1/software/
indexPath=/home/parkdj1/PUA_SPRING_2020/

ml Intel/2017.4.196 IntelMPI Python/3.6.3
source /home/parkdj1/software/trimgalore/bin/activate

for num in $(seq $startNum $numSamp)
do
mkdir ${filePath}$num

file1=${filePath}${start}${num}_1
file2=${filePath}${start}${num}_2

# gunzip -k ${file1}${ext}.gz
# gunzip -k ${file2}${ext}.gz

# ${fqcPath}FastQC/fastqc ${file1}${ext} ${file2}${ext} --outdir=${filePath}$num

#${TGPath}TrimGalore/trim_galore -j 0 -o $filePath --paired --gzip ${file1}${ext} ${file2}${ext}

#${fqcPath}FastQC/fastqc ${file1}_val_1.fq.gz ${file2}_val_2.fq.gz --outdir=${filePath}$num

${salmonPath}salmon-latest_linux_x86_64/bin/salmon quant -i ${indexPath}mt_index  -l A \
        -1 ${file1}.fq.gz\
        -2 ${file2}.fq.gz\
        -p 8 --validateMappings -o ${filePath}$num
done
