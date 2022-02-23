#!/bin/bash
#SBATCH --mail-user=deborah.j.park@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=300



/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s2/NS-2_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s2/NS-2_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/2_quant

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s3/NS-3_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s3/NS-3_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/3_quant

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s4/NS-4_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s4/NS-4_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/4_quant

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s5/NS-5_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s5/NS-5_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/5_quant

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s6/NS-6_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s6/NS-6_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/6_quant

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s8/NS-8_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s8/NS-8_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/8_quant

/home/parkdj1/software/salmon-latest_linux_x86_64/bin/salmon quant -i /home/parkdj1/PUA_SPRING_2020/mt_index  -l A \
        -1 /scratch/pua_lab/s8/NS-8_1_val_1.fq.gz\
        -2 /scratch/pua_lab/s8/NS-8_2_val_2.fq.gz\
        -p 8 --validateMappings -o /scratch/pua_lab/quants/8_quant
