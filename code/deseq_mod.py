#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:21:28 2020

@author: deborah
"""


file_name_before_num = "../FOLDER_NAME/salmon_quant/"
file_name_after_num = "quant"
q_file = ".sf"
o = "_mod.sf"
num_samps = 18


id_file = "Mus_musculus.GRCm38.cdna.all.fa"
keys = ['chromosome', 'gene', 'gene_biotype', 'transcript_biotype', 'gene_symbol', 'description', 'scaffold']
num_keys = 7


new_row = []

for i in range(1,num_samps+1):
    quant_file = file_name_before_num + str(i) + file_name_after_num + q_file
    ofile = file_name_before_num + str(i) + file_name_after_num + o

    f = open(ofile,"w")

    with open(quant_file) as file:
        for line in file:
            split_row = line.split("\t")
            split_row[-1] = split_row[-1].strip('\n')
            if (split_row[-1] != "NumReads"):
                split_row[-1] = int(float(split_row[-1])) + 1
            f.write(str(split_row[0]))
            f.write('\t')
            f.write(str(split_row[-1]))
            f.write('\n')


