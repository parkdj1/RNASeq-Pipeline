#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 15:12:03 2020

@author: deborah
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 23:57:15 2020

@author: deborah
"""


file_name_before_num = "../Serezani/salmon_quant/"
file_name_after_num = "quant"
q_file = ".sf"
o = "_collapsed.sf"
num_samps = 18

id_file = "../Mus_musculus.GRCm38.cdna.all.fa"
keys = ['chromosome', 'gene', 'gene_biotype', 'transcript_biotype', 'gene_symbol', 'description', 'Source']
num_keys = 7


gene_dict = {}
gene_dict[">Name"] = keys

with open(id_file) as f:
    for line in f:
        gene_ids = []
        split_line = line.split(" ", num_keys)
        if (len(split_line) >1):
            temp = split_line[-1].split(" [Source")
            temp[-1] = "Source" + temp[-1].strip(']\n')
            split_line.pop()
            split_line.extend(temp)
            split_line[0] = split_line[0][1:]
            for i in split_line[2:]:
                label, value = i.split(":", 1)
                gene_ids.append(value)
            gene_dict[split_line[0]] = gene_ids


new_row = []

for num in range(1,num_samps+1):
    quant_file = file_name_before_num + str(num) + file_name_after_num + q_file
    ofile = file_name_before_num + str(num) + file_name_after_num + o

    genes = {}

    with open(quant_file) as file:
        for line in file:
            split_row = line.split("\t")
            split_row[-1] = split_row[-1].strip('\n')
            if (split_row[-1] != "NumReads"):
                split_row[-1] = int(float(split_row[-1])) + 1
            if split_row[0] != "Name":
                key = gene_dict[split_row[0]][4]
                if key in genes:
                    genes[key] += split_row[-1]
                else:
                    genes[key] = split_row[-1]
    with open(ofile, mode='w') as f:
        f.write("Name\tNumReads\n")
        for key in genes:
            f.write(str(key))
            f.write('\t')
            f.write(str(genes[key]))
            f.write('\n')