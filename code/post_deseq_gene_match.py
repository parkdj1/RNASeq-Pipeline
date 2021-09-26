#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 23:57:15 2020

@author: deborah
"""
import csv

extension = "../Serezani/"
file_name_1 = "CTL_v_STZ"
file_name_2 = "norm_counts"
q_file = ".csv"
o = "_mod.csv"

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
            for i in split_line[2:]:
                label, value = i.split(":", 1)
                gene_ids.append(value)
            gene_dict[split_line[0]] = gene_ids


new_row = []

quant_file = extension + file_name_1 + q_file
ofile = extension + file_name_1 + o

with open(ofile, mode='w') as f:
    output_writer = csv.writer(f, delimiter=',', quotechar='"')
    with open(quant_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            row[-1] = row[-1].strip('\n')
            new_row.extend(row)
            if (row[0] != ""):
                new_row.extend(gene_dict[">" + row[0]])
            else:
                new_row.extend(gene_dict[">Name"])
            output_writer.writerow(new_row)
            new_row.clear()

new_row = []

quant_file = extension + file_name_2 + q_file
ofile = extension + file_name_2 + o

with open(ofile, mode='w') as f:
    output_writer = csv.writer(f, delimiter=',', quotechar='"')
    with open(quant_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            row[-1] = row[-1].strip('\n')
            new_row.extend(row)
            if (row[0] != ""):
                new_row.extend(gene_dict[">" + row[0]])
            else:
                new_row.extend(gene_dict[">Name"])
            output_writer.writerow(new_row)
            new_row.clear()