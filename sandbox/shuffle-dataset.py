#!/usr/bin/python

import sys
import csv
import random

fh_gl = open("/home/awright/git/PathwayAnalysis/gecco/GECCO.data.release.2.2.for.MP-BioPath.txt", "r")

for _ in range(10):
    next(fh_gl)

r = csv.reader(fh_gl, delimiter="\t")

participant_to_gene_map = {}
all_genes = set()
for row in r:
    participant = row[0]
    genes = row[3].split(",")
    gene = genes[0]
    if participant not in participant_to_gene_map.keys():
        participant_to_gene_map[participant] = set()
    participant_to_gene_map[participant].add(gene)
    all_genes.add(gene)

all_genes_list_orig = list(all_genes)
all_genes_list = list(all_genes)
random.shuffle(all_genes_list)
genes_map = dict(zip(all_genes_list_orig, all_genes_list))

for participant in participant_to_gene_map.keys():
    for gene in participant_to_gene_map[participant]:
        print participant + "\t" + genes_map[gene]
