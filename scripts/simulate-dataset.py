#!/usr/bin/python

import sys
import csv
import random

fh_gl = open("/home/awright/git/PathwayAnalysis/gecco/GECCO.data.release.2.2.for.MP-BioPath.txt", "r")
simulated_number_of_patients = int(sys.argv[1])

for _ in range(10):
    next(fh_gl)

r = csv.reader(fh_gl, delimiter="\t")

count_participants = 0
participant_to_gene_map = {}
gene_count_map = {}
number_of_genes_associated_to_patients = 0
for row in r:
    participant = row[0]
    genes = row[3].split(",")
    gene = genes[0]
    if participant not in participant_to_gene_map.keys():
        participant_to_gene_map[participant] = set()
    participant_to_gene_map[participant].add(gene)
    if gene not in gene_count_map.keys():
        gene_count_map[gene] = 0
    gene_count_map[gene] += 1
    number_of_genes_associated_to_patients += 1

count_participants = number_of_genes_associated_to_patients / len(gene_count_map)
for i in range(simulated_number_of_patients):
    i += 1
    random_number_for_number_of_mutations = random.uniform(0,1)
    simulated_number_of_mutated_genes = int()
    total = 0
    for participant in participant_to_gene_map.keys():
        participant_count = len(participant_to_gene_map[participant])
        window = participant_count / float(number_of_genes_associated_to_patients)
        total += window
        if total > random_number_for_number_of_mutations:
            simulated_number_of_mutated_genes = participant_count
            break
    
    fh_gl.close()
    
    genes = gene_count_map.keys()
    simulated_patient_genes = set()
    remaining_number_of_genes_associated_to_patients = number_of_genes_associated_to_patients
    for j in range(simulated_number_of_mutated_genes):
        random_number = random.uniform(0,1)
        total = 0
        for gene in genes:
            total += gene_count_map[gene] / float(remaining_number_of_genes_associated_to_patients)
            if total > random_number:
                print "patient" + str(i) + "\t" + gene
                simulated_patient_genes.add(gene)
                remaining_number_of_genes_associated_to_patients -= gene_count_map[gene]
                genes.remove(gene)
                break


