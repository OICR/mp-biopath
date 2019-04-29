#!/usr/bin/python

#################################################################################
# Created by: Cadia Chan (Lincoln lab rotation student 2017)
# Date: November 6, 2017
#################################################################################

# Creates heatmap, t-SNE plot and Kaplan Meier survival curve from
# MP-BioPath output file.

# By default, heatmap (-m) is always generated. User can specify to 
# not generate the heatmap visual.

# Both the t-SNE plot (-t) and Kaplan Meier (-k) survival curve will
# be generated if specified by the user.

# A file containing list of pathways and colours associated with them is
# required to run script.

# Ensure that "create_visuals.R" script is in the scripts folder

#################################################################################

# USAGE: (to make heatmap, t-SNE plot, and Kaplan Meier curve)
# python make_output_files.py mpbiopath_results.summary.tsv -p pathway_file.tsv -o outputfile_prefix -t -k

#################################################################################

# OUTPUT: .svg files of output visuals from R script

#################################################################################

import sys, argparse, subprocess

parser = argparse.ArgumentParser(description='Creates output visuals from MP-BioPath output data file')
parser.add_argument("data", help="output file from MP-BioPath")
parser.add_argument("-m", "--heatmap", help="default is 'T' to make heatmap, else specify 'F' for no heatmap",
 default = "T")
parser.add_argument("-t", "--tsne", help="specify to create t-SNE plot visual",
	action = "store_true")
parser.add_argument("-k", "--kaplan", help="specify to create Kaplan Meier curve visual",
	action = "store_true")
parser.add_argument("-p", "--pathway_list", help="list of pathways with colours specified")
parser.add_argument("-o", "--output", help="specify output file name", default="output")
parser.add_argument("-c", "--compound", help="specify compound name")

args = parser.parse_args()

if args.pathway_list is None:
	print("Need pathway list file for heatmap")
else:
	pathway_color_map = {}
	with open(args.pathway_list, 'r') as p:
		for line in p:
			split_line = line.strip().split("\t")
			pathway_color_map[split_line[2]] = split_line[3]
	#creates empty dictionary to contain donor and key output
	key_outputs = list()
	num_of_nodes = -1
	colour_column_string = ""

	with open(args.data, 'r') as d:
		for line in d:
			num_of_nodes += 1
			cols = line.strip().split("\t")
			if line.startswith("Nodes"):
				donors = cols[1:-1]
				continue
			cols[-1] = cols[-1].replace(" ", "_")
			if cols[-1] in pathway_color_map:
				colour_column_string += (pathway_color_map[cols[-1]] + ", ")
				key_outputs.append(cols[-1])

	colour_column_string = colour_column_string[:-2] #equivalent to colour_column_string

	output_string = ""
	colour_string = ""
	for output in set(key_outputs):
		output_string += (output + ", ")
		colour_string += pathway_color_map[output] + ", "
	output_string = output_string[:-2]
	colour_string = colour_string[:-2]

if args.heatmap is "F":
	print("heatmap script is not executed")

command = "Rscript" # --vanilla --slave < "
path2script = "sandbox/create_visuals.R"

#takes into account user's options
#user can specify to create tsne and/or kaplan meier survival curves
if args.tsne:
	if args.kaplan:
		switch = "htk"
	else:
		switch = "ht"
elif args.kaplan:
	switch = "hk"
else:
	switch = "h"
	
#calls and executes Rscript (create_visuals.R) by passing in arguments
cmd = [command, path2script, args.data, colour_column_string, args.output, output_string, colour_string, switch]#, args.compound]
print(cmd)
subprocess.call(cmd, universal_newlines = True)
