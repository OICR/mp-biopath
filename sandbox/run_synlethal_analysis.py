#/usr/bin/python

import sys
from collections import defaultdict

sl_data = sys.argv[1]

syn_lethal = {}
obs_vs_exp = {}
with open(sl_data, 'r') as f:
	next(f)
	for line in f:
		split = line.strip().split("\t")
		if split[7] == "419195":
			node = "p-T19,S20-MRLC"
		if split[7] == "5668934":
			node = "p-T19-MRLC"
		if float(split[3]) > float(split[6]):
			bigname = split[1]
			bigstate = split[2]
			smallname = split[4]
			smallstate = split[5]
			big = float(split[3])
			small = float(split[6])
		else:
			bigname = split[4]
			bigstate = split[5]
			smallname = split[1]
			smallstate = split[2]
			big = float(split[6])
			small = float(split[3])
		diff = big - small
		expected = big * small
		key = node, bigname, bigstate, smallname, smallstate
		# diff = 1 -(( abs(float(split[2]) - float(split[4])))/big)
		# print diff
		syn_eff = float(split[8]) #no corrected value
		if diff < 0.10 and big <= 0.95:
			syn_diff = abs(float(split[8]) - small)
			# syn_eff = float(split[8]) + diff #if corrected value

			# print syn_diff
			if syn_diff > 0.05:
				syn_lethal[key] = split[8], big, small, diff, syn_eff
		obs_diff = float(split[8]) - expected
		obs_vs_exp[key] = split[8], big, small, diff, syn_eff, obs_diff

# print len(syn_lethal)

# for pair in syn_lethal:
# 	# print syn_lethal[pair]
# 	print (pair[1] + "\t" + pair[2] + "\t" + str(syn_lethal[pair][1]) + "\t" +
# 		pair[3] + "\t" + pair[4] + "\t" + str(syn_lethal[pair][2]) + "\t" +
# 		pair[0] + "\t" + str(syn_lethal[pair][0]) +
# 		"\t" + str(syn_lethal[pair][3]) + "\t" + str(syn_lethal[pair][4]))

for pair in obs_vs_exp:
	print (pair[1] + "\t" + pair[2] + "\t" + str(obs_vs_exp[pair][1]) + "\t" +
		pair[3] + "\t" + pair[4] + "\t" + str(obs_vs_exp[pair][2]) + "\t" +
		pair[0] + "\t" + str(obs_vs_exp[pair][0]) +
		"\t" + str(obs_vs_exp[pair][3]) + "\t" + str(obs_vs_exp[pair][4]) +
		"\t" + str(obs_vs_exp[pair][5]))