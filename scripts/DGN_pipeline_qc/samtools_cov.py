# 2026-02-0 JCF

# Purpose: Want a nicely formatted QC output from running the DGN pipeline.
# Depth of coverage from samtools coverage

# This takes in a tsv file from samtools coverage, and returns a tab-sep line
# with input file name, global coverage, and coverage on main arms

###############################################################

import pandas as pd

###############################################################


def global_cov(d):
	'''
	Samtools coverage output is tsv with a meandepth column for each contig,
	so to get whole genome average coverage, take weighted avg of each contig
	with weights equal to their length
	'''
	total_bp     = d['endpos'].sum()
	coverage_out = ((d['endpos'] / total_bp) * d['meandepth']).sum()
	return coverage_out

def main_arm_cov(d):
	# Contigs 3:8 are the major chr arms
	subset_df  = d.iloc[2:8]
	subset_cov = global_cov(subset_df)
	return subset_cov

###############################################################

if __name__ == "__main__":
	import sys
	input_file = sys.argv[1]
	df = pd.read_csv(input_file, sep='\t')
	print(input_file + '\t' + str(global_cov(df)) + '\t' + str(main_arm_cov(df)) )

