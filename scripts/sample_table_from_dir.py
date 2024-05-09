
#!/usr/bin/env python3

#
# 2024-04 JCF
#######################################################################
#
# Purpose:
#
#  Prints sample table .tsv for snakemake input from a directory of fq files.
#   
#   
#
# Input:
# 
# Postional: 1. directory path of fq files 2. outfile name
#   
#
######################################################################

import os
import sys
import glob
import re
import pandas as pd

######################################################################

# Function definitions

def sample_from_fq(s):
	# Split full dir path off of front and then
	#	Our core returns files that end in _S[0-9][0-9][0-9]_R[1-2]_001.fastq.gz
	# 	split the string on _S[0-9][0-9][0-9]_
	s1 = s.split("/")[len(s.split("/"))-1]
	s2 = re.split('_S(?:[0-9]{2,3})_', s1)[0]
	return s2

def get_unit(s):
	# Use the S[0-9][0-9][0-9] string as sample unit ID- unique per sample per lane
	#	(to be truly safe that all units unique could append lane id here too)
	s1 = re.search('_S(?:[0-9]{2,3})_', s).group()
	s2 = s1.split("_")[len(s1.split("_"))-2]
	return s2

def write_sample_table(d):
	all_files = glob.glob(os.path.join(d, "*R1*fastq.gz"))
	all_files.sort()

def get_RG(s):
    # Format RG string from sample info and use function repr to print the tab
    #   escape character
	s1 = "@RG\tID:" + sample_from_fq(s) + "-" + get_unit(s) + "\tSM:" + sample_from_fq(s)
	s2 = "\tPU:" + "none" + "\tLB:" + sample_from_fq(s) + "\tPL:" + "ILLUMINA"
	return repr(s1 + s2)
	
def sample_table(d):
	# Get all R1 fastq files
	all_files = glob.glob(os.path.join(d, "*R1*fastq.gz"))
	all_files.sort()
	
	# Get info for sample table
	sa  = [ sample_from_fq(x) for x in all_files]
	un  = [ get_unit(x) for x in all_files]
	pl  = [ "ILLUMINA" for x in all_files ] 
	fq1 = all_files
	fq2 = re.sub(r'(?is)_R1_', '_R2_', fq1[0])
	RG  = [ get_RG(x) for x in all_files ]
		
	# Shape into dataframe
	df  = pd.DataFrame({'sample': sa, 'unit': un, 'platform': pl, 'fq1': fq1, 'fq2': fq2, 'RG': RG})
	
	return df

######################################################################

if __name__ == "__main__":
	
	where     = sys.argv[1]
	out_name  = sys.argv[2]
	
	out = sample_table(where)
	out.to_csv(out_name, sep="\t", index = False)
	
	
	
	
	
