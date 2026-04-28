#!/usr/bin/env python3

#
# 2026-04 JCF
#######################################################################
#
# Purpose:
#
#  Filters a sample table TSV (produced by sample_table_from_dir.py)
#  to include only sample-unit combinations that exist in a specified directory.
#  (Normally using because I'm running part of a set in batches).    
#
# Input:
# 
#  Positional: 
#    1. Sample table TSV file path (from sample_table_from_dir.py)
#    2. Directory path to search for matching files
#    3. Output TSV file name
#
#  Optional:
#    --pattern: File pattern to match (default: '*.bam')
#
######################################################################

import os
import sys
import glob
import re
import pandas as pd
import argparse

######################################################################

# Function definitions

def extract_sample_unit(row):
	"""
	Extract concatenated sample-unit identifier from table row.
	Format: sample_unit (e.g., '27Mar26-1-ZI-1-3-F_S604')
	"""
	return row['sample'] + '_' + row['unit']

def get_sample_units_from_files(directory, pattern):
	"""
	Extract sample-unit identifiers from files in directory.
	Assumes filenames contain the pattern: ...sample-unitinfo...
	"""
	files = glob.glob(os.path.join(directory, pattern))
	files.sort()
	
	sample_units = set()
	
	for f in files:
		# Get filename without path
		filename = os.path.basename(f)
		
		# Try to extract sample-unit from filename
		# Assumes format like: sample_unit_S###.bam or similar
		# Extract everything up to and including the S number (e.g., 27Mar26-1-ZI-1-3-F_S604)
		match = re.match(r'^(.+_S\d+)', filename)
		if match:
			sample_unit = match.group(1)
			sample_units.add(sample_unit)
	
	return sample_units

def filter_sample_table(table_path, directory, pattern, output_path):
	"""
	Filter sample table to include only entries matching files in directory.
	"""
	# Read the sample table
	df = pd.read_csv(table_path, sep="\t")
	
	# Get sample-unit identifiers from files in directory
	existing_units = get_sample_units_from_files(directory, pattern)
	
	# Create sample-unit column for comparison
	df['sample_unit'] = df.apply(extract_sample_unit, axis=1)
	
	# Filter to only rows where sample-unit exists in directory
	df_filtered = df[df['sample_unit'].isin(existing_units)].copy()
	
	# Remove the temporary sample_unit column
	df_filtered = df_filtered.drop(columns=['sample_unit'])
	
	# Write filtered table
	df_filtered.to_csv(output_path, sep="\t", index=False)
	
	# Print summary
	print(f"Original table: {len(df)} rows")
	print(f"Files in directory: {len(existing_units)} unique sample-units")
	print(f"Filtered table: {len(df_filtered)} rows")
	print(f"Output written to: {output_path}")
	
	return df_filtered

######################################################################

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description="Filter sample table to include only entries matching files in a directory."
	)
	parser.add_argument("table", help="Input sample table TSV file")
	parser.add_argument("directory", help="Directory to search for matching files")
	parser.add_argument("output", help="Output filtered TSV file")
	parser.add_argument(
		"--pattern",
		default="*.bam",
		help="File pattern to match (default: '*.bam')"
	)
	
	args = parser.parse_args()
	
	filter_sample_table(args.table, args.directory, args.pattern, args.output)
