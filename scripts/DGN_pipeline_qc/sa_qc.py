# 2024-11-19 JCF

# Purpose: Want a nicely formatted QC output from running the DGN pipeline.
# Start with library duplication stats + FAS1K stats for now

# Ideal stats would be:
# Duplication
# mapping rate
# Coverage
# Number SNPS passing
# Number indels passing 
# Number of Ns in fas1k
# Number of het sites in fas1k

###############################################################

import glob
import pandas as pd
from dup_qc import DuplicationMetrics
from parse_flagstat import Flagstat

###############################################################

parent_dir = "/home/jamie/DGN_compatible/ZI_N"

###############################################################

def get_dup_file(sample):
	dup_prefix   = parent_dir + '/round2/logs/stampy/'
	dup_suffix   = '_dups.txt'
	return dup_prefix + sample + dup_suffix

def get_stats_file(sample):
	stats_prefix   = '/home/jamie/DGN_compatible/ZI_N/round2/logs/stampy/'
	stats_suffix = '.stats'
	pattern = stats_prefix + sample + '*' + stats_suffix
	stats_files = sorted(glob.glob(pattern))
	return stats_files

def get_fas1k_het_file(sample):
	f1k_het_prefix = '/home/jamie/DGN_compatible/ZI_N/round2/het_stats/'
	f1k_het_suffix = '_round2_het.txt'
	return f1k_het_prefix + sample + f1k_het_suffix

###############################################################

class SampleQC:
	def __init__(metrics, sample):
		metrics.sample                 = sample
		metrics.dup_file               = get_dup_file(metrics.sample)
		metrics.dup_obj                = DuplicationMetrics(metrics.dup_file)
		metrics.stats_files            = get_stats_file(metrics.sample)
		flags = [Flagstat(path) for path in metrics.stats_files]
		metrics.flagstat_obj = flags[0]
		for fs in flags[1:]:
			metrics.flagstat_obj = metrics.flagstat_obj + fs
		metrics.stats_file = metrics.stats_files[0] if len(metrics.stats_files) == 1 else metrics.stats_files
		metrics.percent_duplication    = metrics.dup_obj.percent_duplication
		metrics.estimated_library_size = metrics.dup_obj.estimated_library_size
		metrics.f1k_het_file           = get_fas1k_het_file(metrics.sample)
		metrics.f1k_het_contents       = pd.read_table(metrics.f1k_het_file)
	
	def __stats_out__(metrics):
		d = {"Library": metrics.sample,
	   		"Raw_read_pairs": metrics.flagstat_obj.read1, 
		     "fas1k_Ns": metrics.f1k_het_contents.iloc[5,]["prop_N"],
		     "fas1k_het": metrics.f1k_het_contents.iloc[5,]["prop_het"],
		     "Duplication_rate": metrics.dup_obj.percent_duplication,
		     "Est_library_size_millions" : metrics.dup_obj.estimated_library_size/1000000
		     }
		return d
	
	def f1k_stats(metrics, stat):
		metrics.stat = metrics.f1k_het_contents.iloc[5,][stat]
		return metrics.stat


###############################################################

if __name__ == "__main__":
        import sys
        library_in = sys.argv[1]
        outfile    = sys.argv[2]
        p1         = SampleQC( library_in )
        df         = pd.DataFrame.from_dict( p1.__stats_out__(),orient='index' ).T
        df.to_csv(outfile, sep='\t', index=False)



