
# 2024-11-11 JCF

# Purpose: 
# To pull duplication metrics from a metrics file written by 
#  GATK MarkDuplicates

# Input:
# metrics file as written by GATK markduplicates


###############################################################
# Helper functions
def lines_in(f):
	f        = open(f,'r')
	lines    = f.readlines()
	f.close()
	return lines

def dup_list(lines):
	'''
	Take raw lines of duplication metrics file and parse to DuplicationMetrics object
	'''
	# Line 7 is the data we want
	i = lines[7].strip().split('\t')
	i[0] = '"' + i[0] + '"'
	return i

def lib_name(s):
	return s.split('/')[ len(s.split('/'))-1 ]

###############################################################
# Class definition
class DuplicationMetrics:
  def __init__(metrics, filename):
    metrics.filename                       = filename
    metrics.raw                            = lines_in(filename)
    metrics.list                           = dup_list(metrics.raw)
    metrics.library                        = lib_name(filename)
    metrics.unpaired_reads_examined        = int( metrics.list[1] )
    metrics.read_pairs_examined            = int( metrics.list[2] )
    metrics.secondary_or_supplementary_rds = int( metrics.list[3] )
    metrics.unmapped_reads                 = int( metrics.list[4] )
    metrics.unpaired_read_duplicates       = int( metrics.list[5] )
    metrics.read_pair_duplicates           = int( metrics.list[6] )
    metrics.read_pair_optical_duplicates   = int( metrics.list[7] )
    metrics.percent_duplication            = float( metrics.list[8] )
    metrics.estimated_library_size         = int( metrics.list[9] )

###############################################################

if __name__ == "__main__":
	import sys
	dup_file = sys.argv[1]
	p1       = DuplicationMetrics(dup_file)
	print( p1.library + "\t" + str(p1.percent_duplication) )



