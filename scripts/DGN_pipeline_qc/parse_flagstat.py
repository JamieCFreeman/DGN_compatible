# 2024-11-11 JCF

# Purpose: 
# To pull duplication metrics from a metrics file written by 
#  GATK MarkDuplicates

# Input:
# metrics file as written by GATK markduplicates

###############################################################
# # Example flasgstat output file
# 31194858 + 0 in total (QC-passed reads + QC-failed reads)
# 31194858 + 0 primary
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 0 + 0 primary duplicates
# 30879493 + 0 mapped (98.99% : N/A)
# 30879493 + 0 primary mapped (98.99% : N/A)
# 31194858 + 0 paired in sequencing
# 15597429 + 0 read1
# 15597429 + 0 read2
# 25196938 + 0 properly paired (80.77% : N/A)
# 30780142 + 0 with itself and mate mapped
# 99351 + 0 singletons (0.32% : N/A)
# 1817110 + 0 with mate mapped to a different chr
# 259848 + 0 with mate mapped to a different chr (mapQ>=5)

###############################################################
# Helper functions
def lines_in(f):
	f        = open(f,'r')
	lines    = f.readlines()
	f.close()
	return lines

def transform_lines(lines):
	l = [ x.split(' ')[0] for x in lines ] 
	return l

def lib_name(s, ext):
	s1 = s.split('/')[ len(s.split('/'))-1 ]
	return s1.split(ext)[0]

###############################################################
# Class definition
class Flagstat:
  def __init__(metrics, filename):
    metrics.filename              = filename
    metrics.raw                   = lines_in(filename)
    metrics.list                  = transform_lines(metrics.raw)
    metrics.library               = lib_name(filename, '.stats')
    metrics.total_reads           = int( metrics.list[0] )
    metrics.primary_align         = int( metrics.list[1] )
    metrics.secondary_align       = int( metrics.list[2] )
    metrics.suppl_align           = int( metrics.list[3] )
    metrics.duplicates            = int( metrics.list[4] )
    metrics.primary_duplicates    = int( metrics.list[5] )
    metrics.mapped                = int( metrics.list[6] )
    metrics.primary_mapped        = int( metrics.list[7] )
    metrics.paired_reads          = int( metrics.list[8] )
    metrics.read1                 = int( metrics.list[9] )
    metrics.read2                 = int( metrics.list[10] )
    metrics.properly_paired       = int( metrics.list[11] )
    metrics.both_reads_mapped     = int( metrics.list[12] )
    metrics.singletons            = int( metrics.list[13] )
    metrics.mate_map_diff_chr     = int( metrics.list[14] )
    metrics.mate_map_diff_chr_MQ  = int( metrics.list[15] )

###############################################################
if __name__ == "__main__":
	import sys
	stats_file = sys.argv[1]
	p1 = Flagstat(stats_file)
	s1 ='total_reads' + '\t' + 'perc_mapped_reads' + '\t' + 'properly_paired_reads'
	print( s1 )
	s2 = str(p1.total_reads) + '\t' \
	+ str(round(p1.mapped / p1.total_reads, 4)) + '\t' \
	+ str(round(p1.properly_paired / p1.total_reads, 4))
	print( s2)

