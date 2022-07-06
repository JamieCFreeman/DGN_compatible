# DGN_compatible
Pipeline to analyze DNA sequencing data with legacy compatibility for the Pool lab's Drosophila Genome Nexus data.


To run:
1. Create a sample table in the format of "units.tsv" (Tab-delimited table with columns sample, unit, platform, fq1, fq2, RG)
sample (sample name), unit (could be used to ID samples with multiple sequencing units- not currently implemented), 
platform ("ILLUMINA" or not), fq1 (path to fq1), fq2 (path to fq2), RG (read group string, eg @RG\tID:13Jan22-1_EF5N_S73\tSM:13Jan22-1_EF5N_S73\tPU:none\tLB:13Jan22-1_EF5N_S73\tPL:ILLUMINA)

2. Edit config.yaml for your project:

Methods description from Lack et al. 2015 Genetics 10.1534/genetics.115.174664:


"In brief, we first mapped short read data to the D. melanogaster reference genome (release 5.57; http://flybase.org) using BWA v0.5.9 (Li and Durbin 2010) using default settings, followed by mapping of all unmapped reads using Stampy v1.0.20 (Lunter and Goodson 2010). This approach combines the rapid but strict BWA algorithm to first map the relatively “easy-to-align” reads with the more sensitive but computationally intensive Stampy algorithm, which more effectively and accurately aligns the relatively divergent reads (Lunter and Goodson 2010). All reads with mapping quality scores <20 were excluded. Optical duplicates were then removed using Picard v1.79 (http://picard.sourceforge.net/), and assemblies were improved around indels using the GATK v3.2 Indel Realigner (McKenna et al. 2010; Depristo et al. 2011). The Unified Genotyper (Depristo et al. 2011) was then used to call indels and SNPs for each individual genome. Among the indel calling criteria, >50% of the reads at a given position had to support the existence of that indel, with a minimum of three reads containing the variant. For SNP calling in this first round, we required a minimum phred-scaled quality value of 31 and that >75% of the reads at a given position support the SNP. For the second round of assembly, the SNPs and indels called in the first round were introduced into the D. melanogaster reference, and this modified reference was then used for a sec- ond round of mapping. Following indel realignment, the Unified Genotyper was then used to call all sites in the modified reference genome. To generate reference-numbered consensus sequences, a custom perl script was used to shift all base coor- dinates back to those of the original D. melanogaster reference. Deletions and all sites within 3 bp of a called indel were coded as “N” (based on the error analysis described in the Results), while insertions do not appear in reference-numbered consensus sequences."

![DAG](https://github.com/JamieCFreeman/DGN_compatible/blob/main/README_files/rulegraph.pdf?raw=true)
![DAG2](https://github.com/JamieCFreeman/DGN_compatible/blob/main/README_files/round2_dag.pdf?raw=true)
