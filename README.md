# DGN_compatible
Pipeline to analyze DNA sequencing data with legacy compatibility for the Pool lab's Drosophila Genome Nexus data.
You likely only want to use this pipeline if you have data you want to compare or use in conjunction with the DGN dataset.

To run:
1. Create a sample table in the format of "units.tsv" (Tab-delimited table with columns sample, unit, platform, fq1, fq2, RG)
sample (sample name), unit (I typically use sample number string from raw data file name), 
platform ("ILLUMINA" or not), fq1 (path to fq1), fq2 (path to fq2), RG (read group string, eg @RG\tID:13Jan22-1_EF5N_S73\tSM:13Jan22-1_EF5N_S73\tPU:none\tLB:13Jan22-1_EF5N_S73\tPL:ILLUMINA)


2. Edit config.yaml for your project. Recommended to first run with the ISTESTING variable set to TRUE- this will subset reads, so the whole
workflow will run quickly, and you can make sure it runs cleanly before starting.

3. Use snakemake to run with env managment.
Make sure snakemake is installed in your env. (On Marula, I have a conda env called snakemake). 
Use the dry run function first (snakemake -n) to make sure the job list matches what you expect.
You must run with conda (--use-conda --conda-frontend conda flags). Specify cores with -c, and optionally specify memory usage
with --resources mem_Gb. NOTE: Snakemake doesn't actually monitor resources (eg memory here), but you can add this tag to
resource-intensive steps to limit how many of these type of jobs are run at once. For now mem_Gb is set for the rules 
running MarkDuplicates (3 Gb) and IndelRealigner (2 Gb). <br>
To run rule 'all' for either round 1 or 2 mapping, use:
snakemake -c48 all --use-conda --conda-frontend conda --resources mem_Gb=56

Round 1 and 2 of the pipeline are not connected here because I typically do the mapping on UW's CHTC cluster- 
stampy is very slow and benefits from being highly parallelized.  
See https://github.com/JamieCFreeman/chtc_bwa_stampy for my CHTC pipeline. Once you have completed mapping there, 
you will have a merged bam file for every sample-unit. You simply need to 
place those bam files where Snakemake expects them to be, and as long as the timestamps of the new bams are more 
recent than the input files to create those bams, Snakemake will not rerun the mapping step. Always make sure to try a dry run with the -n 
flag to make sure what Snakemake thinks should be done agrees with your expectation.

After the round 1 mapping is complete:
If transferring .bam files from CHTC, make sure the round 2 indexes are timestamped first, they can be created using:
snakemake -c48 round2_index --use-conda --conda-frontend conda --resources mem_Gb=56

Once the round 2 mapping is complete, to get .fas1k files, run with rule fas1k:
snakemake fas1k --use-conda --conda-frontend conda -c48

On fas1k files can also run rule stats to get heteozygosity stats, genotyping of inversion-associated SNPs (use with caution):
snakemake stats --use-conda --conda-frontend conda -c48 --resources io=30

Methods description from Lack et al. 2015 Genetics 10.1534/genetics.115.174664:

"In brief, we first mapped short read data to the D. melanogaster reference genome (release 5.57; http://flybase.org) using BWA v0.5.9 (Li and Durbin 2010) using default settings, followed by mapping of all unmapped reads using Stampy v1.0.20 (Lunter and Goodson 2010). This approach combines the rapid but strict BWA algorithm to first map the relatively “easy-to-align” reads with the more sensitive but computationally intensive Stampy algorithm, which more effectively and accurately aligns the relatively divergent reads (Lunter and Goodson 2010). All reads with mapping quality scores <20 were excluded. Optical duplicates were then removed using Picard v1.79 (http://picard.sourceforge.net/), and assemblies were improved around indels using the GATK v3.2 Indel Realigner (McKenna et al. 2010; Depristo et al. 2011). The Unified Genotyper (Depristo et al. 2011) was then used to call indels and SNPs for each individual genome. Among the indel calling criteria, >50% of the reads at a given position had to support the existence of that indel, with a minimum of three reads containing the variant. For SNP calling in this first round, we required a minimum phred-scaled quality value of 31 and that >75% of the reads at a given position support the SNP. For the second round of assembly, the SNPs and indels called in the first round were introduced into the D. melanogaster reference, and this modified reference was then used for a sec- ond round of mapping. Following indel realignment, the Unified Genotyper was then used to call all sites in the modified reference genome. To generate reference-numbered consensus sequences, a custom perl script was used to shift all base coor- dinates back to those of the original D. melanogaster reference. Deletions and all sites within 3 bp of a called indel were coded as “N” (based on the error analysis described in the Results), while insertions do not appear in reference-numbered consensus sequences."

![DAG](https://github.com/JamieCFreeman/DGN_compatible/blob/main/README_files/rulegraph.pdf?raw=true)
![DAG2](https://github.com/JamieCFreeman/DGN_compatible/blob/main/README_files/round2_dag.pdf?raw=true)
