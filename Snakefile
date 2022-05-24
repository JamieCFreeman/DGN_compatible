#### 2021-8-1 JCF ####
# Implementing Justin Lack's Drosophila Genome Nexus pipeline (https://github.com/justin-lack/Drosophila-Genome-Nexus) #

# If IS testing is true, go through pipeline with small number of reads for testing the workflow. 
ISTESTING = 'TRUE'

configfile: "config.yaml"

import pandas as pd

from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

#OUTDIR = config["prefix"]
OUTDIR = "round_testing"

samples_table = pd.read_table(config["sample_table"], dtype=str).set_index("sample", drop=False)

# Get sample wildcards as a list
SAMPLES= samples_table['sample'].values.tolist()

ROUND = 2

wildcard_constraints:
	sample="|".join(samples_table['sample'])

def fq1_from_sample(wildcards, ISTESTING):
	if ISTESTING == 'TRUE':
		words = ["test/" , wildcards.sample, "_1.fq"]
		return "".join(words)
	elif ISTESTING == 'FALSE':
		return samples_table.loc[wildcards.sample, "fq1"]

def fq2_from_sample(wildcards, ISTESTING):
	if ISTESTING == 'TRUE':
		words = ["test/" , wildcards.sample, "_2.fq"]
		return "".join(words)
	elif ISTESTING == 'FALSE':
		return samples_table.loc[wildcards.sample, "fq2"]

def get_ref_fa(wildcards, ROUND, OUTDIR):
	if ROUND == 1:
		return config["genome"]
	elif ROUND == 2:
		words = [OUTDIR, '/round1/alt_ref/', wildcards.sample, '_ref.fasta']
		return "".join(words)

def get_ref_idx(wildcards, ROUND, OUTDIR, idxtype):
	if ROUND == 1:
		words = [config["genome"], idxtype]
		return "".join(words)
	if ROUND == 2:
		words = [OUTDIR, '/round1/alt_ref/', wildcards.sample, '_ref.fasta', idxtype]
		return "".join(words)

def RG_from_sample(wildcards):
	return samples_table.loc[wildcards.sample, "RG"]

rule all:
	input:
#		expand(f"{OUTDIR}/logs/bwa_aln/{{sample}}.stats", sample=SAMPLES),
#		expand(f"{OUTDIR}/logs/bwa_mem/{{sample}}.stats", sample=SAMPLES),
#		expand(f"{OUTDIR}/logs/stampy/{{sample}}_dups.txt", sample=SAMPLES),
		expand(f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.stidx", sample=SAMPLES),
		expand(f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.bwt", sample=SAMPLES)

#rule indexes:
#	input:
#		REF_fai = lambda wc: get_ref_idx(wc, ROUND, OUTDIR, "fai") 

rule test:
	input:
		fq1 = lambda wc: fq1_from_sample(wc, 'FALSE'),
		fq2 = lambda wc: fq2_from_sample(wc, 'FALSE')
	output:
		cut1 = "test/{sample}_1.fq",
		cut2 = "test/{sample}_2.fq"
	shell:
		"zcat {input.fq1} | head -n 50000 > {output.cut1}; zcat {input.fq2} | head -n 50000 > {output.cut2}"

#rule ref_fai:
#	input:
#		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
#	output:
#		REF_fai = lambda wc: get_ref_idx(wc, ROUND, OUTDIR, "fai")		
#	shell:
#		"samtools index {input.REF}"

rule bwa_aln_1:
	input:
		cut1 = "test/{sample}_1.fq",
		fq1 = lambda wc: fq1_from_sample(wc, ISTESTING),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		sai1 = temp( f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_1.sai")
#	params: REF = config["genome"]
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}_map1.log"
	threads: 8
	shell:
		"bwa aln -t {threads} {input.REF} {input.fq1} 2> {log} > {output.sai1}"

rule bwa_aln_2:
	input:  
		cut2 = "test/{sample}_2.fq",
		fq2 = lambda wc: fq2_from_sample(wc, ISTESTING),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output: 
		sai2 = temp(f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_2.sai")
#	params: REF = config["genome"]
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}_map2.log"
	threads: 8
	shell:  
		"bwa aln -t {threads} {input.REF} {input.fq2} 2> {log} > {output.sai2}"

rule bwa_sampe:
	input:
		sai1 = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_1.sai",
		sai2 = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_2.sai",
		fq1 = lambda wc: fq1_from_sample(wc, ISTESTING),
		fq2 = lambda wc: fq2_from_sample(wc, ISTESTING),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output: 
		f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}.bam"
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_sampe/{{sample}}_map2.log"
	shell:
		"bwa sampe -P {input.REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} 2> {log} | samtools view -bS - > {output}"

rule aln_flagstat:
	input:
		f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}.bam"
	output:
		f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}.stats"
	shell:
		"samtools flagstat {input} > {output}"
rule bwa_mem:
	input:
		fq1 = "test/{sample}_1.fq",
		fq2 = "test/{sample}_2.fq"
	params: REF = config["genome"]
	log:  f"{OUTDIR}/logs/bwa_mem/{{sample}}_map1.log"
	output:
		f"{OUTDIR}/bwa_mem/{{sample}}.bam"
	threads: 8
	shell:
		"bwa mem -M -t {threads} {params.REF} {input.fq1} {input.fq2} | samtools view -bS - > {output}"

rule flagstat:
	input:
		f"{OUTDIR}/bwa_mem/{{sample}}.bam"
	output:
		f"{OUTDIR}/logs/bwa_mem/{{sample}}.stats"
	shell:
		"samtools flagstat {input} > {output}"


rule stampy_map:
	input:
		bam = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}.bam",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/{{sample}}.sam")
	params: stampy = "/opt/bioscript/stampy/stampy.py"
	log: f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}.log"
	conda: "envs/py2.yaml"
	shell:
		"python {params.stampy} -g {input.REF} -h {input.REF} --bamkeepgoodreads -M {input.bam} -o {output} 2> {log} "
# "python2.6 " . $stampy . "stampy.py -g " . $reference . " -h " . $reference . " --bamkeepgoodreads -M " . $FastqFile[$i] . ".bam -o " . $FastqFile[$i] . "_remapped.sam"; 
# #Stampy MAPPING STEP. THIS IS A RELATIVELY LONG STEP (AS LONG AS 15 HOURS ON SOME OF THE HIGHEST COVERAGE DPGP2 GENOMES)

rule sam2bam:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/{{sample}}.sam"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/{{sample}}.bam")
	shell:
		"samtools view -bS {input} > {output}"

rule stampy_flagstat:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/{{sample}}.bam"
	output:
		f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}.stats"
	shell:
		"samtools flagstat {input} > {output}"
		
rule qfilter_bam:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/{{sample}}.bam",
		stats = f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}.stats"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}.bam")
	shell:
		"samtools view -q 20 -h {input.bam} > {output}"

rule sort_bam:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}.bam"
	output:
		f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_sort.bam"
	shell:
		"samtools sort {input} > {output}"

# Jeremy filters out unmapped reads next with Picard CleanSam.jar
#
#


# Now mark dups java -Xmx" . $mem . "g -jar " . $picard . "MarkDuplicates.jar INPUT=" . $FastqFile[$i] . "sort.bam OUTPUT=" . $FastqFile[$i] . "dups.bam METRICS_FILE=" . $FastqFile[$i] . "dups.metrics"; #IDENTIFIES DUPLICATE READS
#

rule mark_dups:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_sort.bam"
	output:
		bam = temp(f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam"),
		metrics = f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}_dups.txt"
	conda:  "envs/picard.yaml"
	shell:
		"picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"	
#		"java -Xmx4g -jar picard.jar INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"

# Function RG_from_sample provides RG info from sample table
rule add_RG:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam")
	params: RG= lambda wildcards: RG_from_sample(wildcards)
	shell:
		"samtools addreplacerg -r '{params.RG}' -o {output} {input}"

rule dup_bam_bai:
        input:
                f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam"
        output:
                temp(f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam.bai")
        shell:
                "samtools index {input}"

#IDENTIFIES INTERVAL TO BE REALIGNED AROUND INDELS
rule indel_target:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam",
		bai = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam.bai",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.intervals")
	params: 
		#REF = config["genome"],
		GATKjar = config["gatk.jar"]
	conda: "envs/java8.yaml" # GATK needs java 8 (11 default on marula)
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/RTC/{{sample}}.log"		
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T RealignerTargetCreator \
			-R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule indel_realign:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam",
		bai = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/RG_{{sample}}.bam.bai",
		intervals = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.intervals",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam"
	params: 
		#REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/indel_realigner/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx8g -jar {params.GATKjar} -T IndelRealigner \
			-targetIntervals {input.intervals} \
			-R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule gt_snps:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_SNPs.vcf"
	params: 
		#REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/gt_snps/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T UnifiedGenotyper \
			-mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 1 \
			-R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule gt_indels:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_INDELs.vcf"
	params: 
		#REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/gt_indels/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T UnifiedGenotyper \
			-mbq 10 -stand_call_conf 31 -stand_emit_conf 31 \
			-ploidy 1 -minIndelFrac 0.51 -minIndelCnt 3 -glm INDEL \
			-R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule filter_vcf:
	input:
		snp = f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_SNPs.vcf"
	output:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_SNPs_filt.vcf"
#	conda: "envs/snpsift.yaml"
	shell:
		"cp {input.snp} {output}"

rule add_snps_alt_ref:
	input:
		vcf = f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_SNPs_filt.vcf",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.fasta"
	params: 
		#REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/alt_ref_SNPs/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T FastaAlternateReferenceMaker \
			-R {input.REF} -V {input.vcf} -o {output} 2> {log}
		"""

rule alt_ref_index:
	input:
		fa = f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.fasta",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.dict"	
	params: 
		#REF = config["genome"],
		picardjar = config["picard.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/alt_ref_index/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		samtools faidx {input.fa}; \
		java -Xmx4g -jar {params.picardjar}/CreateSequenceDictionary.jar \
			REFERENCE={input.fa} OUTPUT={output} 2> {log}
		"""

rule add_indel_alt_ref:
	input:
		ref = f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.fasta",
		dict = f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.dict",
		indel = f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_INDELs.vcf"
	output:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta"	
	params: GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/alt_ref_indel/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T FastaAlternateReferenceMaker \
			-R {input.ref} -V {input.indel} -o {output} 2> {log} 
		"""

rule bwa_index:
	input:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta"
	output:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.bwt"
	shell:
		"bwa index -a bwtsw {input}"	

rule stampy_index:
	input:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta"
	output:
		stidx = f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.stidx",
		sthash = f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.sthash"
	params: stampy = "/opt/bioscript/stampy/stampy.py"
	log: f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}_index.log" 
	conda: "envs/py2.yaml"
	shell:
		"""
		TEMP_FILE=`echo {output.stidx} | sed 's ^.*alt_ref/  '| sed 's/.stidx//'`; echo $TEMP_FILE #stampy writes stidx file to working dir
		python {params.stampy} -G {wildcards.sample}_ref.fasta {input} # build genome file
		mv $TEMP_FILE.stidx {output.stidx} # move stidx to alt_ref folder
		python {params.stampy} -g {input} -H {wildcards.sample}_ref.fasta # build hash table
		mv $TEMP_FILE.sthash {output.sthash}
		"""
