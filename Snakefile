#### 2021-8-1 JCF ####
# Implementing Justin Lack's Drosophila Genome Nexus pipeline (https://github.com/justin-lack/Drosophila-Genome-Nexus) #

# If IS testing is true, go through pipeline with small number of reads for testing the workflow. 
ISTESTING = 'FALSE'

configfile: "config.yaml"

import pandas as pd

from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

OUTDIR = config["prefix"]

samples_table = pd.read_table(config["sample_table"], dtype=str).set_index("sample", drop=False)

# Get sample wildcards as a list
SAMPLES= samples_table['sample'].values.tolist()

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

def RG_from_sample(wildcards):
	return samples_table.loc[wildcards.sample, "RG"]

rule all:
	input:
#		expand(f"{OUTDIR}/logs/bwa_aln/{{sample}}.stats", sample=SAMPLES),
#		expand(f"{OUTDIR}/logs/bwa_mem/{{sample}}.stats", sample=SAMPLES),
#		expand(f"{OUTDIR}/logs/stampy/{{sample}}_dups.txt", sample=SAMPLES),
		expand(f"{OUTDIR}/alt_ref/{{sample}}_ref.fasta", sample=SAMPLES)

rule test:
	input:
		fq1 = lambda wc: fq1_from_sample(wc, 'FALSE'),
		fq2 = lambda wc: fq2_from_sample(wc, 'FALSE')
	output:
		cut1 = "test/{sample}_1.fq",
		cut2 = "test/{sample}_2.fq"
	shell:
		"zcat {input.fq1} | head -n 50000 > {output.cut1}; zcat {input.fq2} | head -n 50000 > {output.cut2}"

rule bwa_aln_1:
	input:
		cut1 = "test/{sample}_1.fq",
		fq1 = lambda wc: fq1_from_sample(wc, ISTESTING)
	output:
		sai1 = temp( f"{OUTDIR}/bwa_aln/{{sample}}_1.sai")
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/bwa_aln/{{sample}}_map1.log"
	threads: 8
	shell:
		"bwa aln -t {threads} {params.REF} {input.fq1} 2> {log} > {output.sai1}"

rule bwa_aln_2:
	input:  
		cut2 = "test/{sample}_2.fq",
		fq2 = lambda wc: fq2_from_sample(wc, ISTESTING)
	output: 
		sai2 = temp(f"{OUTDIR}/bwa_aln/{{sample}}_2.sai")
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/bwa_aln/{{sample}}_map2.log"
	threads: 8
	shell:  
		"bwa aln -t {threads} {params.REF} {input.fq2} 2> {log} > {output.sai2}"

rule bwa_sampe:
	input:
		sai1 = f"{OUTDIR}/bwa_aln/{{sample}}_1.sai",
		sai2 = f"{OUTDIR}/bwa_aln/{{sample}}_2.sai",
		fq1 = lambda wc: fq1_from_sample(wc, ISTESTING),
		fq2 = lambda wc: fq2_from_sample(wc, ISTESTING)
	params: REF = config["genome"]
	output: 
		f"{OUTDIR}/bwa_aln/{{sample}}.bam"
	log: f"{OUTDIR}/logs/bwa_sampe/{{sample}}_map2.log"
	shell:
		"bwa sampe -P {params.REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} 2> {log} | samtools view -bS - > {output}"

rule aln_flagstat:
	input:
		f"{OUTDIR}/bwa_aln/{{sample}}.bam"
	output:
		f"{OUTDIR}/logs/bwa_aln/{{sample}}.stats"
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
		f"{OUTDIR}/bwa_aln/{{sample}}.bam"
	output:
		temp(f"{OUTDIR}/stampy/{{sample}}.sam")
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/stampy/{{sample}}.log"
	conda: "envs/py2.yaml"
	shell:
		"python /opt/bioscript/stampy/stampy.py -g {params.REF} -h {params.REF} --bamkeepgoodreads -M {input} -o {output} 2> {log} "
# "python2.6 " . $stampy . "stampy.py -g " . $reference . " -h " . $reference . " --bamkeepgoodreads -M " . $FastqFile[$i] . ".bam -o " . $FastqFile[$i] . "_remapped.sam"; 
# #Stampy MAPPING STEP. THIS IS A RELATIVELY LONG STEP (AS LONG AS 15 HOURS ON SOME OF THE HIGHEST COVERAGE DPGP2 GENOMES)

rule sam2bam:
	input:
		f"{OUTDIR}/stampy/{{sample}}.sam"
	output:
		temp(f"{OUTDIR}/stampy/{{sample}}.bam")
	shell:
		"samtools view -bS {input} > {output}"

rule stampy_flagstat:
	input:
		f"{OUTDIR}/stampy/{{sample}}.bam"
	output:
		f"{OUTDIR}/logs/stampy/{{sample}}.stats"
	shell:
		"samtools flagstat {input} > {output}"
		
rule qfilter_bam:
	input:
		bam = f"{OUTDIR}/stampy/{{sample}}.bam",
		stats = f"{OUTDIR}/logs/stampy/{{sample}}.stats"
	output:
		temp(f"{OUTDIR}/stampy/qfilter_{{sample}}.bam")
	shell:
		"samtools view -q 20 -h {input.bam} > {output}"

rule sort_bam:
	input:
		f"{OUTDIR}/stampy/qfilter_{{sample}}.bam"
	output:
		f"{OUTDIR}/stampy/qfilter_{{sample}}_sort.bam"
	shell:
		"samtools sort {input} > {output}"

# Jeremy filters out unmapped reads next with Picard CleanSam.jar
#
#


# Now mark dups java -Xmx" . $mem . "g -jar " . $picard . "MarkDuplicates.jar INPUT=" . $FastqFile[$i] . "sort.bam OUTPUT=" . $FastqFile[$i] . "dups.bam METRICS_FILE=" . $FastqFile[$i] . "dups.metrics"; #IDENTIFIES DUPLICATE READS
#

rule mark_dups:
	input:
		f"{OUTDIR}/stampy/qfilter_{{sample}}_sort.bam"
	output:
		bam = temp(f"{OUTDIR}/stampy/mark_dup/{{sample}}.bam"),
		metrics = f"{OUTDIR}/logs/stampy/{{sample}}_dups.txt"
	conda:  "envs/picard.yaml"
	shell:
		"picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"	
#		"java -Xmx4g -jar picard.jar INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"

# Function RG_from_sample provides RG info from sample table
rule add_RG:
	input:
		bam = f"{OUTDIR}/stampy/mark_dup/{{sample}}.bam"
	output:
		temp(f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam")
	params: RG= lambda wildcards: RG_from_sample(wildcards)
	shell:
		"samtools addreplacerg -r '{params.RG}' -o {output} {input}"

rule dup_bam_bai:
        input:
                f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam"
        output:
                temp(f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam.bai")
        shell:
                "samtools index {input}"

#IDENTIFIES INTERVAL TO BE REALIGNED AROUND INDELS
rule indel_target:
	input:
		bam = f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam",
		bai = f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam.bai"
	output:
		temp(f"{OUTDIR}/stampy/indel_realign/{{sample}}.intervals")
	params: REF = config["genome"],
		GATKjar = config["gatk.jar"]
	conda: "envs/java8.yaml" # GATK needs java 8 (11 default on marula)
	log: f"{OUTDIR}/logs/gatk/RTC/{{sample}}.log"		
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T RealignerTargetCreator \
			-R {params.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule indel_realign:
	input:
		bam = f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam",
		bai = f"{OUTDIR}/stampy/mark_dup/RG_{{sample}}.bam.bai",
		intervals = f"{OUTDIR}/stampy/indel_realign/{{sample}}.intervals"
	output:
		f"{OUTDIR}/stampy/indel_realign/{{sample}}.bam"
	params: REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/logs/gatk/indel_realigner/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx8g -jar {params.GATKjar} -T IndelRealigner \
			-targetIntervals {input.intervals} \
			-R {params.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule gt_snps:
	input:
		bam = f"{OUTDIR}/stampy/indel_realign/{{sample}}.bam"
	output:
		f"{OUTDIR}/vcf/{{sample}}_SNPs.vcf"
	params: REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/logs/gatk/gt_snps/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T UnifiedGenotyper \
			-mbq 10 -stand_call_conf 31 -stand_emit_conf 31 -ploidy 1 \
			-R {params.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule gt_indels:
	input:
		bam = f"{OUTDIR}/stampy/indel_realign/{{sample}}.bam"
	output:
		f"{OUTDIR}/vcf/{{sample}}_INDELs.vcf"
	params: REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/logs/gatk/gt_indels/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T UnifiedGenotyper \
			-mbq 10 -stand_call_conf 31 -stand_emit_conf 31 \
			-ploidy 1 -minIndelFrac 0.51 -minIndelCnt 3 -glm INDEL \
			-R {params.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule filter_vcf:
	input:
		snp = f"{OUTDIR}/vcf/{{sample}}_SNPs.vcf"
	output:
		f"{OUTDIR}/vcf/{{sample}}_SNPs_filt.vcf"
#	conda: "envs/snpsift.yaml"
	shell:
		"cp {input.snp} {output}"

rule add_snps_alt_ref:
	input:
		f"{OUTDIR}/vcf/{{sample}}_SNPs_filt.vcf"
	output:
		f"{OUTDIR}/alt_ref/{{sample}}_SNPs_ref.fasta"
	params: REF = config["genome"],
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/logs/gatk/alt_ref_SNPs/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T FastaAlternateReferenceMaker \
			-R {params.REF} -V {input} -o {output} 2> {log}
		"""

rule alt_ref_index:
	input:
		f"{OUTDIR}/alt_ref/{{sample}}_SNPs_ref.fasta"
	output:
		f"{OUTDIR}/alt_ref/{{sample}}_SNPs_ref.dict"	
	params: REF = config["genome"],
		picardjar = config["picard.jar"]
	log: f"{OUTDIR}/logs/gatk/alt_ref_index/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		samtools faidx {input}; \
		java -Xmx4g -jar {params.picardjar}/CreateSequenceDictionary.jar \
			REFERENCE={input} OUTPUT={output} 2> {log}
		"""

rule add_indel_alt_ref:
	input:
		ref = f"{OUTDIR}/alt_ref/{{sample}}_SNPs_ref.fasta",
		dict = f"{OUTDIR}/alt_ref/{{sample}}_SNPs_ref.dict",
		indel = f"{OUTDIR}/vcf/{{sample}}_INDELs.vcf"
	output:
		f"{OUTDIR}/alt_ref/{{sample}}_ref.fasta"	
	params: GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/logs/gatk/alt_ref_indel/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T FastaAlternateReferenceMaker \
			-R {input.ref} -V {input.indel} -o {output} 2> {log} 
		"""






