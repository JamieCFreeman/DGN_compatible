#### 2021-8-1 JCF ####

configfile: "config.yaml"

import pandas as pd

from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

#SAMPLES = ["6Jul21-10_S111", "6Jul21-11_S112", "6Jul21-12_S113", "6Jul21-13_S114", "6Jul21-14_S115" ]
#SAMPLES = ["6Jul21-10_S111", "6Jul21-11_S112"]
SAMPLES = ["6Jul21-10_S111", "6Jul21-11_S112", "6Jul21-12_S113", "6Jul21-13_S114", "6Jul21-14_S115", "6Jul21-15_S116", "6Jul21-16_S117", "6Jul21-1_S102", "6Jul21-2_S103", "6Jul21-3_S104", "6Jul21-4_S105", "6Jul21-5_S106", "6Jul21-6_S107", "6Jul21-7_S108", "6Jul21-8_S109", "6Jul21-9_S110"]
 
OUTDIR = config["prefix"]

samples_table = pd.read_table("units.tsv", dtype=str).set_index("sample", drop=False)
#samples_table.index = samples_table.index.set_levels( [i.astype(str) for i in samples_table.index.levels] )

def fq1_from_sample(wildcards):
	return samples_table.loc[wildcards.sample, "fq1"]

def fq2_from_sample(wildcards):
        return samples_table.loc[wildcards.sample, "fq2"]

#samples = pd.read_table(config["samples"])
#SAMPLES = samples[["sample"]]
#units = pd.read_table(config["units"], dtype=str).set_index(
#    ["sample", "unit"], drop=False
#)
#units.index = units.index.set_levels(
#    [i.astype(str) for i in units.index.levels]
#)  # enforce str in index

#wildcard_constraints:
#    sample="|".join(samples.index),
#    unit="|".join(units["unit"]),

rule all:
	input:
		expand(f"{OUTDIR}/logs/bwa_aln/{{sample}}.stats", sample=SAMPLES),
		expand(f"{OUTDIR}/logs/bwa_mem/{{sample}}.stats", sample=SAMPLES),
		expand(f"{OUTDIR}/logs/stampy/{{sample}}_dups.txt", sample=SAMPLES)
rule test:
	input:
		fq1 = fq1_from_sample,
		fq2 = fq2_from_sample
#		fq1 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R1_001.fastq.gz",
#		fq2 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R2_001.fastq.gz"
	output:
		cut1 = "test/{sample}_1.fq",
		cut2 = "test/{sample}_2.fq"
	shell:
		"zcat {input.fq1} | head -n 12000 > {output.cut1}; zcat {input.fq2} | head -n 12000 > {output.cut2}"

rule bwa_aln_1:
	input:
		fq1 = "test/{sample}_1.fq"
	output:
		sai1 = temp( f"{OUTDIR}/bwa_aln/{{sample}}_1.sai")
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/bwa_aln/{{sample}}_map1.log"
	threads: 8
	shell:
		"bwa aln -t {threads} {params.REF} {input.fq1} 2> {log} > {output.sai1}"

rule bwa_aln_2:
	input:  
		fq2 = "test/{sample}_2.fq"
		#fq2 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R2_001.fastq.gz"
	output: 
		sai2 = temp(f"{OUTDIR}/bwa_aln/{{sample}}_2.sai")
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/bwa_aln/{{sample}}_map2.log"
	threads: 8
	shell:  
		"bwa aln -t {threads} {params.REF} {input.fq2} {log} > {output.sai2}"

rule bwa_sampe:
	input:
		sai1 = f"{OUTDIR}/bwa_aln/{{sample}}_1.sai",
		sai2 = f"{OUTDIR}/bwa_aln/{{sample}}_2.sai",
		fq1 = "test/{sample}_1.fq",
		fq2 = "test/{sample}_2.fq"
		#fq1 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R1_001.fastq.gz",
                #fq2 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R2_001.fastq.gz"
	params: REF = config["genome"]
	output: 
		f"{OUTDIR}/bwa_aln/{{sample}}.bam"
	shell:
		"bwa sampe -P {params.REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} | samtools view -bS - > {output}"

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
		f"{OUTDIR}/stampy/{{sample}}.bam"
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
		temp(f"{OUTDIR}/stampy/qfilter/{{sample}}.bam")
	shell:
		"samtools view -q 20 -h {input.bam} > {output}"

rule sort_bam:
	input:
		f"{OUTDIR}/stampy/qfilter/{{sample}}.bam"
	output:
		f"{OUTDIR}/stampy/qfilter/{{sample}}_sort.bam"
	shell:
		"samtools sort {input} > {output}"

# Jeremy filters out unmapped reads next with Picard CleanSam.jar
#
#


# Now mark dups java -Xmx" . $mem . "g -jar " . $picard . "MarkDuplicates.jar INPUT=" . $FastqFile[$i] . "sort.bam OUTPUT=" . $FastqFile[$i] . "dups.bam METRICS_FILE=" . $FastqFile[$i] . "dups.metrics"; #IDENTIFIES DUPLICATE READS
#

rule mark_dups:
	input:
		f"{OUTDIR}/stampy/qfilter/{{sample}}_sort.bam"
	output:
		bam = f"{OUTDIR}/stampy/mark_dup/{{sample}}.bam",
		metrics = f"{OUTDIR}/logs/stampy/{{sample}}_dups.txt"
	conda:  "envs/picard.yaml"
	shell:
		"picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"	
#		"java -Xmx4g -jar picard.jar INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"

#rule fastqc:
#	input:
#		get_samples
#	output:
#		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.html",
#		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.zip"
#	threads: 4
#	params: 
#		outdir = f"{OUTDIR}/results/fastqc"
#	log: f"{OUTDIR}/logs/fastqc/{{sample}}.log"
#	conda: "envs/fastqc.yaml"	
#	shell:
#		"fastqc --outdir {params.outdir} --format fastq --threads {threads} {input}"
#    
#    
#rule bwa_index:
#  input:
#  output:
#  conda:
#  shell:
#
