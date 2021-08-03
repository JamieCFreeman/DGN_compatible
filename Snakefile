#### 2021-8-1 JCF ####

configfile: "config.yaml"

import pandas as pd

from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")

SAMPLES = ["6Jul21-10_S111", "6Jul21-11_S112"]
  
OUTDIR = config["prefix"]

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
		expand(f"{OUTDIR}/bwa/{{sample}}.bam", sample=SAMPLES)

rule test:
	input:
		fq1 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R1_001.fastq.gz",
		fq2 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R2_001.fastq.gz"
	output:
		"test/{sample}.txt"
	shell:
		"wc -l {input.fq1} > {output}"

rule bwa_map_1:
	input:
                fq1 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R1_001.fastq.gz"
	output:
		sai1 = f"{OUTDIR}/bwa/{{sample}}_1.sai"
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/bwa_{{sample}}_map1.log"
	shell:
		"bwa aln {params.REF} {input.fq1} > {output.sai1}"

rule bwa_map_2:
        input:  
                fq2 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R2_001.fastq.gz"
        output: 
                sai2 = f"{OUTDIR}/bwa/{{sample}}_2.sai"
	params: REF = config["genome"]
	log: f"{OUTDIR}/logs/bwa_{{sample}}_map1.log"
        shell:  
                "bwa aln {params.REF} {input.fq2} > {output.sai2}"

rule bwa_sampe:
        input:  
                sai1 = f"{OUTDIR}/bwa/{{sample}}_1.sai",
                sai2 = f"{OUTDIR}/bwa/{{sample}}_2.sai",
		fq1 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R1_001.fastq.gz",
                fq2 = "/home/jamie/EF_genomes/210728_AHCHLLDSX2/{sample}_L003_R2_001.fastq.gz"
        params: REF = config["genome"]
	output: 
                f"{OUTDIR}/bwa//{{sample}}.bam"
        shell: "bwa sampe -P {params.REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} | samtools view -bS - > {output}"

#rule stampy_map:
#	input:
#		f"{OUTDIR}/{{sample}}.bam"
#	output:
#		"test/{sample}_remap.sam"
#	params: REF = config["genome"]
#	log: "logs/stampy_{sample}.log"
#	conda: "envs/py2.yaml"
#	shell:
#		"python stampy.py -g {params.REF} -h {params.REF} --bamkeepgoodreads -M {input} -o {output}"
# "python2.6 " . $stampy . "stampy.py -g " . $reference . " -h " . $reference . " --bamkeepgoodreads -M " . $FastqFile[$i] . ".bam -o " . $FastqFile[$i] . "_remapped.sam"; #Stampy MAPPING STEP. THIS IS A RELATIVELY LONG STEP (AS LONG AS 15 HOURS ON SOME OF THE HIGHEST COVERAGE DPGP2 GENOMES)
 
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
#rule bwa_align:
#  input:
#    get_samples
#  output:
#    f"{OUTDIR}/results/bam/{{sample}}.bam"
#  log: f"{OUTDIR}/logs/bwa/{{sample}}.log"
#  shell:
#    """
#    bwa align {input}_1.fq > {output}_1.sai
#    bwa align {input}_2.fq > {output}_2.sai
#    """
