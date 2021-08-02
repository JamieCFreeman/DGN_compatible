#### 2021-8-1 JCF ####

configfile: "config.yaml"
  
OUTDIR = config["prefix"]

def get_samples(wildcards):
    return config["samples"][wildcards.sample]
  
rule fastqc:
	input:
		get_samples
	output:
		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.html",
		f"{OUTDIR}/results/fastqc/{{sample}}_fastqc.zip"
	threads: 4
	params: 
		outdir = f"{OUTDIR}/results/fastqc"
	log: f"{OUTDIR}/logs/fastqc/{{sample}}.log"
	conda: "envs/fastqc.yaml"	
	shell:
		"fastqc --outdir {params.outdir} --format fastq --threads {threads} {input}"
    
    
rule bwa_index:
  input:
  output:
  conda:
  shell:

rule bwa_align:
  input:
    get_samples
  output:
    f"{OUTDIR}/results/bam/{{sample}}.bam"
  log: f"{OUTDIR}/logs/bwa/{{sample}}.log"
  shell:
    """
    bwa align {input}_1.fq > {output}_1.sai
    bwa align {input}_2.fq > {output}_2.sai
    """
