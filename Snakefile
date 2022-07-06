#### 2021-8-1 JCF ####
# Implementing Justin Lack's Drosophila Genome Nexus pipeline (https://github.com/justin-lack/Drosophila-Genome-Nexus) #

import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
min_version("7.0.0")

# Pass through variables from config file
configfile: "config.yaml"
OUTDIR = config["prefix"]
# If IS testing is true, go through pipeline with small number of reads for testing the workflow. 
ISTESTING = config["ISTESTING"]
# Are you doing round 1 mapping or round 2 mapping?
ROUND = config["ROUND"]

# Read sample table
samples_table = pd.read_table(config["sample_table"], dtype=str).set_index("sample", drop=False)

# Get sample wildcards as a list
SAMPLES= samples_table['sample'].values.tolist()

# Constrain sample wildcards to those in sample table
wildcard_constraints:
	sample="|".join(samples_table['sample'])

def fq_from_sample(wildcards, ISTESTING, READ):
	if ISTESTING == 'TRUE':
		words = ["test/" , wildcards.sample, "_", READ, ".fq"]
		return "".join(words)
	elif ISTESTING == 'FALSE':
		return samples_table.loc[wildcards.sample, "fq"+ READ]

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

include: 
	"rules/genome_indexing.smk",
	"rules/qc.smk"

rule all:
	input:
#		expand(f"{OUTDIR}/logs/stampy/{{sample}}_dups.txt", sample=SAMPLES),
		expand(f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.stidx", sample=SAMPLES),
		expand(f"{OUTDIR}/round{ROUND}_index.ok"),
		expand(f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta.bwt", sample=SAMPLES)

# Target rules for round 1 & 2
rule round1:
	# Produces the updated ref fasta
	input:
		expand("test/{sample}_1.fq", sample=SAMPLES),
		f"{OUTDIR}/round1_index.ok",
		expand(f"{OUTDIR}/round1/alt_ref_for_chtc/{{sample}}_ref.fasta.tgz", sample=SAMPLES)

rule round2:
	input:
		f"{OUTDIR}/round2_index.ok",
		expand(f"{OUTDIR}/round2/vcf/{{sample}}_round2_SNPs.vcf.gz", sample=SAMPLES)

# Rule qc doesn't work
rule qc:
	input:
#		expand(f"{OUTDIR}/round{ROUND}/qc/multiqc/{{pre}}_multiqc_report.html", pre=config["prefix"])	
		expand(f"{OUTDIR}/round{ROUND}/qc/bamqc/{{sample}}_stats", sample=SAMPLES)

rule round2_index:
	input:
		expand(f"{OUTDIR}/round2_{{sample}}_index.ok", sample=SAMPLES),
		f"{OUTDIR}/round2_index.ok"

rule mem_stats:
	input:
		expand(f"{OUTDIR}/logs/bwa_mem/{{sample}}.stats", sample=SAMPLES)

rule test:
	input:
		fq1 = lambda wc: fq_from_sample(wc, 'FALSE', '1'),
		fq2 = lambda wc: fq_from_sample(wc, 'FALSE', '2'),
		index_ok = f"{OUTDIR}/round{ROUND}_index.ok"
	output:
		cut1 = "test/{sample}_1.fq",
		cut2 = "test/{sample}_2.fq"
	shell:
		"zcat {input.fq1} | awk '(NR<=50000)' > {output.cut1}; zcat {input.fq2} | awk '(NR<=50000)' > {output.cut2}"

rule bwa_aln_1:
	input:
		fq1 = lambda wc: fq_from_sample(wc, ISTESTING, '1'),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR) #,
#		index = f"{OUTDIR}/round{ROUND}_index.ok"
	output:
		sai1 = temp( f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_1.sai")
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}_map1.log"
	threads: 10
	shell:
		"bwa aln -t {threads} {input.REF} {input.fq1} 2> {log} > {output.sai1}"

rule bwa_aln_2:
	input:  
		fq2 = lambda wc: fq_from_sample(wc, ISTESTING, '2'),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR) #,
#		index = f"{OUTDIR}/round{ROUND}_index.ok"
	output: 
		sai2 = temp(f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_2.sai")
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}_map2.log"
	threads: 10
	shell:  
		"bwa aln -t {threads} {input.REF} {input.fq2} 2> {log} > {output.sai2}"

rule bwa_sampe:
	input:
		sai1 = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_1.sai",
		sai2 = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_2.sai",
		fq1 = lambda wc: fq_from_sample(wc, ISTESTING, '1'),
		fq2 = lambda wc: fq_from_sample(wc, ISTESTING, '2'),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output: 
		temp(f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}.bam")
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
		fq1 = lambda wc: fq_from_sample(wc, ISTESTING, '1'),
		fq2 = lambda wc: fq_from_sample(wc, ISTESTING, '2'),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/bwa_mem/{{sample}}.bam" 
	params: REF = config["genome"]
	log:  f"{OUTDIR}/logs/bwa_mem/{{sample}}.log"
	threads: 8
	shell:
		"bwa mem -M -t {threads} {input.REF} {input.fq1} {input.fq2} 2> {log} | samtools view -bS - > {output}"

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
		temp(f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_sort.bam")
	shell:
		"samtools sort {input} > {output}"

# Jeremy filters out unmapped reads next with Picard CleanSam.jar
#
#

rule mark_dups:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/RG_{{sample}}.bam"
	output:
		bam = temp(f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam"),
		metrics = f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}_dups.txt"
	conda:  "envs/picard.yaml"
	resources:
		mem_Gb = 3
	benchmark:
		f"{OUTDIR}/benchmarks/round{ROUND}/{{sample}}.markdup.benchmark.txt"
	shell:
		"picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"	

# Function RG_from_sample provides RG info from sample table
rule add_RG:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_sort.bam"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/RG_{{sample}}.bam")
	params: RG= lambda wildcards: RG_from_sample(wildcards)
	shell:
		"samtools addreplacerg -r '{params.RG}' -o {output} {input}"

rule dup_bam_bai:
        input:
                f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam"
        output:
                temp(f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam.bai")
        shell:
                "samtools index {input}"

#IDENTIFIES INTERVAL TO BE REALIGNED AROUND INDELS
rule indel_target:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam",
		bai = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam.bai",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.intervals")
	params: 
		GATKjar = config["gatk.jar"]
	conda: "envs/java8.yaml" # GATK needs java 8 (11 default on marula)
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/RTC/{{sample}}.log"		
	benchmark:
		f"{OUTDIR}/benchmarks/round{ROUND}/{{sample}}.RTC.benchmark.txt"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T RealignerTargetCreator \
			-R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule indel_realign:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam",
		bai = f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam.bai",
		intervals = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.intervals",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam"
	params: 
		GATKjar = config["gatk.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/indel_realigner/{{sample}}.log"
	conda: "envs/java8.yaml"
	resources:
		mem_Gb = 2
	benchmark:
		f"{OUTDIR}/benchmarks/round{ROUND}/{{sample}}.indelrealign.benchmark.txt"
	shell:
		"""
		java -Xmx8g -jar {params.GATKjar} -T IndelRealigner \
			-targetIntervals {input.intervals} \
			-R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

def get_UnifGen_params(ROUND):
        if ROUND == 1:
                return ["-ploidy 1"]
        elif ROUND == 2:
                return "".join(["-ploidy ", config["ploidy"], " -out_mode EMIT_ALL_SITES"])
                

rule gt_snps:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_round{ROUND}_SNPs.vcf"
	params: 
		GATKjar = config["gatk.jar"],
		static = "-mbq 10 -stand_call_conf 31 -stand_emit_conf 31",
		round = get_UnifGen_params(ROUND)
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/gt_snps/{{sample}}.log"
	conda: "envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.GATKjar} -T UnifiedGenotyper \
			{params.static} {params.round} -R {input.REF} -I {input.bam} -o {output} 2> {log}
		"""

rule gz_vcf:
	input:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_round{ROUND}_SNPs.vcf"
	output:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_round{ROUND}_SNPs.vcf.gz"
	shell:
		"""
		gzip {input}
		"""	

rule gt_indels:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_round{ROUND}_INDELs.vcf"
	params: 
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
# Jeremy's perl script filters snps supported by <75% (not <=75%!) reads takes all files as
# 	input, must be run in directory, get many div by 0 errors?
	input:
		f"{OUTDIR}/round1/vcf/{{sample}}_round1_SNPs.vcf"
	output:
		f"{OUTDIR}/round1/vcf/{{sample}}_round1_SNPs_filtered.vcf"
	conda: "envs/snpsift.yaml"
	shell:
		"""
		SnpSift filter "( GEN[0].AD[1]/GEN[0].DP > 0.75)" {input} > {output}
		"""

rule add_snps_alt_ref:
	input:
		vcf = f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_round{ROUND}_SNPs_filtered.vcf",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		temp(f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.fasta")
	params: 
		GATKjar = config["gatk.jar"]
	conda: "envs/java8.yaml"
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/alt_ref_snps/{{sample}}.log"
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
		temp(f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_SNPs_ref.dict")	
	params: 
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
		indel = f"{OUTDIR}/round{ROUND}/vcf/{{sample}}_round{ROUND}_INDELs.vcf"
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

rule tgz_alt_ref:
	input:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta"
	output:
		f"{OUTDIR}/round{ROUND}/alt_ref_for_chtc/{{sample}}_ref.fasta.tgz"
	params:
		dir = f"{OUTDIR}/round{ROUND}/alt_ref"
	shell:
		"""
		tar -C {params.dir} -czf {output} $( basename {input} )
		"""
