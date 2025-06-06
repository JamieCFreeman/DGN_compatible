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
CONTAINER_PATH = config["container_path"]

container: CONTAINER_PATH

# Read sample table
samples_table = pd.read_table(config["sample_table"], dtype=str).set_index("sample", drop=False)

# Get sample wildcards as a list
SAMPLES= samples_table['sample'].values.tolist()
# CHR list needed for fas1k wildcard
CHR = ["Chr2L", "Chr2R", "Chr3L", "Chr3R", "Chr4", "ChrX", "mtDNA", "Yhet"]

# Constrain sample wildcards to those in sample table
wildcard_constraints:
	sample="|".join(samples_table['sample']),
	unit="|".join(samples_table['unit'])

def fq_from_sample(wildcards, ISTESTING, read=None):
	# Sometimes we provide read from wildcards, sometimes we want to specify (eg when merging)
	if read is not None:
		R = read
	elif read is None:
		R = wildcards.read 
	if ISTESTING == 'TRUE':
		words = ["test/" , wildcards.sample, "_", wildcards.unit, '_R', R, ".fq"]
		return "".join(words)
	elif ISTESTING == 'FALSE':
		now = samples_table[ (samples_table["sample"] == wildcards.sample) & (samples_table["unit"] == wildcards.unit) ]["fq" + R ].iloc[0]
		return now

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
	now = samples_table[ (samples_table["sample"] == wildcards.sample) & (samples_table["unit"] == wildcards.unit) ]["RG"].iloc[0]
	return now.strip('\'')

def get_sa_units(wildcards, ROUND, OUTDIR):
	now =  list( samples_table[ samples_table["sample"] == wildcards.sample]["unit"] )
	bams = [ OUTDIR + '/round' + str(ROUND) + '/stampy/RG_' + wildcards.sample + '_' +  x + '.bam' for x in now ]
	return bams

include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/genome_indexing.smk",
include: "rules/fas1k.smk"

# Desired output depends on whether we're running round 1 or 2
index_check = f"{OUTDIR}/round{ROUND}_index.ok"
test_fq = expand("test/{sample}_1.fq", sample=SAMPLES)
alt_ref = expand(f"{OUTDIR}/round1/alt_ref_for_chtc/{{sample}}_ref.fasta.tgz", sample=SAMPLES)
allsites_vcf = expand(f"{OUTDIR}/round2/vcf/{{sample}}_round2_SNPs.vcf.gz", sample=SAMPLES)

rule_all_input_list=[index_check]
round1_input_list=[test_fq, alt_ref]
round1_input_list=[alt_ref]
round2_input_list=[allsites_vcf]

if ROUND == 1:
	rule_all_input_list.extend(round1_input_list)
	print("Doing round 1 mapping")
	print(rule_all_input_list)
elif ROUND == 2:
	rule_all_input_list.extend(round2_input_list)
	print("Doing round 2 mapping")

#rule_all_input_list=['FR_N/round1_index.ok', 'test/6Jul21-1_EF2N_S102_1.fq', 'test/6Jul21-2_EF6N_S103_1.fq', 'FR_N/round1/alt_ref_for_chtc/6Jul21-1_EF2N_S102_ref.fasta.tgz', 'FR_N/round1/alt_ref_for_chtc/6Jul21-2_EF6N_S103_ref.fasta.tgz']


rule all:
	input:
		rule_all_input_list

rule fas1k:
	input:
		expand(f"{OUTDIR}/round2/fas1k/{{sample}}_round2_{{chr}}_diploid.fas1k", sample=SAMPLES, chr=CHR)	

rule stats:
	input:
		expand(f"{OUTDIR}/round{ROUND}/het_stats/{{sample}}_round{ROUND}_het.txt",sample=SAMPLES),
		f"{OUTDIR}/round{ROUND}/inv_calls/inv_snp_freq_table.tsv"

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
		fq = lambda wc: fq_from_sample(wc, 'FALSE'),
		index_ok = f"{OUTDIR}/round{ROUND}_index.ok"
	output:
		cut = "test/{sample}_{unit}_R{read}.fq"
	shell:
		"zcat {input.fq} | awk '(NR<=100000)' > {output.cut}"

rule bwa_aln:
	input:
		fq = lambda wc: fq_from_sample(wc, ISTESTING),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR),
		index = ancient( f"{OUTDIR}/round{ROUND}_index.ok" )
	output:
		sai = temp( f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_{{unit}}_R{{read}}.sai")
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}_{{unit}}_map{{read}}.log"
	threads: 10
	shell:
		"bwa aln -t {threads} {input.REF} {input.fq} 2> {log} > {output.sai}"

rule bwa_sampe:
	input:
		sai1 = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_{{unit}}_R1.sai",
		sai2 = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_{{unit}}_R2.sai",
		fq1 = lambda wc: fq_from_sample(wc, ISTESTING, '1'),
		fq2 = lambda wc: fq_from_sample(wc, ISTESTING, '2'),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output: 
		temp(f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_{{unit}}.bam")
	log: f"{OUTDIR}/round{ROUND}/logs/bwa_sampe/{{sample}}_{{unit}}_sampe.log"
	shell:
		"bwa sampe -P {input.REF} {input.sai1} {input.sai2} {input.fq1} {input.fq2} 2> {log} | samtools view -bS - > {output}"

rule aln_flagstat:
	input:
		f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_{{unit}}.bam"
	output:
		f"{OUTDIR}/round{ROUND}/logs/bwa_aln/{{sample}}_{{unit}}.stats"
	shell:
		"samtools flagstat {input} > {output}"
rule bwa_mem:
	input:
		fq1 = lambda wc: fq_from_sample(wc, ISTESTING, '1'),
		fq2 = lambda wc: fq_from_sample(wc, ISTESTING, '2'),
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/bwa_mem/{{sample}}_{{unit}}.bam" 
	params: REF = config["genome"]
	log:  f"{OUTDIR}/logs/bwa_mem/{{sample}}_{{unit}}.log"
	threads: 8
	shell:
		"bwa mem -M -t {threads} {input.REF} {input.fq1} {input.fq2} 2> {log} | samtools view -bS - > {output}"

rule flagstat:
	input:
		f"{OUTDIR}/bwa_mem/{{sample}}_{{unit}}.bam"
	output:
		f"{OUTDIR}/logs/bwa_mem/{{sample}}_{{unit}}.stats"
	shell:
		"samtools flagstat {input} > {output}"

rule stampy_map:
	input:
		bam = f"{OUTDIR}/round{ROUND}/bwa_aln/{{sample}}_{{unit}}.bam",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round{ROUND}/stampy/{{sample}}_{{unit}}.sam"
#		temp(f"{OUTDIR}/round{ROUND}/stampy/{{sample}}_{{unit}}.sam")
	params: stampy = "/opt/bioscript/stampy/stampy.py"
	log: f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}_{{unit}}.log"
	conda: "envs/py2.yaml"
	shell:
		"python {params.stampy} -g {input.REF} -h {input.REF} --bamkeepgoodreads -M {input.bam} -o {output} 2> {log} "
# "python2.6 " . $stampy . "stampy.py -g " . $reference . " -h " . $reference . " --bamkeepgoodreads -M " . $FastqFile[$i] . ".bam -o " . $FastqFile[$i] . "_remapped.sam"; 
# #Stampy MAPPING STEP. THIS IS A RELATIVELY LONG STEP (AS LONG AS 15 HOURS ON SOME OF THE HIGHEST COVERAGE DPGP2 GENOMES)

rule sam2bam:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/{{sample}}_{{unit}}.sam"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/{{sample}}_{{unit}}.bam")
	shell:
		"samtools view -bS {input} > {output}"

rule stampy_flagstat:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/{{sample}}_{{unit}}.bam"
	output:
		f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}_{{unit}}.stats"
	shell:
		"samtools flagstat {input} > {output}"

rule qfilter_bam:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/{{sample}}_{{unit}}.bam",
		stats = f"{OUTDIR}/round{ROUND}/logs/stampy/{{sample}}_{{unit}}.stats"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_{{unit}}.bam")
	shell:
		"samtools view -q 20 -h {input.bam} > {output}"

rule sort_bam:
	input:
		f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_{{unit}}.bam"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_{{unit}}_sort.bam")
	conda:  "envs/samtools.yaml"
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
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/mark_dup/{{sample}}.log"		
	benchmark:
		f"{OUTDIR}/benchmarks/round{ROUND}/{{sample}}.markdup.benchmark.txt"
	shell:
		"picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} >& {log}"	

# Function RG_from_sample provides RG info from sample table
rule add_RG:
	input:
		bam = f"{OUTDIR}/round{ROUND}/stampy/qfilter_{{sample}}_{{unit}}_sort.bam"
	output:
		temp(f"{OUTDIR}/round{ROUND}/stampy/RG_{{sample}}_{{unit}}.bam")
	params: RG= lambda wildcards: RG_from_sample(wildcards)
	conda:  "envs/samtools.yaml"
	shell:
		"samtools addreplacerg -r '{params.RG}' -o {output} {input}"

rule merge_sample_bams:
	input:
		lambda wc: get_sa_units(wc, ROUND, OUTDIR)
	output:
		bam = f"{OUTDIR}/round{ROUND}/stampy/RG_{{sample}}.bam"
	conda:  "envs/samtools.yaml"
	shell:
		"samtools merge -o {output.bam} {input}"
	

rule dup_bam_bai:
        input:
                f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam"
        output:
                temp(f"{OUTDIR}/round{ROUND}/stampy/mark_dup/{{sample}}.bam.bai")
	conda:  "envs/samtools.yaml"
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
	retries: 2 #sometimes fails on marula, with success on rerun
	resources:
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

