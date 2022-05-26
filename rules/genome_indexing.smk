
rule round1_index_ok:
# Check for presence of all index files for Dmel ref genome
	input:
		fai = lambda wc: get_ref_idx(wc, 1, OUTDIR, ".fai"),
		stidx = lambda wc: get_ref_idx(wc, 1, OUTDIR, ".stidx"),
		sthash = lambda wc: get_ref_idx(wc, 1, OUTDIR, ".sthash"),
		bwt = lambda wc: get_ref_idx(wc, 1, OUTDIR, ".bwt"),
		dict = lambda wc: get_ref_idx(wc, 1, OUTDIR, ".dict")
	output:
		f"{OUTDIR}/round1_index.ok"
	shell:
		"touch {output}"	
		
rule round2_index_ok:
# Check for presence of all sample-specific genome index files
	input:
		fai = f"{OUTDIR}/round1/alt_ref/{{sample}}_ref.fasta.fai",
		stidx =  f"{OUTDIR}/round1/alt_ref/{{sample}}_ref.fasta.stidx",
		sthash = f"{OUTDIR}/round1/alt_ref/{{sample}}_ref.fasta.sthash",
		bwt = f"{OUTDIR}/round1/alt_ref/{{sample}}_ref.fasta.bwt",
		dict = f"{OUTDIR}/round1/alt_ref/{{sample}}_ref.dict"
	output:
		temp(f"{OUTDIR}/round2_{{sample}}_index.ok")
	shell:
		"touch {output}"	

rule all_round2_index_ok:
# Input all sample-specific ref, output 1 file to indicate all ref good
	input:
		expand(f"{OUTDIR}/round2_{{sample}}_index.ok", sample=SAMPLES)
	output:
		f"{OUTDIR}/round2_index.ok"
	shell:
		"touch {output}"
		
rule bwa_index:
	input:
		f"{OUTDIR}/round{ROUND}/alt_ref/{{sample}}_ref.fasta",
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR),
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
	conda: "../envs/py2.yaml"
	shell:
		"""
		TEMP_FILE=`echo {output.stidx} | sed 's ^.*alt_ref/  '| sed 's/.stidx//'`; echo $TEMP_FILE #stampy writes stidx file to working dir
		python {params.stampy} -G {wildcards.sample}_ref.fasta {input} # build genome file
		mv $TEMP_FILE.stidx {output.stidx} # move stidx to alt_ref folder
		python {params.stampy} -g {input} -H {wildcards.sample}_ref.fasta # build hash table
		mv $TEMP_FILE.sthash {output.sthash}
		"""


rule gatk_seq_dict:
	input:
		REF = lambda wc: get_ref_fa(wc, ROUND, OUTDIR)
	output:
		f"{OUTDIR}/round1/alt_ref/{{sample}}_ref.dict"	
	params: 
		picardjar = config["picard.jar"]
	log: f"{OUTDIR}/round{ROUND}/logs/gatk/alt_ref_index/{{sample}}.log"
	conda: "../envs/java8.yaml"
	shell:
		"""
		java -Xmx4g -jar {params.picardjar}/CreateSequenceDictionary.jar \
			REFERENCE={input.REF} OUTPUT={output} 2> {log}
		"""
