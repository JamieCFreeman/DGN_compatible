rule bamqc:
	input:
		bai = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bai",
		bam = f"{OUTDIR}/round{ROUND}/stampy/indel_realign/{{sample}}.bam"
	output:
		dir = directory(f"{OUTDIR}/round{ROUND}/qc/bamqc/{{sample}}_stats"),
		file = f"{OUTDIR}/round{ROUND}/qc/bamqc/{{sample}}_stats/qualimapReport.html"
	threads: 6
	log:
		f"{OUTDIR}/round{ROUND}/logs/bamqc/{{sample}}.log"
	conda: "../envs/qualimap.yaml"
	shell:
		"qualimap bamqc -bam {input.bam} -nt {threads} -outdir {output.dir} --collect-overlap-pairs --java-mem-size=2G &> {log} "

#rule multiqc:
#	input:
#		[f"{OUTDIR}/round{ROUND}/qc/bamqc/{sample}_stats/qualimapReport.html" for sample in SAMPLES]
#	output:
#		expand(f"{OUTDIR}/round{ROUND}/qc/multiqc/{{pre}}_multiqc_report.html",
#			pre=config["prefix"])
#	params:
#		indir = f"{OUTDIR}/round{ROUND}",
#		name = config["prefix"],
#		outdir =  f"{OUTDIR}/round{ROUND}/qc/multiqc"
#	conda: "envs/multiqc.yaml"
#	shell:
#		"multiqc {params.indir} --outdir {params.outdir} --title {params.name}"
                                                                                           
