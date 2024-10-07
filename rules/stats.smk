rule calc_fas1k_het:
	input:
		fas1k = f"{OUTDIR}/round{ROUND}/fas1k/{{sample}}_round{ROUND}_Chr2L_diploid.fas1k"
	output:
		f"{OUTDIR}/round{ROUND}/het_stats/{{sample}}_round{ROUND}_het.txt"
	params:
		dir = f"{OUTDIR}/round{ROUND}/het_stats"
	shell:
		"python /home/jamie/FAS1K_utils/run_het_stats.py {input.fas1k} {params.dir}"

rule call_inv_snps:
	input:
		fas1k = f"{OUTDIR}/round{ROUND}/fas1k/{{sample}}_round{ROUND}_Chr2L_diploid.fas1k"
	output:
		f"{OUTDIR}/round{ROUND}/inv_calls/{{sample}}_round{ROUND}_inv_SNP_counts.tsv"
	params:
		dir = f"{OUTDIR}/round{ROUND}/inv_calls"
	shell:
		"python /home/jamie/FAS1K_utils/run_inv_snp_freq.py --fas1k {input.fas1k} --out {output}"

rule proc_inv_snps:
	input:
		expand(f"{OUTDIR}/round{ROUND}/inv_calls/{{sample}}_round{ROUND}_inv_SNP_counts.tsv", sample=SAMPLES)
	output:
		f"{OUTDIR}/round{ROUND}/inv_calls/inv_snp_freq_table.tsv",
		f"{OUTDIR}/round{ROUND}/inv_calls/inv_snp_freq_counts_ge0.05.tsv"
	params:
		dir = f"{OUTDIR}/round{ROUND}/inv_calls"
	shell:
		"python /home/jamie/FAS1K_utils/process_inv_freq.py {params.dir}"
