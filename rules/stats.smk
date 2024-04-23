rule calc_fas1k_het:
	input:
		fas1k = f"{OUTDIR}/round{ROUND}/fas1k/{{sample}}_round{ROUND}_Chr2L_diploid.fas1k"
	output:
		f"{OUTDIR}/round{ROUND}/het_stats/{{sample}}_round{ROUND}_het.txt"
	params:
		dir = f"{OUTDIR}/round{ROUND}/het_stats"
	shell:
		"python /home/jamie/FAS1K_utils/run_het_stats.py {input.fas1k} {params.dir}"

