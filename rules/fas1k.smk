rule prep_indel_shift:
	input:
		SNPs = f"{OUTDIR}/round2/vcf/{{sample}}_round2_SNPs.vcf.gz",
		indel = f"{OUTDIR}/round1/vcf/{{sample}}_round1_INDELs.vcf"
	output:
		SNPs = temp( f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_sites.vcf.gz"),
		indel = temp( f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_INDELS.vcf")
	shell:
		"""
		cp {input.SNPs} {output.SNPs}; cp {input.indel} {output.indel}
		"""

rule run_indel_shift:
	input:
		sites = f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_sites.vcf.gz",
		indel = f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_INDELS.vcf"
	output:
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_shifted.vcf.gz"
#		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2IndelCoordinates.out"
	params:
		dir = f"{OUTDIR}/round2/shifted_vcf/"
	shell:
		"""
		cp scripts/IndelShift_parallel_JCF.pl {params.dir}; cd {params.dir};
		perl IndelShift_parallel_JCF.pl $( echo {input.sites} | sed 's ^.*/  ' )
		"""

def get_fas_list(wildcards, OUTDIR):
		s = [ OUTIDIR + "round2/shifted_vcf/" + wildcards.sample + "round2" + x + "_diploid.fas" for x in CHR ]
		return s

rule fas:
	input:
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_shifted.vcf.gz"
	output:
#		lambda wc: get_fas_list(wc, OUTDIR)
#		expand(f"{OUTDIR}/round2/fas1k/{{sample}}_round2_{chr}_diploid.fas1k", chr=CHR, allow_missing=True)
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_Chr2L_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_Chr2R_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_Chr3L_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_Chr3R_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_Chr4_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_ChrX_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_mtDNA_diploid.fas",
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_Yhet_diploid.fas"
	params:
		dir = f"{OUTDIR}/round2/shifted_vcf/"
	resources: io=1 #zipping and unzipping is i/o intensive, want to limit (eg --resources io=30)
	shell:
		"""
		perl scripts/VCF_to_Seq_diploid_ambiguities.pl {input}
		"""

rule indel_mask_bed:
	input:
		indel = f"{OUTDIR}/round1/vcf/{{sample}}_round1_INDELs.vcf"
	output:
		temp( f"{OUTDIR}/round1/vcf/{{sample}}_round1_INDELs.bed" )
	params:
		indel_mask_buffer = 3 # how many basepairs around indels to mask?
	conda: "../envs/bed.yaml"
	shell:
		"""
		bcftools query -f '%CHROM\t%POS0\t%END\n' {input.indel} | \
		bedops --range {params.indel_mask_buffer} -u - > {output}
		"""

# Need fasta header line for bedtools
rule add_fas_header:
	input:
		fas = f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid.fas"	
	output:
		temp( f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid.fas.tmp" )
	run:
		import sys
		sys.path.append(r'/home/jamie/FAS1K_utils')
		import fas1k_utils as f1k
		arm = f1k.arm_to_int(wildcards.chr)
		nt  = f1k.read_in(input.fas)
		with open(output[0], 'w') as f:
			f.write( '>' + str(arm) + '\n')
			f.write(nt[0] + '\n')

rule mask_around_indels:
	input:
		bed = f"{OUTDIR}/round1/vcf/{{sample}}_round1_INDELs.bed",
		fas = f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid.fas.tmp"	
	output:
		tmp  = temp( f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid_indelMASK.fas.tmp" ),
		mask = temp( f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid_indelMASK.fas" )	
	shell:
		"""
		bedtools maskfasta -fi {input.fas} -bed {input.bed} -fo {output.tmp}
		tail -n1 {output.tmp} > {output.mask} # remove fasta header
		"""

rule wrap_fas1k:
	input:
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid_indelMASK.fas"
	output:
		temp( f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid.fas1k" )
	shell:
		"""
		perl scripts/wrap_fasta_JCF.pl {input}
		mv {input}1k {output} # wrap script just globs filename
		"""

rule mv_fas1k:
	input:
		f"{OUTDIR}/round2/shifted_vcf/{{sample}}_round2_{{chr}}_diploid.fas1k"
	output:
		f"{OUTDIR}/round2/fas1k/{{sample}}_round2_{{chr}}_diploid.fas1k"
	params:
		dir = f"{OUTDIR}/round2/fas1k/"	
	shell:
		"""
		mv {input} {params.dir}
		"""
