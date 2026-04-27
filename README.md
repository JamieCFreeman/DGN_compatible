# DGN_compatible

Pipeline to analyze DNA sequencing data with legacy compatibility for the Pool lab's Drosophila Genome Nexus data.

> **Note:** You likely only want to use this pipeline if you have data you want to compare or use in conjunction with the DGN dataset.

## Getting Started

### Step 1: Prepare Sample Table

Create a sample table in the format of `units.tsv` (tab-delimited table).

**Required columns:**
- `sample`: Sample name
- `unit`: I typically use sample number string from raw data file name
- `platform`: "ILLUMINA" or other sequencing platform
- `fq1`: Path to forward reads file
- `fq2`: Path to reverse reads file
- `RG`: Read group string

**Example RG format:**
```
@RG\tID:13Jan22-1_EF5N_S73\tSM:13Jan22-1_EF5N_S73\tPU:none\tLB:13Jan22-1_EF5N_S73\tPL:ILLUMINA
```

### Step 2: Configure Your Project

Edit `config.yaml` for your project settings.

**Recommendation:** First run with the `ISTESTING` variable set to `TRUE`. This will subset reads so the entire workflow runs quickly, allowing you to verify it runs cleanly before starting the full analysis.

### Step 3: Run the Pipeline with Snakemake

**Requirements:**
- Snakemake must be installed in your environment (on Marula, use the `snakemake` conda environment)
- Run with conda: `--use-conda --conda-frontend conda`

**Always start with a dry run:**
```bash
snakemake -n
```

This verifies the job list matches your expectations.

**Command examples:**

For Round 1 or 2 mapping (all rule):
```bash
snakemake -c48 all --use-conda --conda-frontend conda --resources mem_Gb=56
```

For Round 2 indexes:
```bash
snakemake -c48 round2_index --use-conda --conda-frontend conda --resources mem_Gb=56
```

For fas1k files:
```bash
snakemake fas1k --use-conda --conda-frontend conda -c48
```

For statistics (heterozygosity and inversion SNP genotyping):
```bash
snakemake stats --use-conda --conda-frontend conda -c48 --resources io=30
```

**Resource flags:**
- `-c`: Number of cores to use
- `--resources mem_Gb`: Memory cap (useful for memory-intensive steps like MarkDuplicates at 3 Gb and IndelRealigner at 2 Gb)

> **Note:** Snakemake doesn't actively monitor resources, but these tags limit how many resource-intensive jobs run simultaneously.

## Pipeline Workflow

### Two-Round Mapping Strategy

Rounds 1 and 2 are not connected in this pipeline because mapping is typically run on UW's CHTC cluster (Stampy is computationally intensive and benefits from high parallelization).

See the [CHTC pipeline](https://github.com/JamieCFreeman/chtc_bwa_stampy) for parallel mapping.

**Workflow after CHTC mapping:**

1. Obtain merged `.bam` files for every sample-unit from CHTC
2. Place `.bam` files where Snakemake expects them
3. Ensure the timestamps of new `.bam` files are more recent than their input files (Snakemake won't rerun mapping otherwise)
4. **Always verify** with a dry run (`-n` flag) before executing

**After Round 1 mapping:**
- If transferring `.bam` files from CHTC, first create round 2 indexes (they should be timestamped):
```bash
snakemake -c48 round2_index --use-conda --conda-frontend conda --resources mem_Gb=56
```

## Methods

Reference mapping and variant calling approach from Lack et al. (2015) *Genetics* [10.1534/genetics.115.174664](https://doi.org/10.1534/genetics.115.174664):

> In brief, we first mapped short read data to the *D. melanogaster* reference genome (release 5.57; http://flybase.org) using BWA v0.5.9 (Li and Durbin 2010) using default settings, followed by mapping of all unmapped reads using Stampy v1.0.20 (Lunter and Goodson 2010). This approach combines the rapid but strict BWA algorithm to first map the relatively “easy-to-align” reads with the more sensitive but computationally intensive Stampy algorithm, which more effectively and accurately aligns the relatively divergent reads (Lunter and Goodson 2010). All reads with mapping quality scores <20 were excluded. Optical duplicates were then removed using Picard v1.79 (http://picard.sourceforge.net/), and assemblies were improved around indels using the GATK v3.2 Indel Realigner (McKenna et al. 2010; Depristo et al. 2011). The Unified Genotyper (Depristo et al. 2011) was then used to call indels and SNPs for each individual genome. Among the indel calling criteria, >50% of the reads at a given position had to support the existence of that indel, with a minimum of three reads containing the variant. For SNP calling in this first round, we required a minimum phred-scaled quality value of 31 and that >75% of the reads at a given position support the SNP. For the second round of assembly, the SNPs and indels called in the first round were introduced into the *D. melanogaster* reference, and this modified reference was then used for a second round of mapping. Following indel realignment, the Unified Genotyper was then used to call all sites in the modified reference genome. To generate reference-numbered consensus sequences, a custom perl script was used to shift all base coordinates back to those of the original *D. melanogaster* reference. Deletions and all sites within 3 bp of a called indel were coded as “N” (based on the error analysis described in the Results), while insertions do not appear in reference-numbered consensus sequences."

![DAG - Round 1](https://github.com/JamieCFreeman/DGN_compatible/blob/main/README_files/rulegraph.pdf?raw=true)

![DAG - Round 2](https://github.com/JamieCFreeman/DGN_compatible/blob/main/README_files/round2_dag.pdf?raw=true)
