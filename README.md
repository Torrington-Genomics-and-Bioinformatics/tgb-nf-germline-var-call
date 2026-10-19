# tgb-nf-germline-var-call

Nextflow (DSL2) GATK germline short-variant pipeline for whole-exome sequencing (Celemics WES, GRCh38),
written by Savindu (original repository: savindu93/tgb-nf-germline-var-call) and maintained by
Torrington Genomics and Bioinformatics.

Per sample: fastp → BWA-MEM (read groups per flowcell-lane) → samtools sort → GATK MarkDuplicates →
VerifyBamID2 contamination + CrosscheckFingerprints → BaseRecalibrator/ApplyBQSR (+ post-BQSR table) →
HaplotypeCaller `-ERC GVCF` on the Celemics targets (±100 bp) with the VerifyBamID2 FREEMIX as
`-contamination`. Joint genotyping (GenomicsDBImport + GenotypeGVCFs) is run separately on the cohort.

```
nextflow run . --samplesheet samplesheet.csv --outdir results
```
Sample sheet columns: `sample_id,fastq_1,fastq_2`. Reference, truth sets, targets and interval paths are
set in `nextflow.config` (`params`).

Every version used for production samples is tagged; see [CHANGELOG.md](CHANGELOG.md) and
[ENVIRONMENT.md](ENVIRONMENT.md).
