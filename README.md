# tgb-nf-germline-var-call

Nextflow (DSL2) GATK germline short-variant pipeline for whole-exome sequencing (Celemics WES, GRCh38),
written by Savindu (original repository: savindu93/tgb-nf-germline-var-call) and maintained by
Torrington Genomics and Bioinformatics.

Per sample: fastp → BWA-MEM (read groups per flowcell-lane) → samtools sort → GATK MarkDuplicates →
VerifyBamID2 contamination + CrosscheckFingerprints → BaseRecalibrator/ApplyBQSR (+ post-BQSR table) →
HaplotypeCaller `-ERC GVCF` on the Celemics targets (±100 bp) with the VerifyBamID2 FREEMIX as
`-contamination`. Joint genotyping (GenomicsDBImport + GenotypeGVCFs) is run separately on the cohort.

From v3.1-qc each run also writes alignment QC per sample (CollectHsMetrics on the targets, samtools stats)
and one QC table and report for the run: `qc/batch_qc.tsv`, `qc/batch_qc_report.html`, `qc/multiqc_report.html`.

```
nextflow run . --samplesheet samplesheet.csv --outdir results
```
Sample sheet columns: `sample_id,fastq_1,fastq_2`, plus an optional `fastq_s` (v3.2): unpaired reads, e.g. the
singletons left by a BAM→FASTQ conversion, trimmed single-end and aligned into the same BAM. Reference, truth sets, targets and interval paths are
set in `nextflow.config` (`params`).

Every version used for production samples is tagged; see [CHANGELOG.md](CHANGELOG.md) and
[ENVIRONMENT.md](ENVIRONMENT.md).
