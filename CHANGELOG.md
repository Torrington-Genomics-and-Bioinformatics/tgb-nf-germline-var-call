# Versions

Each tag is the exact code used for a set of production samples. Changes that can affect variant calls
are marked **(calls)**; the rest only change resources, file handling or which outputs are kept.

| Tag | Where it ran | Samples |
|---|---|---|
| `v0-imac` | lab iMac (macOS), external WD drive | ~99 |
| `v1-vm` | GCP VM `tgb-nf-germline-var-call` | 22 |
| `v2.0-batch` | GCP Spot batch workers | benchmark run + batch b01 |
| `v2.1-batch-qc` | GCP Spot batch workers | batches b02–b19 |
| `v3.0-multilane` | GCP Spot batch workers | batches b20–b49 |

## v0-imac
Original pipeline. macOS `gzcat`/`gcut`; VerifyBamID built from source; GenomicsDBImport creating one
workspace per interval shard; conda env in `env/conda_env.imac.yml` (GATK 4.6.2.0, bwa 0.7.19,
samtools 1.21, fastp 1.1.0).

## v1-vm
- Linux `zcat`/`cut` (same function).
- **(calls, small)** contamination from `verifybamid2` 2.0.3 (bioconda) instead of the source build;
  FREEMIX is passed to HaplotypeCaller `-contamination`.
- Known-sites VCFs selected by file name instead of the first element of the channel.
- More CPUs/memory; GenomicsDBImport updates an existing workspace.
- **(calls, small)** environment: fastp 1.3.7 and samtools/htslib 1.24 (see ENVIRONMENT.md).

## v2.0-batch
- GenomicsDB import removed from the per-batch run (cohort database built once, later).
- `bwa mem -t ${task.cpus}`: alignments unchanged, because `-K 10000000` makes BWA output independent
  of the thread count.
- GATK tools get `-Xmx` = task memory − 1 GB; retries on exit 137/140/247.
- Results copied (`publishDir mode: 'copy'`) instead of symlinked.

## v2.1-batch-qc
- fastp reports and MarkDuplicates metrics published to `<outdir>/qc`.

## v3.0-multilane
- **(calls, multi-lane samples only)** one read group per flowcell-lane: each read gets the read group
  of its own lane (ID = PU = flowcell.lane, same SM and LB). Single-lane samples run the unchanged
  command. BQSR then models each lane; CrosscheckFingerprints (BY READGROUP) checks the lanes agree.
