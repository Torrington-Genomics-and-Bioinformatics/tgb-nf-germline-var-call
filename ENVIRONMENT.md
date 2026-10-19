# Software environment

| Tool | v0-imac (iMac) | v1–v3 (GCP VM / batch workers) |
|---|---|---|
| GATK | 4.6.2.0 | 4.6.2.0 |
| BWA | 0.7.19 | 0.7.19 |
| samtools / htslib | 1.21 / 1.21 | 1.24 / 1.24 |
| fastp | 1.1.0 | 1.3.7 |
| VerifyBamID2 | built from source | 2.0.3 (bioconda) |
| Nextflow | 26.04.6 | 26.04.6 |
| OpenJDK | 17.0.18 | 17.0.14 |

- `env/conda_env.imac.yml`: full pinned export of the iMac conda environment.
- `env/ngs-tools.vm.yml`: the VM environment definition (unpinned).
- `env/ngs-tools.vm-batch.lock.tsv`: exact installed package, version and build of the environment used
  on the VM and the batch workers (the workers were cloned from the VM; read from the worker image).

Reference: GRCh38 primary assembly (GENCODE `GRCh38.primary_assembly.genome.fa`). Targets: Celemics WES
V2 (CLMXPD037-V2) lifted to GRCh38, 216,330 targets, 36.8 Mb; HaplotypeCaller `-L` targets `-ip 100`.
Known sites: Mills and 1000G gold standard indels, dbSNP 138, 1000G phase 1 high-confidence SNPs.
