process BWA_MEM {

	tag "$sample_id"

	input:
	tuple val(sample_id), path(fastq1), path(fastq2)
	path ref
	path ref_idx
	output:
	tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

	script:
	"""
	set -o pipefail

	# Prepare RG info
	LB="${sample_id}_lib1"
	PL="ILLUMINA"

	# Flowcell:lane of every read (Illumina read name fields 3 and 4). Celemics FASTQs can hold reads
	# from several lanes/flowcells (top-up runs), and each lane needs its own read group for
	# MarkDuplicates and BQSR.
	zcat ${fastq1} | awk 'NR%4==1' | cut -d: -f3,4 --output-delimiter="." | LC_ALL=C sort -u > lanes.txt
	echo "${sample_id}: \$(wc -l < lanes.txt) read group(s): \$(tr '\\n' ' ' < lanes.txt)"

	if [ "\$(wc -l < lanes.txt)" -eq 1 ]; then
		FLOWCELL_LANE=\$(cat lanes.txt)

		# Single lane (unchanged): run alignment while adding read groups
		bwa mem -t ${task.cpus} \\
		-K 10000000 -Y \\
		-R "@RG\\tID:\${FLOWCELL_LANE}\\tSM:${sample_id}\\tLB:\${LB}\\tPL:\${PL}\\tPU:\${FLOWCELL_LANE}" \\
		${ref} \\
		${fastq1} ${fastq2} \\
		| samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
	else
		# Several lanes: one @RG per flowcell.lane (same SM and LB); each read is tagged with the
		# read group of its own lane, taken from its read name. Same result as aligning each lane
		# separately with bwa mem -R and merging.
		awk -v sm="${sample_id}" -v lb="\${LB}" -v pl="\${PL}" \\
			'{printf "@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:%s\\tPU:%s\\n", \$1, sm, lb, pl, \$1}' lanes.txt > rg_header.txt

		bwa mem -t ${task.cpus} \\
		-K 10000000 -Y \\
		${ref} \\
		${fastq1} ${fastq2} \\
		| LC_ALL=C awk -v hdr=rg_header.txt '
			BEGIN { OFS = "\\t" }
			/^@/ { if (\$1 == "@PG" && !done) { while ((getline l < hdr) > 0) print l; done = 1 } print; next }
			{ if (!done) { while ((getline l < hdr) > 0) print l; done = 1 }
			  split(\$1, a, ":"); print \$0, "RG:Z:" a[3] "." a[4] }' \\
		| samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -

		# the BAM header must declare one @RG per lane found
		samtools view -H ${sample_id}.sorted.bam | grep -c '^@RG' | grep -qx "\$(wc -l < lanes.txt)"
	fi

	# Index sorted bam
	samtools index ${sample_id}.sorted.bam



	"""




}
