process MARK_DUPLICATES {
	
	tag "$sample_id"
	
	input:
	tuple val(sample_id), path(bam), path(bai)

	output:
	tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

	script:
	"""
	
	gatk MarkDuplicates -I ${bam} \\
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
	-O ${sample_id}.dedup.bam \\
	-M ${sample_id}.metrics.txt

	samtools index ${sample_id}.dedup.bam
	
	
	"""

}

process BQSR {
	
	tag "$sample_id"
	debug true
	publishDir ("${params.outdir}/bqsr"), mode: 'symlink'

	input:
	tuple val(sample_id), path(bam), path(bai)
	path ref
	path ref_index
	path ref_dict
	path dat
	path map
	path kg
	path mills
	path dbsnp

	output:
	tuple val(sample_id), path("${sample_id}.recal.bam"), path("${sample_id}.recal.bam.bai"), path("fm.txt"), emit: bam
	path "${sample_id}.selfSM" 
	path "${sample_id}.Ancestry"
	path "${sample_id}.crosscheck_metrics"
	path "${sample_id}.recal.table"
	path "${sample_id}.after.recal.table"
	

	script:
	"""
	
	VerifyBamID \\
	--SVDPrefix 1000g.phase3.100k.b38.vcf.gz.dat \\
	--Reference ${ref} \\
	--BamFile ${bam} \\
	--DisableSanityCheck \\
	--NumPC 4 \\
	--Output ${sample_id}
	
	awk 'NR==2 {print \$7}' ${sample_id}.selfSM > fm.txt
	
	fm=\$(awk 'NR==2 {print \$7}' ${sample_id}.selfSM)
	echo "FREEMIX"
	echo "##########"
	echo "FM value: ###  \${fm}"

	gatk CrosscheckFingerprints \\
	--INPUT ${bam} \\
	--HAPLOTYPE_MAP ${map} \\
	--LOD_THRESHOLD -20.0 \\
	--CROSSCHECK_BY READGROUP \\
	--OUTPUT ${sample_id}.crosscheck_metrics
		
	gatk BaseRecalibrator -I ${bam} -R ${ref} \\
	--known-sites ${mills[0]} \\
	--known-sites ${dbsnp[0]} \\
	--known-sites ${kg[0]} \\
	-O ${sample_id}.recal.table

	gatk ApplyBQSR -I ${bam} -R ${ref} \\
	--bqsr-recal-file ${sample_id}.recal.table \\
	--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \\
	-O ${sample_id}.recal.bam

	samtools index ${sample_id}.recal.bam

	## Generate the covariate table after applying BQSR
	gatk BaseRecalibrator -I ${sample_id}.recal.bam -R ${ref} \\
	--known-sites ${mills[0]} \\
	--known-sites ${dbsnp[0]} \\
	--known-sites ${kg[0]} \\
	-O ${sample_id}.after.recal.table


	"""



}
