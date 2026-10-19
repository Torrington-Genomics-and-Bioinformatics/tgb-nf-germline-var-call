process MARK_DUPLICATES {
	
	tag "$sample_id"
	
	input:
	tuple val(sample_id), path(bam), path(bai)

	output:
	tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")

	script:
	"""
	
	gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" MarkDuplicates -I ${bam} \\
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
	-O ${sample_id}.dedup.bam \\
	-M ${sample_id}.metrics.txt

	samtools index ${sample_id}.dedup.bam
	
	
	"""

}

process BQSR {
	
	tag "$sample_id"
	debug true
	publishDir ("${params.outdir}/bqsr"), mode: 'copy'

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
	
	def mills_vcf = mills.find { it.name =~ /\.(vcf|vcf\.gz)$/ }
	def dbsnp_vcf = dbsnp.find { it.name =~ /\.(vcf|vcf\.gz)$/ }
	def kg_vcf = kg.find { it.name =~ /\.(vcf|vcf\.gz)$/ }
	
	"""
	verifybamid2 \\
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

	gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" CrosscheckFingerprints \\
	--INPUT ${bam} \\
	--HAPLOTYPE_MAP ${map} \\
	--LOD_THRESHOLD -20.0 \\
	--CROSSCHECK_BY READGROUP \\
	--OUTPUT ${sample_id}.crosscheck_metrics
		
	gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" BaseRecalibrator -I ${bam} -R ${ref} \\
	--known-sites ${mills_vcf} \\
	--known-sites ${dbsnp_vcf} \\
	--known-sites ${kg_vcf} \\
	-O ${sample_id}.recal.table

	gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" ApplyBQSR -I ${bam} -R ${ref} \\
	--bqsr-recal-file ${sample_id}.recal.table \\
	--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \\
	-O ${sample_id}.recal.bam

	samtools index ${sample_id}.recal.bam

	## Generate the covariate table after applying BQSR
	gatk --java-options "-Xmx${task.memory.toGiga() - 1}g" BaseRecalibrator -I ${sample_id}.recal.bam -R ${ref} \\
	--known-sites ${mills_vcf} \\
	--known-sites ${dbsnp_vcf} \\
	--known-sites ${kg_vcf} \\
	-O ${sample_id}.after.recal.table


	"""



}
