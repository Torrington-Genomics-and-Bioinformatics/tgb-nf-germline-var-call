process HAPLOTYPE_CALLER {

	tag "$sample_id"
	publishDir "${params.outdir}/gvcfs", mode: 'symlink'

	input:
	tuple val(sample_id), path(bam), path(bai), path(fmt)
	path ref
	path ref_index
	path ref_dict	
	path intervals


	output:
	tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi")

	script:
	"""
	
	fm=\$(cat ${fmt})
	
	gatk HaplotypeCaller -I ${bam} -R ${ref} \\
	-L ${intervals} \\
	-ip 100 \\
	-contamination \${fm} \\
	-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \\
	-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
	-ERC GVCF \\
	-O ${sample_id}.g.vcf.gz


	"""
	


}
