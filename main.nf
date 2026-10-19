include { QC } from './modules/001_fastqc.nf'
include { BWA_MEM } from './modules/002_align.nf'
include { MARK_DUPLICATES; BQSR} from './modules/003_markdup_bqsr.nf'
include { HAPLOTYPE_CALLER } from './modules/004_var_calling.nf'
include { GENOMICSDB_IMPORT } from './modules/005_gendbimp.nf'
//include { GENOTYPE_GVCFS; MERGE } from './modules/006_joint_genotyping.nf'


workflow {

        // Create channels to stream input files from their respective folders

        ref_ch = Channel.fromPath(params.ref).collect()
        ref_index_ch = Channel.fromPath("${params.ref}.{fai,amb,ann,bwt,pac,sa,dict}").collect()        
	ref_dict_ch = Channel.fromPath(params.ref.replaceAll(/\.fa$/, '.dict')).collect()
	kg_ch = Channel.fromPath("${params.kg}{,.idx}").collect()
	indel_ch = Channel.fromPath("${params.indel}{,.idx}").collect()
	omni_ch = Channel.fromPath("${params.omni}{,.tbi}").collect()        
	mills_ch = Channel.fromPath("${params.mills}{,.tbi}").collect()        
	dbsnp_ch = Channel.fromPath("${params.dbsnp}{,.idx}").collect()
	int_ch = Channel.fromPath(params.intervals).collect()
	bed_ch = Channel.fromPath(params.bed).collect()
	dat_ch = Channel.fromPath("${params.dat}.*").collect()
	map_ch = Channel.fromPath(params.map).collect()
	scatter_ch = Channel.fromPath("${params.interval_dir}/*.interval_list")

        reads_ch = Channel
                .fromPath(params.samplesheet)
                .splitCsv(header: true)
                .map { row -> tuple(row.sample_id, file(row.fastq_1), file(row.fastq_2)) }

        
        // 3. Wire processes using channels
	qc_ch = QC(reads_ch)
        aligned_ch = BWA_MEM(qc_ch, ref_ch, ref_index_ch)
        dedup_ch = MARK_DUPLICATES(aligned_ch)
        recal_ch = BQSR(dedup_ch, ref_ch, ref_index_ch, ref_dict_ch, dat_ch, map_ch, kg_ch, mills_ch, dbsnp_ch)
        gvcf_ch = HAPLOTYPE_CALLER(recal_ch.bam, ref_ch, ref_index_ch, ref_dict_ch, int_ch)
        
	// [benchmark clone] GenomicsDB import disabled; gVCFs get imported into the main DB later
//	// Getting the sample map for joint genotyping 
//	sample_map_ch = gvcf_ch
//                .map { id, gvcf, tbi -> "${id}\t${gvcf}" }                
//		.collectFile(name: 'cohort.sample_map', newLine: true)
//		.collect()
//
//	// Getting existing genomicsdb
//	gendb_ch = Channel
//		.fromPath("${params.gendb}/genomicsdb_*", type: 'dir')
//		.map { d -> 
//			def real = d.toRealPath()
//			tuple(d.name.replaceFirst(/^genomicsdb_/, ''), real.toString())}
//
//	        
//	db_ch = GENOMICSDB_IMPORT(sample_map_ch, gendb_ch)
        // geno_ch = GENOTYPE_GVCFS(db_ch, ref_ch, ref_index_ch, ref_dict_ch)

}
