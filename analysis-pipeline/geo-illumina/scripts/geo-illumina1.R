#################################################################
#################################################################
############### GEO Illumina RNA-Seq - R Support #################
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tibble))
suppressPackageStartupMessages(require(glue))
suppressPackageStartupMessages(require(tidyr))

##### 2. Other libraries #####

#######################################################
#######################################################
########## S3. Align
#######################################################
#######################################################

#############################################
########## 4. Junction counts
#############################################

get_junction_counts <- function(infiles, outfile, junction_file) {

    # Read junctions
    junction_dataframe <- fread(junction_file) %>% mutate(junction_start=junction_start+1, junction_end=junction_end-1) %>% rename('chrom'='seqid') %>% select(-exon_ids)

    # Read SJ counts
    count_dataframe <- lapply(infiles, function(x) {
        
        # Read STAR junctions
        sj_dataframe <- fread(x, col.names = c('chrom', 'junction_start', 'junction_end', 'strand', 'intron_motif', 'annotated', 'uniquely_mapped_reads', 'multimapping_reads', 'max_spliced_alignment_overhang')) %>% 
            filter(strand !=0 ) %>% mutate(strand=recode(strand, '1'='+', '2'='-')) %>% select(-intron_motif, -annotated, -max_spliced_alignment_overhang)
        
        # Merge to introns
        merged_dataframe <- junction_dataframe %>% left_join(sj_dataframe, by=c('chrom', 'junction_start', 'junction_end', 'strand')) %>% replace_na(list(uniquely_mapped_reads=0, multimapping_reads=0)) %>% mutate(sample=gsub('.*/(.*)/.*', '\\1', x))

    }) %>% bind_rows

    # Count
    summary_dataframe <- count_dataframe %>% group_by(gene_id, transcript_id, sample) %>% summarize(
        min_unique_sj=min(uniquely_mapped_reads),
        mean_unique_sj=mean(uniquely_mapped_reads),
        max_unique_sj=max(uniquely_mapped_reads),
        min_multi_sj=min(multimapping_reads),
        mean_multi_sj=mean(multimapping_reads),
        max_multi_sj=max(multimapping_reads)
    )# %>% drop_na

    # Write
    fwrite(summary_dataframe, file=outfile, sep='\t')

}

#######################################################
#######################################################
########## S4. Expression
#######################################################
#######################################################

#############################################
########## 2. Aggregate
#############################################

aggregate_counts <- function(infiles, outfile) {
    
    # Load
    suppressPackageStartupMessages(require(tximeta)) #tximeta_1.8.2
    suppressPackageStartupMessages(require(SummarizedExperiment))
    suppressPackageStartupMessages(require(DESeq2))

    # Get sample dataframe
    sample_dataframe <- fread(infiles[1])

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[2]) %>% select(-source) %>% filter(type=='transcript')

    # Get transcript information
    transcript_dataframe <- gtf %>% select(gene_id, gene_name, gene_status, transcript_id, transcript_name, transcript_status, transcript_biotype, NNC_transcript, NIC_transcript, intergenic_transcript, antisense_transcript) %>% distinct
    rownames(transcript_dataframe) <- transcript_dataframe$transcript_id

    # Get gene information
    gene_dataframe <- transcript_dataframe %>% group_by(gene_id, gene_name, gene_status) %>% summarize(
        nr_transcripts=length(unique(transcript_id)),
        novel_transcripts=sum(transcript_status=='NOVEL'),
        NNC_transcripts=sum(NNC_transcript=='TRUE', na.rm=TRUE),
        NIC_transcripts=sum(NIC_transcript=='TRUE', na.rm=TRUE),
        intergenic_transcripts=sum(intergenic_transcript=='TRUE', na.rm=TRUE),
        antisense_transcripts=sum(antisense_transcript=='TRUE', na.rm=TRUE)
    ) %>% as.data.frame
    rownames(gene_dataframe) <- gene_dataframe$gene_id

    # Read isoform counts
    se <- tximeta(sample_dataframe, type='rsem', txIn=TRUE, txOut=TRUE)
    rowData(se) <- transcript_dataframe[rownames(rowData(se)),]

    # Read gene counts
    gse <- tximeta(sample_dataframe, type='rsem', txIn=TRUE, txOut=FALSE, tx2gene=transcript_dataframe %>% select(transcript_id, gene_id))
    rowData(gse) <- gene_dataframe[rownames(rowData(gse)),]

    # Create DESeq dataset
    dds_list <- list(
        'transcript' = DESeqDataSet(se, design=~cell_type),
        'gene' = DESeqDataSet(gse, design=~cell_type)
    )

    # Save
    save(se, gse, dds_list, file=outfile)

}

#############################################
########## 3. Get size factors
#############################################

get_size_factors <- function(infile, outfile) {
    
    # Library
    suppressPackageStartupMessages(require(DESeq2))

    # Load
    load(infile)

    # Get counts
    count_matrix <- counts(dds_list[['gene']])

    # Get factors
    normalization_factors <- edgeR::calcNormFactors(object = count_matrix, method = "TMM")

    # Get library sizes
    library_sizes <- colSums(count_matrix)

    # Get dataframe
    normalization_dataframe <- data.frame(normalization_factor=normalization_factors, library_size=library_sizes) %>% rownames_to_column('sample_name') %>% mutate(size_factor=normalization_factor*library_size/1000000, size_factor_reciprocal=1/size_factor)

    # Write
    fwrite(normalization_dataframe, file=outfile, sep='\t')

}

#######################################################
#######################################################
########## S5. Primate analysis
#######################################################
#######################################################

#############################################
########## 1. Copy lifted GTF
#############################################

copy_lifted_gtf <- function(infile, outfile) {
    
    # Read gtf
    gtf <- rtracklayer::import(infile)

    # Replace source
    gtf$source <- 'liftOver'

    # Add chromosome
    if (grepl('macaque', outfile)) {
        gtf <- diffloop::rmchr(gtf)
    }

    # Wrote
    rtracklayer::export(gtf, outfile, format='gtf')

}

#######################################################
#######################################################
########## S5. Sashimi
#######################################################
#######################################################

#############################################
########## 3. Novel primate genes
#############################################

get_novel_primate_genes <- function(infiles, outfile) {

    # Library
    require(DESeq2)
    
    # Get average expression
    expression_dataframe <- lapply(infiles, function(x) {
        load(x)
        dds <- estimateSizeFactors(dds_list[['gene']])
        average_dataframe <- counts(dds, normalized=TRUE) %>% rowMeans %>% as.data.frame %>% rownames_to_column('gene_id') %>% dplyr::rename('expression'='.') %>% mutate(organism=gsub('.*/(.*)-counts.rda', '\\1', x))
    }) %>% bind_rows %>% pivot_wider(id_cols = gene_id, names_from = organism, values_from = expression) %>% filter(grepl('TALON', gene_id))

    # Write
    fwrite(expression_dataframe, file=outfile, sep='\t')

}

#######################################################
###############################################z########
########## S10. Datasets
#######################################################
#######################################################

#############################################
########## 1. GTF
#############################################

create_gtf <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::import(infiles[1])

    # Get filtered transcript IDs
    novel_transcript_ids <- unique(rtracklayer::import(infiles[2])$transcript_id)

    # Filter
    gtf_filtered <- gtf[grepl('ENS', gtf$transcript_id) | gtf$transcript_id %in% novel_transcript_ids,]

    # Write
    rtracklayer::export(gtf_filtered, outfile)
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################