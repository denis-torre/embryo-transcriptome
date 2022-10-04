#################################################################
#################################################################
############### Embryo ATAC-Seq - R Support #################
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
########## Plots
#######################################################
#######################################################

#######################################################
#######################################################
########## S5. Peak Counts
#######################################################
#######################################################

#############################################
########## 2. Get counts
#############################################

get_peak_counts <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(DiffBind))
    
    # Read data
    atac <- dba(sampleSheet=infile)

    # Add blacklist
    atac <- dba.blacklist(atac, blacklist=DBA_BLACKLIST_HG38, greylist=FALSE)

    # Get consensus peakset maybe - or not

    # Count
    atac <- dba.count(atac, summits=500, bUseSummarizeOverlaps=FALSE)

    # Save
    save(atac, file=outfile)

}

#############################################
########## 3. Get size factors
#############################################

get_size_factors <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(DiffBind))

    # Load
    load(infile)

    # Get count matrix
    count_matrix <- dba.peakset(dba.count(atac, peaks=NULL, score=DBA_SCORE_READS), bRetrieve=TRUE, DataType=DBA_DATA_FRAME) %>% mutate(CHR=glue('{CHR}:{START}-{END}')) %>% select(-START, -END) %>% column_to_rownames('CHR') %>% as.matrix

    # Get normalization factors (if you prefer to use the DESeq2 strategy use method="RLE" instead)
    normalization_factors <- edgeR::calcNormFactors(object = count_matrix, method = "TMM")

    # Library size
    library_sizes <- colSums(count_matrix)

    # Create dataframe
    normalization_dataframe <- data.frame(normalization_factor=normalization_factors, library_size=library_sizes) %>% rownames_to_column('sample_name') %>% mutate(size_factor=normalization_factor*library_size/1000000, size_factor_reciprocal=1/size_factor)

    # Write
    fwrite(normalization_dataframe, file=outfile, sep='\t')


}

#############################################
########## 6. Consensus peaks
#############################################

get_consensus_peaks <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(DiffBind))

    # Load
    load(infile)

    # Get peaks
    peak_dataframe <- dba.peakset(atac, bRetrieve=TRUE) %>% as.data.frame %>% mutate(id=paste0('atac_peak_', 1:n()), score=100) %>% select(seqnames, start, end, id, score, strand) %>% mutate(seqnames=gsub('chr', '', seqnames), start=start-1)
    # peak_dataframe <- dba.peakset(atac, bRetrieve=TRUE) %>% as.data.frame %>% mutate(id=paste0('atac_peak_', 1:n()), score=100) %>% select(seqnames, start, end, id, score, strand) %>% mutate(seqnames=gsub('chr', '', seqnames), start=start-1)# %>% mutate(start=start-1000, end=end+1000)

    # Write
    fwrite(peak_dataframe, file=outfile, sep='\t', col.names=FALSE)

}

#############################################
########## 7. Differential peaks
#############################################

get_differential_peaks <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(DiffBind))

    # Load
    load(infile)

    # Normalize
    atac <- dba.normalize(atac)

    # Get contrasts
    atac <- dba.contrast(atac, minMembers = 2)
    atac$contrasts <- atac$contrasts[c(1,7,12,18,20,21)]

    # Analyze
    atac <- dba.analyze(atac)

    # Save
    save(atac, file=outfile)

}

#######################################################
#######################################################
########## S5. TSS coverage
#######################################################
#######################################################

#############################################
########## 1. Get TSS BED
#############################################

get_tss_bed <- function(infile, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infile)

    # Extract TSSs
    tss_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end))) %>% distinct %>% mutate(start=tss-1, end=tss, score=100) %>% select(chr, start, end, transcript_id, score, strand)

    # Write
    fwrite(tss_dataframe, file=outfile, sep='\t', col.names=FALSE)

    # Genomic window
    window_length <- 500
    window_dataframe <- tss_dataframe %>% mutate(start=start-window_length/2, end=end+window_length/2)# %>% select(chr, start, end, transcript_id, score, strand)

    # Write
    fwrite(window_dataframe, file=gsub('.bed', glue('_{window_length}bp.bed'), outfile), sep='\t', col.names=FALSE)

}

#############################################
########## 3. Get counts
#############################################

get_tss_counts <- function(infile, outfile) {

    # Library
    suppressPackageStartupMessages(require(DiffBind))
    
    # Read data
    atac <- dba(sampleSheet=infile)

    # Count
    atac <- dba.count(atac, summits=FALSE, bUseSummarizeOverlaps=FALSE)

    # Save
    save(atac, file=outfile)

}

#############################################
#############################################
########## S6. TSS scores
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

split_tss_types <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read classification
    classification_dataframe <- fread(infiles[2]) %>% select(transcript_id, Transcript_novelty)
    head(classification_dataframe)

    # Extract TSSs
    transcript_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end)), tss_coordinates=paste0('chr', chr, ':', tss)) %>% distinct %>% select(transcript_id, tss_coordinates, strand)
    head(transcript_dataframe)

    # Merge
    merged_dataframe <- transcript_dataframe %>% left_join(classification_dataframe, by='transcript_id')
    head(merged_dataframe)

    # Pivot
    tss_dataframe <- merged_dataframe %>% group_by(tss_coordinates, strand) %>% summarize(transcript_types=paste(unique(Transcript_novelty), collapse=',')) %>% 
        mutate(tss_category=ifelse(grepl('Known', transcript_types), 
                                'Known_TSS',
                                ifelse(grepl('Intergenic', transcript_types), 
                                        'Intergenic_TSS', 
                                        ifelse(grepl('Antisense', transcript_types), 
                                                'Antisense_TSS',
                                                'Novel_TSS')))) %>% group_by(tss_category)

    # Split
    dataframes <- setNames(tss_dataframe %>% group_split, tss_dataframe %>% group_keys %>% pull(tss_category))

    # Get directory
    outdir <- dirname(outfile)

    # Loop through TSS
    for (tss_class in names(dataframes)) {

        # Loop through window
        for (window_size in c(0, as.numeric(gsub('.*_(.*)bp.bed', '\\1', outfile)))) {
            
            # Get bed (500bp)
            bed_dataframe <- dataframes[[tss_class]] %>% mutate(
                chr=gsub('chr(.*):.*', '\\1', tss_coordinates),
                start_int=as.numeric(gsub('.*:(.*)', '\\1', tss_coordinates)),
                start=format(start_int-1-(window_size/2), scientific=FALSE, trim=TRUE),
                end=format(start_int+(window_size/2), scientific=FALSE, trim=TRUE),
                score=0) %>% select(chr, start, end, transcript_types, score, strand)
        
            # Add TSS number
            bed_dataframe$transcript_types <- paste0(tss_class, '_', 1:nrow(bed_dataframe))

            # Get outfile
            bed_outfile <- ifelse(window_size>0, glue('{outdir}/{tss_class}_{window_size}bp.bed'), glue('{outdir}/{tss_class}.bed'))

            # Export
            fwrite(bed_dataframe, bed_outfile, sep='\t', col.names = FALSE)

        }

    }

}

#############################################
########## 3. Get promoter BED
#############################################

get_promoter_bed <- function(infile, outfile) {

    # Make TxDB
    txdb <- GenomicFeatures::makeTxDbFromGFF(infile)

    # Get promoters
    promoter_ranges <- GenomicFeatures::promoters(txdb, upstream=3000, downstream=500)

    # Convert to dataframe
    promoter_dataframe <- promoter_ranges %>% as.data.frame %>% mutate(start=start-1, score=100) %>% select(seqnames, start, end, tx_name, score, strand) %>% filter(start > 0 & end > 0)

    # Write
    fwrite(promoter_dataframe, file=outfile, sep='\t', col.names=FALSE)

}


#######################################################
#######################################################
########## S8. Gene coverage
#######################################################
#######################################################

#############################################
########## 1. Get gene BED
#############################################

get_gene_bed <- function(infile, outfile) {
    
}

#############################################
#############################################
########## 1. Isoform heatmap
#############################################
#############################################

#############################################
########## 1.1 Split GTF by isoform class
#############################################

split_gtf <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read classification
    classification_dataframe <- fread(infiles[2]) %>% select(transcript_id, Transcript_novelty)

    # Merge
    merged_gtf <- gtf %>% left_join(classification_dataframe, by='transcript_id')

    # Split
    gtf_split <- split(merged_gtf, merged_gtf$Transcript_novelty)

    # Get directory
    outdir <- dirname(outfile)

    # Loop
    for (transcript_class in names(gtf_split)) {
        
        # Get outfile
        gtf_outfile <- glue('{outdir}/{transcript_class}.gtf')
        
        # Export
        rtracklayer::export(gtf_split[[transcript_class]], gtf_outfile, format='gtf') # %>% head(5000)
        
    }

}

#############################################
########## 1.1.1 Merge shuffled bed
#############################################

create_shuffled_gtf <- function(infiles, outfile) {

    # Get filtered transcripts
    gtf <- rtracklayer::readGFF(infiles[1])

    # Get filtered transcripts
    transcript_ids <- gtf %>% pull(transcript_id) %>% unique

    # Read bed
    exon_gtf <- lapply(infiles[2:length(infiles)], function(x) {
        fread(x, col.names = c('seqid', 'start', 'end', 'exon_id', 'score', 'strand')) %>% mutate(start=start+1, transcript_id=gsub('(.*)_exon.*?-(.*)', '\\1-\\2', exon_id), seqid=gsub('chr', '', seqid)) %>% select(-exon_id) %>% mutate(type='exon')
    }) %>% bind_rows

    # Get transcript gtf
    transcript_gtf <- exon_gtf %>% group_by(seqid, score, strand, transcript_id) %>% summarize(start=min(start), end=max(end)) %>% select(seqid, start, end, score, strand, transcript_id) %>% mutate(type='transcript')

    # Merge
    shuffled_gtf <- rbind(transcript_gtf, exon_gtf) %>% arrange(seqid, transcript_id, start, strand) %>% mutate(source_transcript_id=gsub('(.*?)-.*', '\\1', transcript_id)) %>% filter(source_transcript_id %in% transcript_ids)# %>% select(-transcript_id)

    # Result
    result_gtf <- gtf %>% rename('source_transcript_id'='transcript_id') %>% select(gene_id, source_transcript_id, gene_name, gene_status, talon_gene, transcript_status, transcript_name, talon_transcript) %>% distinct %>%
        inner_join(shuffled_gtf, by='source_transcript_id') %>% select(-source_transcript_id)

    # Export
    rtracklayer::export(result_gtf, outfile, format='gtf')

}

#######################################################
#######################################################
########## S8. HOMER on novel gene clusters
#######################################################
#######################################################

#############################################
########## 8.1 Get cluster TSSs
#############################################

get_novel_gene_cluster_tss <- function(infiles, outfileRoot) {

    # Load clusters
    load(infiles[1])

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[2])

    # Get cluster dataframe
    cluster_dataframe <- data.frame(cluster=cluster_results$cluster) %>% rownames_to_column('gene_id') %>% mutate(cluster=paste0('cluster_', cluster))

    # Get TSS locations
    # tss_dataframe <- gtf %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% select(gene_id, transcript_id, seqid, tss, strand) %>% mutate(start=tss-1, end=tss, score=1) %>% inner_join(cluster_dataframe, by='gene_id') %>% select(seqid, start, end, transcript_id, score, strand, cluster)
    tss_dataframe <- gtf %>% filter(type=='transcript') %>% mutate(tss=ifelse(strand=='+', start, end)) %>% select(gene_id, transcript_id, seqid, tss, strand) %>% mutate(start=tss-1, end=tss, score=1) %>% inner_join(cluster_dataframe, by='gene_id') %>% select(seqid, start, end, transcript_id, score, strand, cluster)

    # Split
    tss_dataframes <- split(tss_dataframe, tss_dataframe$cluster)

    # Write
    for (cluster_nr in names(tss_dataframes)) {
        
        # Get outfile
        outfile <- glue(outfileRoot)
        
        # Write
        fwrite(tss_dataframes[[cluster_nr]] %>% select(-cluster), file=outfile, sep='\t', col.names=FALSE)
        
    }

    # Get all clusters
    cluster_nr <- 'all'
    outfile <- glue(outfileRoot)
    
    # Write
    fwrite(tss_dataframe %>% select(-cluster), file=outfile, sep='\t', col.names=FALSE)

}

#############################################
#############################################
########## 2. TSS heatmap
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

split_tss <- function(infiles, outfile) {

    # # Read GTF
    # gtf <- rtracklayer::readGFF(infiles[1])

    # # Read abundance
    # abundance_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='annot_transcript_id') %>% mutate(transcript_novelty_v2=ifelse(gene_novelty=='Known', transcript_novelty, gene_novelty)) %>% select(transcript_id, gene_novelty, transcript_novelty_v2)

    # # Extract TSSs
    # transcript_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end)), tss_coordinates=paste0('chr', chr, ':', tss)) %>% distinct %>% select(transcript_id, tss_coordinates, strand)

    # # Merge
    # merged_dataframe <- transcript_dataframe %>% left_join(abundance_dataframe, by='transcript_id') %>% replace_na(list(gene_novelty='Known', transcript_novelty='Known')) %>% replace_na(list(transcript_novelty_v2='Known'))
    # head(merged_dataframe)

    # # Pivot
    # tss_dataframe <- merged_dataframe %>% group_by(tss_coordinates, strand) %>% summarize(transcript_types=paste(unique(transcript_novelty_v2), collapse=',')) %>% 
    #     mutate(tss_category=ifelse(grepl('Known', transcript_types), 
    #                             'Known_TSS',
    #                             ifelse(grepl('Intergenic', transcript_types), 
    #                                     'Intergenic_TSS', 
    #                                     ifelse(grepl('Antisense', transcript_types), 
    #                                             'Antisense_TSS',
    #                                             'Novel_TSS')))) %>% group_by(tss_category)

    # # Split
    # dataframes <- setNames(tss_dataframe %>% group_split, tss_dataframe %>% group_keys %>% pull(tss_category))

    # # Get directory
    # outdir <- dirname(outfile)

    # # Loop
    # for (tss_class in names(dataframes)) {
        
    #     # Get bed
    #     bed_dataframe <- dataframes[[tss_class]] %>% mutate(
    #         chr=gsub('chr(.*):.*', '\\1', tss_coordinates),
    #         start=as.numeric(gsub('.*:(.*)', '\\1', tss_coordinates))-1, #-1 is actually unnecessary, because BED is 0-based
    #         end=start+2, score=0) %>% select(chr, start, end, score, transcript_types, strand)

    #     # Subset
    #     # nr_rows <- min(nrow(bed_dataframe), 500)
    #     # row_idx <- sample(1:nrow(bed_dataframe), nr_rows)
    #     # bed_dataframe <- as.data.frame(bed_dataframe)[row_idx,]
        
    #     # Get outfile
    #     bed_outfile <- glue('{outdir}/{tss_class}.bed')

    #     # Export
    #     fwrite(bed_dataframe, bed_outfile, sep='\t', col.names = FALSE)

    # }

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################