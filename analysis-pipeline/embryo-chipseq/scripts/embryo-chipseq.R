#################################################################
#################################################################
############### Embryo ChIP-Seq - R Support #################
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
    chip <- dba(sampleSheet=infile)

    # Add blacklist
    chip <- dba.blacklist(chip, blacklist=DBA_BLACKLIST_HG38, greylist=FALSE)

    # Count
    chip <- dba.count(chip, summits=1000, bUseSummarizeOverlaps=FALSE)

    # Save
    save(chip, file=outfile)

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
    count_matrix <- dba.peakset(dba.count(chip, peaks=NULL, score=DBA_SCORE_READS), bRetrieve=TRUE, DataType=DBA_DATA_FRAME) %>% mutate(CHR=glue('{CHR}:{START}-{END}')) %>% select(-START, -END) %>% column_to_rownames('CHR') %>% as.matrix

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
    peak_dataframe <- dba.peakset(chip, bRetrieve=TRUE) %>% as.data.frame %>% mutate(id=paste0('chip_peak_', 1:n()), score=100) %>% select(seqnames, start, end, id, score, strand) %>% mutate(seqnames=gsub('chr', '', seqnames), start=start-1)

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
    chip <- dba.normalize(chip)

    # Get contrasts
    chip <- dba.contrast(chip, minMembers = 2)

    # Filter
    if (grepl('H3K4me3', infile)) {
        chip$contrasts <- chip$contrasts[c(1,3)]
    }

    # Analyze
    chip <- dba.analyze(chip)

    # Save
    save(chip, file=outfile)

}

#######################################################
#######################################################
########## S6. TSS coverage
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
########## S6. TSS coverage
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

split_tss_types <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read abundance
    abundance_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='annot_transcript_id') %>% mutate(transcript_novelty_v2=ifelse(gene_novelty=='Known', transcript_novelty, gene_novelty)) %>% select(transcript_id, gene_novelty, transcript_novelty_v2)

    # Extract TSSs
    transcript_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end)), tss_coordinates=paste0('chr', chr, ':', tss)) %>% distinct %>% select(transcript_id, tss_coordinates, strand)

    # Merge
    merged_dataframe <- transcript_dataframe %>% left_join(abundance_dataframe, by='transcript_id') %>% replace_na(list(gene_novelty='Known', transcript_novelty='Known')) %>% replace_na(list(transcript_novelty_v2='Known'))
    head(merged_dataframe)

    # Pivot
    tss_dataframe <- merged_dataframe %>% group_by(tss_coordinates, strand) %>% summarize(transcript_types=paste(unique(transcript_novelty_v2), collapse=',')) %>% 
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

    # Loop
    for (tss_class in names(dataframes)) {
        
        # Get bed
        bed_dataframe <- dataframes[[tss_class]] %>% mutate(
            chr=gsub('chr(.*):.*', '\\1', tss_coordinates),
            start_int=as.numeric(gsub('.*:(.*)', '\\1', tss_coordinates)),
            start=format(start_int-251, scientific=FALSE, trim=TRUE), #-1 is actually unnecessary, because BED is 0-based
            end=format(start_int+250, scientific=FALSE, trim=TRUE),
            score=0) %>% select(chr, start, end, transcript_types, score, strand)
    
        # Add TSS number
        bed_dataframe$transcript_types <- paste0(tss_class, '_', 1:nrow(bed_dataframe))

        # Subset
        # nr_rows <- min(nrow(bed_dataframe), 500)
        # row_idx <- sample(1:nrow(bed_dataframe), nr_rows)
        # bed_dataframe <- as.data.frame(bed_dataframe)[row_idx,]
        
        # Get outfile
        bed_outfile <- glue('{outdir}/{tss_class}.bed')

        # Export
        fwrite(bed_dataframe, bed_outfile, sep='\t', col.names = FALSE)

    }

}

#######################################################
#######################################################
########## Plots
#######################################################
#######################################################

#############################################
#############################################
########## 2. TSS heatmap
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

split_tss <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read abundance
    abundance_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='annot_transcript_id') %>% mutate(transcript_novelty_v2=ifelse(gene_novelty=='Known', transcript_novelty, gene_novelty)) %>% select(transcript_id, gene_novelty, transcript_novelty_v2)

    # Extract TSSs
    transcript_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id, strand) %>% summarize(chr=seqid, strand=strand, tss=ifelse(strand=='+', min(start), max(end)), tss_coordinates=paste0('chr', chr, ':', tss)) %>% distinct %>% select(transcript_id, tss_coordinates, strand)

    # Merge
    merged_dataframe <- transcript_dataframe %>% left_join(abundance_dataframe, by='transcript_id') %>% replace_na(list(gene_novelty='Known', transcript_novelty='Known')) %>% replace_na(list(transcript_novelty_v2='Known'))
    head(merged_dataframe)

    # Pivot
    tss_dataframe <- merged_dataframe %>% group_by(tss_coordinates, strand) %>% summarize(transcript_types=paste(unique(transcript_novelty_v2), collapse=',')) %>% 
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

    # Loop
    for (tss_class in names(dataframes)) {
        
        # Get bed
        bed_dataframe <- dataframes[[tss_class]] %>% mutate(
            chr=gsub('chr(.*):.*', '\\1', tss_coordinates),
            start=as.numeric(gsub('.*:(.*)', '\\1', tss_coordinates))-1, #-1 is actually unnecessary, because BED is 0-based
            end=start+2, score=0) %>% select(chr, start, end, score, transcript_types, strand)

        # Subset
        # nr_rows <- min(nrow(bed_dataframe), 500)
        # row_idx <- sample(1:nrow(bed_dataframe), nr_rows)
        # bed_dataframe <- as.data.frame(bed_dataframe)[row_idx,]
        
        # Get outfile
        bed_outfile <- glue('{outdir}/{tss_class}.bed')

        # Export
        fwrite(bed_dataframe, bed_outfile, sep='\t', col.names = FALSE)

    }

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################