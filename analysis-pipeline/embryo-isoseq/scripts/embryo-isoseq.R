#################################################################
#################################################################
############### Embryo IsoSeq - R Support #################
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
########## S2. Align
#######################################################
#######################################################

#############################################
########## 1. Copy FASTQ
#############################################

copy_fastq <- function(infile, outfile) {
    
    # Library
    suppressPackageStartupMessages(require(ShortRead))
    
    # Read
    message(glue('Reading {infile}...'))
    fastq <- readFastq(infile)

    # Find duplicate names
    duplicate_bool <- duplicated(as.character(id(fastq)))

    # Check duplicate names
    if (any(duplicate_bool)) {
        
        # Get indices
        duplicate_idx <- which(duplicate_bool)
        message('Following duplicate sequence names found:')
        message(paste(unique((as.character(id(fastq)[duplicate_bool]))), collapse=', '))
        
        # Get name
        name_dataframe <- id(fastq)[duplicate_idx] %>% as.data.frame %>% group_by(x) %>% mutate(count=row_number()) %>% ungroup %>% mutate(index=duplicate_idx) %>% mutate(new_name=paste0(x, '_', count))    
        
        # Loop
        for (i in 1:nrow(name_dataframe)) {
        
            # Get sequence index
            sequence_index <- name_dataframe[i,'index'][[1]]
            
            # Replace sequence with new name
            fastq[sequence_index] <- ShortReadQ(sread=sread(fastq[sequence_index]), quality=quality(fastq[sequence_index]), id=BStringSet(name_dataframe[i,'new_name'][[1]]))
            
        }

    } else {
        message('No duplicate sequence names found.')
    }

    # Write
    writeFastq(fastq, outfile, compress=FALSE)
}

#######################################################
#######################################################
########## S3. Illumina alignment
#######################################################
#######################################################

#############################################
########## 4. Merge SJ files
#############################################

merge_sjs <- function(infiles, outfile) {
    
    # Split infiles
    infiles_split <- sapply(c('pass1', 'pass2'), function(x) infiles[grepl(x, infiles)], simplify=FALSE)

    # Get colnames
    sj_colnames <- c('chr', 'intron_start', 'intron_end', 'strand', 'intron_motif', 'annotated', 'unique_junction_reads', 'multimapping_junction_reads', 'max_spliced_alignment_overhang')

    # Get annotations
    annotation_dataframe <- lapply(infiles_split[['pass1']], function(x) {
        df <- fread(x);
        colnames(df) <- sj_colnames;
        df <- df %>% select(chr, intron_start, intron_end, strand, intron_motif, annotated)
        return(df)
    }) %>% bind_rows %>% distinct %>% rename('annotated_ensembl'='annotated')

    # Read
    sj_dataframe <- lapply(infiles_split[['pass2']], function(x) {
        df <- fread(x);
        colnames(df) <- sj_colnames;
        return(df)
    }) %>% bind_rows %>% group_by(chr, intron_start, intron_end, strand, intron_motif, annotated) %>% summarize(tot_unique_junction_reads=sum(unique_junction_reads), tot_multimapping_junction_reads=sum(multimapping_junction_reads), max_spliced_alignment_overhang=max(max_spliced_alignment_overhang))
#  %>% filter(strand != 0)

    ### NEW
    # Fix strands
    strand_dataframe <- lapply(c(1, 2), function(x) {
        sj_dataframe %>% filter(strand==0) %>% mutate(strand=x)
    }) %>% bind_rows

    # Concatenate
    result_dataframe <- rbind(sj_dataframe %>% filter(strand != 0), strand_dataframe)
    
    # Merge
    merged_dataframe <- result_dataframe %>% left_join(annotation_dataframe, by=c('chr', 'intron_start', 'intron_end', 'strand', 'intron_motif')) %>% select(-annotated_ensembl)
    ###

    # Merge
    # merged_dataframe <- sj_dataframe %>% left_join(annotation_dataframe, by=c('chr', 'intron_start', 'intron_end', 'strand', 'intron_motif')) %>% select(-annotated_ensembl)

    # Write
    fwrite(merged_dataframe, file=outfile, sep='\t', col.names=FALSE)

}

#######################################################
#######################################################
########## S6. Illumina alignment
#######################################################
#######################################################

#############################################
########## 4. Aggregate counts
#############################################

# aggregate_salmon_counts <- function(infiles, outfile) {

#     # Load
#     require(tximport)

#     # Fix names
#     names(infiles) <- gsub('.*/(.*)/quant.sf', '\\1', infiles)

#     # Import
#     txi <- tximport(infiles, type = "salmon", txOut = TRUE, countsFromAbundance='scaledTPM')

#     # Get dataframe
#     tpm_dataframe <- txi$abundance %>% as.data.frame %>% rownames_to_column('ID')

#     # Write
#     fwrite(tpm_dataframe, file=outfile, sep='\t')

# }

#######################################################
#######################################################
########## S3. Collapse
#######################################################
#######################################################

#############################################
########## 3. Statistics
#############################################

# get_collapse_statistics <- function(infiles, outfile) {
    
#     # Get parameters
#     sample <- gsub('.*/(.*)-5p.*_corrected.gtf', '\\1', infiles[1])
#     max_5_diff <- gsub('.*-5p(.*)-J(.*)-3p(.*)_corrected.gtf', '\\1', infiles[1])
#     max_fuzzy_junction <- gsub('.*-5p(.*)-J(.*)-3p(.*)_corrected.gtf', '\\2', infiles[1])
#     max_3_diff <- gsub('.*-5p(.*)-J(.*)-3p(.*)_corrected.gtf', '\\3', infiles[1])

#     # Read GTF
#     gtf <- rtracklayer::readGFF(infiles[1])

#     # Read SQANTI
#     sqanti_dataframe <- fread(infiles[2]) %>% select(isoform, structural_category, associated_gene) %>% rename('transcript_id'='isoform') %>% mutate(gene_id=gsub('(.*)\\..*', '\\1', transcript_id))

#     # Get genes
#     pb_gene_dataframe <- sqanti_dataframe %>% group_by(gene_id) %>% summarize(ensembl_gene_id=paste(unique(associated_gene), collapse=','))

#     ### Statistics
#     # Number of transcripts per gene
#     transcript_dataframe <- gtf %>% filter(type=='transcript') %>% group_by(gene_id) %>% tally(name = 'nr_transcripts') %>% mutate(sample=sample, max_5_diff=max_5_diff, max_fuzzy_junction=max_fuzzy_junction, max_3_diff=max_3_diff) %>% merge(pb_gene_dataframe, by='gene_id')

#     # Number of exons per transcript
#     exon_dataframe <- gtf %>% filter(type=='exon') %>% group_by(transcript_id) %>% tally(name = 'nr_exons') %>% merge(sqanti_dataframe, by='transcript_id') %>% mutate(sample=sample, max_5_diff=max_5_diff, max_fuzzy_junction=max_fuzzy_junction, max_3_diff=max_3_diff)

#     # Number of junction-matched transcripts
#     matching_dataframe <- gtf %>% filter(type=='exon') %>% mutate(coordinates=paste0(start, '-', end)) %>% group_by(seqid, gene_id, transcript_id) %>% 
#         summarize(nr_exons=length(coordinates), coordinates_merged=paste0(sort(coordinates), collapse='-')) %>% filter(nr_exons > 1) %>%
#         mutate(junctions=paste0(glue('chr{seqid}:'), '5p-', gsub('.*?-(.*)-.*', '\\1', coordinates_merged), '-3p')) %>% group_by(gene_id, junctions) %>% tally(name='nr_matched_transcripts') %>%
#         merge(pb_gene_dataframe, by='gene_id') %>% arrange(-nr_matched_transcripts) %>% mutate(sample=sample, max_5_diff=max_5_diff, max_fuzzy_junction=max_fuzzy_junction, max_3_diff=max_3_diff)

#     ### Write
#     fwrite(transcript_dataframe, file=outfile, sep='\t', row.names=FALSE)
#     fwrite(exon_dataframe, file=gsub('transcript', 'exon', outfile), sep='\t', row.names=FALSE)
#     fwrite(matching_dataframe, file=gsub('transcript', 'matching_transcript', outfile), sep='\t', row.names=FALSE)

# }

#############################################
########## 4. Make report
#############################################

# make_collapse_report <- function(infiles, outfile) {

#     # Library
#     suppressPackageStartupMessages(require(ggplot2))
#     suppressPackageStartupMessages(require(ggrepel))
    
#     # Transcripts
#     transcript_dataframe <- lapply(infiles[grepl('\\.transcript_stats.tsv', infiles)], fread) %>% bind_rows %>% group_by(max_5_diff, max_3_diff) %>% summarize(avg_transcripts=mean(nr_transcripts))

#     # Exons
#     exon_dataframe <- lapply(infiles[grepl('exon_stats.tsv', infiles)], fread) %>% bind_rows %>% group_by(max_5_diff, max_3_diff) %>% summarize(avg_exons=mean(nr_exons), nr_monoexons=sum(nr_exons==1))

#     # Matching transcripts
#     matching_dataframe <- lapply(infiles[grepl('matching_transcript_stats.tsv', infiles)], fread) %>% bind_rows %>% filter(nr_matched_transcripts > 1) %>% mutate(nr_matched_transcripts=nr_matched_transcripts-1) %>% group_by(max_5_diff, max_3_diff) %>% summarize(avg_matching_transcripts=mean(nr_matched_transcripts), nr_matching_transcripts=sum(nr_matched_transcripts))

#     # Merge
#     merged_dataframe <- transcript_dataframe %>% merge(exon_dataframe, by=c('max_5_diff', 'max_3_diff')) %>% merge(matching_dataframe, by=c('max_5_diff', 'max_3_diff'))

#     # Plot
#     gp <- ggplot(merged_dataframe, aes(x=nr_monoexons, y=nr_matching_transcripts, color=as.factor(max_5_diff), group=max_3_diff)) +
#         geom_point() +
#         geom_line() +
#         geom_label_repel(data=merged_dataframe %>% filter(max_5_diff==max(merged_dataframe$max_5_diff)) %>% mutate(label=glue('{formatC(max_3_diff, big.mark=",")}bp 3\'')), aes(label=label), nudge_y=10000, size=3, color='black') +
#         scale_color_brewer(type='div', palette=7, direction=-1) +
#         scale_x_continuous(labels = scales::comma) +
#         scale_y_continuous(labels = scales::comma) +
#         labs(x='Mono-exonic transcripts', y='Transcripts with alternate start/end', color='5\' collapse\nthreshold (bp)', title=glue('Transcript collapsing parameters for {strsplit(basename(outfile), "-")[[1]][1]}')) +
#         theme_classic() + theme(plot.title = element_text(hjust = 0.5))

#     # Save
#     ggsave(outfile, plot=gp, height=5, width=9)

# }

#######################################################
#######################################################
########## S4. Merged
#######################################################
#######################################################

#############################################
########## 4. Fix gene IDs
#############################################

# fix_chained_gff <- function(infile, outfile) {

#     # Library
#     require(rtracklayer)

#     # Read
#     gtf <- import(infile)

#     # Get transcript coordinates
#     coordinate_dataframe <- gtf %>% as.data.frame %>% filter(type=='exon') %>% mutate(boundaries=paste0(start, '-', end)) %>% group_by(gene_id, transcript_id) %>% summarize(coordinates=paste0(boundaries, collapse='-'))

#     # Duplicated transcripts
#     duplicate_transcripts <- coordinate_dataframe$transcript_id[which(duplicated(coordinate_dataframe$coordinates))]

#     # Filter
#     gtf_filtered <- gtf[!gtf$transcript_id %in% duplicate_transcripts,]

#     # Write
#     export(gtf_filtered, outfile, format='gtf')
# }

#######################################################
#######################################################
########## S8. SQANTI Second Pass
#######################################################
#######################################################

#############################################
########## 3. ISM filter
#############################################

# filter_isms <- function(infiles, outfile) {

#     # Library
#     require(rtracklayer)

#     # Read GTF
#     gtf <- readGFF(infiles[1])

#     # Get ISMs
#     ism_ids <- fread(infiles[2]) %>% filter(structural_category=='incomplete-splice_match') %>% pull(isoform)

#     # Remove
#     gtf_filtered <- gtf %>% filter(!transcript_id %in% ism_ids)

#     # Write
#     export(gtf_filtered, outfile)
# }

#######################################################
#######################################################
########## S5. TALON
#######################################################
#######################################################

#############################################
########## 8. Get SJs
#############################################

get_junctions <- function(infile, outfile) {
    
    # Library
    require(parallel)

    # Read GTF
    gtf <- rtracklayer::readGFF(infile) %>% select(-source) %>% filter(type=='exon')# %>% head(500)# %>% select(gene_id, transcript_id, seqid, start, end)

    # Get multiexonic transcripts
    exon_counts <- table(gtf$transcript_id)
    multiexon_transcripts <- names(exon_counts)[exon_counts > 1]

    # Filter
    gtf_filtered <- gtf %>% filter(transcript_id %in% multiexon_transcripts)

    # Split
    gtf_split <- split(gtf_filtered, gtf_filtered$transcript_id)

    # Get cores
    cores <- 10

    # Make cluster
    cluster <- makeCluster(cores)

    # Load libraries
    clusterEvalQ(cluster, library("dplyr"));

    # Loop
    junction_dataframe <- parSapply(cluster, gtf_split, function(x) {
        
        # Get GTF
        transcript_gtf <- x %>% arrange(start) %>% select(gene_id, transcript_id, exon_id, seqid, start, end, strand)
        
        # Get strand
        strand <- unique(transcript_gtf$strand)
        
        # Get junctions
        junctions <- transcript_gtf %>% group_by(gene_id, transcript_id, seqid, strand) %>% mutate(junction_start=end, junction_end=lead(start), exon_ids=paste0(exon_id, '-', lead(exon_id))) %>% ungroup %>% head(-1) %>% select(gene_id, transcript_id, exon_ids, seqid, junction_start, junction_end, strand)
        
    }, simplify=FALSE) %>% bind_rows

    # Write
    fwrite(junction_dataframe, outfile, sep='\t')#, compress='gzip')

}

#############################################
########## 8. Get SJs
#############################################

annotate_splice_sites <- function(infile, outfile) {

    # Library
    require(GenomicRanges)
    
    # Read GTF
    gtf_dataframe <- rtracklayer::readGFF(infile) %>% filter(type=='exon')
    
    # Get exons per gene
    exon_dataframe <- gtf_dataframe %>% select(seqid:exon_id) %>% group_by(transcript_id) %>% mutate(nr_exons=n(), exon_position=ifelse(exon_number==1, 'first_exon', ifelse(exon_number==nr_exons, 'last_exon', 'intermediate_exon'))) %>% filter(nr_exons>1) %>% mutate(exon_5p=paste0('chr', seqid, ':', ifelse(strand=='+', start, end)), exon_3p=paste0('chr', seqid, ':', ifelse(strand=='+', end, start))) %>%  
        select(gene_id, gene_name, transcript_id, exon_id, exon_number, strand, exon_position, exon_5p, exon_3p) %>% pivot_longer(-c(gene_id, gene_name, transcript_id, exon_id, exon_number, strand, exon_position), names_to = 'splice_site', values_to='coordinates')

    # Get splice sites
    ss_dataframe <- exon_dataframe %>% filter(!((splice_site=='exon_5p' & exon_position=='first_exon') | (splice_site=='exon_3p' & exon_position=='last_exon'))) %>% group_by(gene_id, gene_name, coordinates) %>% summarize(nr_transcripts=n(), splice_site_novelty=ifelse(any(grepl('ENS', exon_id)), 'known', 'novel'))

    # Get ranges
    range_dataframe <- ss_dataframe %>% arrange(coordinates) %>% group_by(gene_id, gene_name, splice_site_novelty) %>% mutate(splice_site_nr=1:n(), splice_site_name=paste0(gene_id, '_', splice_site_novelty, '_', splice_site_nr)) %>% arrange(gene_id) %>% separate(coordinates, into=c('seqid', 'splice_site_coordinates'), convert=TRUE) %>% mutate(start=splice_site_coordinates-10, end=splice_site_coordinates+10) %>% filter(start>0)

    # Convert to ranges
    site_ranges <- makeGRangesFromDataFrame(range_dataframe)
    names(site_ranges) <- range_dataframe$splice_site_name

    # Write
    fwrite(ss_dataframe, file=outfile, sep='\t')

    # Export bed
    rtracklayer::export(site_ranges, gsub('s.tsv', '_ranges_10bp.bed', outfile))

}

#############################################
########## 7. Get SJ support
#############################################

get_sj_counts <- function(infile, outfile) { # rewrite to count number of samples with coverage ≥1

    # Read
    junction_dataframe <- fread(infile) %>% rename('transcript_id'='isoform')

    # Get columns
    unique_sj_cols <- colnames(junction_dataframe)[grepl('SJ.out_unique', colnames(junction_dataframe))]

    # Count
    filtered_dataframe <- lapply(seq(1, 5), function(min_sj_reads) {
        junction_dataframe %>% group_by(transcript_id) %>% summarize_at(all_of(unique_sj_cols), function(x) all(x>=min_sj_reads)) %>% mutate(samples_with_unique_sj_support=rowSums(across(all_of(unique_sj_cols)))) %>% select(transcript_id, samples_with_unique_sj_support) %>% mutate(min_sj_reads=min_sj_reads)
    }) %>% bind_rows

    # Write
    fwrite(filtered_dataframe, file=outfile, sep='\t')

}

#############################################
########## 8. Filter GTF
#############################################

filter_gtf <- function(infiles, outfile) {

    # Read GTF
    gtf <- rtracklayer::import(infiles[1])

    # Read SJ counts
    sj_dataframe <- fread(infiles[2])

    # Read FL counts
    abundance_dataframe <- fread(infiles[3])
    abundance_dataframe$fl_counts <- apply(abundance_dataframe[,12:ncol(abundance_dataframe)], 1, sum)

    # Get SJ-supported multiexonic transcripts
    min_sj_samples <- ifelse(grepl('human', outfile), 3, 2)
    sj_multiexon_transcripts <- sj_dataframe %>% filter(min_sj_reads == 1 & samples_with_unique_sj_support >= min_sj_samples) %>% pull(transcript_id)
    # length(sj_multiexon_transcripts)

    # Get PB-supported monoexonic transcripts
    pb_monoexon_transcripts <- abundance_dataframe %>% filter(transcript_novelty=='Known' & n_exons==1 & fl_counts >=1) %>% pull(annot_transcript_id)
    # length(pb_monoexon_transcripts)

    # Get ISMs
    ism_transcripts <- abundance_dataframe %>% filter(transcript_novelty == 'ISM') %>% pull(annot_transcript_id)
    # length(ism_transcripts)

    # Merge
    transcript_ids <- setdiff(c(sj_multiexon_transcripts, pb_monoexon_transcripts), ism_transcripts)
    # length(transcript_ids)

    # Filter GTF
    gtf_filtered <- gtf[gtf$type != 'gene' & gtf$transcript_id %in% transcript_ids]

    # Write
    rtracklayer::export(gtf_filtered, outfile)

}

#############################################
########## 10. Classify transcripts
#############################################

create_transcript_classification <- function(infiles, outfile) {

    # Read GTF
    if (grepl('SJ', outfile)==TRUE) {
        transcript_dataframe <- rtracklayer::readGFF(infiles[1]) %>% filter(type=='exon') %>% mutate(exon_length=end-start+1) %>% group_by(gene_id, transcript_id) %>% summarize(transcript_length=sum(exon_length), nr_exons=length(exon_length))
    } else {
        transcript_dataframe <- rtracklayer::readGFF(infiles[1]) %>% select(-source) %>% filter(type=='exon') %>% mutate(exon_length=end-start+1) %>% group_by(gene_id, transcript_id) %>% summarize(transcript_length=sum(exon_length), nr_exons=length(exon_length))
    }

    # Read classification
    abundance_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='annot_transcript_id') %>% mutate(fl_counts=rowSums(across(starts_with('human_')))) %>% select(transcript_id, gene_novelty, transcript_novelty, fl_counts)

    # Read SQANTI
    sqanti_dataframe <- fread(infiles[3]) %>% rename('transcript_id'='isoform') %>% select(transcript_id, structural_category, subcategory)

    # Merge
    merged_dataframe <- transcript_dataframe %>% left_join(abundance_dataframe, by='transcript_id') %>% left_join(sqanti_dataframe, by='transcript_id')  %>% replace_na(list(gene_novelty='Known', transcript_novelty='Known', 'fl_counts'=0))

    # Novel gene dataframe
    gene_dataframe <- merged_dataframe %>% filter(gene_novelty != 'Known') %>% group_by(gene_id, gene_novelty) %>% summarize(sqanti_any_antisense=any(structural_category=='antisense'), sqanti_all_intergenic=all(structural_category=='intergenic'))

    # Fix novelty of Antisense and Intergenic (rename antisense of intergenic genes to intergenic, rename intergenic genes called as antisense by SQANTI to antisense)
    gene_dataframe$gene_novelty_merged <- apply(gene_dataframe, 1, function(x) {
        if (x['gene_novelty']=='Antisense') {
            if (x['sqanti_all_intergenic']) {
                result <- 'Intergenic'
            } else {
                result <- 'Antisense'
            }
        } else if (x['gene_novelty'] == 'Intergenic') {
            if (x['sqanti_any_antisense']) {
                result <- 'Antisense'
            } else {
                result <- 'Intergenic'
            }
        }
    })

    # Merge
    result_dataframe <- merged_dataframe %>% left_join(gene_dataframe %>% select(gene_id, gene_novelty_merged), by='gene_id') %>% replace_na(list(gene_novelty_merged='Known'))

    # Fix TALON NIC issue to NNC as SQANTI, rename genomic to NNC
    result_dataframe$Transcript_novelty <- apply(result_dataframe, 1, function(x) {
        if (x['gene_novelty_merged']=='Known') {
            if ((x['transcript_novelty'] == 'NIC') && (x['structural_category'] == 'novel_not_in_catalog')) {
                result <- 'NNC'
            } else if (x['transcript_novelty'] == 'Genomic') {
                result <- 'NNC'
            } else {
                result <- x['transcript_novelty']
            }
        } else {
            result <- x['gene_novelty_merged']
        }
    })

    # Write
    fwrite(result_dataframe, file=outfile, sep='\t')

}

#############################################
########## 10. Split transcript GTF
#############################################

split_transcript_gtf <- function(infiles, outfileRoot) {

    # Read GTF
    gtf <- rtracklayer::import(infiles[1])
    gtf <- gtf[!gtf$type =='gene',]

    # Read abundance
    abundance_dataframe <- fread(infiles[2])
    # abundance_dataframe$fl_counts <- rowSums(abundance_dataframe[,12:ncol(abundance_dataframe)])
    # abundance_dataframe <- abundance_dataframe %>% filter(fl_counts >= 3)

    # Get transcript classes
    transcript_classes <- split(abundance_dataframe$transcript_id, abundance_dataframe$Transcript_novelty) # merged classification file
    # transcript_classes <- split(abundance_dataframe$annot_transcript_id, abundance_dataframe$transcript_novelty) # TALON abundance file

    # Fix known
    all_transcripts <- unique(gtf$transcript_id)
    transcript_classes$Known <- all_transcripts[grepl('ENS', all_transcripts)]

    # Write transcript only
    gtf_transcript <- gtf[gtf$type == 'transcript',]
    rtracklayer::export(gtf_transcript, gsub('\\{transcript_class\\}', 'transcript', outfileRoot))

    # Loop
    for (transcript_class in names(transcript_classes)) {
        
        # Get transcripts
        transcript_ids <- transcript_classes[[transcript_class]]
        
        # Filter
        gtf_filtered <- gtf[gtf$transcript_id %in% transcript_ids,]
        
        # Get outfile
        outfile <- glue(outfileRoot)
        
        # Write
        if (length(gtf_filtered) > 0) {
            rtracklayer::export(gtf_filtered, outfile)
        }

    }
}

#######################################################
#######################################################
########## S6. CPAT
#######################################################
#######################################################

#############################################
########## 3. Split GTF
#############################################

split_gtf <- function(infile, outfileRoot) {
    
    # Library
    suppressPackageStartupMessages(require(GenomicFeatures))
    suppressPackageStartupMessages(require(rtracklayer))

    # Read GTF
    gtf <- import(infile)

    # Get gene chunks
    gene_ids <- unique(gtf$gene_id)
    gene_chunks <- split(gene_ids, cut(seq_along(1:length(gene_ids)), 50, labels=FALSE))

    # Loop
    for (i in names(gene_chunks)) {

        # Get subset
        message(glue('Doing chunk {i}...'))
        gtf_subset <- gtf[gtf$gene_id %in% gene_chunks[[i]]]
        
        # Get outfile
        chunk_nr <- stringr::str_pad(i, width=2, pad='0')
        outfile <- glue(outfileRoot)
        
        # Write
        export(gtf_subset, outfile, format='gtf')
        
    }
}

#############################################
########## 4. Add CDS
#############################################

add_cds <- function(infiles, outfile, coding_cutoff) {
    
    # Library
    suppressPackageStartupMessages(require(ensembldb))
    suppressPackageStartupMessages(require(rtracklayer))

    ### 1. Create EnsemblDB
    # Read GTF
    gtf <- import(infiles[1])

    # Remove gene boundaries (some novel transcripts are outside of Ensembl-defined ones)
    gtf <- gtf[!gtf$type == 'gene',]

    # Add fake CDS from first exon of first gene
    exon <- gtf[gtf$type=='exon'][1,]
    gtf_cds <- c(gtf, GRanges(seqnames(exon), ranges(exon), strand=strand(exon), type='CDS', source='PacBio', gene_id=exon$gene_id, transcript_id=exon$transcript_id))

    # Get genome info
    organism <- ifelse(grepl('mouse', outfile), 'Mus_musculus', 'Homo_sapiens')
    genomeVersion <- ifelse(grepl('mouse', outfile), 'GRCm38', 'GRCh38')

    # Create database
    db_outfile <- gsub('.gtf', '.sqlite', outfile, fixed=TRUE)
    db_path <- ensDbFromGRanges(gtf_cds, organism = organism, genomeVersion = genomeVersion, version = 102, outfile=db_outfile)
    ensdb <- EnsDb(db_path)

    ### 2. Read CPAT results
    # Read CPAT
    cpat_dataframe <- fread(infiles[2]) %>% dplyr::filter(seq_ID %in% unique(gtf$transcript_id) & Coding_prob >= as.numeric(coding_cutoff))

    # ORF ranges
    orf_ranges <- IRanges(start=cpat_dataframe$ORF_start, width=cpat_dataframe$ORF, names=cpat_dataframe$seq_ID)

    ### 3. Add CDS coordinates
    # Get CDS genomic coordinates
    message('Mapping ORFS to genomic coordinates...')
    cds_ranges <- transcriptToGenome(orf_ranges, ensdb)
    cds_ranges

    # Add phase
    message('Adding phase info...')
    cds_ranges_phase <- lapply(cds_ranges, function(x) {

        # Get CDS mod
        cds_phase <- x %>% as.data.frame %>% arrange(ifelse(strand=='+', start, -start)) %>% mutate(
            cds_mod = (end-start+1)%%3,
            phase = 0,
        )

        # Add phase
        if (nrow(cds_phase) > 1) {
            for (i in 2:nrow(cds_phase)) {
                cds_phase[i,'phase'] <- (3-cds_phase[i-1, 'cds_mod']+cds_phase[i-1, 'phase']) %% 3
            }
        }
        
        # Convert
        result <- makeGRangesFromDataFrame(cds_phase, keep.extra.columns=TRUE);

        # Return
        return(result)

    })

    # Split GTF by transcript to add CDS
    gtf_split <- split(gtf, gtf$transcript_id)

    # Add CDS to GTF if available, merging CDS coordinates by exon ID to preserve exon metadata
    message('Merging...')
    result_gtf <- do.call('c', lapply(names(gtf_split), function(x) {
        result <- gtf_split[[x]];
        if (x %in% names(cds_ranges_phase)) {
            cds_ranges_annotated <- as.data.frame(cds_ranges_phase[[x]]) %>% dplyr::select(start, end, width, phase, exon_id) %>% left_join(as.data.frame(result) %>% dplyr::select(-start, -end, -width, -phase), by='exon_id') %>% mutate(type='CDS') %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE);
            result <- c(result, cds_ranges_annotated)
        }
        return(result)
    }))

    # Subset columns
    result_gtf_subset <- result_gtf %>% as.data.frame %>% dplyr::select(seqnames, start, end, width, strand, source, type, score, phase, gene_id, transcript_id, exon_id, gene_name, transcript_name) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    # Export
    export(result_gtf_subset, outfile, format='gtf')
    
}

#############################################
########## 5. Merge
#############################################

merge_gtf <- function(infiles, outfile) {
    
    # Load
    suppressPackageStartupMessages(require(rtracklayer))

    # Read, concatenate and sort
    merged_gtf <- do.call('c', lapply(infiles, import)) %>% sortSeqlevels %>% sort

    # Write
    export(merged_gtf, outfile, format='gtf')

}

#############################################
########## 6. Filter
#############################################

# filter_gtf <- function(infiles, outfile, transcript_type) {
    
#     # Read GTF
#     gtf <- rtracklayer::import(infiles[1])

#     # Get transcript IDs
#     transcript_ids <- rtracklayer::import(infiles[2])$transcript_id %>% unique

#     # Filter
#     if (transcript_type == 'all') {
#         filtered_transcript_ids <- transcript_ids
#     } else if (transcript_type == 'novel') {
#         filtered_transcript_ids <- transcript_ids[grepl('TALON', transcript_ids)]
#     }

#     # Filter GTF
#     gtf_filtered <- gtf[gtf$transcript_id %in% filtered_transcript_ids,]
    
#     # Write
#     rtracklayer::export(gtf_filtered, outfile, format='gtf')

# }

#######################################################
#######################################################
########## S7. Pfam
#######################################################
#######################################################

#############################################
########## 1. Translate
#############################################

translate_orfs <- function(infiles, outfile) {
    
    # Library
    suppressPackageStartupMessages(require(Biostrings))

    # Read ORFs
    orf_dataframe <- fread(infiles[1])

    # Read FASTA
    nucleotide_fasta <- readDNAStringSet(infiles[2], format="fasta")

    # Fix names
    names(nucleotide_fasta) <- gsub('(.*ORF_.).*', '\\1', names(nucleotide_fasta))

    # Subset
    nucleotide_fasta <- nucleotide_fasta[orf_dataframe$ID]

    # Translate
    aa_fasta <- translate(nucleotide_fasta, if.fuzzy.codon='solve')

    # Write
    writeXStringSet(aa_fasta, file=outfile)

}

#######################################################
#######################################################
########## S8. RepeatMasker
#######################################################
#######################################################

#############################################
########## 3. Merge
#############################################

merge_repeatmasker <- function(infiles, outfile) {

    # Read
    repeatmasker_dataframe <- lapply(infiles, function(x) {

        # Read
        df_string <- readChar(x, file.info(x)$size)

        # Replace asterisk
        df_string <- gsub(' *', '*', df_string, fixed=TRUE)

        # Read
        con <- textConnection(df_string)
        df <- read.table(con, fill=TRUE, skip=2)
        close(con)

        # Add column names
        colnames(df) <- c('sw_score', 'pct_div', 'pct_del', 'pct_ins', 'query_sequence', 'query_begin', 'query_end', 'query_left', 'match', 'matching_repeat', 'repeat_class', 'repeat_begin', 'repeat_end', 'repeat_left', 'ID')

        # Return
        return(df)

    }) %>% bind_rows

    # Write
    fwrite(repeatmasker_dataframe, file=outfile, sep='\t', row.names=FALSE, quote=FALSE)

}

#######################################################
#######################################################
########## S9. Evolutionary conservation
#######################################################
#######################################################

#############################################
########## 1. Convert GTF
#############################################

gtf_to_bed <- function(infile, outfile, feature_type) {
    
    # Read
    gtf <- rtracklayer::readGFF(infile)

    # Exon
    if (feature_type == 'exon') {
        bed_dataframe <- gtf %>% filter(type=='exon') %>% mutate(seqid=paste0('chr', seqid), start=format(start-1, scientific=FALSE, trim=TRUE), end=format(end, scientific=FALSE, trim=TRUE), score=0, transcript_exon_id=paste0(transcript_id, '_exon', exon_number)) %>% select(seqid, start, end, transcript_exon_id, score, strand)
    } else if (feature_type == 'transcript') {
        bed_dataframe <- gtf %>% filter(type=='transcript') %>% mutate(seqid=paste0('chr', seqid), start=format(start-1, scientific=FALSE, trim=TRUE), end=format(end, scientific=FALSE, trim=TRUE), score=0) %>% select(seqid, start, end, transcript_id, score, strand)
    }

    # Write
    fwrite(bed_dataframe, file=outfile, sep='\t', col.names=FALSE)

}

#############################################
########## 3. Get excluded regions
#############################################

get_excluded_regions <- function(infiles, outfile) {

    # Read
    gtf <- rtracklayer::readGFF(infiles[1])

    # Read exons
    exon_dataframe <- gtf %>% filter(type=='exon') %>% mutate(seqid=paste0('chr', seqid), start=format(start-1, scientific=FALSE, trim=TRUE), end=format(end, scientific=FALSE, trim=TRUE), score=0, id=paste0(transcript_id, '_exon', exon_number), strand='*') %>% select(seqid, start, end, id, score, strand)

    # Read gaps
    gap_dataframe <- rtracklayer::import(infiles[2]) %>% as.data.frame %>% select(-width) %>% rename('seqid'='seqnames') %>% group_by(seqid) %>% mutate(score=0, nr=1:n(), id=glue('chr{seqid}_gap{nr}')) %>% select(seqid, start, end, id, score, strand)

    # Concatenate
    result_dataframe <- rbind(exon_dataframe, gap_dataframe)

    # Write
    fwrite(result_dataframe, file=outfile, sep='\t', col.names=FALSE)

}


#############################################
########## 5. Get random background
#############################################

get_shuffled_exons <- function(infiles, outfile) {
    
    # Read shuffled bed
    shuffled_transcript_bed <- fread(infiles[1], col.names = c('chr', 'transcript_start_shuffled', 'transcript_end_shuffled', 'transcript_id', 'score', 'strand'))

    # Read exon bed
    exon_bed <- fread(infiles[2], col.names = c('chr', 'start', 'end', 'transcript_exon_id', 'score', 'strand')) %>% mutate(transcript_id=gsub('(.*)_.*', '\\1', transcript_exon_id)) %>% group_by(transcript_id) %>% mutate(start_shift=start-min(start), end_shift=end-min(start))

    # Merge
    shuffle_nr <- gsub('.*shuffled(.*).bed', '\\1', infiles[1])
    shuffled_exon_bed <- shuffled_transcript_bed %>% left_join(exon_bed %>% select(-chr, -score, -strand), by='transcript_id') %>% mutate(exon_start_shuffled=transcript_start_shuffled+start_shift, exon_end_shuffled=transcript_start_shuffled+end_shift, transcript_exon_id=paste0(transcript_exon_id, '-shuffle', shuffle_nr)) %>% select(chr, exon_start_shuffled, exon_end_shuffled, transcript_exon_id, score, strand)

    # Write
    fwrite(shuffled_exon_bed, file=outfile, sep='\t', col.names=FALSE)

}

#############################################
########## 7. Merge
#############################################

merge_conservation_scores <- function(infiles, outfile) {
    
    # Read results
    result_dataframe <- lapply(infiles, function(x) 
        fread(x) %>% mutate(transcript_id=gsub('(.*?)_.*', '\\1', V1)) %>% group_by(transcript_id) %>% summarize(mean_score=mean(V6)) %>% mutate(method=gsub('.*hg38.(.*way).*', '\\1', x), shuffled=grepl('shuffled', x))
    ) %>% bind_rows

    # Write
    fwrite(result_dataframe, file=outfile, sep='\t')

}

#############################################
########## 8. Get shuffled N content
#############################################

# get_shuffled_n <- function(infiles, outfile) {
    
#     # Library
#     require(bedr)

#     # Read ranges
#     bed_ranges <- rtracklayer::import(infiles[1])

#     # Convert to dataframe
#     bed_dataframe <- bed_ranges %>% as.data.frame %>% rename('chr'='seqnames') %>% mutate(chr=gsub('chr', '', chr)) %>% arrange(chr, start)

#     # Get fasta
#     range_fasta <- get.fasta(bed_dataframe, fasta=infiles[2], check.chr=FALSE)

#     # Get counts
#     count_dataframe <- apply(range_fasta, 1, function(x) {
#         nucleotide_counts <- table(strsplit(toupper(x['sequence']), "")[[1]])
#         suppressWarnings(data.frame(coordinates=x['index'], nucleotide_counts))
#     }) %>% bind_rows %>% rename('nucleotide'='Var1', 'count'='Freq')

#     # Get fraction per sequence
#     fraction_dataframe <- bed_dataframe %>% mutate(coordinates=glue('{chr}:{start}-{end}'), transcript_shuffle=gsub('(.*)_exon.*?-(.*)', '\\1_\\2', name)) %>% 
#         select(coordinates, transcript_shuffle) %>% inner_join(count_dataframe, by='coordinates') %>% group_by(transcript_shuffle, nucleotide) %>% summarize(nucleotide_count=sum(count)) %>%
#         mutate(nucleotide_percent=prop.table(nucleotide_count)) %>% pivot_wider(id_cols = transcript_shuffle, names_from = nucleotide, values_from = nucleotide_percent, values_fill = 0)

#     # Write
#     fwrite(fraction_dataframe, file=outfile, sep='\t')

# }

#######################################################
#######################################################
########## S10. liftOver
#######################################################
#######################################################

#############################################
########## 4. Filter
#############################################

filter_genepred <- function(infiles, outfile) {

    # Read GenePred
    genepred_dataframe <- fread(infiles[1]) %>% mutate(V2=gsub('chr', '', V2))

    # Read transcripts
    transcript_ids <- rtracklayer::readGFF(infiles[2]) %>% pull(transcript_id) %>% unique

    # Filter
    result_dataframe <- genepred_dataframe %>% filter(V1 %in% transcript_ids)

    # Write
    fwrite(result_dataframe, file=outfile, sep='\t', col.names = FALSE)
}

#############################################
########## 6. Add gene ID
#############################################

add_gene_id <- function(infiles, outfile) {

    # Read GTF
    lifted_gtf <- rtracklayer::readGFF(infiles[1]) %>% mutate(source=basename(as.character(source)))

    # Transcript dataframe
    transcript_dataframe <- rtracklayer::readGFF(infiles[2]) %>% select(gene_id, transcript_id) %>% distinct

    # Merge
    merged_gtf <- lifted_gtf %>% select(-gene_id) %>% left_join(transcript_dataframe, by='transcript_id')

    # Write
    rtracklayer::export(merged_gtf, outfile, format='gtf')
}

#######################################################
#######################################################
########## S10. BLAST
#######################################################
#######################################################

#############################################
########## 2. Merge genomes
#############################################

merge_blast_genomes <- function(infile, outfile) {

    # Library
    require(seqinr)

    # Read FASTA
    fasta <- read.fasta(infile, whole.header=TRUE, forceDNAtolower=FALSE)

    # Rename
    genome <- gsub('.*/(.*).fa.gz', '\\1', infile)
    sequence_names <- paste0(genome, '_', names(fasta))

    # Write FASTA
    write.fasta(fasta, sequence_names, outfile)   

}

#############################################
########## 6. Filter BLAST
#############################################

filter_blast <- function(infile, outfile) {

    # Read
    blast_dataframe <- fread(infile, col.names = c('query_acc_ver', 'subject_acc_ver', 'pct_identity', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'))

    # # Filter (v1)
    # filtered_dataframe <- blast_dataframe %>% group_by(query_acc_ver, subject_acc_ver) %>% slice_max(order_by = alignment_length, n = 1) %>% 
    #     mutate(genome=gsub('(.*?)_.*', '\\1', subject_acc_ver)) %>% group_by(query_acc_ver, genome) %>% 
    #     slice_max(order_by = alignment_length, n = 1) %>% slice_max(order_by = pct_identity, n=1) %>% 
    #     slice_max(order_by = bit_score, n=1) %>% slice_min(order_by = evalue, n=1) %>% 
    #     slice_min(order_by = q_start, n=1) %>%
    #     arrange(query_acc_ver, -pct_identity) %>% filter(evalue < 0.05)

    # Filter (v2)
    filtered_dataframe <- blast_dataframe %>% filter(pct_identity > 75 & alignment_length > 100 & evalue < 0.05) %>% 
        mutate(genome=gsub('(.*?)_.*', '\\1', subject_acc_ver)) %>% select(query_acc_ver, genome) %>% distinct
    
    # Write
    fwrite(filtered_dataframe, file=outfile, sep='\t')

}

#############################################
########## 6. Filter BLAST
#############################################

count_filtered_blast <- function(infile, outfile) {

    # Read
    blast_dataframe <- fread(infile, col.names = c('query_acc_ver', 'subject_acc_ver', 'pct_identity', 'alignment_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score'))

    # # Filter (v1)
    # filtered_dataframe <- blast_dataframe %>% group_by(query_acc_ver, subject_acc_ver) %>% slice_max(order_by = alignment_length, n = 1) %>% 
    #     mutate(genome=gsub('(.*?)_.*', '\\1', subject_acc_ver)) %>% group_by(query_acc_ver, genome) %>% 
    #     slice_max(order_by = alignment_length, n = 1) %>% slice_max(order_by = pct_identity, n=1) %>% 
    #     slice_max(order_by = bit_score, n=1) %>% slice_min(order_by = evalue, n=1) %>% 
    #     slice_min(order_by = q_start, n=1) %>%
    #     arrange(query_acc_ver, -pct_identity) %>% filter(evalue < 0.05)

    # Filter (v2)
    filtered_dataframe <- lapply(seq(75, 100, by=5), function(x) {
        message(x)
        blast_dataframe %>% filter(pct_identity > x & alignment_length > 100 & evalue < 0.05) %>% mutate(genome=gsub('(.*?)_.*', '\\1', subject_acc_ver)) %>% select(query_acc_ver, genome) %>% distinct %>% mutate(pct_identity=x)
    }) %>% bind_rows
    
    # Write
    fwrite(filtered_dataframe, file=outfile, sep='\t')

}

#######################################################
#######################################################
########## Summary
#######################################################
#######################################################

#############################################
########## 1. Create
#############################################

get_transcript_summary <- function(infiles, outfile) {

    # Read TALON summary
    talon_dataframe <- fread(infiles[1])

    # Read CPAT
    cpat_dataframe <- fread(infiles[2]) %>% rename('transcript_id'='Sequence Name', 'ORF_length'='ORF size', 'coding_probability'='Coding Probability', 'coding'='Coding Label') %>% select(transcript_id, ORF_length, coding_probability, coding)

    # Read Pfam
    pfam_dataframe <- fread(infiles[3]) %>% rename('transcript_id'='seq_id') %>% group_by(transcript_id) %>% summarize(nr_domains=length(hmm_name), domains=paste0(unique(hmm_name), collapse=','))

    # Read RepeatMasker
    repeat_dataframe <- fread(infiles[4]) %>% rename('transcript_id'='query_sequence') %>% group_by(transcript_id) %>% summarize(nr_repeats=length(matching_repeat), repeats=paste0(unique(matching_repeat), collapse=','))

    # Merge
    merged_dataframe <- list(talon_dataframe, cpat_dataframe, pfam_dataframe, repeat_dataframe) %>% purrr::reduce(left_join, by = "transcript_id") %>% replace_na(list(nr_repeats=0, nr_domains=0, coding='no'))

    # Write
    fwrite(merged_dataframe, file=outfile, sep='\t')
    
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################