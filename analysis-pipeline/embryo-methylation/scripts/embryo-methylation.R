#################################################################
#################################################################
############### Embryo methylation - R Support #################
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
########## S4. Methylation estimates
#######################################################
#######################################################

#############################################
########## 1. Tile
#############################################

tile_methylation <- function(infile, outfile) {

    # Library
    require(methylKit)
    
    # Read results
    bismark_results <- methRead(infile, sample.id=gsub('.*/results/(.*)/methylation/.*', '\\1', infile), assembly='hg38', pipeline='bismarkCytosineReport', mincov = 1, context='CpG')
    
    # Tile
    bismark_results_tiled <- tileMethylCounts(bismark_results, win.size=100, step.size=100, cov.bases = 1, mc.cores=3)

    # Write
    fwrite(bismark_results_tiled, file=outfile, sep='\t') 

}

#############################################
########## 2. Average
#############################################

average_tiled_methylation <- function(infiles, outfile, chromsize_file) {

    print(infiles)
    print(chromsize_file)

    # Library
    require(GenomeInfoDb)
    require(GenomicRanges)
    
    # Average methylation across replicates
    methylation_dataframe <- lapply(infiles, function(x) {
        fread(x) %>% mutate(methylation_level=numCs/(numCs+numTs), sample=gsub('.*/(.*)-.*', '\\1', x))
    }) %>% bind_rows %>% group_by(chr, start, end, strand) %>% summarize(average_methylation_level=mean(methylation_level, na.rm=TRUE))

    # Read chromosome sizes
    chromsize_dataframe <- fread(chromsize_file, col.names = c('seqnames', 'seqlengths'))

    # Convert to seqinfo
    chromsizes <- Seqinfo(seqnames=chromsize_dataframe$seqnames, seqlengths = chromsize_dataframe$seqlengths, genome='hg38')

    # Convert to ranges
    methylation_ranges <- methylation_dataframe %>% dplyr::rename('score'='average_methylation_level') %>% makeGRangesFromDataFrame(seqinfo=chromsizes, keep.extra.columns = TRUE) # %>% mutate(start=start+50, end=start)

    # Export
    rtracklayer::export(methylation_ranges, outfile, format='BigWig')

}

#######################################################
#######################################################
########## S5. TSS Average
#######################################################
#######################################################

#############################################
########## 1. TSS
#############################################

get_tss_methylation <- function(infiles, outfile) {
    
    # Library
    require(methylKit)

    # Read results
    bismark_results <- methRead(infiles[1], sample.id=gsub('.*/results/(.*)/methylation/.*', '\\1', infiles[1]), assembly='hg38', pipeline='bismarkCytosineReport', mincov = 1, context='CpG')

    # Read TSS ranges
    tss_ranges <- rtracklayer::import(infiles[2])# %>% sort %>% head(100)
    strand(tss_ranges) <- '*'

    # Get counts
    region_counts <- regionCounts(bismark_results, tss_ranges)

    # Merge
    merged_dataframe <- tss_ranges %>% as.data.frame %>% dplyr::rename('chr'='seqnames') %>% left_join(region_counts, by=c('chr', 'start', 'end', 'strand')) %>% mutate(percent_methylation=numCs/(numCs+numTs)) %>% dplyr::select(-score)

    # Write
    fwrite(merged_dataframe, file=outfile, sep='\t')

}

# #############################################
# ########## 1. Smoothen
# #############################################

# smoothen_methylation <- function(infiles, outfile) {

#     # Library
#     require(bsseq)
    
#     # Create metadata
#     sample_names <- gsub('.*results/(.*?)/.*', '\\1', infiles)
#     metadata_dataframe <- data.frame(row.names=sample_names, sample_name=sample_names, cell_type=gsub('(human_.*?)_.*', '\\1', sample_names))

#     # Read bismark results
#     bismark_results <- read.bismark(infiles, colData=metadata_dataframe)

#     # Smooth
#     bs_smooth <- BSmooth(orderBSseq(bismark_results), BPPARAM = BiocParallel::MulticoreParam(workers = 1))
    
#     # Save
#     save(bs_smooth, file=outfile)

# }

#######################################################
#######################################################
########## S3. Bismark
#######################################################
#######################################################

#############################################
########## 1. Average replicates
#############################################

get_average_methylation <- function(infiles, outfile) {

    # Read data
    bedgraph_dataframe <- lapply(infiles, function(x) {
        fread(x, skip = 1, col.names = c('chr', 'start', 'end', 'pct')) %>% mutate(sample=basename(x))
    }) %>% bind_rows %>% pivot_wider(id_cols = c(chr, start, end), names_from = sample, values_from = pct, values_fill = 0)

    # Get average
    bedgraph_dataframe$average <- rowMeans(bedgraph_dataframe[,4:ncol(bedgraph_dataframe)])

    # Subset
    average_dataframe <- bedgraph_dataframe[,c('chr', 'start', 'end', 'average')]

    # Write
    fwrite(average_dataframe, file=outfile, sep='\t', col.names = FALSE)

}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################