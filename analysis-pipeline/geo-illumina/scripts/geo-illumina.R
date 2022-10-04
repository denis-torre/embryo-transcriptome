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
#######################################################
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
########## S7. Single-cell RNA-Seq data
#######################################################
#######################################################

#############################################
########## 2. Filter GTF
#############################################

filter_gtf <- function(infiles, outfile) {

    # Library
    require(GenomicRanges)

    # Read GTFs
    ensembl_gtf <- rtracklayer::import(infiles[1])
    talon_gtf <- rtracklayer::import(infiles[2])

    # Read SQANTI classification
    classification_dataframe <- fread(infiles[3])

    # Get biotypes
    gene_biotypes <- c('protein_coding', 'lncRNA', 'antisense', 'IG_LV_gene', 'IG_V_gene', 'IG_V_pseudogene', 'IG_D_gene', 'IG_J_gene', 'IG_J_pseudogene', 'IG_C_gene', 'IG_C_pseudogene', 'TR_V_gene', 'TR_V_pseudogene', 'TR_D_gene', 'TR_J_gene', 'TR_J_pseudogene', 'TR_C_gene')

    # Get genes
    gene_ids <- ensembl_gtf %>% as.data.frame %>% filter(type=='gene') %>% select(gene_id, gene_biotype) %>% distinct %>% filter(gene_biotype %in% gene_biotypes) %>% pull(gene_id)

    # Get fusion transcripts
    fusion_transcripts <- classification_dataframe %>% filter(structural_category=='fusion' & grepl('TALON', isoform)) %>% pull(isoform)

    # Get novel transcripts to add
    novel_gtf <- talon_gtf[grepl('TALON', talon_gtf$transcript_id) & !talon_gtf$transcript_id %in% fusion_transcripts,]

    # Merge and sort
    merged_gtf <- c(ensembl_gtf, novel_gtf) %>% sortSeqlevels %>% sort

    # Filter
    filtered_gtf <- merged_gtf[merged_gtf$gene_id %in% gene_ids | grepl('TALON', merged_gtf$gene_id)]

    # Export
    rtracklayer::export(filtered_gtf, outfile)

}

#############################################
########## 7. Create Seurat object
#############################################

create_seurat_object <- function(infile, outfile) {

    # Source
    require(Seurat)
    require(ggplot2)
    source('/hpc/users/torred23/pipelines/projects/early-embryo/geo-illumina/pipeline/scripts/utils.R')

    # Read results
    starsolo_matrix <- ReadSTARsolo(infile)

    # Get project name
    project_name <- gsub('.*/(.*?)-Solo.out/.*', '\\1', infile)

    # Create object
    seurat_object <- CreateSeuratObject(counts = starsolo_matrix, project = project_name, min.cells = 0, min.features = 0)

    # Add mitochondrial
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

    # Get dataframe
    metadata_dataframe <- data.frame(nCount_RNA=seurat_object$nCount_RNA, nFeature_RNA=seurat_object$nFeature_RNA, percent.mt=seurat_object$percent.mt)

    # Pivot
    value_dataframe <- metadata_dataframe %>% rownames_to_column('barcode') %>% pivot_longer(-barcode)

    # Get colors
    feature_colors <- setNames(RColorBrewer::brewer.pal(name='Set1', n = 3), unique(value_dataframe$name))

    # Plot
    plot_file <- gsub('rda', 'pdf', outfile)
    pdf(plot_file, height=5, width=15)
    gp <- ggplot(value_dataframe, aes(x=value, fill=name)) +
        geom_histogram(bins=300) +
        scale_fill_manual(values=feature_colors) +
        guides(fill=FALSE) +
        facet_wrap(~name, nrow=1, scales='free') +
        labs(x='Value', y='Count', title=glue('QC values for {project_name}'))
    print(gp)

    gp <-ggplot(metadata_dataframe, aes(x=nFeature_RNA)) +
        geom_histogram(bins=500, fill=feature_colors['nFeature_RNA']) +
        scale_fill_manual(values=feature_colors) +
        scale_x_continuous(breaks=seq(0, max(metadata_dataframe$nFeature_RNA), by=200), labels=formatC) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        labs(y='Count')
    print(gp)

    gp <- ggplot(metadata_dataframe, aes(x=percent.mt)) +
        geom_histogram(bins=500, fill=feature_colors['percent.mt']) +
        scale_x_continuous(lim=c(0, 25), breaks=seq(0, 25, by=1))+
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        labs(y='Count')
    print(gp)

    gp1 <- ggplot(metadata_dataframe, aes(x=nCount_RNA, y=percent.mt)) +
    geom_point(size=1, alpha=0.1) #+
    gp2 <- ggplot(metadata_dataframe, aes(x=nCount_RNA, y=nFeature_RNA)) +
        geom_point(size=1, alpha=0.1)
    p <- cowplot::plot_grid(gp1, gp2, nrow=1)
    print(p)
    dev.off()

    # Save
    save(seurat_object, file=outfile)
    
}

#############################################
########## 7. Create Seurat object
#############################################

filter_seurat_object <- function(infile, outfile) {

    # Library
    require(Seurat)

    # Load
    load(infile)

    # Filter 
    if (grepl('yu', infile)) {
        seurat_object <- subset(seurat_object, subset = nFeature_RNA > 1300 & percent.mt < 15)
    } else if (grepl('taub', infile)) {

        # Get percent largest gene
        seurat_object$percent.largest.gene <- apply(seurat_object@assays$RNA@counts, 2, function(x)(100*max(x))/sum(x))

        # Filter
        seurat_object <- subset(seurat_object, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 10 & percent.largest.gene < 20)
        
        # Normalize
        seurat_object <- NormalizeData(seurat_object, normalization.method = "CLR")

        # Find variable features
        seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 500)

        # Scale
        seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

        # Run PCA
        seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose=FALSE)

        # Find neighbors and clusters
        dims <- 1:30
        seurat_object <- FindNeighbors(seurat_object, dims = dims)
        seurat_object <- FindClusters(seurat_object, resolution = 0.5, verbose=FALSE)
        seurat_object <- RunUMAP(seurat_object, dims = dims, verbose=FALSE)

        # Load clusters
        load('arion/illumina/s04-expression.dir/human/human-novel_gene_clusters.rda')

        # Get cluster names
        cluster_dataframe <- cluster_results$cluster %>% as.data.frame %>% mutate(cluster=paste0('cluster_', .)) %>% select(-.) %>% rownames_to_column('gene_id')  %>% mutate(Cluster=recode_factor(cluster, 'cluster_3'='early_cluster', 'cluster_2'='2C_cluster', 'cluster_5'='4C_cluster', 'cluster_4'='8C_cluster', 'cluster_1'='late_cluster'))

        # Convert to list
        novel_gene_clusters <- split(cluster_dataframe$gene_id, cluster_dataframe$Cluster)

        # Get module scores
        seurat_object <- AddModuleScore(seurat_object, features=novel_gene_clusters, name = names(novel_gene_clusters))

        # Save
        save(seurat_object, file=outfile)

        # # Normalize
        # seurat_object <- SCTransform(seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
        
        # # Run PCA
        # seurat_object <- RunPCA(seurat_object, verbose=FALSE)


    } else if (grepl('iblastoid', infile)) {
        
        # Read cell types
        metadata_dataframe <- fread('arion/datasets/liu_iblastoid/liu_iblastoid-supplementary_table_1.csv') %>% rename('barcode'='V1') %>% mutate(barcode=gsub('(.*)-.', '\\1', barcode), group=gsub('(.*?)_.*', '\\1', original_id))

        # Filter
        # seurat_object <- subset(seurat_object, subset = nFeature_RNA > 1300 & percent.mt < 10)
        seurat_object <- subset(seurat_object, cells=metadata_dataframe$barcode)

        # Add clusters
        seurat_object <- AddMetaData(seurat_object, metadata=metadata_dataframe %>% pull(group, barcode), col.name='liu_cluster')

        # Normalize
        seurat_object <- SCTransform(seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
        
        # Run PCA
        seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), verbose=FALSE)

        # UMAP
        seurat_object <- RunUMAP(seurat_object, dims = 1:15, n.neighbors = 30, min.dist=.3, verbose=FALSE)
        
        # Plot
        umap_dataframe <- Embeddings(seurat_object[['umap']]) %>% as.data.frame %>% rownames_to_column('barcode') %>% inner_join(metadata_dataframe, by='barcode')

        # Save
        save(seurat_object, metadata_dataframe, umap_dataframe, file=outfile)

    } else if (grepl('kagawa', infile)) {

        # Filter
        seurat_object <- subset(seurat_object, subset = nFeature_RNA > 2000 & percent.mt < 15)

        # Normalize
        seurat_object <- NormalizeData(seurat_object)

        # Find variable features
        seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")

        # Scale
        seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

        # Run PCA
        seurat_object <- RunPCA(seurat_object, verbose=FALSE)

        # Find neighbors and clusters
        dims <- 1:30
        seurat_object <- FindNeighbors(seurat_object, dims = dims)
        seurat_object <- FindClusters(seurat_object, resolution = 0.05, verbose=FALSE)

        # Get clusters
        cluster_dataframe <- seurat_object[[]] %>% rownames_to_column('barcode') %>% mutate(cluster=paste0('cluster_', seurat_clusters)) %>% select(barcode, cluster)

        # Save
        save(seurat_object, file=outfile)

        # Write clusters
        fwrite(cluster_dataframe, file=gsub('.rda', '_clusters.tsv', outfile), sep='\t')

    }
    #add petropoulos, mazid, yan

    # # Normalize
    # seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

    # # Scale
    # seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

    # # Normalize
    # seurat_object <- SCTransform(seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

    # # Run PCA
    # seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

    # # # Find variable features
    # # seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")

    # # Cluster
    # dims <- 1:15
    # seurat_object <- FindNeighbors(seurat_object, dims = dims, k.param=20)
    # # seurat_object <- FindClusters(seurat_object, resolution = 0.5)

    # # Plot
    # umap_dataframe <- Embeddings(seurat_object[['umap']]) %>% as.data.frame %>% rownames_to_column('barcode') %>% inner_join(metadata_dataframe, by='barcode') %>% mutate(group=gsub('(.*?)_.*', '\\1', original_id))

    # # Save
    # save(seurat_object, file=outfile)

    # # Get gene set
    # gene_set <- list(
    #     'EPI'=c('POU5F1B', 'HIST2H2BF', 'KLF5', 'TBX3', 'ZEB2', 'WNT1', 'GFAP', 'FGF17', 'LAMC1', 'FNIP1', 'BMP2', 'KLHL2', 'PRDM14', 'BHLHA15', 'SMAD9', 'SETD1B', 'TGIF1', 'ARFRP1', 'METRN', 'FZD7', 'CER1', 'FN1', 'DIDO1', 'RLIM', 'ASCL2', 'EPHB3', 'FST', 'DPPA3', 'EOMES', 'BMP4', 'DTX4', 'POU3F1', 'CDH1', 'TBX4', 'PIGB', 'TJP3', 'KLF4', 'RASGRF1', 'PRICKLE1', 'TFCP2L1', 'FGF2', 'FZD5', 'PCGF2', 'HIST1H2BC', 'HIST1H2BD', 'FMN1', 'EED', 'FERMT2', 'CHRD', 'HS3ST5', 'AMOT', 'FGF3', 'INA', 'TRO', 'PRDM1', 'GATA4', 'B3GAT1', 'UTF1', 'MIR203A', 'TET1', 'VPS52', 'PCSK5', 'NOTCH1', 'GRIP1', 'GBX2', 'PGK1', 'SFMBT2', 'BMP8B', 'HIST1H2BA', 'SLAIN1', 'CDH2', 'OTX1', 'GDF3', 'KMT2B', 'SFRP5', 'NR5A2', 'SNAI2', 'SMAD2', 'MYOD1', 'GDF11', 'POMT1', 'TBX20', 'KDM7A', 'SMAD5', 'FGF1', 'H3F3A', 'PAK4', 'NODAL', 'ICAM5', 'POU5F1', 'TCF7L1', 'PCSK6', 'ECSIT', 'AKAP17A', 'INHBA', 'NOG', 'HOXB7', 'POU3F2', 'CDX2', 'PAX3', 'NUP133', 'DNMT3B', 'ACVR1B', 'IGF2R', 'GMNN', 'PINK1', 'FGFR2', 'BPTF', 'DLK1', 'SMAD4', 'ZSCAN4', 'VANGL2', 'TCF4', 'ALPPL2', 'CNBP', 'FOXA2', 'MYF5', 'DNMT3A', 'CTBP2', 'NANOG', 'PROM1', 'ALPP', 'LEFTY1', 'CXCR4', 'CBX3', 'TDGF1', 'RBPJL', 'E2F6', 'FZD10', 'CUL4A', 'CDX4', 'SOX1', 'CTBP1', 'PRTG', 'HOXD8', 'FOXRED1', 'DES', 'HIST1H2BF', 'ADD3', 'CDC45', 'DAB2', 'CCR4', 'FAM3B', 'IFITM3', 'OTX2', 'HOXB6', 'TGIF2', 'STAT3', 'CYTH1', 'CRB2', 'HDAC9', 'HMX2', 'DAG1', 'CHD3', 'FNIP2', 'TBX6', 'HIST3H2BB', 'SHH', 'LIFR', 'HOXB1', 'SBDS', 'LEFTY2', 'CHURC1', 'CUBN', 'CUL4B', 'RANGAP1', 'FUT4', 'DDX4', 'LAMB1', 'SOX3', 'CBX1', 'TBX5', 'TFAP2C', 'FZD8', 'KLF2', 'SOX10', 'TRIM47', 'SOX2', 'NES', 'DVL2', 'LAMA1', 'PDGFRA', 'KIF16B', 'CYP26A1', 'PLET1', 'CDT1', 'POMT2', 'PORCN', 'ROCK2', 'GINS4', 'WNT9A', 'HIST1H2BH', 'HOXA1', 'SALL4', 'RYBP', 'SHC4', 'FGF8', 'KDM5B', 'MAGEB16', 'ASTL', 'H3F3B', 'WWTR1', 'BMPR1A', 'WNT3A', 'SETDB1', 'HESX1', 'EGOT', 'PRMT5', 'TRIM28', 'ABR', 'WLS', 'FGF4', 'TENM1', 'FGFR1', 'PIK3C3', 'BIK', 'CNOT1', 'DOCK1', 'DTX1', 'HAND1', 'INHBB', 'ESRRB', 'TRIM33', 'QRSL1', 'L1TD1', 'GJB5', 'ELF5', 'NID1', 'SOX15', 'STAG3', 'DNMT3L', 'TET2', 'LHX1', 'SOX17', 'DNMT1', 'PITX2', 'GATA6', 'MNX1', 'FOXD3', 'MESDC2', 'KDM1B', 'FOXH1', 'EZH2', 'GSK3B', 'BYSL', 'GLB1', 'DKK1', 'MEOX2', 'ETS2', 'ZIC3', 'ACVR1', 'FURIN', 'CBX5', 'ALPL', 'CITED2', 'POFUT2', 'TGFB1', 'WNT3', 'HIST2H3PS2', 'ACTC1', 'FGF5', 'COPS3', 'CITED1', 'KRT8', 'TLN1', 'NANOS3', 'SIN3A'),
    #     'TE'=c('MYL6', 'IFIH1', 'SLC2A6', 'SLC1A5', 'TDGF1', 'CSH1', 'SLC2A8', 'TROAP', 'SRY', 'ESX1', 'FGFBP1', 'SLC38A6', 'SETD7', 'VCL', 'HIST2H3A', 'ISG20', 'DNAH8', 'H3F3A', 'EIF4EBP1', 'GTF2IRD1', 'CST3', 'IFI27', 'CITED1', 'GATA6', 'DNTT', 'SMARCA5', 'CDC25A', 'STAT3', 'HDAC1', 'NR0B1', 'FRAT2', 'PAG1', 'GCH1', 'CBX3', 'EZH2', 'FST', 'NRG4', 'IQGAP2', 'KLF4', 'MAP3K10', 'EFNA1', 'CTNNB1', 'CXXC1', 'STAT1', 'PRICKLE1', 'BAX', 'IGF2', 'ZFP42', 'RYBP', 'MSI2', 'DICER1', 'INS', 'RPL4', 'ASXL1', 'LPAR1', 'COX7B', 'FGF1', 'MUC15', 'LRP2', 'TGFB3', 'IGF2BP1', 'NRG2', 'ISM1', 'OTP', 'LATS1', 'DSC3', 'ST5', 'FUT4', 'ASAP1', 'SCRIB', 'TUBB2A', 'FGFR2', 'PGR', 'IQGAP1', 'NLRP4', 'YY2', 'HOXD8', 'IFNG', 'ITGA8', 'MIR93', 'ITGB1', 'KDM6B', 'ATP1B3', 'CTBS', 'SLC25A1', 'DSP', 'TBCE', 'ZP3', 'SPP1', 'ESR1', 'ALPPL2', 'GNB1L', 'RPS6', 'RSAD2', 'CUZD1', 'TGFB1', 'PXK', 'IFNAR1', 'NODAL', 'SEBOX', 'SUDS3', 'EOMES', 'CNOT3', 'SALL4', 'WDR74', 'SLC7A7', 'PRKCD', 'RLIM', 'EREG', 'ERBB4', 'FGF7', 'GATA4', 'SLC43A2', 'SLC2A12', 'AQP8', 'OCM2', 'CUL4A', 'SLC2A4', 'ROCK2', 'ACTG1', 'AHNAK', 'ITGB5', 'ZNF35', 'TUBB', 'OCLN', 'RBM19', 'DTX1', 'GCM1', 'ECM2', 'HOXB7', 'STPG1', 'SLC7A2', 'ITGB3', 'KLF8', 'H2AFZ', 'CGN', 'SNAI2', 'CDH2', 'ISG15', 'KLF5', 'RFPL1', 'RAB13', 'CLDN14', 'ATP11A', 'BYSL', 'GDF9', 'TRIM27', 'TET2', 'BMP10', 'NLRP5', 'UTF1', 'CTNNA1', 'FOLH1', 'HIST2H3C', 'EED', 'CLDN1', 'NMI', 'SIN3A', 'LGALS1', 'MAPKAP1', 'GBP2', 'RPTOR', 'SLC7A6', 'MTNR1A', 'PRDM14', 'DPPA2', 'DLL3', 'POU5F1B', 'PGK1', 'MTOR', 'FGF5', 'IGFBP1', 'SLC6A14', 'FGFR1', 'OVOL2', 'SRP14', 'JUP', 'GJD4', 'CTSL', 'TET1', 'SYNE2', 'EGFL7', 'EZR', 'NOBOX', 'ELF5', 'LEFTY1', 'TGFB2', 'LPAR3', 'DPPA4', 'TTK', 'CARM1', 'ODC1', 'DPPA5', 'ESRRB', 'DNMT3A', 'SLC7A8', 'DSG1', 'LAMB2', 'HLA', 'TBP', 'SCYL1', 'SLC38A4', 'HAND1', 'SLC16A3', 'TUBA1C', 'SETDB1', 'EPHA1', 'GJB3', 'NRG3', 'AHCTF1', 'DNMT1', 'PPP1R13L', 'BMPR1A', 'UQCR10', 'GATA3', 'IRF2', 'TCF4', 'SLC5A11', 'F11R', 'CSF2RB', 'YAF2', 'KRT18', 'SDC2', 'CDX4', 'DNMT3L', 'SOX17', 'SLC1A1', 'ITGAV', 'COX6C', 'ETS2', 'IRF1', 'ATP6V0B', 'MTNR1B', 'PARD6B', 'CDH1', 'PRR15', 'TEAD2', 'DVL3', 'TBX3', 'POU5F1', 'SFN', 'IGF1R', 'KRT8', 'OXTR', 'WNT3A', 'TFAP2C', 'PGF', 'SLC25A5', 'NANOG', 'CSH2', 'SLC7A1', 'FGF4', 'CXCL14', 'SLC16A8', 'SLC2A5', 'PSTPIP2', 'FOXD3', 'NDUFA4', 'ACVR2B', 'RPL7L1', 'ITGB8', 'TRO', 'HBEGF', 'RRP7A', 'CTNNA3', 'CHRFAM7A', 'DKKL1', 'TFAP2A', 'TJP1', 'GJA1', 'LLGL1', 'DKK3', 'AGMAT', 'AMOTL2', 'CDX1', 'PRAMEF10', 'ASH1L', 'TINAGL1', 'GLB1', 'PRSS12', 'PATZ1', 'ARID3A', 'SLC2A3', 'FN1', 'SLC6A19', 'DAB2', 'PRKCZ', 'IRF9', 'PORCN', 'MXD1', 'AMOT', 'SLC2A9', 'ATP1A4', 'GPI', 'ITGA7', 'ISM2', 'GTF2I', 'ASCL2', 'RHOU', 'SLC2A1', 'SLC7A3', 'DKK1', 'FGF10', 'HOXA10', 'MUC1', 'CFC1', 'GABRB3', 'NDUFA13', 'ALPP', 'POU2F1', 'SOX2', 'MBD3L2', 'ZSCAN4', 'ECEL1', 'HOXB6', 'HIST2H3D', 'DNMT3B', 'TEAD1', 'CDX2', 'FKBP4', 'OOEP', 'HMX2', 'GAPDH', 'ALB', 'LIFR', 'TEAD4', 'UQCRFS1', 'SLC1A4', 'SLC2A2', 'MAGT1', 'PADI6', 'TEX11', 'FGF2', 'WNT7A', 'H3F3B', 'IGF1', 'FOXJ2', 'PIKFYVE', 'SOX15', 'PRMT2', 'SLC5A1', 'TUBA4A', 'ACTB', 'IFNA1', 'AQP3', 'TJP2', 'SMARCA4', 'ETF1', 'HIST2H3PS2', 'ATP1A1', 'CUL4B', 'ROCK1P1', 'LAMA1', 'KDM8', 'GHSR', 'GFPT1', 'GATA2', 'TSPO', 'ROCK1', 'TLN1', 'CSF2', 'BMP4', 'ALPL', 'H1FOO', 'SLC7A9', 'OCM', 'DSC2', 'IL6ST', 'IGF2BP3', 'CD44', 'CTSZ'),
    #     'PE'=c('GCNT2', 'NID1', 'FGF4', 'FOXD3', 'WNT6', 'PLAUR', 'MIR17', 'TEAD4', 'DICER1', 'ZFP42', 'TBX3', 'GLB1', 'LEFTY1', 'PLS1', 'HNF4A', 'TTR', 'ASXL1', 'ELF5', 'HHEX', 'CNOT3', 'WNT9A', 'FOXH1', 'HESX1', 'FN1', 'CBX3', 'NANOG', 'EED', 'KIF16B', 'GNAI2', 'GATA6', 'EZR', 'FOXA1', 'GCNT1', 'SALL4', 'TAB1', 'GATA4', 'SETD2', 'HIST1H2BC', 'GATA5', 'LAMC1', 'KLF5', 'LAMB1', 'DSP', 'ARHGEF1', 'DAB2', 'ACTC1', 'PRR15', 'SUDS3', 'FGFR2', 'POFUT2', 'OTX2', 'CFC1', 'ROCK1P1', 'ISM1', 'PCSK6', 'ISM2', 'SOX2', 'TCF12', 'FGFR3', 'PARP2', 'CDH1', 'THBS4', 'RDX', 'AKTIP', 'FOXA2', 'GSK3B', 'FGF1', 'MYC', 'MSN', 'SUZ12', 'LIF', 'CBX5', 'POU5F1B', 'CNOT1', 'CDX2', 'ECEL1', 'ACTG1', 'FUT4', 'RSU1', 'HIST1H2BA', 'UTF1', 'NXN', 'GPI', 'GCNT4', 'S1PR2', 'OSTM1', 'PGK1', 'HIST1H2BH', 'THBD', 'FGFR1', 'NDST2', 'ADCY2', 'ACTB', 'ROCK2', 'PINX1', 'PDGFRA', 'DIDO1', 'ACVR1C', 'SRY', 'ACVR1B', 'ARID1A', 'WNT3A', 'PLCB1', 'ITGB1', 'BMP4', 'GBX2', 'KLF2', 'DVL2', 'HMX2', 'FZD2', 'EOMES', 'RGS19', 'HIST1H2BF', 'CHD1', 'IGHV1', 'FGFR4', 'TET1', 'FGF5', 'GRB2', 'SOS1', 'SOX17', 'NODAL', 'PLAT', 'LIFR', 'FOXM1', 'HIST2H2BF', 'POU5F1', 'MNX1', 'NDST1', 'UGCG', 'AFP', 'LRP2', 'TDGF1', 'DVL3', 'HIST1H2BD', 'FZD1', 'HIST3H2BB', 'MSC', 'TBCE', 'IPMK', 'TGFB1', 'HNF1B', 'EZH2', 'STAT3')
    # )

    # # Get module scores
    # seurat_object <- AddModuleScore(seurat_object, features=gene_set, name = names(gene_set))
    # UMAP - default n.neighbors = 30, min.dist=.3
    # seurat_object <- RunUMAP(seurat_object, dims = dims, n.neighbors = 30, min.dist=.3)

}


#############################################
########## 8. UMAP
#############################################

run_umap <- function(infile, outfile, additional_params) {

    # Library
    require(Seurat)

    # Load
    load(infile)

    # Get parameters
    nr_neighbors <- as.numeric(additional_params[1])
    min_dist <- as.numeric(additional_params[2])

    # Find neighbors
    dims <- 1:30
    seurat_object <- FindNeighbors(seurat_object, k.param=nr_neighbors, dims = dims, verbose = FALSE)

    # Run UMAP
    seurat_object <- RunUMAP(seurat_object, dims = dims, n.neighbors=nr_neighbors, min.dist=min_dist, verbose = FALSE)

    # Get coordinates
    umap_dataframe <- Embeddings(seurat_object[['umap']]) %>% as.data.frame %>% rownames_to_column('barcode') %>% mutate(nr_neighbors=nr_neighbors, min_dist=min_dist)

    # Write
    fwrite(umap_dataframe, outfile, sep='\t')

}

#############################################
########## 7. Merge RSEM expression
#############################################

merge_rsem <- function(infile, outfile) {

    # Read metadata
    petropoulos_dataframe <- fread('arion/datasets/petropoulos/petropoulos-samples.tsv') %>% select('Comment[ENA_RUN]', 'Characteristics[individual]') %>% rename('barcode'='Comment[ENA_RUN]', 'group'='Characteristics[individual]') %>% mutate(group=gsub('(.*?)\\..*', '\\1', group), dataset='petropoulos')

    # Read mazid samples
    mazid_sample_dataframe <- fread('arion/datasets/mazid/mazid-smartseq_samples.tsv') %>% rename('accession_id'='Experiment ID', 'barcode'='Run ID') %>% select(accession_id, barcode)

    # Read mazid
    mazid_dataframe <- fread('arion/datasets/mazid/mazid-smartseq_sample_metadata.tsv') %>% left_join(mazid_sample_dataframe, by='accession_id') %>% mutate(group=gsub('(.*?)_.*', '\\1', title)) %>% select(barcode, group) %>% mutate(dataset='mazid')

    # Read kagawa
    kagawa_dataframe <- fread('arion/geo_illumina/s07-scrna.dir/seurat/kagawa/normalized/kagawa-seurat_filtered_normalized_clusters.tsv') %>% rename('group'='cluster') %>% mutate(dataset='kagawa')#, group=recode_factor(group, 'cluster_0'='epiblast', 'cluster_1'='TE', 'cluster_2'='PrE', 'cluster_3'='TSCs', 'cluster_4'='primed PSCs'), group=as.character(group))

    # Merge
    metadata_dataframe <- rbind(petropoulos_dataframe, mazid_dataframe, kagawa_dataframe)

    # Classification
    classification_dataframe <- fread('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered-transcript_classification.tsv')
    classification_dataframe$Transcript_novelty <- factor(classification_dataframe$Transcript_novelty, levels=c('Known', 'NIC', 'NNC', 'Antisense', 'Intergenic'))

    # Read expression
    expression_dataframe <- lapply(Sys.glob(paste0(infile, '/*/*.isoforms.results')), function(x) {
        fread(x) %>% mutate(barcode=gsub('.*/(.*).isoforms.results', '\\1', x))
    }) %>% bind_rows %>% left_join(metadata_dataframe, by='barcode') %>% group_by(transcript_id, group) %>% summarize(mean_TPM=mean(TPM))

    # Result
    merged_dataframe <- classification_dataframe %>% select(gene_id, transcript_id, Transcript_novelty) %>% left_join(expression_dataframe, by='transcript_id')

    # Write
    fwrite(merged_dataframe, file=outfile, sep='\t')


}

# #############################################
# ########## 7. Create Seurat object
# #############################################

# filter_seurat_object <- function(infiles, outfile) {

#     # Library
#     require(Seurat)

#     # Add names
#     names(infiles) <- gsub('.*/(.*)-seurat.rda', '\\1', infiles)

#     # Merge liu
#     if (grepl('liu', infiles[1])) {
        
#         # Read
#         seurat_objects <- sapply(infiles, function(x) {
#             load(x)
#             seurat_object
#         }, simplify=FALSE)
        
#         # Merge
#         seurat_object <- merge(seurat_objects[[1]], y=seurat_objects[2:length(seurat_objects)], add.cell.ids=names(seurat_objects))

#         # Filter
#         seurat_object <- subset(seurat_object, subset = nFeature_RNA > 1300 & percent.mt < 15)
        
#     } else if (grepl('taub', infiles[1])) {
        
#         # Load
#         load(infiles[1])
        
#         # Filter
#         seurat_object <- subset(seurat_object, subset = nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 10)
#     }

#     # Normalize
#     seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

#     # Find variable features
#     seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst")

#     # Scale
#     seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))

#     # Run PCA
#     seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

#     # Cluster
#     dims <- 1:15
#     seurat_object <- FindNeighbors(seurat_object, dims = dims)
#     seurat_object <- FindClusters(seurat_object, resolution = 0.5)

#     # Get gene set
#     gene_set <- list(
#         'EPI'=c('POU5F1B', 'HIST2H2BF', 'KLF5', 'TBX3', 'ZEB2', 'WNT1', 'GFAP', 'FGF17', 'LAMC1', 'FNIP1', 'BMP2', 'KLHL2', 'PRDM14', 'BHLHA15', 'SMAD9', 'SETD1B', 'TGIF1', 'ARFRP1', 'METRN', 'FZD7', 'CER1', 'FN1', 'DIDO1', 'RLIM', 'ASCL2', 'EPHB3', 'FST', 'DPPA3', 'EOMES', 'BMP4', 'DTX4', 'POU3F1', 'CDH1', 'TBX4', 'PIGB', 'TJP3', 'KLF4', 'RASGRF1', 'PRICKLE1', 'TFCP2L1', 'FGF2', 'FZD5', 'PCGF2', 'HIST1H2BC', 'HIST1H2BD', 'FMN1', 'EED', 'FERMT2', 'CHRD', 'HS3ST5', 'AMOT', 'FGF3', 'INA', 'TRO', 'PRDM1', 'GATA4', 'B3GAT1', 'UTF1', 'MIR203A', 'TET1', 'VPS52', 'PCSK5', 'NOTCH1', 'GRIP1', 'GBX2', 'PGK1', 'SFMBT2', 'BMP8B', 'HIST1H2BA', 'SLAIN1', 'CDH2', 'OTX1', 'GDF3', 'KMT2B', 'SFRP5', 'NR5A2', 'SNAI2', 'SMAD2', 'MYOD1', 'GDF11', 'POMT1', 'TBX20', 'KDM7A', 'SMAD5', 'FGF1', 'H3F3A', 'PAK4', 'NODAL', 'ICAM5', 'POU5F1', 'TCF7L1', 'PCSK6', 'ECSIT', 'AKAP17A', 'INHBA', 'NOG', 'HOXB7', 'POU3F2', 'CDX2', 'PAX3', 'NUP133', 'DNMT3B', 'ACVR1B', 'IGF2R', 'GMNN', 'PINK1', 'FGFR2', 'BPTF', 'DLK1', 'SMAD4', 'ZSCAN4', 'VANGL2', 'TCF4', 'ALPPL2', 'CNBP', 'FOXA2', 'MYF5', 'DNMT3A', 'CTBP2', 'NANOG', 'PROM1', 'ALPP', 'LEFTY1', 'CXCR4', 'CBX3', 'TDGF1', 'RBPJL', 'E2F6', 'FZD10', 'CUL4A', 'CDX4', 'SOX1', 'CTBP1', 'PRTG', 'HOXD8', 'FOXRED1', 'DES', 'HIST1H2BF', 'ADD3', 'CDC45', 'DAB2', 'CCR4', 'FAM3B', 'IFITM3', 'OTX2', 'HOXB6', 'TGIF2', 'STAT3', 'CYTH1', 'CRB2', 'HDAC9', 'HMX2', 'DAG1', 'CHD3', 'FNIP2', 'TBX6', 'HIST3H2BB', 'SHH', 'LIFR', 'HOXB1', 'SBDS', 'LEFTY2', 'CHURC1', 'CUBN', 'CUL4B', 'RANGAP1', 'FUT4', 'DDX4', 'LAMB1', 'SOX3', 'CBX1', 'TBX5', 'TFAP2C', 'FZD8', 'KLF2', 'SOX10', 'TRIM47', 'SOX2', 'NES', 'DVL2', 'LAMA1', 'PDGFRA', 'KIF16B', 'CYP26A1', 'PLET1', 'CDT1', 'POMT2', 'PORCN', 'ROCK2', 'GINS4', 'WNT9A', 'HIST1H2BH', 'HOXA1', 'SALL4', 'RYBP', 'SHC4', 'FGF8', 'KDM5B', 'MAGEB16', 'ASTL', 'H3F3B', 'WWTR1', 'BMPR1A', 'WNT3A', 'SETDB1', 'HESX1', 'EGOT', 'PRMT5', 'TRIM28', 'ABR', 'WLS', 'FGF4', 'TENM1', 'FGFR1', 'PIK3C3', 'BIK', 'CNOT1', 'DOCK1', 'DTX1', 'HAND1', 'INHBB', 'ESRRB', 'TRIM33', 'QRSL1', 'L1TD1', 'GJB5', 'ELF5', 'NID1', 'SOX15', 'STAG3', 'DNMT3L', 'TET2', 'LHX1', 'SOX17', 'DNMT1', 'PITX2', 'GATA6', 'MNX1', 'FOXD3', 'MESDC2', 'KDM1B', 'FOXH1', 'EZH2', 'GSK3B', 'BYSL', 'GLB1', 'DKK1', 'MEOX2', 'ETS2', 'ZIC3', 'ACVR1', 'FURIN', 'CBX5', 'ALPL', 'CITED2', 'POFUT2', 'TGFB1', 'WNT3', 'HIST2H3PS2', 'ACTC1', 'FGF5', 'COPS3', 'CITED1', 'KRT8', 'TLN1', 'NANOS3', 'SIN3A'),
#         'TE'=c('MYL6', 'IFIH1', 'SLC2A6', 'SLC1A5', 'TDGF1', 'CSH1', 'SLC2A8', 'TROAP', 'SRY', 'ESX1', 'FGFBP1', 'SLC38A6', 'SETD7', 'VCL', 'HIST2H3A', 'ISG20', 'DNAH8', 'H3F3A', 'EIF4EBP1', 'GTF2IRD1', 'CST3', 'IFI27', 'CITED1', 'GATA6', 'DNTT', 'SMARCA5', 'CDC25A', 'STAT3', 'HDAC1', 'NR0B1', 'FRAT2', 'PAG1', 'GCH1', 'CBX3', 'EZH2', 'FST', 'NRG4', 'IQGAP2', 'KLF4', 'MAP3K10', 'EFNA1', 'CTNNB1', 'CXXC1', 'STAT1', 'PRICKLE1', 'BAX', 'IGF2', 'ZFP42', 'RYBP', 'MSI2', 'DICER1', 'INS', 'RPL4', 'ASXL1', 'LPAR1', 'COX7B', 'FGF1', 'MUC15', 'LRP2', 'TGFB3', 'IGF2BP1', 'NRG2', 'ISM1', 'OTP', 'LATS1', 'DSC3', 'ST5', 'FUT4', 'ASAP1', 'SCRIB', 'TUBB2A', 'FGFR2', 'PGR', 'IQGAP1', 'NLRP4', 'YY2', 'HOXD8', 'IFNG', 'ITGA8', 'MIR93', 'ITGB1', 'KDM6B', 'ATP1B3', 'CTBS', 'SLC25A1', 'DSP', 'TBCE', 'ZP3', 'SPP1', 'ESR1', 'ALPPL2', 'GNB1L', 'RPS6', 'RSAD2', 'CUZD1', 'TGFB1', 'PXK', 'IFNAR1', 'NODAL', 'SEBOX', 'SUDS3', 'EOMES', 'CNOT3', 'SALL4', 'WDR74', 'SLC7A7', 'PRKCD', 'RLIM', 'EREG', 'ERBB4', 'FGF7', 'GATA4', 'SLC43A2', 'SLC2A12', 'AQP8', 'OCM2', 'CUL4A', 'SLC2A4', 'ROCK2', 'ACTG1', 'AHNAK', 'ITGB5', 'ZNF35', 'TUBB', 'OCLN', 'RBM19', 'DTX1', 'GCM1', 'ECM2', 'HOXB7', 'STPG1', 'SLC7A2', 'ITGB3', 'KLF8', 'H2AFZ', 'CGN', 'SNAI2', 'CDH2', 'ISG15', 'KLF5', 'RFPL1', 'RAB13', 'CLDN14', 'ATP11A', 'BYSL', 'GDF9', 'TRIM27', 'TET2', 'BMP10', 'NLRP5', 'UTF1', 'CTNNA1', 'FOLH1', 'HIST2H3C', 'EED', 'CLDN1', 'NMI', 'SIN3A', 'LGALS1', 'MAPKAP1', 'GBP2', 'RPTOR', 'SLC7A6', 'MTNR1A', 'PRDM14', 'DPPA2', 'DLL3', 'POU5F1B', 'PGK1', 'MTOR', 'FGF5', 'IGFBP1', 'SLC6A14', 'FGFR1', 'OVOL2', 'SRP14', 'JUP', 'GJD4', 'CTSL', 'TET1', 'SYNE2', 'EGFL7', 'EZR', 'NOBOX', 'ELF5', 'LEFTY1', 'TGFB2', 'LPAR3', 'DPPA4', 'TTK', 'CARM1', 'ODC1', 'DPPA5', 'ESRRB', 'DNMT3A', 'SLC7A8', 'DSG1', 'LAMB2', 'HLA', 'TBP', 'SCYL1', 'SLC38A4', 'HAND1', 'SLC16A3', 'TUBA1C', 'SETDB1', 'EPHA1', 'GJB3', 'NRG3', 'AHCTF1', 'DNMT1', 'PPP1R13L', 'BMPR1A', 'UQCR10', 'GATA3', 'IRF2', 'TCF4', 'SLC5A11', 'F11R', 'CSF2RB', 'YAF2', 'KRT18', 'SDC2', 'CDX4', 'DNMT3L', 'SOX17', 'SLC1A1', 'ITGAV', 'COX6C', 'ETS2', 'IRF1', 'ATP6V0B', 'MTNR1B', 'PARD6B', 'CDH1', 'PRR15', 'TEAD2', 'DVL3', 'TBX3', 'POU5F1', 'SFN', 'IGF1R', 'KRT8', 'OXTR', 'WNT3A', 'TFAP2C', 'PGF', 'SLC25A5', 'NANOG', 'CSH2', 'SLC7A1', 'FGF4', 'CXCL14', 'SLC16A8', 'SLC2A5', 'PSTPIP2', 'FOXD3', 'NDUFA4', 'ACVR2B', 'RPL7L1', 'ITGB8', 'TRO', 'HBEGF', 'RRP7A', 'CTNNA3', 'CHRFAM7A', 'DKKL1', 'TFAP2A', 'TJP1', 'GJA1', 'LLGL1', 'DKK3', 'AGMAT', 'AMOTL2', 'CDX1', 'PRAMEF10', 'ASH1L', 'TINAGL1', 'GLB1', 'PRSS12', 'PATZ1', 'ARID3A', 'SLC2A3', 'FN1', 'SLC6A19', 'DAB2', 'PRKCZ', 'IRF9', 'PORCN', 'MXD1', 'AMOT', 'SLC2A9', 'ATP1A4', 'GPI', 'ITGA7', 'ISM2', 'GTF2I', 'ASCL2', 'RHOU', 'SLC2A1', 'SLC7A3', 'DKK1', 'FGF10', 'HOXA10', 'MUC1', 'CFC1', 'GABRB3', 'NDUFA13', 'ALPP', 'POU2F1', 'SOX2', 'MBD3L2', 'ZSCAN4', 'ECEL1', 'HOXB6', 'HIST2H3D', 'DNMT3B', 'TEAD1', 'CDX2', 'FKBP4', 'OOEP', 'HMX2', 'GAPDH', 'ALB', 'LIFR', 'TEAD4', 'UQCRFS1', 'SLC1A4', 'SLC2A2', 'MAGT1', 'PADI6', 'TEX11', 'FGF2', 'WNT7A', 'H3F3B', 'IGF1', 'FOXJ2', 'PIKFYVE', 'SOX15', 'PRMT2', 'SLC5A1', 'TUBA4A', 'ACTB', 'IFNA1', 'AQP3', 'TJP2', 'SMARCA4', 'ETF1', 'HIST2H3PS2', 'ATP1A1', 'CUL4B', 'ROCK1P1', 'LAMA1', 'KDM8', 'GHSR', 'GFPT1', 'GATA2', 'TSPO', 'ROCK1', 'TLN1', 'CSF2', 'BMP4', 'ALPL', 'H1FOO', 'SLC7A9', 'OCM', 'DSC2', 'IL6ST', 'IGF2BP3', 'CD44', 'CTSZ'),
#         'PE'=c('GCNT2', 'NID1', 'FGF4', 'FOXD3', 'WNT6', 'PLAUR', 'MIR17', 'TEAD4', 'DICER1', 'ZFP42', 'TBX3', 'GLB1', 'LEFTY1', 'PLS1', 'HNF4A', 'TTR', 'ASXL1', 'ELF5', 'HHEX', 'CNOT3', 'WNT9A', 'FOXH1', 'HESX1', 'FN1', 'CBX3', 'NANOG', 'EED', 'KIF16B', 'GNAI2', 'GATA6', 'EZR', 'FOXA1', 'GCNT1', 'SALL4', 'TAB1', 'GATA4', 'SETD2', 'HIST1H2BC', 'GATA5', 'LAMC1', 'KLF5', 'LAMB1', 'DSP', 'ARHGEF1', 'DAB2', 'ACTC1', 'PRR15', 'SUDS3', 'FGFR2', 'POFUT2', 'OTX2', 'CFC1', 'ROCK1P1', 'ISM1', 'PCSK6', 'ISM2', 'SOX2', 'TCF12', 'FGFR3', 'PARP2', 'CDH1', 'THBS4', 'RDX', 'AKTIP', 'FOXA2', 'GSK3B', 'FGF1', 'MYC', 'MSN', 'SUZ12', 'LIF', 'CBX5', 'POU5F1B', 'CNOT1', 'CDX2', 'ECEL1', 'ACTG1', 'FUT4', 'RSU1', 'HIST1H2BA', 'UTF1', 'NXN', 'GPI', 'GCNT4', 'S1PR2', 'OSTM1', 'PGK1', 'HIST1H2BH', 'THBD', 'FGFR1', 'NDST2', 'ADCY2', 'ACTB', 'ROCK2', 'PINX1', 'PDGFRA', 'DIDO1', 'ACVR1C', 'SRY', 'ACVR1B', 'ARID1A', 'WNT3A', 'PLCB1', 'ITGB1', 'BMP4', 'GBX2', 'KLF2', 'DVL2', 'HMX2', 'FZD2', 'EOMES', 'RGS19', 'HIST1H2BF', 'CHD1', 'IGHV1', 'FGFR4', 'TET1', 'FGF5', 'GRB2', 'SOS1', 'SOX17', 'NODAL', 'PLAT', 'LIFR', 'FOXM1', 'HIST2H2BF', 'POU5F1', 'MNX1', 'NDST1', 'UGCG', 'AFP', 'LRP2', 'TDGF1', 'DVL3', 'HIST1H2BD', 'FZD1', 'HIST3H2BB', 'MSC', 'TBCE', 'IPMK', 'TGFB1', 'HNF1B', 'EZH2', 'STAT3')
#     )

#     # Get module scores
#     seurat_object <- AddModuleScore(seurat_object, features=gene_set, name = names(gene_set))

#     # Save
#     save(seurat_object, file=outfile)

#     # UMAP - default n.neighbors = 30, min.dist=.3
#     # seurat_object <- RunUMAP(seurat_object, dims = dims, n.neighbors = 30, min.dist=.3)

# }

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################