# Analysis code for isoform-resolved human preimplantation embryo transcriptome

This repository contains the R and Python code used to process the short- and long-read RNA-Seq data in order to generate the isoform-resolved transcriptome, and integrate published multi-omics datasets for orthogonal isoform validation.

## Overview
The code consists of five separate modules, each containing pipelines built for a specific analysis on distinct data types and sources used in the manuscript:

- **[embryo-isoseq](embryo-isoseq)** contains code to pre-process the PacBio IsoSeq data, integrated with Illumina short-read RNA-Seq data for splice junction validation, to generate the isoform-resolved transcriptome and predict biological properties of the transcripts contained in it.
- **[embryo-illumina](embryo-illumina)** contains code to quantify the expression of genes and isoforms in the novel transcriptome using Illumina short-read RNA-Seq data, and leverage these to compute differential expression, alternative splicing, gene co-expression networks and more, throughout human preimplantation development.
- **[geo-illumina](geo-illumina)** contains code to process publicly available short-read RNA-Seq datasets from  human and non-human primate embryo studies to validate expression of novel isoforms and genes in the transcriptome.
- **[embryo-atacseq](embryo-atacseq)** contains code to process publicly available ATAC-Seq data ([Liu et al., Nature Communications 2019](https://www.nature.com/articles/s41467-018-08244-0)) to investigate chromatin accessibility at novel isoform and gene TSSs.
- **[embryo-chipseq](embryo-chipseq)** contains code to process publicly available CU&RUN data ([Xia et al., Science 2019](https://www.science.org/doi/10.1126/science.aaw5118) to investigate H3K4me3 and H3K27ac marks at novel isoform and gene TSSs.
- **[embryo-methylation](embryo-methylation)** contains code to process publicly available RRBS data ([Guo et al., Nature 2014](https://www.nature.com/articles/nature13544) to investigate H3K4me3 and H3K27ac marks at novel isoform and gene TSSs.