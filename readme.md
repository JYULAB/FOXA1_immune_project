# Code for FOXA1 Paper
Citation:

The analyses were conducted in a High-Performance Computing (HPC) environment using a combination of Bash scripting, R (v4.0.3), and Python (v3.6 or later). Detailed software versions are described in the Methods section of the paper.

The scripts are organized by technology or analysis type. Each script's expected output corresponds to specific figures in the paper.
    
### Spatial_Transcriptomics

#### Visium
- [**TESLA_analyses.py**](/Spatial_Transcriptomics/Visium_analyses/TESLA_analyses.py): Tumor Edge Structure and Lymphocyte multi-level Annotation Analysis, Fig 5

### scRNA-seq

- [**scRNA-seq_processing.R**](scRNA-seq/scRNA-seq_processing.R): Process scRNA-Seq samples, Fig 2, 3
- [**scRNA-seq_integration.R**](scRNA-seq/scRNA-seq_integration.R): Integrate and annotate scRNA-Seq samples, Fig 2, 3
- [**10x_scRNAseq_analysis_whole_prostate.R**](scRNA-seq/10x_scRNAseq_analysis_whole_prostate.R): Analyses of the 10x scRNA-seq data with seurat, Fig 2, 3
- [**trajectory_analysis.R**](scRNA-seq/trajectory_analysis.R): Trajectory analysis of epithelial populations Fig 2
- [**CellChat_analysis.R**](scRNA-seq/CellChat_analysis.R): Analysis of 10x whole tumor scRNA-seq data with CellChat, Fig 6
- [**Pip-seq_processing_processing.R**](scRNA-seq/Pip-seq_processing_processing.R): Processing, integrating and annotation of PIP-Seq data, Fig 3
- [**pipseq_analysis_CD45.R**](scRNA-seq/pipseq_analysis_CD45.R): Analyses of the CD45+ sorted PIP-Seq data, Fig 4
- [**scRNAseq_analysis_human.R**](scRNA-seq/scRNAseq_analysis_human.R): Analyse of published human PCa GSE181294 scRNA-seq data, Fig 7

### ChIP-seq

- [**mapping.sh**](mapping.sh): map ChIP-seq samples, Fig 3
- [**spike_in_downsampling.sh**](spike_in_downsampling.sh): Process spike-in ChIP-seq samples, Fig 3
- [**peakcalling.sh**](peakcalling.sh): call peaks of ChIP-seq samples, Fig 3
- [**heatmap.sh**](heatmap.sh): draw heatmap of ChIP-seq samples, Fig 3
- [**venn.R**](venn.R): R code for venn diagram, Fig 3
- [**figure3.R**](figure_3.R): R code for figure, Fig 3
