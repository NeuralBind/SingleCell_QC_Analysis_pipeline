# Single-Cell RNA-seq Analysis Pipeline

🧬 A comprehensive end-to-end pipeline for single-cell RNA sequencing data analysis, from quality control to cell-cell communication inference.

## 🔬 Overview

This repository provides a complete workflow for analyzing single-cell RNA-seq data using state-of-the-art computational methods. The pipeline integrates both R/Seurat and Python/scVI frameworks to deliver robust results from raw data to biological insights.

**Key Features:**
- Quality control and doublet detection
- Multiple normalization strategies (log1p, SCTransform)
- Data integration using Seurat V5/Harmony and scVI
- Automated and manual cell type annotation
- Differential abundance testing with MILOR
- Multiple cell-cell communication methods
- Differential expression analysis
- Comprehensive visualization tools

## 📁 Pipeline Structure

### 1. Quality Control & Preprocessing
- **`1.QC_DoubletF_Log1p.qmd`** - QC pipeline with doublet filtering and log1p normalization
- **`1.QC_DoubletF_SCT.qmd`** - Alternative QC using SCTransform normalization and doubletFinder
- **`basic_analysis_steps_MISC.Rmd`** - Additional QC utilities and helper functions

### 2. Data Integration
- **`2.Data intergration.Rmd`** - Batch correction and sample integration (Seurat V5)
- **`Single_cell_scVI.ipynb`** - Deep learning-based integration using scVI
- **`Subclustering_SCVI.ipynb`** - Subclustering analysis with scVI
- **`change_scvi_continuous.ipynb`** - Handling continuous covariates in scVI

### 3. Cell Type Annotation
- **`3.Annotation.Rmd`** - Automated and manual cell type annotation
- **`High_annotation_multicontrast.Rmd`** - High-resolution annotation across conditions
- **`High_level_multicontrast.Rmd`** - Broad cell type classification
- **`Low_annotation_multiniche_multicontrast.Rmd`** - Fine-grained annotation analysis

### 4. Differential Abundance Analysis
- **`4. MILOR_dif_abud.Rmd`** - MILOR-based differential abundance testing

### 5. Subclustering & Advanced Analysis
- **`5.Subclustering.Rmd`** - Detailed subclustering of cell populations
- **`Subclustering_SCVI.ipynb`** - Python-based subclustering with scVI

### 6. Cell-Cell Communication Analysis
- **`5.cc_interactions_niche.Rmd`** - Niche-based cell communication analysis
- **`5.cellchat.Rmd`** - CellChat communication inference
- **`6.cc_interactions_niche.Rmd`** - Extended niche interaction analysis
- **`6.cellchat.Rmd`** - Advanced CellChat workflows
- **`6.cellphonedb.Rmd`** - CellPhoneDB communication analysis
- **`6.Multinichenet.Rmd`** - MultiNicheNet analysis
- **`6.1.Multinichenet_output_analysis.Rmd`** - MultiNicheNet results interpretation

### 7. Differential Expression Analysis
- **`7.DEG_conditions_subtypes.Rmd`** - Differential gene expression across conditions and cell types

### 8. Metabolic Analysis
- **`mebocost_rasV.ipynb`** - Metabolic cost analysis using MEBOCOST

### 9. Visualization & Utilities
- **`tSNE.Rmd`** - t-SNE dimensionality reduction and plotting
- **`UMAP_color_Change.Rmd`** - UMAP visualization customization
- **`relative_percentages_plots.Rmd`** - Cell proportion visualization
- **`object_format_convert.Rmd`** - Format conversion utilities

## 🚀 Quick Start

### Prerequisites
```r
# R packages (Seurat V5 required)
install.packages("Seurat") # Version 5.0+
install.packages(c("SingleCellExperiment", "scater"))
BiocManager::install(c("miloR", "CellChat", "MultiNicheNet"))

# Python packages
pip install scvi-tools scanpy pandas numpy matplotlib seaborn
```

### Basic Usage
1. **Start with QC**: Run `1.QC_DoubletF_Log1p.qmd` or `1.QC_DoubletF_SCT.qmd`
2. **Integrate data**: Use `2.Data intergration.Rmd` or `Single_cell_scVI.ipynb`
3. **Annotate cells**: Apply `3.Annotation.Rmd`
4. **Analyze communications**: Choose from CellChat, CellPhoneDB, or MultiNicheNet scripts
5. **Find DEGs**: Run `7.DEG_conditions_subtypes.Rmd`

## 📊 Methods Implemented

**Quality Control:**
- Doublet detection and filtering
- Mitochondrial gene filtering
- Low-quality cell removal

**Normalization:**
- Log1p normalization
- SCTransform
- scVI normalization

**Integration Methods:**
- Seurat V5 CCA/RPCA/Harmony
- scVI deep learning integration

**Cell-Cell Communication:**
- CellChat
- CellPhoneDB  
- MultiNicheNet
- Custom niche analysis
- Mebocost

**Differential Analysis:**
- MILOR (differential abundance)
- Seurat V5 FindMarkers/FindAllMarkers
- edgeR/DESeq2 integration

## 📈 Expected Outputs

- Quality control reports and plots
- Integrated single-cell object
- Cell type annotations
- UMAP/t-SNE visualizations
- Cell-cell communication networks
- Differential expression results
- Abundance analysis results

## 🛠️ File Formats

- **Input**: 10X Genomics H5/MTX, CSV, H5AD, RDS
- **Output**: Seurat objects (RDS), AnnData (H5AD), CSV results, HTML reports

---
**Note**: This pipeline was developed for a private murine dataset. Adjust parameters as needed for your dataset.
