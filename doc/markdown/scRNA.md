---
html:
    toc: true
---
# CDesk: 2. scRNA pipeline{ignore}
Our CDesk scRNA module comprises of 8 function submodules. Here we present you the CDesk scRNA working pipeline and how to use it to analyze your scRNA data. 
[toc]

## 2.1 scRNA: Preprocess


## 2.2 scRNA: cluster
CDesk scRNA cluster module automates scRNA-seq clustering analysis by providing two analytical modes: Seurat (R) and Scanpy (Python). The workflow includes standard scRNA-seq analysis steps—data loading, quality control (filtering of genes and cells, mitochondrial gene content assessment), normalization, feature selection, dimensionality reduction (PCA, UMAP, t-SNE), and clustering. It supports multiple input formats, including H5, H5AD, RDS, and 10x Genomics data, and generates comprehensive PDF reports containing clustering results and visualizations.

Here is an example about how to use the CDesk scRNA clustering module. 
```shell
CDesk scRNA cluster \
-i /.../input \
-o /.../output_directory \
-n test --mode seurat(scanpy)
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|The input scRNA data (.h5,.txt/csv/tsv.gz,10x input directory,.rds for seurat mode,.h5ad for scanpy mode)|
|-o,--output*|The output directory|
|-n,--name*|The output prefix name|
|--mode*|Analysis mode|{seurat,scanpy}
|--min_cells|minimum cells threshold for data filtration preprocess|3
|--min_features|minimum features threshold data filtration preprocess|200
|--nFeature_RNA_min|clustering feature minimum|200
|--nFeature_RNA_max|clustering feature maximum|2500
|--mt_percent|clustering mitochondria percent threshold|5
|--variable_features|clustering variable features|2000
|--dim_prefer|clustering dimesions|30
|--res|clustering resolution|1
|--width|Plot width|10
|--height|Plot height|8

If the pipeline runs successfully, there would be clustering plots and h5ad output file in the output directory for scanpy mode. There would be elbow plots, cluste clustering plots, cluster tree plots and rds output file in the output directory for seurat mode.
<div align="center">
<img src="./images/scRNA_cluster_scanpy1.png" width="330"><img src="./images/scRNA_cluster_scanpy2.png" width="330">
</div>
<center><i>CDesk scRNA scanpy cluster example result</i></center>

<div align="center">
<img src="./images/scRNA_cluster_seurat1.png" width="330"><img src="./images/scRNA_cluster_seurat2.png" width="330">
<img src="./images/scRNA_cluster_seurat3.png" width="330"><img src="./images/scRNA_cluster_seurat4.png" width="330">
</div>
<center><i>CDesk scRNA seurat cluster example result</i></center>

<details>
<summary>A successful CDesk scRNA cluster running process</summary>
<pre><blockcode>

</blockcode></pre>
</details>

## 2.3 scRNA: annotation
CDesk scRNA annotation module realizes the automation of annotation of single-cell data. It offers two annotation strategies: marker gene-based annotation (marker mode) and reference-guided label transfer (transfer mode). In marker mode, cell types are assigned by comparing the expression levels and detection rates of known marker genes across clusters. In transfer mode, cell type labels from a reference dataset are transferred to the query dataset using Seurat’s integration framework. Both methods produce annotated cell type labels and automatically generate visualizations of the annotations projected onto low-dimensional embeddings.

1. marker mode
   
Cell type annotation is performed using user-provided marker genes that are specifically associated with defined cell types. Three key thresholds are applied to determine cell type identity:

1. The marker gene must have an average expression level in the cluster exceeding the threshold expression_thre (log1p normalized expression). 
2. The gene must be expressed in a proportion of cells within the cluster greater than percentage_thre.
3. For a candidate cell type to be considered, a minimum fraction (marker_percentage) of its marker genes must simultaneously meet both of the above criteria.
   
Finally, the average proportion score is calculated for each qualifying cell type, and the cell type with the highest score is assigned as the final annotation for each cluster.

Here is an example about how to use the CDesk scRNA annotation marker module. 

```shell
CDesk scRNA annotation marker \
-i /.../input.rds \
--marker /.../marker.csv \
--meta meta \
-o /.../output_directory
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|The input .rds file to annotate|
|-o,--output*|The output directory|
|--marker*|The marker file with cell marker information|
|--meta*|The meta colname of reference data used|
|--cluster|Plot clustering|umap
|--expression_thre|The average expression threshold|1.5
|--percentage_thre|The express > 0 percentage threshold|0.7
|--marker_thre|The proportion threshold for passing the screening|0.7
|--width|Plot width|10
|--height|Plot height|8

If the pipeline runs successfully, there would be a cell annotation file and a plot of the annotations projected onto low-dimensional embeddings.

<details>
<summary>
What should the input file look like?</summary>
<pre>
<blockcode>
=== marker.csv ===
Celltype,Marker
C2,Zscan4a
C2,Zscan4b
C2,Eif1a
C2,Dux
C4,Obox4
C4,Khdc1b
C4,Usp17l1
C4,Zfp352
C8,Cdx2
C8,Pou5f1
C8,Nanog
C8,Eomes
C8,Tead4
TE,Cdx2
TE,Eomes
TE,Gata3
TE,Krt8
TE,Krt18
TE,Elf5
PrE,Gata6
PrE,Gata4
PrE,Sox17
PrE,Pdgfra
PrE,Dab2
PrE,Lrp2
</blockcode>
- Celltype: The celltypes to specify
- Marker: The corresponding markers
</pre>
</details>
<br>

2. transfer mode
   
It enables reference-based cell type label transfer. It first loads the reference and query datasets, identifies shared genes, and performs independent preprocessing—including normalization, variable feature selection, scaling, and PCA. Anchors between datasets are then identified using Seurat’s FindTransferAnchors in PCA space, followed by label transfer via TransferData. Cell type annotations from the reference (specified by ref_meta) are predicted and added to the query dataset’s metadata.

Here is an example about how to use the CDesk scRNA  annotation transfer module. 

```shell
CDesk scRNA annotation transfer \
--query /.../qry.rds \
--ref /.../ref.rds \
--meta meta \
-o /.../output_directory
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|--query*|The query .rds data to annotate|
|--ref*|The reference .rds data|
|-o,--output*|The output directory|
|--meta*|The meta colname of reference data used|
|--cluster|Plot type|umap
|--nfeatures|The number of variable features|5000
|--dims|The number of PC dimensions used|30
|--width|Plot width|10
|--height|Plot height|8

If the pipeline runs successfully, there would be a cell annotation file and a plot of the annotations projected onto low-dimensional embeddings.

<div align="center">
<img src="./images/scRNA_annotation_marker.png" height="300"><img src="./images/scRNA_annotation_transfer.png" height="300">
</div>
<center><i>CDesk scRNA annotation example result (left: marker, right: transfer)</i></center>

## 2.4 scRNA: marker
CDesk marker module supports two modes of differential expression analysis for single-cell RNA-seq data:
1. All-clusters mode (all): Uses FindAllMarkers to automatically identify differentially expressed genes (DEGs) for each cluster compared to all others, and generates a heatmap visualizing the top 10 marker genes per cluster.
2. Two-group mode (2group): Uses FindMarkers to perform pairwise comparison between two user-specified cell types or clusters, identifying their specific DEGs.
   
Both modes allow customization of key parameters, including log-fold change threshold, adjusted p-value cutoff, and minimum expression percentage. Significantly differentially expressed genes are exported to CSV files for downstream analysis.

1. all

Here is an example about how to use the CDesk scRNA marker all module.

```shell
CDesk scRNA marker all \
-i /.../input.rds -o /.../output_directory \
--meta meta -m 0.25 -p 0.05 -f 0.5
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|The input Seurat object rds file|
|-o,--output*|The output directory|
|--meta*|The meta colname of reference data used|
|-fc|Log Fold Change threshold|0.5
|-p|Adjusted p-value threshold|0.05
|-m|Minimum percentage of expressed cells|0.25
|--width|Plot width|16
|--height|Plot height|14

If the pipeline runs successfully, there would be a cell type marker file and a heatmap plot showing the top differential marker expression in all cell types.

1. 2 group
   
Here is an example about how to use the CDesk scRNA marker 2group module.

```shell
CDesk scRNA marker 2group \
-i /.../input.rds -o /.../output_directory \
--meta meta --type1 type1 --type2 type2 \
-f 0.5 -p 0.05 -m 0.25
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|The input Seurat object rds file|
|-o,--output*|The output directory|
|-m,--meta*|The meta colname of reference data used|
|--type1*|Specify the first cell type|
|--type2*|Specify the second cell type|
|-fc|Log Fold Change threshold|0.5
|-p|Adjusted p-value threshold|0.05
|-m|Minimum percentage of expressed cells|0.25
|--width|Plot width|16
|--height|Plot height|14

If the pipeline runs successfully, there would be a cell type marker file and a heatmap plot showing the marker expression in two cell types.

<div align="center">
<img src="./images/scRNA_marker_all.png" height="300"><img src="./images/scRNA_marker_2.png" height="300">
</div>
<center><i>CDesk scRNA marker example result (left: all, right: 2 group)</i></center>

## 2.5 scRNA: trajectory
CDesk trajectory module supports four single-cell transcriptomic analysis functions for trajectory inference and RNA velocity:

1. cstreet: Performs cell trajectory inference using the CStreet tool.
2. diffusion: Constructs developmental trajectories via diffusion map embedding and computes pseudotime.
3. monocle: Implements advanced trajectory inference with Monocle3, enabling differential gene analysis and co-expression module detection.
4. velocity: Conducts RNA velocity analysis to predict cell fate transitions by integrating spliced and unspliced mRNA dynamics, revealing the directionality of cellular state changes.

All modules support parameterized input and generate visual outputs, providing a comprehensive suite for transitioning from static clustering to dynamic developmental analysis of single-cell data.

1. cstreet

Here is an example about how to use the CDesk trajectory cstreet module.

```shell
CDesk scRNA trajectory cstreet \
-i /.../ExpressionMatrix_list.txt \
-s /.../CellStates_list.txt \
-n test -o /.../output_directory
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|Input expression matrixes list file|
|-s,--state*|Input cell state list file|
|-o,--output*|The output directory|
|-n,--name*|The meta colname of reference data used|

If the pipeline runs successfully, there would be CStreet result in the output directory.

<div align="center">
<img src="./images/scRNA_trajectory_cstreet.png" height="400"></div>
<center><i>CDesk scRNA trajectory cstreet example result</i></center>

2. diffusion

Here is an example about how to use the CDesk trajectory diffusion module.

```shell
CDesk scRNA trajectory diffusion \
-i /.../input.rds \
-o /.../output_directory \
-m meta -s celltype
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|Input seurat rds file|
|-s,--start*|Cell type start point for trajectory analysis|
|-o,--output*|The output directory|
|-m,--meta*|Meta column for trajectory analysis|
|--width|Plot width|8
|--height|Plot height|8

If the pipeline runs successfully, it would generate a PDF file containing all diffusion map embeddings and rds object with diffusion map.

<div align="center">
<img src="./images/scRNA_trajectory_density1.png" height="300"><img src="./images/scRNA_trajectory_density2.png" height="300"></div>
<div align="center">
<img src="./images/scRNA_trajectory_density3.png" height="300"><img src="./images/scRNA_trajectory_density4.png" height="300"></div>
<center><i>CDesk scRNA trajectory diffusion example result</i></center>

1. monocle

Here is an example about how to use the CDesk trajectory monocle module.

```shell
CDesk scRNA trajectory monocle \
i /.../input.rds \
-m meta -s celltype \
-o /.../output_directory
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|-i,--input*|Input seurat rds file|
|-s,--start*|Cell type start point for trajectory analysis|
|-o,--output*|The output directory|
|-m,--meta*|Meta column for trajectory analysis|
|--width|Plot width|8
|--height|Plot height|8

If the pipeline runs successfully, it generates the following outputs:
- Trajectory plot annotated with cell type labels
- Pseudotime trajectory visualization
- Expression trend plots for the top 10 differentially expressed genes
- Feature plots for the top 10 DEGs
- Co-expression module heatmap
- Seurat object with pseudotime annotations
- Gene module information

<div align="center">
<img src="./images/scRNA_trajectory_monocle_pseudotime0.png" width="350"></div>
<div align="center">
<img src="./images/scRNA_trajectory_monocle_pseudotime1.png" width="350"><img src="./images/scRNA_trajectory_monocle_pseudotime2.png" width="350"></div>
<center><i>CDesk scRNA trajectory monocle pseudotime trajectory visualization example result</i></center>
<div align="center">
<img src="./images/scRNA_trajectory_monocle_topDEGs1.png" width="350"><img src="./images/scRNA_trajectory_monocle_topDEGs2.png" width="350"></div>
<center><i>CDesk scRNA trajectory monocle feature plots for the top 10 DEGs example result</i></center>

<div align="center">
<img src="./images/scRNA_trajectory_monocle_coexpression.png" width="400"></div>
<center><i>CDesk scRNA trajectory monocle feature co-expression module heatmap example result</i></center>

4. velocity

Here is an example about how to use the CDesk trajectory velocity module.

```shell
CDesk scRNA trajectory velocity \
--cellranger /.../cellranger \
--rds /.../input.rds \
--genes /.../genes.txt \
-o /.../output_directory \
--meta meta \
--gtf /.../gtf.gtf
```
|Parameters^(*necessary)^|Description|Default value|
|----|----|----|
|--cellranger*|Input cellranger output directory|
|-o,--output*|The output directory|
|--rds*|Input seurat rds file|
|-m,--meta*|Meta column for trajectory analysis|
|--gtf|GTF file for velocity analysis if velocyto result in the cellranger output|
|--genes|Interested genes file|
|-t,--thread|Number of threads|10
|--width|Plot width|10
|--height|Plot height|10

If the pipeline runs successfully, it would generate:
- Loom files generated by velocyto, containing spliced and unspliced count matrices
- Cell IDs, UMAP coordinates, and cluster annotations extracted from the Seurat object
- Bar plots showing the spliced vs. unspliced read ratios for each cluster
- Visualizations of RNA velocity results, including grid velocity and manifold (embedding) plots
- List of identified key genes with significant velocity signals
- Velocity stream plots for user-specified genes of interest (if provided)
- RNA velocity pseudotime plot
- PAGA (Partition-based Graph Abstraction) graph illustrating the developmental relationships between clusters

<div align="center">
<img src="./images/scRNA_trajectory_velocity_proportion.png" height="200"></div>

<div align="center">
<img src="./images/scRNA_trajectory_velocity1.png" height="400"><img src="./images/scRNA_trajectory_velocity2.png" height="400"></div>

<div align="center">
<img src="./images/scRNA_trajectory_velocity_genes.png" height="250"></div>

<div align="center">
<img src="./images/scRNA_trajectory_pseudotime1.png" height="400"><img src="./images/scRNA_trajectory_pseudotime2.png" height="400"></div>
<center><i>CDesk scRNA trajectory velocity example result</i></center>

## 2.6 scRNA: similarity


## 2.7 scRNA: interaction

## 2.8 scRNA: integrate