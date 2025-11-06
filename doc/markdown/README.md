# CDesk pipeline handbook{ignore}

CDesk is an integrated multi-omics analysis pipeline designed for processing data from various sequencing-based assays, including RNA-seq, scRNA-seq, ATAC-seq, CUT&Tag, ChIP-seq, and Hi-C. It comprises multiple subcommands that cover a comprehensive range of analysis tasks, from raw sequencing data process to downstream various advanced functions. Dedicated conda environment YAML files are supplied. To install, simply create the Conda environment from thess files (some functions may require additional software) and prepare the necessary species-specific data. Once configured, users can perform the desired analyses by entering the corresponding command on the command line.

[0. Installation ](BulkRNA.md)

[1. BulkRNA ](BulkRNA.md)

1.1 Preprocss

1.2 Quality Control

[2. scRNA ](scRNA.md)

[3. ATAC ](ATAC.md)

[4. ChIPseq&CUTTag ](ChIPseq&CUTTag.md)

[5. HiC ](HiC.md)

5.1 HiC: Preprocess

5.2 HiC: Sample Correlation

5.3 HiC: Matrix balancing and Format transformation

5.4 HiC: TAD

5.5 HiC: Compartment

5.6 HiC: Loop

5.7 HiC: 3D reconstruction

5.8 HiC: Distance-contact

5.9 HiC: Contact compare