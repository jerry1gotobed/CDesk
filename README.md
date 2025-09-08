# CDesk {ignore}
[TOC]
## Install
### 1. Donwload the scripts
```
git clone https://github.com/jerry1gotobed/CDesk_develop.git
```
### 2. Prepare the conda environments
```
mamba env create -f CDesk.yml
mamba env create -f CDesk_py3.7.yml
mamba env create -f CDesk_py2.7.yml
mamba env create -f CDesk_R.yml
```
We prepare the conda environments for the softwares, R and python environment needed in the scripts (some require additional manual installation). We recommend you use mamba instead of conda to install or the environment process might be too long. Every time running the CDesk, it would first search whether there are CDesk conda or mamba environments. You can also set the env path in the configuration file, it would search the path in the configuration file if there is no CDesk conda env. If no environment detected, it would use your own environment.

### 3. Prepare the data and write the configuration file
You need to prepare the data of different species and write the path in the config.json configuration file that stores the data path and additional softwares path. Below is an example, you can customize your configuration file, prepare for other species, for example. Not all is necessary, it depends on the function you run, it would exit and reports error if you don't have necessary data or software for the specific task.
```
{
  "software":{
    "cellranger":"/home/sugon/Install/cellranger-7.1.0/cellranger",
    "DrSeq":"/usr/local/bin/DrSeq",
    "dnbc4tools":"/mnt/linzejie/software/dnbc4tools2.1.3/dnbc4tools",
    "juicer_tools_jar":"/mnt/kuangjunqi/kuangjunqi/tools/hic/juicer_tools_1.22.01.jar",
    "dipc":"/mnt/zhaochengchen/Work/OTHER/4.KJQ/2.HiC-Pro/4.hickit/dip-c",
    "hickit":"/mnt/linzejie/software/hickit-0.1.1_x64-linux/hickit",
    "hickit_js":"/mnt/linzejie/software/hickit-0.1.1_x64-linux/hickit.js"
  },
  "conda_env":{
    "CDesk":"/mnt/linzejie/miniconda3/envs/CDesk",
    "CDesk_R":"/mnt/linzejie/miniconda3/envs/CDesk_R",
    "CDesk_py3.7":"/mnt/linzejie/miniconda3/envs/CDesk_py3.7",
    "CDesk_py2.7":"/mnt/linzejie/miniconda3/envs/CDesk_py2.7"
  },
  "data":{
    "mm10": {
      "mapping_index": "/mnt/zhaochengchen/Data/mm10/mm10",
      "refseq_gtf": "/mnt/zhaochengchen/Data/mm10/mm10.ncbiRefSeq.WithUCSC.gtf",
      "refseq_bed": "/mnt/zhaochengchen/Data/mm10/mm10.refseq.bed",
      "chromInfo": "/mnt/zhaochengchen/Data/mm10/mm10.len",
      "fasta": "/mnt/zhaochengchen/Data/mm10/mm10.fa",
      "tf_file": "/mnt/zhaochengchen/Data/mm10/Mus_musculus_TF.txt",
      "promoter_file": "/mnt/liudong/data/Genome/mm10/mm10.promoter.ncbiRefSeq.WithUCSC.bed",
      "TE_idx": "/mnt/liudong/data/Genome/mm10/mm10.exclusive.idx",
      "scRef10x": "/mnt/liudong/data/Genome/mm10Self",
      "effective_genome_size": "mm",
      "gff3":"/mnt/linzejie/data/mm10.gff3",
      "refgenes":"/mnt/linzejie/data/Drseq_data/refgenes/mm10.refgenes.txt",
      "bowtie2_mapindex":"/mnt/linzejie/data/Drseq_data/bowtie2/mm10/mm10",
      "singleron_mapindex":"/mnt/linzejie/data/celescope_data/mm10",
      "dnbc_mapindex":"/mnt/linzejie/data/dnbc4tools_data/mm10"
    },
    "rn7": {
      "mapping_index": "/mnt/liudong/data/Genome/rn7/rn7",
      "refseq_gtf": "/mnt/liudong/data/Genome/rn7/rn7.refGene.gtf",
      "refseq_bed": "/mnt/liudong/data/Genome/rn7/rn7.refGene.fix.bed",
      "chromInfo": "/mnt/liudong/data/Genome/rn7/rn7.len",
      "fasta": "/mnt/linzejie/data/fasta/rn7.fa",
      "tf_file": "/mnt/liudong/data/Genome/rn7/Rattus_norvegicus_TF.txt",
      "promoter_file": "/mnt/liudong/data/Genome/rn7/rn7.refGene.promoter.bed",
      "TE_idx": "/mnt/linzejie/data/rn7.exclusive.idx",
      "scRef10x": "/mnt/linzejie/data/cellranger_data/rn7",
      "effective_genome_size": "2.37e9",
      "gff3":"/mnt/linzejie/data/rn7.gff3",
      "refgenes":"/mnt/linzejie/data/Drseq_data/refgenes/rn7.refgenes.txt",
      "bowtie2_mapindex":"/mnt/linzejie/data/Drseq_data/bowtie2/rn7/rn7",
      "singleron_mapindex":"/mnt/linzejie/data/celescope_data/rn7",
      "dnbc_mapindex":"/mnt/linzejie/data/dnbc4tools_data/rn7"
    },
    "hg38": {
      "mapping_index": "/mnt/zhaochengchen/Data/hg38/hg38",
      "refseq_gtf": "/mnt/zhaochengchen/Data/hg38/hg38.ncbiRefSeq.WithUCSC.gtf",
      "refseq_bed": "/mnt/zhaochengchen/Data/hg38/hg38.refseq.bed",
      "chromInfo": "/mnt/zhaochengchen/Data/hg38/hg38.len",
      "fasta": "/mnt/zhaochengchen/Data/hg38/hg38.fa",
      "tf_file": "/mnt/zhaochengchen/Data/hg38/Homo_sapiens_TF.txt",
      "promoter_file": "/mnt/zhaochengchen/Data/hg38/hg38.Promoter.bed",
      "TE_idx": "/mnt/liudong/data/Genome/hg38/hg38.exclusive.idx",
      "scRef10x": "/mnt/linzejie/data/cellranger_data/hg38",
      "effective_genome_size": "hs",
      "gff3":"/mnt/linzejie/data/hg38.gff3",
      "refgenes":"/mnt/linzejie/data/Drseq_data/refgenes/hg38.refgenes.txt",
      "bowtie2_mapindex":"/mnt/linzejie/data/Drseq_data/bowtie2/hg38/hg38",
      "singleron_mapindex":"/mnt/linzejie/data/celescope_data/hg38",
      "dnbc_mapindex":"/mnt/linzejie/data/dnbc4tools_data/hg38"
    },
    "susScr11": {
      "mapping_index": "/mnt/zhaochengchen/Data/susScr11/susScr11",
      "refseq_gtf": "/mnt/liudong/data/Genome/susScr11/susScr11.ncbiRefSeq.gtf",
      "refseq_bed": "/mnt/liudong/data/Genome/susScr11/susScr11.refseq.bed",
      "chromInfo": "/mnt/zhaochengchen/Data/susScr11/susScr11.chrom.sizes",
      "fasta": "/mnt/linzejie/data/fasta/susScr11.fa",
      "tf_file": "/mnt/liudong/data/Genome/susScr11/Sus_scrofa_TF.txt",
      "promoter_file": "/mnt/liudong/data/Genome/susScr11/susScr11.Promoter.bed",
      "TE_idx": "/mnt/liudong/data/Genome/susScr11/susScr11.exclusive.idx",
      "scRef10x": "/mnt/liudong/data/Genome/susScr11Self",
      "effective_genome_size": "2.36e9",
      "gff3":"/mnt/linzejie/data/susScr11.gff3",
      "refgenes":"/mnt/linzejie/data/Drseq_data/refgenes/susScr11.refgenes.txt",
      "bowtie2_mapindex":"/mnt/linzejie/data/Drseq_data/bowtie2/susScr11/susScr11",
      "singleron_mapindex":"/mnt/linzejie/data/celescope_data/susScr11",
      "dnbc_mapindex":"/mnt/linzejie/data/dnbc4tools_data/susScr11"
    },
    "galGal6": {
      "mapping_index": "/mnt/zhaochengchen/Data/Gal6/galGal6",
      "refseq_gtf": "/mnt/zhaochengchen/Data/Gal6/galGal6.ncbiRefSeq.WithUCSC.gtf",
      "refseq_bed": "/mnt/zhaochengchen/Data/Gal6/galGal6.refseq.bed",
      "chromInfo": "/mnt/zhaochengchen/Data/Gal6/galGal6.chrom.sizes",
      "fasta": "/mnt/linzejie/data/fasta/galGal6.fa",
      "TE_idx": "/mnt/linzejie/data/galGal6.exclusive.idx",
      "scRef10x": "/mnt/linzejie/data/cellranger_data/galGal6",
      "effective_genome_size": "1.0e9",
      "tf_file": "/mnt/zhaochengchen/Data/Gal6/Gallus_gallus_TF.txt",
      "promoter_file": "/mnt/linzejie/data/galGal6.Promoter.bed",
      "gff3":"/mnt/linzejie/data/galGal6.gff3",
      "refgenes":"/mnt/linzejie/data/Drseq_data/refgenes/galGal6.refgenes.txt",
      "bowtie2_mapindex":"/mnt/linzejie/data/Drseq_data/bowtie2/galGal6/galGal6",
      "singleron_mapindex":"/mnt/linzejie/data/celescope_data/galGal6",
      "dnbc_mapindex":"/mnt/linzejie/data/dnbc4tools_data/galGal6"
    },
    "rmats_test":{
    	"gff3":"/mnt/linzejie/CDesk_test/data/1.RNA/7.splice/annotation.gff3",
	"refseq_gtf": "/mnt/linzejie/CDesk_test/data/1.RNA/7.splice/annotation.gtf"
    },
    "mm10_plusmajsat":{
        "fasta":"/mnt/zhangzheting/reference/mm10_majsat/mm10.plusmajsat.fa",
        "chromInfo": "/mnt/zhangzheting/reference/mm10_majsat/mm10.plusmajsat.chrom.size"
    }
  }
}

```
## Hints