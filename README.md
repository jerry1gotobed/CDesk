# CDesk_develop
This is the developing CDesk

How to use?
```
git clone https://github.com/jerry1gotobed/CDesk_develop.git
```

Modify the config.json to change the software and data path

It is OK on pei2: /mnt/linzejie/CDesk

Test command(pei2):
```
# RNA:
#1. preprocess
python CDesk.py bulkRNA preprocess -i /mnt/zhaochengchen/CDesk_TestData/RNAseq -o /mnt/linzejie/CDesk_test/result/1.RNA/1.preprocess/test -s mm10 -t 200 -l 2 -bw -te --te_gene /mnt/linzejie/CDesk_test/data/1.RNA/1.preprocess/TE_name.txt

python CDesk.py bulkRNA preprocess -i /mnt/linzejie/CDesk_test/data/1.RNA/1.preprocess/bam_fq -o /mnt/linzejie/CDesk_test/result/1.RNA/1.preprocess/bam_fq_test -s mm10 -t 100 -l 2 -bw -te

python CDesk.py bulkRNA preprocess -i /mnt/linzejie/CDesk_test/data/1.RNA/1.preprocess/bam -o /mnt/linzejie/CDesk_test/result/1.RNA/1.preprocess/bam_test -s mm10 -t 100 -l 2 -bw -te --te_gene /mnt/linzejie/CDesk_test/data/1.RNA/1.preprocess/TE_name.txt

# 2. QC
python CDesk.py bulkRNA QC -i /mnt/linzejie/CDesk_test/data/1.RNA/2.QC/fpkm.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/2.QC --group /mnt/linzejie/CDesk_test/data/1.RNA/2.QC/group.xlsx

# 3. DEGD
python CDesk.py bulkRNA DEG deseq2 -i /mnt/linzejie/CDesk_test/data/1.RNA/3.DEGD/merged_count.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/3.DEGD/deseq2 -m /mnt/linzejie/CDesk_test/data/1.RNA/3.DEGD/DEG_yap.meta.xlsx -p -g Yap1,Trp53,Nelfa,Hdac11,Mtf2,Jade1,Dppa2,Dppa4,Nanog,Pou5f1,Sox2,Zscan4c,Zscan4d,Zscan4f,Duxf3,Tcstv1,Tcstv3,Gm5662 -t 30

python CDesk.py bulkRNA DEG adjusted_t -i /mnt/linzejie/CDesk_test/data/1.RNA/3.DEGD/merged_count.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/3.DEGD/adjusted_t -m /mnt/linzejie/CDesk_test/data/1.RNA/3.DEGD/DEG_yap.meta.xlsx -p -g Yap1,Trp53,Nelfa,Hdac11,Mtf2,Jade1,Dppa2,Dppa4,Nanog,Pou5f1,Sox2,Zscan4c,Zscan4d,Zscan4f,Duxf3,Tcstv1,Tcstv3,Gm5662 -t 30

python CDesk.py bulkRNA DEG gfold -i /mnt/linzejie/CDesk_test/result/1.RNA/1.preprocess/test/Bam -o /mnt/linzejie/CDesk_test/result/1.RNA/3.DEGD/gfold -m /mnt/linzejie/CDesk_test/data/1.RNA/3.DEGD/meta2.xlsx -s mm10 -p -g LEUTX,POU5F1 -t 100

# 4. enrichment
python CDesk.py bulkRNA enrichment custom -i /mnt/linzejie/CDesk_test/data/1.RNA/4.enrichment/custom.txt -c /mnt/linzejie/CDesk_test/data/1.RNA/4.enrichment/custom_ref.txt -o /mnt/linzejie/CDesk_test/result/1.RNA/4.enrichment/custom

python CDesk.py bulkRNA enrichment GO -i /mnt/linzejie/CDesk_test/data/1.RNA/4.enrichment/gene_list.txt -s pig -t 0.5 -o /mnt/linzejie/CDesk_test/result/1.RNA/4.enrichment/GO

python CDesk.py bulkRNA enrichment KEGG -i /mnt/linzejie/CDesk_test/data/1.RNA/4.enrichment/test.txt -s pig -o /mnt/linzejie/CDesk_test/result/1.RNA/4.enrichment/KEGG

# 5. GSEA
python CDesk.py bulkRNA GSEA analyze --rnk /mnt/linzejie/CDesk_test/result/1.RNA/3.DEGD/deseq2/type.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/5.GSEA/CDesk -s human --prepreocess
python CDesk.py bulkRNA GSEA analyze --rnk /mnt/linzejie/CDesk_test/data/1.RNA/5.GSEA/gsea_test.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/5.GSEA/test -s human
python CDesk.py bulkRNA GSEA plot -i /mnt/linzejie/CDesk_test/result/1.RNA/5.GSEA/CDesk --term HALLMARK_MYC_TARGETS_V2,HALLMARK_MYC_TARGETS_V1

# 6. Similarity
python CDesk.py bulkRNA similarity --sample1 /mnt/linzejie/CDesk_test/data/1.RNA/6.Similarity/sample1_fpkm.csv --sample2 /mnt/linzejie/CDesk_test/data/1.RNA/6.Similarity/sample2_fpkm.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/6.Similarity/group --group /mnt/linzejie/CDesk_test/data/1.RNA/6.Similarity/meta.xlsx --gene /mnt/linzejie/CDesk_test/data/1.RNA/6.Similarity/test_gene_list.txt

python CDesk.py bulkRNA similarity --sample1 /mnt/linzejie/CDesk_test/data/1.RNA/6.Similarity/sample1_fpkm.csv --sample2 /mnt/linzejie/CDesk_test/data/1.RNA/6.Similarity/sample2_fpkm.csv -o /mnt/linzejie/CDesk_test/result/1.RNA/6.Similarity/no_group

# 7. WGCNA
python CDesk.py bulkRNA WGCNA --count_file /mnt/linzejie/CDesk_test/data/1.RNA/7.WGCNA/test.txt --pheno_file /mnt/linzejie/CDesk_test/data/1.RNA/7.WGCNA/pheno.txt --trait 3D -o /mnt/linzejie/CDesk_test/result/1.RNA/7.wgcna

# 8. Splice
python CDesk.py bulkRNA splice detect --b1 /mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789776.bam,/mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789777.bam --b2 /mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789778.bam,/mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789779.bam --species mm10 -o /mnt/linzejie/CDesk_test/result/1.RNA/8.splice/detect -t 100

python CDesk.py bulkRNA splice draw --b1 /mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789776.bam,/mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789777.bam --b2 /mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789778.bam,/mnt/linzejie/CDesk_test/data/1.RNA/8.splice/GSM7789779.bam --species mm10 -o /mnt/linzejie/CDesk_test/result/1.RNA/8.splice/draw --interval chr1:+:3000000:3200000 --group /mnt/linzejie/CDesk_test/data/1.RNA/8.splice/grouping.gf
```