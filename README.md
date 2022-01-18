# ALTSS (Alternative Transcription Start Site Analysis Tool)
Using single-nuclei ATAC-seq data to study alternative transcription start site (TSS) usage.

### Introduction
A transcription start site (TSS) is the location where the first DNA nucleotide is transcribed into RNA. This position is very critical for learning the regulation of the gene. For example, in research works, the core promoter is usually defined as a region around the TSS (a few hundred/thousand bps). Also, cis-regulatory elements such as enhancers and repressors are usually thought to be relatively close to the TSS, although in a much larger region (a few hundred kbps) compared with promoter.

In some genes, researchers have observed that they have multiple TSSs (this phenomenon is also called alternative transcription initiation). The alternative TSS usage would lead to different 5'UTR of genes, and some of them would also cause alterations in coding region. Given the importance of the TSS and the core promoter in gene regulation, it is very likely that the alternative TSS usage would affect the gene expression. It would be interesting to study the alternative TSS events at single-cell level to learn the complexity in gene regulation at single cell level.

Here, we develop a pipeline named ALTSS, aiming at using single-nuclei ATAC-seq data as the input and compute a TSS enrichment score for all the possible alternative TSSs, for each single cell. Currently, most single-cell RNA-seq data are 3'-end sequencing, which are not able to resolve the 5'-end in most cases. We reason that snATAC-seq would have a better resolution and coverage close to the TSS region of many genes, and thus it could be the right data type for this task. Of course, this pipeline works with 10x Multiome data with the "atac_fragments.tsv.gz".

The way we compute the TSS enrichment is very simple. For a certain TSS in a given gene, the enrichment is calculated:

![image](https://user-images.githubusercontent.com/53788946/149866315-bf715d79-4546-4eeb-b048-9ddcef83ad57.png)

- **ES** is enrichment score; 
- **T-norm** is normalized TSS region reads; 
- **P-norm** is normalized promoter region reads, and:

![image](https://user-images.githubusercontent.com/53788946/149866527-60aacc47-b2bf-475d-916b-a1356532e896.png)

![image](https://user-images.githubusercontent.com/53788946/149866571-9a3902fa-e777-4135-9805-adca1f0bd6ee.png)

- **T-reads** is the number of fragments (in the fragments.tsv.gz) overlapping with TSS region (user defined **T-length** bp region upstream of the tss; recommended 100 bp);
- **T-length** is the defined length of TSS region, as mentioned above;
- **P-reads** and **P-length** are the same thing for a user defined promoter region. The recommended length is 2000 bp;
- **scale** is a scale factor, recommended to use 10000.

In this method, we solely rely on the gene annotation (`refGene file`) for filtering out the genes that do not have alternative TSS. In other words, we only calculate the TSS enrichment score for genes with alternative TSS annotated. The way to calculate this score is simply inspired by the 'routine' quality control method for snATAC-seq data, which is to calculated the overall TSS enrichment score for a cell. The cells very low scores are recommended to be filtered out from further analysis. Here we expand this to each TSS.

### System requirement and how to run this pipeline
To run this pipeline, the followings are required:
- UNIX system
- Basic UNIX tools: gunzip, zcat, awk, grep, uniq
- Python3, with pandas and numpy installed (python 3.6 tested; should work for other python3 versions)
- bedtools (v2.25.0 tested)

**First**, clone this repository to somewhere local:
```
git clone https://github.com/qingnanl/ALTSS
```

**Second**, collect the following data to one folder:
- `refGene file`; e.g., `hg38.refGene.txt.gz`; could be found http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/; it is strongly recommended to use the refGene files from this `'goldenPath'` resource, and others may not work due to their differences in organizing the data (basically, this pipeline is sensitive to the order of columns; however, if others have to be used and have issues, we can change the code in a reasonable time :grinning:)
- `chrom.size file`; e.g., `hg38.chrom.sizes.txt`; available from https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
- `black list`; e.g., `hg38-blacklist.v2.bed.gz`; available here https://github.com/Boyle-Lab/Blacklist/tree/master/lists; it could be gzipped
- `fragment file`; e.g., `fragments.tsv.gz`; this pipeline is currently only for 10x CellRanger-ATAC output

If you work with human hg38 data, the first three files could be found in the `~/ALTSS/data/` folder.

**Third**, go to the `~/ALTSS/code/sample.sh` file, and modify the parameters accordingly (it has instruction on how to set the parameters). Then run:
```
sh ./sample.sh
```

### Output of this pipeline
By default, this pipeline will give one output named 'dup_enrichment.csv'

It could be used for downstream analysis in R or python, for example in R:
```
dup <- read.csv('dup_enrichment.csv', row.names = 1)
head(dup[1:10, 1:20])
```
![image](https://user-images.githubusercontent.com/53788946/149873050-7e2b4467-c55c-4d14-9dab-5e5283af06e5.png)

Each row is a cell barcode; each column is a gene-name_tss.
