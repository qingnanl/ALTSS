# ALTSS (Alternative Transcription Start Site Analysis Tool)
Using single-nuclei ATAC-seq data to study alternative transcription start site (TSS) usage.

### Introduction
A transcription start site (TSS) is the location where the first DNA nucleotide is transcribed into RNA. This position is very critical for learning the regulation of the gene. For example, in research works, the core promoter is usually defined as a region around the TSS (a few hundred/thousand bps). Also, cis-regulatory elements such as enhancers and repressors are usually thought to be relatively close to the TSS, although in a much larger region (a few hundred kbps) compared with promoter.

In some genes, researchers have observed that they have multiple TSSs (this phenomenon is also called alternative transcription initiation). The alternative TSS usage would lead to different 5'UTR of genes, and some of them would also cause alterations in coding region. Given the importance of the TSS and the core promoter in gene regulation, it is very likely that the alternative TSS usage would affect the gene expression. It would be interesting to study the alternative TSS events at single-cell level to learn the complexity in gene regulation at single cell level.

Here, we develop a pipeline named ALTSS, aiming at using single-nuclei ATAC-seq data as the input and compute a TSS enrichment score for all the possible alternative TSSs, for each single cell. Currently, most single-cell RNA-seq data are 3'-end sequencing, which are not able to resolve the 5'-end in most cases. We reason that snATAC-seq would have a better resolution and coverage close to the TSS region of many genes, and thus it could be the right data type for this task. Of course, this pipeline works with 10x Multiome data with the "atac_fragments.tsv.gz".

The way we compute the TSS enrichment 

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
- `fragment file`; e.g., `fragments.tsv.gz`; this pipeline is only for 10x CellRanger-ATAC output

If you work with human hg38 data, the first three files could be found in the `~/ALTSS/data/` folder.

**Third**, go to the `~/ALTSS/code/sample.sh` file, and modify the parameters accordingly (it has instruction on how to set the parameters). Then run:
```
sh ./sample.sh
```

### Output of this pipeline
By default, this pipeline will give one output named 'dup_enrichment.csv'
