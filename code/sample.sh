#!bin/sh

# -s script directory; it will be /where/this/repository/is/cloned/code/; please note the directory variables should have the last '/'.
# -i input directory; it must contain four files: refGene file; chrom.size file; tsv file; and black list file
# -o output_directory; just set an output_directory; will create one if it does not exist; make sure you have permission to it
# -r refGene file; something like hg38.refGene.txt.gz; curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
# -c chrom.size file; available from https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
# -t tss size; can consider something like 50 or 100
# -p promoter size; can consider something like 2000 or so
# -n python; which python interpreter to use; it must have pandas and numpy installed there
# -v tsv file; usually fragments.tsv.gz file from 10x CellRanger pipeline
# -b black list; e.g., hg38-blacklist.v2.bed.gz; available here https://github.com/Boyle-Lab/Blacklist/tree/master/lists; it could be gzipped

sh ALTSS.sh -s /storage/chenlab/Users/qingnanl/multiome/1_10x_data/ALTSS/code/ \
            -i /storage/chenlab/Users/qingnanl/multiome/1_10x_data/ALTSS_test/ \
            -o /storage/chenlab/Users/qingnanl/multiome/1_10x_data/ALTSS_test/output/ \
            -r hg38.refGene.txt.gz \
            -c hg38.chrom.sizes.txt \
            -t 100 \
            -p 2000 \
            -n /storage/chen/home/qingnanl/anaconda3/envs/scanpy/bin/python\
            -v pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz\
            -b hg38-blacklist.v2.bed.gz\


#
#
