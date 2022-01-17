#!bin/sh
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