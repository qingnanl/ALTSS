#! /bin/bash
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qingnan.liang@bcm.edu
python /storage/chenlab/Users/qingnanl/multiome/1_10x_data/bin_by_gene/split_bam.py \
/storage/chenlab/Users/qingnanl/multiome/1_10x_data/bin_by_gene/input/pbmc_granulocyte_sorted_3k_atac_possorted_bam.bam \
/storage/chenlab/Users/qingnanl/multiome/1_10x_data/bin_by_gene/input/cell_id_AB1.txt \
"HSP90" \
/storage/chenlab/Users/qingnanl/multiome/1_10x_data/bin_by_gene/output/ \

