# This is a pipeline to generate a matrix containing the information of TSS usage for single-cell ATAC-seq data.
# Currently this is designed for the output data from the 10x Genomics Cell Ranger (atac) pipeline.
# There are several input files for this pipeline:
# 
# First, it takes a gene annotation file, e.g., hg38.refGene.txt.gz as the input. 
# This file have the coordinates for genes, and thus allows for getting the TSS information
# For each individual TSS, two window files are created: one for TSS (recommended 100 bp region upstream of the TSS region), 
# the other for 'promoter' (recommended 2000 bp region upstream of the TSS region). The genes were only considered if they
# have alternative tss annotated (if only one tss, the gene will be filtered out)
# 
# Second, the window files (after filtering out blacklist regions) were used to intersect with the xxx.fragments.tsv.gz file,
# to count reads in each window (for each cell). The count data for tss and promoter would be summarized into a matrix format (cell by tss/promoter count).
# If a cell has very low total tss or promoter count, it will be filtered out (just like how most snATAC-seq pipeline do)
#
# Third, the tss enrichment is calculated simply by (tss-read-count/tss-length)/(promoter-read-count/promoter-length) for each
# tss in each cell. The final output will be a cell by tss-enrichment matrix.
# 

# get opts 
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
# -a sample name; e.g., GSMxxxxx

while getopts s:i:o:r:c:t:p:n:v:b:a: flag
do
    case "${flag}" in
        s) script_dir=${OPTARG};;
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
        r) refGene=${OPTARG};;
        c) chrom_size=${OPTARG};;
        t) tss_size=${OPTARG};;
        p) promoter_size=${OPTARG};;
        n) python=${OPTARG};;
        v) tsv=${OPTARG};;
        b) black_list=${OPTARG};;
        a) sample_name=${OPTARG};;
    esac
done

# echo for directories

echo "script directory: $script_dir"
echo "input directory: $input_dir"
echo "output directory: $output_dir"
echo "sample: $sample_name"

# create output dir

mkdir -p $output_dir

# step 1: create tss/promoter window file
echo "using refGene file: $refGene"
echo "creating a bed file with gene annotation"

gunzip -c ${input_dir}$refGene | awk 'BEGIN {FS="\t"; OFS="\t"} {print $3, $5, $6, $13, $13, $4}' | uniq | grep -v "_" > ${output_dir}genes.bed
echo "done"

echo "TSS size is defined as $tss_size bp and Promoter size is defined as $promoter_size bp"
echo "create initial window files with $tss_size bp and $promoter_size bp"
echo "using chomosome size file: $chrom_size"

bedtools flank -i ${output_dir}genes.bed -g ${input_dir}$chrom_size -l $tss_size -r 0 -s | uniq | bedtools sort -i | uniq > ${output_dir}genes.tss.bed
bedtools flank -i ${output_dir}genes.bed -g ${input_dir}$chrom_size -l $promoter_size -r 0 -s | uniq | bedtools sort -i | uniq > ${output_dir}genes.promoter.bed

# remove the cases that multiple genes use the same tss (usually the same gene (some are non-coding RNAs) get duplicated annotations); leave only one of them
# remove the genes with only one tss (no alternative)

echo "using python: $python"
echo "get genes with more than one TSS in $refGene"
$python ${script_dir}FilterBed.py -d $output_dir -i "genes.tss.bed" -n "tss" 
$python ${script_dir}FilterBed.py -d $output_dir -i "genes.promoter.bed" -n "promoter"

# step 2: intersect the tss/promoter windows with ATAC fragment data

unzipped=${output_dir}fragment.tsv
zcat ${input_dir}$tsv > $unzipped

# filter window files to remove regions with blacklist
# then create the tss/promoter matrix for the genes that have more than one tss
for item in tss promoter
do
echo "process sample $item"
window_file=${output_dir}${item}.dup.bed
bedtools intersect -a $window_file -b ${input_dir}$black_list -v > ${output_dir}${item}.dup.blackfilter.bed
bedtools intersect -a $unzipped \
                   -b ${output_dir}${item}.dup.blackfilter.bed \
                   -wa -wb | \
                   cut -f4,13,5 | \
                   awk '{ k = $1 OFS $3 } { sum[k] += $2; count[k]++ } END{ for (i in sum) if (count[i] >= 1) print i, sum[i] }' > ${output_dir}${sample_name}.${item}.dup.mtx.tsv
done

rm $unzipped
rm ${output_dir}promoter.dup.blackfilter.bed
rm ${output_dir}tss.dup.blackfilter.bed

# step 3: calculate tss enrichment

echo "using python: $python"
echo "get genes with more than one TSS in $refGene"
$python ${script_dir}ConvertMtx.py \
    -d $output_dir \
    -t ${sample_name}.tss.dup.mtx.tsv \
    -p ${sample_name}.promoter.dup.mtx.tsv \
    -a $tss_size \
    -b $promoter_size \
    -s 500 \
    -r 1000 \
    -c 10000\
    -n ${sample_name}_dup_tss_enrichment.csv\

# remove intermediate files
rm ${output_dir}genes.bed
rm ${output_dir}genes.promoter.bed
rm ${output_dir}genes.tss.bed
rm ${output_dir}promoter.dup.bed
rm ${output_dir}promoter.unique.bed
rm ${output_dir}tss.dup.bed
rm ${output_dir}tss.unique.bed
echo "Finished"
#
#

