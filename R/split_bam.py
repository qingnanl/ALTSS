import pysam
import csv
import sys
#python split_bam.py a b c d 
args = sys.argv
bam_file = args[1]
cluster_file = args[2]
sample = args[3]
outdir = args[4]

cluster_dict = {}
with open(cluster_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    #skip header
    header = next(csv_reader)
    for row in csv_reader:
        cluster_dict[row[0]] = row[1]
clusters = set(x for x in cluster_dict.values())

fin = pysam.AlignmentFile(bam_file, "rb")

n = 0
m = 0
# open the number of bam files as the same number of clusters, and map the out file handler to the cluster id, write to a bam with wb
fouts_dict = {}
for cluster in clusters:
    fout = pysam.AlignmentFile(outdir + cluster + "-"+ sample + ".bam", "wb", template = fin)
    fouts_dict[cluster] = fout

for read in fin:
    m += 1
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    # the bam files may contain reads not in the final clustered barcodes
    # will be None if the barcode is not in the clusters.csv file
    else:
        continue
    cluster_id = cluster_dict.get(cell_barcode)
    if cluster_id:
        fouts_dict[cluster_id].write(read)
        n += 1
## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()

print(m,n)

