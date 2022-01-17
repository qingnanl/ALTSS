import pandas as pd
import argparse

# parse cmdline arguments
parser = argparse.ArgumentParser(description='get arguments')
parser.add_argument('-d', '--directory', dest = "directory", help='directory that stores the input bed file')
parser.add_argument('-i', '--input_bed', dest = "bed", help='input bed file')
parser.add_argument('-n', '--name', dest = "name", help='name identifier for output file (e.g. tss or promoter)')
args = parser.parse_args()

def filterbed(df, file_name):
    df['position'] = df.Chr.astype(str) + ":" + df.Start.astype(str) + "-" + df.End.astype(str) + "_" + df.Strand
    df = df.drop_duplicates(subset='position', keep='first')
    df = df.reset_index()
    lst = []
    for i in range(len(df.index)):
        strand = df['Strand'].iloc[i]
        if strand == '-':
            tss = df.Gene1.iloc[i] + '_' + df.Start.iloc[i].astype(str)
        if strand == '+':
            tss = df.Gene1.iloc[i] + '_' + df.End.iloc[i].astype(str)
        lst.append(tss)
        if i % 5000 == 0:
            print("progress report: ", i)
    df['tss'] = lst
    dupmark = df.duplicated(subset=['Gene1'], keep=False)
    dupmark.value_counts()
    df_unique = df[~dupmark]
    df_dup = df[dupmark]
    df_unique.iloc[:, 1:].to_csv(file_name + '.unique.bed', sep='\t', index=False, header=False)
    df_dup.iloc[:, 1:].to_csv(file_name + '.dup.bed', sep='\t', index=False, header=False)
#
rawdir = args.directory
bed = args.bed
name = args.name

df = pd.read_csv(rawdir + bed, sep="\t", 
                 names=["Chr", "Start", "End", "Gene1", "Gene2", "Strand"])
filterbed(df = df, file_name= rawdir + name)
