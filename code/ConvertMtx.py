# conda activate scanpy
import pandas as pd
import numpy as np
import argparse

# parse cmdline arguments
parser = argparse.ArgumentParser(description='get arguments')
parser.add_argument('-d', '--directory', dest = "directory", help='directory that stores the input mtx file')
parser.add_argument('-t', '--input_tss_mtx', dest = "tmtx", help='input tss mtx file')
parser.add_argument('-p', '--input_promoter_mtx', dest = "pmtx", help='input promoter mtx file')
parser.add_argument('-a', '--tss_length', dest = "trange", help='tss length(50/100 bp?)')
parser.add_argument('-b', '--promoter_length', dest = "prange", help='promoter length (2000 bp?)')
parser.add_argument('-s', '--tss_thresh', dest = "t_thresh", help='threshold for the number of gene (with at least one tss fragment) found in a cell')
parser.add_argument('-r', '--promoter_thresh', dest = "p_thresh", help='threshold for the number of gene (with at least one promoter fragment) found in a cell')
parser.add_argument('-c', '--scale', dest = "scale", help='scaling factor for tss enrichment')
parser.add_argument('-n', '--name', dest = "name", help='name for output file (e.g. dup_tss_enrichment.csv)')
args = parser.parse_args()

rawdir = args.directory
tmtx = args.tmtx
pmtx = args.pmtx
tthresh = args.t_thresh
pthresh = args.p_thresh
trange = args.trange
prange = args.prange
scale = args.scale
name = args.name
# read bed files
tss = pd.read_csv(rawdir + tmtx, sep=" ",
                       names=["Cell", "Tss", "Count"])
promoter = pd.read_csv(rawdir + pmtx, sep=" ",
                       names=["Cell", "Tss", "Count"])
                       
def ConvertMat(data, thresh):
    cell_cover = data.Cell.value_counts()
    cell_retain = cell_cover[cell_cover > thresh].index
    data = data.loc[data['Cell'].isin(cell_retain)]
    mtx = data.pivot(index='Cell', columns='Tss', values='Count')
    mtx.replace(np.nan,0,inplace=True)
    return mtx

tmat = ConvertMat(data=tss, thresh=int(tthresh))
pmat = ConvertMat(data=promoter, thresh=int(pthresh))


def ComputeTssEnrichment(tss, promoter, scale, trange, prange):
    common_cells = tss.index.intersection(promoter.index)
    common_tss = tss.columns.intersection(promoter.columns)
    tss_filter = tss.loc[common_cells, common_tss]
    tss_trans = np.log1p((tss_filter*scale)/trange)
    promoter_filter = promoter.loc[common_cells, common_tss]
    promoter_trans = np.log1p((tss_filter*scale)/prange)
    return tss_trans.subtract(promoter_trans)

dup_enrichment = ComputeTssEnrichment(tss=tmat, promoter=pmat, trange=int(trange), prange=int(prange), scale=int(scale))
dup_enrichment.to_csv(rawdir+ name)

