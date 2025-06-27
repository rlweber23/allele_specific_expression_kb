import pandas as pd
import scanpy as sc
import argparse
import anndata as ad

parser = argparse.ArgumentParser(description='kb out directory')
parser.add_argument('--directory', '-d', default= None,required = True,  help='kb out directory')
parser.add_argument('--output', '-o', default= None,required = True,  help='outname')

args = parser.parse_args()
kb_path = args.directory
outputname = args.output

def make_adata(path):
    
    mtx = path+'cells_x_genes.total.mtx'
    cgb = path+'cells_x_genes.barcodes.txt'
    cgg = path+'cells_x_genes.genes.txt'
    cggn = path+'cells_x_genes.genes.names.txt'
    data = sc.read_mtx(mtx)
    obs = pd.read_csv(cgb, header = None, sep="\t")
    obs.columns = ["bc"]
    var = pd.read_csv(cgg, header = None, sep="\t")
    var.columns = ["gene_id"]
    genes = pd.read_csv(cggn, header = None, sep="\t")
    genes.columns = ["gene_name"]
    var['gene_name'] = genes['gene_name']

    X = data.X

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs.index = obs["bc"]
    adata.var_names = adata.var['gene_name']

    return adata

cunfmod = kb_path + '/counts_unfiltered_modified/'

adata = make_adata(cunfmod)
adata.layers['raw_counts'] = adata.X


adata.write_h5ad(outputname)