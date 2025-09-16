import argparse
from cellbender.remove_background.downstream import anndata_from_h5
import scanpy as sc

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input cellbender filtered .h5 file')
    parser.add_argument('--output', required=True, help='Output .h5ad file')
    args = parser.parse_args()

    adata = anndata_from_h5(args.input)
    adata.write_h5ad(args.output)

if __name__ == "__main__":
    main()
