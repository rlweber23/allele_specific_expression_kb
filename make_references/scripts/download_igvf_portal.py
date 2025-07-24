import os
import argparse
import pandas as pd
import igvf_utils
from igvf_utils.connection import Connection

# --- Parse command-line arguments ---
parser = argparse.ArgumentParser(description="Download a file from IGVF portal using accession ID.")
parser.add_argument("--acc", required=True, help="Accession ID to download (e.g., IGVFFI9282QLXO)")
parser.add_argument("--outdir", required=True, help="Output directory to save the downloaded file")
args = parser.parse_args()

# --- Environment credentials ---
api_key = os.environ['IGVF_API_KEY']
secret_key = os.environ['IGVF_SECRET_KEY']
conn = Connection("prod")

# --- Download ---
conn.download(rec_id=args.acc, directory=args.outdir)
