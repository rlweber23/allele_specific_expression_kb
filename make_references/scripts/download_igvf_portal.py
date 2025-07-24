import os
import argparse
import pandas as pd
import igvf_utils
from igvf_utils.connection import Connection




api_key = os.environ['IGVF_API_KEY']
secret_key = os.environ['IGVF_SECRET_KEY']
conn = Connection("prod")

conn.download("IGVFFI9282QLXO")