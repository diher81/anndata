
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import scanpy as sc
import pandas as pd
import pyranges as pr
from scipy import sparse

print("Loading dataset...")

adata = sc.read_h5ad(
    "/home/dirk/NOBINFBACKUP/data/gastruloid_tchic_counts/gast_H3K27me3_2000_bins_ctype.h5ad"
)

print("\nDataset loaded successfully.")

print("object overview:\n")
print(adata)

print("\nMatrix information:")
print("Shape:", adata.X.shape)
print("Storage type:", type(adata.X))
print("Sparse matrix:", sparse.issparse(adata.X))
print("Data type:", adata.X.dtype)

print("\nObservation metadata:")
print(adata.obs.head())

print("\nVariable metadata:")
print(adata.var.head())

print("\nFirst 5 observation names:")
print(adata.obs_names[:5])

print("\nFirst 5 variable names:")
print(adata.var_names[:5])

print("\nCopy var to csv.")
var_df = adata.var.copy()

# Save the dataframe
var_df.to_csv(
    "/home/dirk/NOBINFBACKUP/data/gastruloid_tchic_counts/var.csv",
    index=True,
    index_label="var_index"
)

print("var_df.head:\n", var_df.head())


