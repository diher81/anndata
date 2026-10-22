
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import scanpy as sc
import pandas as pd
import pyranges as pr
from scipy import sparse

base_dir = "/home/dirk/NOBINFBACKUP/data/gastruloid_tchic_counts"
gtf_path = "/home/dirk/NOBINFBACKUP/references/mm39/gencode.vM39.annotation.gtf.gz"

# --- LOADING DATA ---
adata = sc.read_h5ad(
    f"{base_dir}/gast_H3K27me3_2000_bins_ctype.h5ad"
)

# --- EXPLORING ANNDATA OBJECT ---

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

# Convert adata.var into a dataframe with explicit columns
var_df = adata.var.copy()

# --- SAVE DATAFRAME AS CSV (successfully created once)---

# var_df.to_csv(
#     "/home/dirk/NOBINFBACKUP/data/gastruloid_tchic_counts/var.csv",
#     index=True,
#     index_label="var_index"
# )
#
# print("var_df.head:\n", var_df.head())

var_df.index.name = "var_index"
var_df = var_df.reset_index()
# Standardize coordinate column names for PyRanges
intervals_df = var_df.rename(columns={
    "chrom": "Chromosome",
    "start": "Start",
    "end": "End"
})

intervals = pr.PyRanges(intervals_df)

# --- LOAD MOUSE REFERENCE GENOME ---
# Load the GTF annotation
gtf_columns = [
    "Chromosome", "Feature", "Start", "End",
    "Score", "Strand", "Attribute"
]

gene_chunks = []

for chunk in pd.read_csv(
    gtf_path,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "Chromosome", "Source", "Feature", "Start", "End",
        "Score", "Strand", "Frame", "Attribute"
    ],
    usecols=[0, 2, 3, 4, 6, 8],
    compression="infer",
    chunksize=50000
):
    # Discard transcripts and exons immediately
    chunk = chunk.loc[chunk["Feature"] == "gene"].copy()

    if chunk.empty:
        continue

    # Extract gene identifiers and names from GTF attributes
    chunk["gene_id"] = chunk["Attribute"].str.extract(
        r'gene_id "([^"]+)"'
    )
    chunk["gene_name"] = chunk["Attribute"].str.extract(
        r'gene_name "([^"]+)"'
    )

    # GTF uses 1-based inclusive coordinates.
    # Convert Start to PyRanges' 0-based convention.
    chunk["Start"] = chunk["Start"] - 1

    gene_chunks.append(
        chunk[[
            "Chromosome", "Start", "End", "Strand",
            "gene_id", "gene_name"
        ]]
    )

# Combine only the gene records
genes_df = pd.concat(gene_chunks, ignore_index=True)
genes = pr.PyRanges(genes_df)

print("Annotated genes:", len(genes))
print(genes.df.head())

print("Input intervals:", len(intervals))
print("Annotated genes:", len(genes))
print(intervals.df.head())
print(genes.df.head())

# --- FIND OVERLAPPING GENES ---

overlaps = intervals.join(genes)

overlaps_df = overlaps.df

print("Number of interval-gene overlap pairs:", len(overlaps_df))
print(overlaps_df.head())

overlaps_df.to_csv(
    f"{base_dir}/interval_gene_overlaps_mm39.csv",
    index=False
)

# --- Count how many genes overlap each interval ---

interval_gene_counts = intervals.count_overlaps(genes).df

interval_gene_counts.to_csv(
    f"{base_dir}/interval_gene_counts_mm39.csv",
    index=False
)
