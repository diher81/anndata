
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
adata_H3K27me3 = sc.read_h5ad(
    f"{base_dir}/gast_H3K27me3_2000_bins_ctype.h5ad"
)
adata_H3K4me3 = sc.read_h5ad(
    f"{base_dir}/gast_H3K4me3_2000_bins_ctype.h5ad"
)

# --- EXPLORING ANNDATA OBJECT FOR H3K27me3 ---
print("object overview:\n")
print(adata_H3K27me3)

print("\nMatrix information:")
print("Shape:", adata_H3K27me3.X.shape)
print("Storage type:", type(adata_H3K27me3.X))
print("Sparse matrix:", sparse.issparse(adata_H3K27me3.X))
print("Data type:", adata_H3K27me3.X.dtype)

print("\nObservation metadata:")
print(adata_H3K27me3.obs.head())

print("\nVariable metadata:")
print(adata_H3K27me3.var.head())

print("\nFirst 5 observation names:")
print(adata_H3K27me3.obs_names[:5])

print("\nFirst 5 variable names:")
print(adata_H3K27me3.var_names[:5])

# Convert adata.var into a dataframe with explicit columns
var_H3K27me3_df = adata_H3K27me3.var.copy()
var_H3K4me3_df = adata_H3K27me3.var.copy()

# --- SAVE DATAFRAME AS CSV (already successfully created once)---
var_H3K27me3_df.to_csv(
    "/home/dirk/NOBINFBACKUP/data/gastruloid_tchic_counts/var_H3K27me3.csv",
    index=True,
    index_label="var_index"
)
print("var_df.head:\n", var_H3K27me3_df.head())

var_H3K4me3_df.to_csv(
    "/home/dirk/NOBINFBACKUP/data/gastruloid_tchic_counts/var_H3K4me3.csv",
    index=True,
    index_label="var_index"
)
print("var_df.head:\n", var_H3K4me3_df.head())


var_H3K27me3_df.index.name = "var_index"
var_H3K27me3_df = var_H3K27me3_df.reset_index()
# Standardize coordinate column names for PyRanges
intervals_H3K27me3_df = var_H3K27me3_df.rename(columns={
    "chrom": "Chromosome",
    "start": "Start",
    "end": "End"
})
intervals_H3K27me3 = pr.PyRanges(intervals_H3K27me3_df)

var_H3K4me3_df.index.name = "var_index"
var_H3K4me3_df = var_H3K4me3_df.reset_index()
# Standardize coordinate column names for PyRanges
intervals_H3K4me3_df = var_H3K4me3_df.rename(columns={
    "chrom": "Chromosome",
    "start": "Start",
    "end": "End"
})
intervals_H3K4me3 = pr.PyRanges(intervals_H3K4me3_df)

# --- LOAD MOUSE REFERENCE GENOME ---
# Load the GTF annotation
gtf_columns = [
    "Chromosome", "Feature", "Start", "End",
    "Score", "Strand", "Attribute"
]

# Reading reference genome genes in chunks:
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

print("Annotated genes:", len(genes))                               ### 78289
print(genes.df.head())

print("Input intervals H3K27me3:", len(intervals_H3K27me3))         ### 1364141
print(intervals_H3K27me3.df.head())
print("Input intervals H3K4me3:", len(intervals_H3K4me3))           ### 1364141
print(intervals_H3K4me3.df.head())

print(genes.df.head())

# --- FIND OVERLAPPING GENES ---
overlaps_H3K27me3 = intervals_H3K27me3.join(genes)
overlaps_H3K27me3_df = overlaps_H3K27me3.df

overlaps_H3K4me3 = intervals_H3K4me3.join(genes)
overlaps_H3K4me3_df = overlaps_H3K4me3.df

print("Number of interval-gene overlap pairs (H3K27me3):", len(overlaps_H3K27me3_df))   ### 938147
print(overlaps_H3K27me3_df.head())
print("Number of interval-gene overlap pairs (H3K4me3):", len(overlaps_H3K4me3_df))     ### 938147
print(overlaps_H3K4me3_df.head())

overlaps_H3K27me3_df.to_csv(
    f"{base_dir}/interval_gene_overlaps_mm39_H3K27me3.csv",
    index=False
)
overlaps_H3K4me3_df.to_csv(
    f"{base_dir}/interval_gene_overlaps_mm39_H3K4me3.csv",
    index=False
)

# --- Count how many genes overlap each interval ---
interval_gene_counts = intervals_H3K27me3.count_overlaps(genes).df
interval_gene_counts.to_csv(
    f"{base_dir}/interval_gene_counts_mm39_H3K27me3.csv",
    index=False
)
interval_gene_counts = intervals_H3K4me3.count_overlaps(genes).df
interval_gene_counts.to_csv(
    f"{base_dir}/interval_gene_counts_mm39_H3K4me3.csv",
    index=False
)
