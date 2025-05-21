import scanpy as sc
import pandas as pd
import os
import argparse

def create_pyscenic_input_from_10x(
    mtx_dir: str,
    output_tsv_path: str,
    gene_column_name: str = 'gene_symbols', # 'gene_symbols' or 'gene_ids' from features.tsv
    transpose_to_genes_x_cells: bool = False
):
    """
    Reads 10x Genomics MTX data, processes it, and saves the expression matrix
    as a TSV file (genes x cells) suitable for pySCENIC input.

    Args:
        mtx_dir (str): Path to the directory containing 10x matrix.mtx.gz,
                       features.tsv.gz (or genes.tsv.gz), and barcodes.tsv.gz.
        output_tsv_path (str): Path to save the output TSV file.
        gene_column_name (str): Column name in features.tsv.gz to use for gene names.
                                 Typically 'gene_symbols' (common) or 'gene_ids'.
                                 Scanpy's read_10x_mtx default is 'gene_symbols'.
        transpose_to_genes_x_cells (bool): If True (default), transposes the matrix
                                           to genes as rows and cells as columns.
                                           This is the format pySCENIC CLI expects.
    """
    print(f"Reading 10x data from: {mtx_dir}")
    try:
        # var_names defines which column from features.tsv to use for gene names
        adata = sc.read_10x_mtx(mtx_dir, var_names=gene_column_name, cache=True)
        adata.var_names_make_unique() # Important for unique gene names
        adata.obs_names_make_unique() # Important for unique cell barcodes
    except Exception as e:
        print(f"Error reading 10x data: {e}")
        print("Please ensure the path is correct and files (matrix.mtx.gz, features.tsv.gz/genes.tsv.gz, barcodes.tsv.gz) exist.")
        return

    print(f"Successfully loaded AnnData object with shape: {adata.shape} (cells x genes)")

    # --- Optional: Basic Filtering (Uncomment and adapt if needed) ---
    print("Applying basic filtering (min_genes=200 per cell, min_cells=3 per gene)...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"Shape after basic filtering: {adata.shape}")
    if adata.n_obs == 0 or adata.n_vars == 0:
         print("Error: All cells or genes were filtered out. Check filtering parameters and data quality.")
         return

    # pySCENIC typically works well with raw counts.
    # If adata.X is not raw counts (e.g., if it was normalized), ensure you use the correct layer.
    # For this script, we assume adata.X contains the desired counts (usually raw).

    # Convert to DataFrame
    # AnnData stores expression data as cells x genes.
    if isinstance(adata.X, pd.DataFrame): # Should not happen with read_10x_mtx
        expr_df = adata.X
    elif hasattr(adata.X, "toarray"): # Sparse matrix
        expr_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
    else: # Dense numpy array
        expr_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    print(f"Expression matrix DataFrame shape: {expr_df.shape} (cells x genes)")

    if transpose_to_genes_x_cells:
        print("Transposing matrix to genes x cells format...")
        expr_df_transposed = expr_df.T
        print(f"Transposed matrix shape: {expr_df_transposed.shape} (genes x cells)")
    else:
        expr_df_transposed = expr_df # Keep as cells x genes if specified

    # Ensure output directory exists
    output_dir = os.path.dirname(output_tsv_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")

    print(f"Saving expression matrix to: {output_tsv_path}")
    try:
        expr_df_transposed.to_csv(output_tsv_path, sep='\t', header=True, index=True)
        print("Successfully saved the expression matrix.")
    except Exception as e:
        print(f"Error saving TSV file: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare an expression matrix TSV file (genes x cells) from 10x Genomics data for pySCENIC."
    )
    parser.add_argument(
        "mtx_dir",
        type=str,
        help="Path to the 10x Genomics mtx directory (containing matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)."
    )
    parser.add_argument(
        "output_tsv",
        type=str,
        help="Path for the output TSV file (e.g., expr_matrix_for_scenic.tsv)."
    )
    parser.add_argument(
        "--gene_names",
        type=str,
        default="gene_symbols", # scanpy's default
        help="Column from features.tsv.gz to use as gene names (e.g., 'gene_symbols' or 'gene_ids'). Default: 'gene_symbols'."
    )

    args = parser.parse_args()

    create_pyscenic_input_from_10x(
        mtx_dir=args.mtx_dir,
        output_tsv_path=args.output_tsv,
        gene_column_name=args.gene_names
    )

    # --- Example Usage (if not running from CLI) ---
    # Replace with your actual paths if you run this part directly in a Python interpreter/notebook
    # example_mtx_dir = "path/to/your/10x_filtered_feature_bc_matrix_directory"
    # example_output_file = "path/to/your/output_directory/expr_matrix_for_scenic.tsv"
    #
    # if not os.path.exists(example_mtx_dir) and __name__ != "__main__": # Only run example if dir exists and not via CLI
    #     print(f"\nSkipping example usage: Directory '{example_mtx_dir}' not found.")
    # elif __name__ != "__main__":
    #     print("\n--- Running example usage ---")
    #     create_pyscenic_input_from_10x(example_mtx_dir, example_output_file)
    #     print("--- Example usage finished ---")
