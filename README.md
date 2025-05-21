# Single-Cell Gene Regulatory Network Analysis with pySCENIC

This repository contains a suite of scripts to perform and analyze gene regulatory networks from single-cell RNA-sequencing data using pySCENIC. The workflow is divided into three main parts:
1.  **Expression Matrix Preparation**: Convert 10x Genomics output into a format suitable for pySCENIC.
2.  **pySCENIC Workflow Execution**: Run the core pySCENIC steps (GRN inference, regulon prediction, and AUCell scoring) using a Docker container for reproducibility.
3.  **Downstream Analysis and Visualization**: Analyze and visualize the pySCENIC results, including regulon activity scores and integration with single-cell metadata.

---

## Table of Contents

* [Prerequisites](#prerequisites)
* [Setup](#setup)
* [Workflow & Scripts](#workflow--scripts)
    * [Step 1: Prepare Expression Matrix (`1.expression_to_csv.py`)](#step-1-prepare-expression-matrix-1expression_to_csvpy)
    * [Step 2: Run pySCENIC Workflow (`2.pySENIC.docker.sh`)](#step-2-run-pyscenic-workflow-2pysenicdockersh)
    * [Step 3: Downstream Analysis (`3.post-Scenic.ipynb`)](#step-3-downstream-analysis-3post-scenicipynb)
* [Input Data Requirements](#input-data-requirements)
* [Suggested Directory Structure](#suggested-directory-structure)
* [Contributing](#contributing)
* [License](#license)

---

## Prerequisites

Before you begin, ensure you have the following installed:

* **Python 3.7+**: For running the expression matrix preparation script and the Jupyter Notebook.
* **Docker**: For running the pySCENIC workflow in a containerized environment. Ensure the Docker daemon is running.
* **Git**: For cloning this repository.
* **Jupyter Notebook or JupyterLab**: For running the `3.post-Scenic.ipynb` notebook.

---

## Setup

1.  **Clone the Repository:**
    ```bash
    git clone <your-repository-url>
    cd <your-repository-name>
    ```

2.  **Python Packages:**
    Install the required Python packages. You can install them using pip:
    ```bash
    pip install scanpy pandas numpy matplotlib seaborn anndata loompy jupyterlab
    ```
    It's recommended to use a virtual environment.

3.  **pySCENIC Docker Image:**
    Pull the official pySCENIC Docker image (or your preferred version):
    ```bash
    docker pull aertslab/pyscenic:latest
    ```
    You can replace `:latest` with a specific version if needed. Note the image name/tag, as you'll need it for the Docker script.

---

## Workflow & Scripts

This repository provides scripts to streamline the pySCENIC analysis.

### Step 1: Prepare Expression Matrix (`1.expression_to_csv.py`) 🧬➡️📄

This Python script converts 10x Genomics sparse matrix output (MTX format) into a tab-separated value (TSV) file where rows are genes and columns are cells. This format is required for pySCENIC.

* **Purpose**: To generate a gene-by-cell expression matrix from 10x Genomics `filtered_feature_bc_matrix` output.
* **Dependencies**: `scanpy`, `pandas`, `argparse`
* **Usage**:
    ```bash
    python 1.expression_to_csv.py <mtx_dir> <output_tsv> [--gene_names <column_name>]
    ```
    **Arguments**:
    * `mtx_dir`: Path to the 10x Genomics mtx directory (e.g., `filtered_feature_bc_matrix/` which contains `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz`).
    * `output_tsv`: Path for the output TSV file (e.g., `expr_matrix_for_scenic.tsv`).
    * `--gene_names` (optional): Column from `features.tsv.gz` to use as gene names. Common choices are `gene_symbols` (default) or `gene_ids`.

    **Example**:
    ```bash
    python 1.expression_to_csv.py \
        /path/to/your/10x_filtered_feature_bc_matrix/ \
        /path/to/your/output_directory/expr_matrix_genes_x_cells.tsv \
        --gene_names gene_symbols
    ```
* **Input**: A directory containing:
    * `matrix.mtx.gz`
    * `features.tsv.gz` (or `genes.tsv.gz`)
    * `barcodes.tsv.gz`
* **Output**: A single TSV file (`output_tsv`) with genes as rows and cells as columns.

---

### Step 2: Run pySCENIC Workflow (`2.pySENIC.docker.sh`) 🐳⚙️

This shell script automates the three main steps of the pySCENIC workflow (GRN, ctx, AUCell) using a Docker container. This ensures a consistent and reproducible environment.

* **Purpose**: To execute the core pySCENIC pipeline for GRN inference, regulon prediction, and cell enrichment scoring (AUCell).
* **Prerequisites**:
    * Docker installed and running.
    * pySCENIC Docker image pulled (e.g., `aertslab/pyscenic`).
    * **Expression Matrix**: The genes-by-cells TSV file generated from Step 1.
    * **Transcription Factor (TF) List**: A plain text file listing TFs, one per line. (See [Input Data Requirements](#input-data-requirements)).
    * **Motif Annotation Files**: Database files (e.g., `.tbl` files from cisTarget) for motif-to-gene annotation. (See [Input Data Requirements](#input-data-requirements)).
* **Configuration**:
    ⚠️ **Before running, you MUST configure the following variables at the top of the `2.pySENIC.docker.sh` script:**
    * `PYSCENIC_IMAGE_NAME`: Name of the pySCENIC Docker image you pulled (e.g., `"aertslab/pyscenic:latest"`).
    * `HOST_EXPRESSION_MATRIX_FNAME`: Absolute path to your input expression matrix TSV file (output from Step 1).
    * `HOST_TF_LIST_FNAME`: Absolute path to your transcription factor list file.
    * `HOST_MOTIF_ANNOTATIONS_DIR`: Absolute path to the directory containing your motif annotation files (`.tbl` files).
    * `HOST_OUTPUT_DIR`: Absolute path to the directory where pySCENIC outputs will be saved. **Ensure this directory exists.**
    * `NUM_WORKERS`: Number of CPU cores to use for parallel processing steps.

    The script also defines output filenames (defaults are `adj.tsv`, `reg.csv`, `auc_mtx.csv`). You can change these within the script if needed.

* **Usage**:
    ```bash
    bash 2.pySENIC.docker.sh <step>
    ```
    **Arguments**:
    * `<step>`: Specifies which part(s) of the pySCENIC pipeline to run.
        * `grn`: Run only GRN inference.
        * `ctx`: Run only Regulon prediction (requires GRN output).
        * `aucell`: Run only AUCell scoring (requires Regulon output).
        * `grn_ctx`: Run GRN inference then Regulon prediction.
        * `ctx_aucell`: Run Regulon prediction then AUCell scoring (requires GRN output for `ctx`).
        * `all`: Run all steps sequentially (GRN → ctx → AUCell). This is the default if no step is provided.

    **Example**:
    ```bash
    # First, make the script executable:
    chmod +x 2.pySENIC.docker.sh

    # Run all steps:
    ./2.pySENIC.docker.sh all
    # Or, run only the GRN step:
    # ./2.pySENIC.docker.sh grn
    ```

* **Input**:
    * Configured paths in the script pointing to:
        * Expression matrix TSV.
        * TF list file.
        * Motif annotation files directory.
* **Output**:
    Files will be saved in the `HOST_OUTPUT_DIR` you specified. Key outputs include:
    * `adj.tsv` (or `${ADJACENCIES_FNAME}`): Adjacency matrix from the GRN inference step.
    * `reg.csv` (or `${REGULONS_FNAME}`): Discovered regulons from the `ctx` (motif enrichment) step.
    * `auc_mtx.csv` (or `${AUCELL_MATRIX_FNAME}`): AUCell matrix with enrichment scores for each regulon in each cell.
    * Other intermediate files and logs may also be present.

---

### Step 3: Downstream Analysis (`3.post-Scenic.ipynb`) 📊🔬

This Jupyter Notebook provides a framework for loading, analyzing, and visualizing the outputs from the pySCENIC workflow, particularly the AUCell matrix, in conjunction with the original single-cell expression data.

* **Purpose**: To integrate pySCENIC results (AUCell scores, regulons) with single-cell analysis (e.g., UMAP, clustering) for biological interpretation.
* **Dependencies**: `os`, `pandas`, `scanpy`, `numpy`, `matplotlib`, `seaborn`, `anndata`, `loompy` (if loading `.loom` AUCell files directly, though the notebook primarily uses CSV).
* **Setup**:
    ⚠️ **Before running the cells, you MUST configure the following path variables and parameters in the initial cells of the notebook:**
    * `EXPRESSION_CSV_PATH`: Path to your cell-by-gene (or gene-by-cell, the notebook can transpose) expression data in CSV format. This could be the raw counts matrix used for pySCENIC input, or a normalized one.
    * `AUCELL_PATH`: Path to the AUCell matrix CSV file generated by `2.pySENIC.docker.sh` (e.g., `auc_mtx.csv`).
    * `REGULONS_PATH` (optional): Path to the regulons CSV file (e.g., `reg.csv`).
    * `CELL_ANNOTATIONS_PATH` (optional): Path to a CSV file containing cell annotations or metadata (e.g., cell types, conditions). The CSV should have a column that can be used as an index to match cells (e.g., cell barcodes).
    * `CLUSTER_KEY`: The key in `adata.obs` that stores your cell cluster information (e.g., `'leiden'`, `'louvain'`). This is used for identifying top regulons per cluster.
    * `REGULON_KEY_IN_ADATA`: The key in `adata.obsm` where AUCell scores will be stored (default: `'X_AUCell'`).
    * `FIGURES_PATH`: Directory to save output figures.

* **Usage**:
    1.  Start Jupyter Notebook or JupyterLab:
        ```bash
        jupyter lab
        # or
        # jupyter notebook
        ```
    2.  Open `3.post-Scenic.ipynb`.
    3.  Configure the paths and parameters in the "1. Setup" section.
    4.  Run the cells sequentially.

* **Input**:
    * Expression data (CSV).
    * AUCell matrix (CSV from pySCENIC output).
    * Optionally, regulons data (CSV) and cell metadata/annotations (CSV).
* **Key Analyses Performed**:
    * Loading and preprocessing of expression data (normalization, scaling, PCA, neighbors, UMAP, clustering using Scanpy).
    * Loading AUCell scores and regulon data.
    * Aligning cells between expression data and AUCell scores.
    * Adding AUCell scores to the `AnnData` object (`adata.obsm`).
    * Visualization of AUCell scores for specific regulons on UMAP plots.
    * Heatmap visualization of regulon activities across cell clusters.
    * Identification of top active regulons per cell cluster.

---

## Input Data Requirements

* **10x Genomics Data**: Standard output from Cell Ranger (`filtered_feature_bc_matrix` directory).
* **Transcription Factor (TF) Lists**:
    * A plain text file with one gene symbol (must match those in your expression matrix) per line.
    * TFs lists can be obtained from resources like [AnimalTFDB](http://bioinfo.life.hust.edu.cn/AnimalTFDB/) or custom curated.
    * Ensure the gene identifiers match those used in your expression matrix and motif databases.
* **Motif Annotation Databases (`.tbl` files)**:
    * These are required for the `pyscenic ctx` step to link TFs to target genes via cis-regulatory motifs.
    * You can download pre-compiled cisTarget databases from resources like [resources.aertslab.org/cistarget/](https://resources.aertslab.org/cistarget/) (e.g., for human or mouse).
    * Select databases appropriate for your species and the type of gene IDs used (e.g., gene symbols or Ensembl IDs). The `2.pySENIC.docker.sh` script expects these to be in the directory specified by `HOST_MOTIF_ANNOTATIONS_DIR`.

---

## Suggested Directory Structure

Organizing your project files can be helpful. Here's a suggestion:

pyscenic_project/
│
├── data/                          # Raw and processed data
│   ├── 10x_matrix/                # Input 10x genomics data
│   │   ├── filtered_feature_bc_matrix/
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   └── ...
│   ├── expression_matrix_for_scenic.tsv # Output from 1.expression_to_csv.py
│   ├── cell_annotations.csv       # Optional cell metadata
│   └── ...
│
├── scenic_resources/              # Resources for pySCENIC
│   ├── tf_lists/
│   │   └── human_tfs.txt
│   ├── motif_databases/           # cisTarget .tbl files
│   │   ├── hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.tbl
│   │   └── ...
│   └── ...
│
├── scenic_analysis/               # Output from pySCENIC (2.pySENIC.docker.sh)
│   ├── adj.tsv
│   ├── reg.csv
│   ├── auc_mtx.csv
│   └── ...
│
├── notebooks/                     # Jupyter notebooks
│   └── 3.post-Scenic.ipynb
│
├── scripts/                       # Custom scripts
│   ├── 1.expression_to_csv.py
│   └── 2.pySENIC.docker.sh
│
├── figures/                       # Output figures from the notebook
│   └── ...
│
└── README.md                      # This file