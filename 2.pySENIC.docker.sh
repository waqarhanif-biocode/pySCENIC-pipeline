#!/bin/bash

# ##############################################################################
# Shell Script to Run pySCENIC using Docker (Modular)
# ##############################################################################
#
# This script automates the three main steps of the pySCENIC workflow:
# 1. GRN inference (pyscenic grn)
# 2. Regulon prediction (pyscenic ctx)
# 3. AUCell scoring (pyscenic aucell)
#
# It uses a Docker container to ensure a consistent environment and can run
# all steps, individual steps, or specific sequences of steps.
#
# USAGE:
#   ./pyscenic_docker.sh <step>
#
#   <step> can be one of:
#     grn         - Run only GRN inference
#     ctx         - Run only Regulon prediction (requires GRN output)
#     aucell      - Run only AUCell scoring (requires Regulon output)
#     grn_ctx     - Run GRN inference then Regulon prediction
#     ctx_aucell  - Run Regulon prediction then AUCell scoring (requires GRN output for ctx part)
#     all         - Run all steps sequentially (grn -> ctx -> aucell) (default if no step is provided)
#
# IMPORTANT:
# - Ensure Docker is installed and running.
# - Pull the pySCENIC Docker image (e.g., aertslab/pyscenic).
# - Update the HOST_* path variables below to match your system.
# - Ensure your input expression matrix is in genes x cells TSV format.
#
# ##############################################################################

# --- Configuration ---

# Name or ID of the pySCENIC Docker image
DOCKER_IMAGE="aertslab/pyscenic:0.12.1" # User-specified version

# Number of CPUs to use within the Docker container for pySCENIC
N_CPU=16 # User-specified

# --- Host Machine Paths ---
HOST_DATA_DIR="/home/waqar/Documents/pySCENIC/raw_feature_bc_matrix/expressionData"
HOST_RESOURCES_DIR="/home/waqar/Documents/pySCENIC/raw_feature_bc_matrix/auxiliary_files"
HOST_OUTPUT_DIR="/home/waqar/Documents/pySCENIC/raw_feature_bc_matrix/results"

# --- File Names (relative to the HOST directories defined above) ---
EXPRESSION_MTX_FNAME="expr_matrix_for_scenic.tsv"
TF_LIST_FNAME="allTFs_hg38.txt"
MOTIF_ANNOTATIONS_FNAME="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
DATABASE_GLOB_PATTERN="hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"

# --- Output File Names (will be created in HOST_OUTPUT_DIR) ---
ADJACENCIES_FNAME="adjacencies.tsv"
REGULONS_FNAME="regulons.csv" # User-specified CSV output for regulons
AUCELL_MTX_FNAME="aucell_matrix.csv"

# --- Docker Mount Paths ---
DOCKER_DATA_MOUNT="/data"
DOCKER_RESOURCES_MOUNT="/resources"
DOCKER_OUTPUT_MOUNT="/output"

# ##############################################################################
# Helper Functions
# ##############################################################################
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

run_docker_pyscenic() {
    local pyscenic_command_args=("$@")
    log "Executing Docker command: pyscenic ${pyscenic_command_args[*]}"
    docker run --rm \
        -v "${HOST_DATA_DIR}:${DOCKER_DATA_MOUNT}:ro" \
        -v "${HOST_RESOURCES_DIR}:${DOCKER_RESOURCES_MOUNT}:ro" \
        -v "${HOST_OUTPUT_DIR}:${DOCKER_OUTPUT_MOUNT}:rw" \
        --user "$(id -u):$(id -g)" \
        "${DOCKER_IMAGE}" \
        pyscenic "${pyscenic_command_args[@]}"
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        log "ERROR: Docker command failed with exit code $exit_code."
        log "Command was: pyscenic ${pyscenic_command_args[*]}"
        return $exit_code # Return error code
    fi
    log "Docker command finished successfully."
    return 0 # Return success
}

# ##############################################################################
# Script Arguments for Modularity
# ##############################################################################
STEP_TO_RUN=${1:-all} # Default to 'all' if no argument is provided

# ##############################################################################
# Main Script Logic
# ##############################################################################

log "Starting pySCENIC Docker pipeline for step: ${STEP_TO_RUN}"

# --- Preparations (common to all steps that write output) ---
if [[ "$STEP_TO_RUN" == "all" || "$STEP_TO_RUN" == "grn" || "$STEP_TO_RUN" == "ctx" || "$STEP_TO_RUN" == "aucell" || "$STEP_TO_RUN" == "grn_ctx" || "$STEP_TO_RUN" == "ctx_aucell" ]]; then
    mkdir -p "${HOST_OUTPUT_DIR}"
    log "Host output directory: ${HOST_OUTPUT_DIR}"

    # Check if essential host directories exist
    if [ ! -d "${HOST_DATA_DIR}" ]; then
        log "ERROR: Host data directory not found: ${HOST_DATA_DIR}"
        exit 1
    fi
    if [ ! -d "${HOST_RESOURCES_DIR}" ]; then
        log "ERROR: Host resources directory not found: ${HOST_RESOURCES_DIR}"
        exit 1
    fi
fi

# --- Step 1: GRN Inference (pyscenic grn) ---
run_grn() {
    log "--- Starting Step 1: GRN Inference ---"
    # Check for input files specific to GRN
    if [ ! -f "${HOST_DATA_DIR}/${EXPRESSION_MTX_FNAME}" ]; then
        log "ERROR: Expression matrix not found for GRN step: ${HOST_DATA_DIR}/${EXPRESSION_MTX_FNAME}"
        return 1
    fi
    if [ ! -f "${HOST_RESOURCES_DIR}/${TF_LIST_FNAME}" ]; then
        log "ERROR: TF list file not found for GRN step: ${HOST_RESOURCES_DIR}/${TF_LIST_FNAME}"
        return 1
    fi

    run_docker_pyscenic grn \
        "${DOCKER_DATA_MOUNT}/${EXPRESSION_MTX_FNAME}" \
        "${DOCKER_RESOURCES_MOUNT}/${TF_LIST_FNAME}" \
        -o "${DOCKER_OUTPUT_MOUNT}/${ADJACENCIES_FNAME}" \
        --num_workers "${N_CPU}"
    local status=$?
    if [ $status -eq 0 ]; then
        log "GRN inference complete. Adjacencies file: ${HOST_OUTPUT_DIR}/${ADJACENCIES_FNAME}"
    fi
    return $status
}

# --- Step 2: Regulon Prediction (pyscenic ctx) ---
run_ctx() {
    log "--- Starting Step 2: Regulon Prediction (cisTarget) ---"
    # Check for input files specific to CTX
    if [ ! -f "${HOST_OUTPUT_DIR}/${ADJACENCIES_FNAME}" ]; then
        log "ERROR: Adjacencies file not found for CTX step: ${HOST_OUTPUT_DIR}/${ADJACENCIES_FNAME}. Run GRN step first or ensure it exists."
        return 1
    fi
    if [ ! -f "${HOST_RESOURCES_DIR}/${MOTIF_ANNOTATIONS_FNAME}" ]; then
        log "ERROR: Motif annotations file not found for CTX step: ${HOST_RESOURCES_DIR}/${MOTIF_ANNOTATIONS_FNAME}"
        return 1
    fi
    if [ ! -f "${HOST_DATA_DIR}/${EXPRESSION_MTX_FNAME}" ]; then # Also needs expression matrix
        log "ERROR: Expression matrix not found for CTX step: ${HOST_DATA_DIR}/${EXPRESSION_MTX_FNAME}"
        return 1
    fi

    DB_FILES_HOST_PATHS=(${HOST_RESOURCES_DIR}/${DATABASE_GLOB_PATTERN})
    if [ ${#DB_FILES_HOST_PATHS[@]} -eq 0 ] || [ ! -e "${DB_FILES_HOST_PATHS[0]}" ]; then
        log "ERROR: No cisTarget database files found in ${HOST_RESOURCES_DIR} matching pattern ${DATABASE_GLOB_PATTERN}"
        return 1
    fi
    DB_FILES_CONTAINER_PATHS=()
    for db_file_host_path in "${DB_FILES_HOST_PATHS[@]}"; do
        db_filename=$(basename "${db_file_host_path}")
        DB_FILES_CONTAINER_PATHS+=("${DOCKER_RESOURCES_MOUNT}/${db_filename}")
    done
    log "Found cisTarget databases (host): ${DB_FILES_HOST_PATHS[*]}"
    log "cisTarget databases (container paths for ctx): ${DB_FILES_CONTAINER_PATHS[*]}"

    run_docker_pyscenic ctx \
        "${DOCKER_OUTPUT_MOUNT}/${ADJACENCIES_FNAME}" \
        "${DB_FILES_CONTAINER_PATHS[@]}" \
        --annotations_fname "${DOCKER_RESOURCES_MOUNT}/${MOTIF_ANNOTATIONS_FNAME}" \
        --expression_mtx_fname "${DOCKER_DATA_MOUNT}/${EXPRESSION_MTX_FNAME}" \
        -o "${DOCKER_OUTPUT_MOUNT}/${REGULONS_FNAME}" \
        --mode "dask_multiprocessing" \
        --num_workers "${N_CPU}"
    local status=$?
    if [ $status -eq 0 ]; then
        log "Regulon prediction complete. Regulons file: ${HOST_OUTPUT_DIR}/${REGULONS_FNAME}"
    fi
    return $status
}

# --- Step 3: AUCell Scoring (pyscenic aucell) ---
run_aucell() {
    log "--- Starting Step 3: AUCell Scoring ---"
    # Check for input files specific to AUCell
    if [ ! -f "${HOST_DATA_DIR}/${EXPRESSION_MTX_FNAME}" ]; then
        log "ERROR: Expression matrix not found for AUCell step: ${HOST_DATA_DIR}/${EXPRESSION_MTX_FNAME}"
        return 1
    fi
    if [ ! -f "${HOST_OUTPUT_DIR}/${REGULONS_FNAME}" ]; then
        log "ERROR: Regulons file ('${REGULONS_FNAME}') not found for AUCell step: ${HOST_OUTPUT_DIR}/${REGULONS_FNAME}. Run CTX step first or ensure it exists."
        return 1
    fi

    run_docker_pyscenic aucell \
        "${DOCKER_DATA_MOUNT}/${EXPRESSION_MTX_FNAME}" \
        "${DOCKER_OUTPUT_MOUNT}/${REGULONS_FNAME}" \
        -o "${DOCKER_OUTPUT_MOUNT}/${AUCELL_MTX_FNAME}" \
        --num_workers "${N_CPU}"
    local status=$?
    if [ $status -eq 0 ]; then
        log "AUCell scoring complete. AUCell matrix: ${HOST_OUTPUT_DIR}/${AUCELL_MTX_FNAME}"
    fi
    return $status
}


# --- Execute selected step(s) ---
EXIT_CODE=0
case "$STEP_TO_RUN" in
    grn)
        run_grn
        EXIT_CODE=$?
        ;;
    ctx)
        run_ctx
        EXIT_CODE=$?
        ;;
    aucell)
        run_aucell
        EXIT_CODE=$?
        ;;
    grn_ctx)
        run_grn && run_ctx
        EXIT_CODE=$?
        ;;
    ctx_aucell)
        run_ctx && run_aucell
        EXIT_CODE=$?
        ;;
    all)
        run_grn && run_ctx && run_aucell
        EXIT_CODE=$? # Will be the exit code of the first command that fails, or 0 if all succeed
        ;;
    *)
        log "ERROR: Invalid step '${STEP_TO_RUN}'."
        log "Usage: $0 [grn|ctx|aucell|grn_ctx|ctx_aucell|all]"
        EXIT_CODE=1
        ;;
esac

if [ $EXIT_CODE -eq 0 ]; then
    log "pySCENIC Docker pipeline (step: ${STEP_TO_RUN}) finished successfully!"
    log "All outputs are in: ${HOST_OUTPUT_DIR}"
    # List key outputs, acknowledging some might not exist if not all steps were run
    echo "---"
    echo "Key output files (will exist if corresponding step was run successfully):"
    [ -f "${HOST_OUTPUT_DIR}/${ADJACENCIES_FNAME}" ] && echo "  Adjacencies:   ${HOST_OUTPUT_DIR}/${ADJACENCIES_FNAME}"
    [ -f "${HOST_OUTPUT_DIR}/${REGULONS_FNAME}" ]    && echo "  Regulons:      ${HOST_OUTPUT_DIR}/${REGULONS_FNAME}"
    [ -f "${HOST_OUTPUT_DIR}/${AUCELL_MTX_FNAME}" ]  && echo "  AUCell Scores: ${HOST_OUTPUT_DIR}/${AUCELL_MTX_FNAME}"
    echo "---"
else
    log "pySCENIC Docker pipeline (step: ${STEP_TO_RUN}) failed with exit code ${EXIT_CODE}."
fi

exit $EXIT_CODE
