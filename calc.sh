#!/bin/bash

# STELAR-X Triplet Score Calculator
# Usage: ./calc.sh <gene_trees_file> <species_tree_file> [computation_mode]
#
# Arguments:
#   gene_trees_file   - File containing gene trees in Newick format
#   species_tree_file - File containing the species tree to score
#   computation_mode  - CPU_SINGLE, CPU_PARALLEL (default), GPU_PARALLEL
#
# Examples:
#   ./calc.sh all_gt_angio.tre s_angio.tre
#   ./calc.sh gene_trees.tre species.tre CPU_PARALLEL

# Check for help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "STELAR-X Triplet Score Calculator"
    echo ""
    echo "Calculates the triplet score between gene trees and a given species tree."
    echo "This is a 'score-only' mode - no inference/DP is performed."
    echo ""
    echo "Usage: $0 <gene_trees_file> <species_tree_file> [computation_mode]"
    echo ""
    echo "Arguments:"
    echo "  gene_trees_file   File containing gene trees in Newick format"
    echo "  species_tree_file File containing the species tree to score"
    echo "  computation_mode  CPU_SINGLE, CPU_PARALLEL (default), GPU_PARALLEL"
    echo ""
    echo "Examples:"
    echo "  $0 all_gt_angio.tre s_angio.tre"
    echo "  $0 gene_trees.tre species.tre CPU_PARALLEL"
    echo "  $0 gene_trees.tre inferred_tree.tre GPU_PARALLEL"
    exit 0
fi

# Default configuration
DEFAULT_COMPUTATION_MODE="CPU_PARALLEL"

# Get arguments
GENE_TREES_FILE=${1:-}
SPECIES_TREE_FILE=${2:-}
COMPUTATION_MODE=${3:-$DEFAULT_COMPUTATION_MODE}

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Validate arguments
if [[ -z "$GENE_TREES_FILE" || -z "$SPECIES_TREE_FILE" ]]; then
    echo -e "${RED}Error: Both gene trees file and species tree file are required.${NC}"
    echo ""
    echo "Usage: $0 <gene_trees_file> <species_tree_file> [computation_mode]"
    echo "Run '$0 --help' for more information."
    exit 1
fi

# Check if files exist
if [[ ! -f "$GENE_TREES_FILE" ]]; then
    echo -e "${RED}Error: Gene trees file '$GENE_TREES_FILE' does not exist.${NC}"
    exit 1
fi

if [[ ! -f "$SPECIES_TREE_FILE" ]]; then
    echo -e "${RED}Error: Species tree file '$SPECIES_TREE_FILE' does not exist.${NC}"
    exit 1
fi

# Validate computation mode
VALID_MODES=("CPU_SINGLE" "CPU_PARALLEL" "GPU_PARALLEL")
if [[ ! " ${VALID_MODES[@]} " =~ " ${COMPUTATION_MODE} " ]]; then
    echo -e "${RED}Error: Invalid computation mode '$COMPUTATION_MODE'${NC}"
    echo "Valid modes: ${VALID_MODES[*]}"
    exit 1
fi

# Check if binaries exist
if [[ ! -f "target/stelar-x-1.0.0-SNAPSHOT.jar" ]]; then
    echo -e "${RED}Error: JAR file not found. Please run build.sh first.${NC}"
    exit 1
fi

echo "=== STELAR-X Triplet Score Calculator ==="
echo "Gene trees file:   $GENE_TREES_FILE"
echo "Species tree file: $SPECIES_TREE_FILE"
echo "Computation mode:  $COMPUTATION_MODE"
echo

# Run the score calculation
echo -e "${YELLOW}Calculating triplet score...${NC}"
echo

java -Xms4g -Xmx128g \
    -Djava.library.path="$(pwd)/cuda" \
    -Djna.debug_load=false \
    -Djna.platform.library.path="$(pwd)/cuda" \
    -cp target/stelar-x-1.0.0-SNAPSHOT.jar \
    Main -i "$GENE_TREES_FILE" -c "$SPECIES_TREE_FILE" -m "$COMPUTATION_MODE"

EXIT_CODE=$?

if [[ $EXIT_CODE -ne 0 ]]; then
    echo -e "${RED}Score calculation failed!${NC}"
    exit $EXIT_CODE
fi

echo -e "${GREEN}Score calculation completed successfully!${NC}"

