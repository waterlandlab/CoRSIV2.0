#!/bin/bash

# Directory containing the original BED files
DIR="query"
# Reference genome sizes file
GENOME_FILE="$DIR/hg38.chrom.sizes.txt"

# Number of random intervals to generate
N=4000

echo "Generating random background regions..."
echo "-----------------------------------"

# File names
black_bed="$DIR/Black_CpGge5_Rangegt20.bed"
white_bed="$DIR/White_CpGge5_Rangegt20.bed"

if [[ ! -f "$GENOME_FILE" ]]; then
    echo "Error: Genome file '$GENOME_FILE' not found."
    exit 1
fi

# Function to calculate average length and generate random bed file
generate_random_bed() {
    local input_bed="$1"
    local prefix="$2"
    
    if [[ ! -f "$input_bed" ]]; then
        echo "Error: Input file '$input_bed' not found."
        return
    fi
    
    # Calculate average region size and round it to nearest integer
    # Rounding is necessary because bedtools random -l requires an integer
    local avg_len=$(awk '{sum += $3 - $2} END {if (NR > 0) printf "%.0f", sum / NR; else print 0}' "$input_bed")
    
    local output_bed="$DIR/${prefix}_random_background.bed"
    
    echo "Processing $prefix..."
    echo "  - Average region length: $avg_len bp"
    echo "  - Generating $N random regions..."
    
    # Generate random bed using bedtools
    # We use awk to keep only the first 3 columns (chr, start, end) matching typical BED format
    bedtools random -l "$avg_len" -n "$N" -g "$GENOME_FILE" | awk '{print $1"\t"$2"\t"$3}' > "$output_bed"
    
    echo "  - Saved to: $output_bed"
}

# Generate random backgrounds for Black regions
generate_random_bed "$black_bed" "Black"

# Generate random backgrounds for White regions
generate_random_bed "$white_bed" "White"

echo "-----------------------------------"
echo "Random background generation complete."
