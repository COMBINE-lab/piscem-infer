#!/usr/bin/env bash
#
# End-to-end simulation pipeline for validating piscem-infer multi-quant.
#
# Prerequisites:
#   - piscem (on PATH or set PISCEM below)
#   - piscem-infer (cargo build --release, or on PATH)
#   - R with polyester installed (for simulate_reads.R)
#
# Usage:
#   bash simulation/run_pipeline.sh [sim_data_dir]
#
# The script expects simulate_reads.R to have been run first,
# producing the reference FASTA, read files, and manifest.

set -euo pipefail

SIMDIR="${1:-sim_data}"
THREADS="${THREADS:-8}"
KLEN="${KLEN:-31}"
MLEN="${MLEN:-19}"

# Tool paths (override with env vars if needed)
PISCEM="${PISCEM:-piscem}"
PISCEM_INFER="${PISCEM_INFER:-cargo run --release --}"

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== Simulation Pipeline ==="
echo "  Data directory: $SIMDIR"
echo "  Threads: $THREADS"
echo "  piscem: $PISCEM"
echo ""

# ---- Step 0: Validate inputs ----

REF_FASTA="$SIMDIR/reference.fa"
MANIFEST="$SIMDIR/manifest.csv"
SAMPLE_INFO="$SIMDIR/sample_info.csv"

if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA"
    echo "Run 'Rscript simulation/simulate_reads.R $SIMDIR' first."
    exit 1
fi

if [ ! -f "$SAMPLE_INFO" ]; then
    echo "ERROR: Sample info not found: $SAMPLE_INFO"
    echo "Run 'Rscript simulation/simulate_reads.R $SIMDIR' first."
    exit 1
fi

# ---- Step 1: Build piscem index ----

INDEX_DIR="$SIMDIR/index"
INDEX_PREFIX="$INDEX_DIR/ref"

if [ -f "${INDEX_PREFIX}.sshash" ] || [ -f "${INDEX_PREFIX}.piscem_idx" ]; then
    echo "Step 1: Piscem index already exists, skipping."
else
    echo "Step 1: Building piscem index..."
    mkdir -p "$INDEX_DIR"
    ulimit -n 2048 2>/dev/null || true
    $PISCEM build \
        -s "$REF_FASTA" \
        -k "$KLEN" \
        -m "$MLEN" \
        -t "$THREADS" \
        -o "$INDEX_PREFIX"
    echo "  Index built at $INDEX_PREFIX"
fi
echo ""

# ---- Step 2: Map reads with piscem ----

MAPPED_DIR="$SIMDIR/mapped"
mkdir -p "$MAPPED_DIR"

echo "Step 2: Mapping reads with piscem..."

# Read sample_info.csv (skip header), extract sample_name, r1_path, r2_path
tail -n +2 "$SAMPLE_INFO" | while IFS=, read -r sample_name condition r1_path r2_path; do
    # Strip quotes that R's write.csv may add
    sample_name=$(echo "$sample_name" | tr -d '"')
    condition=$(echo "$condition" | tr -d '"')
    r1_path=$(echo "$r1_path" | tr -d '"')
    r2_path=$(echo "$r2_path" | tr -d '"')

    SAMPLE_MAP_DIR="$MAPPED_DIR/$sample_name"
    SAMPLE_MAP_PREFIX="$SAMPLE_MAP_DIR/$sample_name"

    if [ -f "${SAMPLE_MAP_PREFIX}.rad" ]; then
        echo "  $sample_name: already mapped, skipping."
        continue
    fi

    echo "  Mapping $sample_name ($condition)..."
    mkdir -p "$SAMPLE_MAP_DIR"

    $PISCEM map-bulk \
        -i "$INDEX_PREFIX" \
        -1 "$r1_path" \
        -2 "$r2_path" \
        -t "$THREADS" \
        -o "$SAMPLE_MAP_PREFIX"
done

echo "  Mapping complete."
echo ""

# ---- Step 3: Run piscem-infer single-sample quant (baseline) ----

SINGLE_DIR="$SIMDIR/quant_single"
mkdir -p "$SINGLE_DIR"

echo "Step 3: Running single-sample quantification (baseline)..."

tail -n +2 "$SAMPLE_INFO" | while IFS=, read -r sample_name condition r1_path r2_path; do
    sample_name=$(echo "$sample_name" | tr -d '"')
    SAMPLE_MAP_PREFIX="$MAPPED_DIR/$sample_name/$sample_name"
    SAMPLE_OUT="$SINGLE_DIR/$sample_name/$sample_name"

    if [ -f "${SAMPLE_OUT}.quant" ]; then
        echo "  $sample_name: already quantified, skipping."
        continue
    fi

    echo "  Quantifying $sample_name (single-sample EM)..."
    mkdir -p "$SINGLE_DIR/$sample_name"

    $PISCEM_INFER quant \
        -i "$SAMPLE_MAP_PREFIX" \
        -o "$SAMPLE_OUT" \
        -l auto \
        --num-threads "$THREADS" \
        2>&1 | tail -3
done

echo "  Single-sample quantification complete."
echo ""

# ---- Step 4: Run piscem-infer multi-quant (hierarchical) ----

MULTI_DIR="$SIMDIR/quant_multi"
JOINT_DIR="$SIMDIR/quant_multi/joint"

echo "Step 4: Running multi-sample hierarchical quantification..."

# The manifest.csv from simulate_reads.R has rad_path and output_dir set
# We need to regenerate it with correct paths now that mapping is done
MULTI_MANIFEST="$SIMDIR/manifest_multi.csv"
{
    echo "sample_name,condition,rad_path,output_dir"
    tail -n +2 "$SAMPLE_INFO" | while IFS=, read -r sample_name condition r1_path r2_path; do
        sample_name=$(echo "$sample_name" | tr -d '"')
        condition=$(echo "$condition" | tr -d '"')
        echo "$sample_name,$condition,$MAPPED_DIR/$sample_name/$sample_name,$MULTI_DIR/$sample_name"
    done
} > "$MULTI_MANIFEST"

echo "  Manifest: $MULTI_MANIFEST"

$PISCEM_INFER multi-quant \
    -m "$MULTI_MANIFEST" \
    -o "$JOINT_DIR" \
    -l auto \
    --num-threads "$THREADS" \
    --num-outer-iters 7 \
    --prior-weight 0.25 \
    2>&1 | tail -10

echo "  Multi-sample quantification complete."
echo ""

# ---- Step 5: Evaluate ----

echo "Step 5: Evaluating results..."
echo "  Run: Rscript simulation/evaluate.R $SIMDIR"
echo ""

echo "=== Pipeline complete ==="
echo "Results:"
echo "  Single-sample: $SINGLE_DIR/"
echo "  Multi-sample:  $MULTI_DIR/"
echo "  Joint params:  $JOINT_DIR/"
echo "  Ground truth:  $SIMDIR/ground_truth.csv"
