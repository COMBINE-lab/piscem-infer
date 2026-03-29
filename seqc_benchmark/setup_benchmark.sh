#!/usr/bin/env bash
set -euo pipefail

# SEQC/MAQC-III Benchmark Pipeline for piscem-infer
# BGI site, Illumina HiSeq 2000, 2x50bp PE
#
# Samples:
#   A = UHRR (Universal Human Reference RNA)
#   B = HBRR (Human Brain Reference RNA)
#   C = 75% A + 25% B
#   D = 25% A + 75% B
#
# Using 4 replicates per sample, first lane of each replicate.
# Each lane gives ~5-10M read pairs — sufficient for transcript quant.

BENCHDIR="seqc_benchmark"
CONDA_ENV="seqc-bench"
THREADS=8
GENCODE_VERSION="v49"
GENCODE_FA="gencode.${GENCODE_VERSION}.pc_transcripts.fa"

cd "$(dirname "$0")/.."

# ============================================================
# Step 0: Sample metadata
# ============================================================

# First SRR accession per replicate (lane 1 of flowcell 1)
# 4 replicates × 4 samples = 16 libraries
declare -A SAMPLES
SAMPLES[A_1]=SRR896663; SAMPLES[A_2]=SRR896679; SAMPLES[A_3]=SRR896695; SAMPLES[A_4]=SRR896711
SAMPLES[B_1]=SRR896743; SAMPLES[B_2]=SRR896759; SAMPLES[B_3]=SRR896775; SAMPLES[B_4]=SRR896791
SAMPLES[C_1]=SRR896823; SAMPLES[C_2]=SRR896839; SAMPLES[C_3]=SRR896855; SAMPLES[C_4]=SRR896871
SAMPLES[D_1]=SRR896903; SAMPLES[D_2]=SRR896919; SAMPLES[D_3]=SRR896935; SAMPLES[D_4]=SRR896951

# ============================================================
# Step 1: Download reference transcriptome
# ============================================================

echo "=== Step 1: Reference transcriptome ==="
mkdir -p "${BENCHDIR}/ref"
if [ ! -f "${BENCHDIR}/ref/${GENCODE_FA}" ]; then
    echo "Downloading GENCODE ${GENCODE_VERSION} protein-coding transcripts..."
    curl -L "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.${GENCODE_VERSION}.pc_transcripts.fa.gz" \
        -o "${BENCHDIR}/ref/${GENCODE_FA}.gz"
    gunzip "${BENCHDIR}/ref/${GENCODE_FA}.gz"
else
    echo "Reference already exists."
fi

# ============================================================
# Step 2: Build piscem index
# ============================================================

echo "=== Step 2: Build piscem index ==="
if [ ! -d "${BENCHDIR}/index" ]; then
    conda run -n ${CONDA_ENV} piscem build \
        -s "${BENCHDIR}/ref/${GENCODE_FA}" \
        -k 31 -m 19 --threads ${THREADS} \
        -o "${BENCHDIR}/index/gencode_${GENCODE_VERSION}"
else
    echo "Index already exists."
fi

# ============================================================
# Step 3: Download FASTQ files
# ============================================================

echo "=== Step 3: Download FASTQ files ==="
mkdir -p "${BENCHDIR}/fastq"

for key in "${!SAMPLES[@]}"; do
    srr="${SAMPLES[$key]}"
    outdir="${BENCHDIR}/fastq/${key}"
    if [ ! -f "${outdir}/seg_0.fq.gz" ]; then
        echo "  Downloading ${key} (${srr})..."
        mkdir -p "${outdir}"
        xsra dump -s -c g -o "${outdir}" "${srr}"
    else
        echo "  ${key} already downloaded."
    fi
done

# ============================================================
# Step 4: Map with piscem
# ============================================================

echo "=== Step 4: Map with piscem ==="
mkdir -p "${BENCHDIR}/mapped"

for key in "${!SAMPLES[@]}"; do
    srr="${SAMPLES[$key]}"
    r1="${BENCHDIR}/fastq/${key}/seg_0.fq.gz"
    r2="${BENCHDIR}/fastq/${key}/seg_1.fq.gz"
    outdir="${BENCHDIR}/mapped/${key}"
    if [ ! -f "${outdir}/map.rad" ]; then
        echo "  Mapping ${key}..."
        mkdir -p "${outdir}"
        conda run -n ${CONDA_ENV} piscem map-bulk \
            -i "${BENCHDIR}/index/gencode_${GENCODE_VERSION}" \
            -1 "${r1}" -2 "${r2}" \
            -t ${THREADS} \
            -o "${outdir}/map"
    else
        echo "  ${key} already mapped."
    fi
done

# ============================================================
# Step 5: Build salmon index
# ============================================================

echo "=== Step 5: Build salmon index ==="
if [ ! -d "${BENCHDIR}/salmon_index" ]; then
    conda run -n ${CONDA_ENV} salmon index \
        -t "${BENCHDIR}/ref/${GENCODE_FA}" \
        -i "${BENCHDIR}/salmon_index" \
        -p ${THREADS}
else
    echo "Salmon index already exists."
fi

# ============================================================
# Step 6: Quantify with salmon
# ============================================================

echo "=== Step 6: Quantify with salmon ==="
mkdir -p "${BENCHDIR}/quant_salmon"

for key in "${!SAMPLES[@]}"; do
    srr="${SAMPLES[$key]}"
    r1="${BENCHDIR}/fastq/${key}/seg_0.fq.gz"
    r2="${BENCHDIR}/fastq/${key}/seg_1.fq.gz"
    outdir="${BENCHDIR}/quant_salmon/${key}"
    if [ ! -f "${outdir}/quant.sf" ]; then
        echo "  Salmon quant ${key}..."
        conda run -n ${CONDA_ENV} salmon quant \
            -i "${BENCHDIR}/salmon_index" \
            -l A -1 "${r1}" -2 "${r2}" \
            -p ${THREADS} \
            --validateMappings \
            -o "${outdir}"
    else
        echo "  ${key} already quantified."
    fi
done

# ============================================================
# Step 7: Build kallisto index
# ============================================================

echo "=== Step 7: Build kallisto index ==="
if [ ! -f "${BENCHDIR}/kallisto_index/gencode.idx" ]; then
    mkdir -p "${BENCHDIR}/kallisto_index"
    conda run -n ${CONDA_ENV} kallisto index \
        -i "${BENCHDIR}/kallisto_index/gencode.idx" \
        "${BENCHDIR}/ref/${GENCODE_FA}"
else
    echo "Kallisto index already exists."
fi

# ============================================================
# Step 8: Quantify with kallisto
# ============================================================

echo "=== Step 8: Quantify with kallisto ==="
mkdir -p "${BENCHDIR}/quant_kallisto"

for key in "${!SAMPLES[@]}"; do
    srr="${SAMPLES[$key]}"
    r1="${BENCHDIR}/fastq/${key}/seg_0.fq.gz"
    r2="${BENCHDIR}/fastq/${key}/seg_1.fq.gz"
    outdir="${BENCHDIR}/quant_kallisto/${key}"
    if [ ! -f "${outdir}/abundance.tsv" ]; then
        echo "  Kallisto quant ${key}..."
        mkdir -p "${outdir}"
        conda run -n ${CONDA_ENV} kallisto quant \
            -i "${BENCHDIR}/kallisto_index/gencode.idx" \
            -o "${outdir}" \
            -t ${THREADS} \
            "${r1}" "${r2}"
    else
        echo "  ${key} already quantified."
    fi
done

# ============================================================
# Step 9: Quantify with piscem-infer (plain EM)
# ============================================================

echo "=== Step 9: Quantify with piscem-infer (plain EM) ==="
mkdir -p "${BENCHDIR}/quant_em"

for key in "${!SAMPLES[@]}"; do
    outdir="${BENCHDIR}/quant_em/${key}"
    if [ ! -f "${outdir}.quant" ]; then
        echo "  piscem-infer EM ${key}..."
        ./target/release/piscem-infer quant \
            -i "${BENCHDIR}/mapped/${key}/map" \
            -l auto \
            -o "${outdir}" \
            --num-threads ${THREADS} \
            --no-squarem
    else
        echo "  ${key} already quantified."
    fi
done

# ============================================================
# Step 10: Quantify with piscem-infer consensus-quant variants
# ============================================================

echo "=== Step 10: Consensus-quant variants ==="

# Generate manifests
for variant in cons_support sel_support sel_adapt sel_adapt_rescue; do
    manifest="${BENCHDIR}/manifest_${variant}.csv"
    if [ ! -f "${manifest}" ]; then
        echo "sample_name,condition,rad_path,output_dir" > "${manifest}"
        for key in "${!SAMPLES[@]}"; do
            sample_type="${key%%_*}"  # A, B, C, or D
            echo "${key},${sample_type},${BENCHDIR}/mapped/${key}/map,${BENCHDIR}/quant_${variant}/${key}" >> "${manifest}"
        done
    fi
done

# Sel+Support (base consensus)
manifest="${BENCHDIR}/manifest_sel_support.csv"
outcheck="${BENCHDIR}/quant_sel_support/A_1/A_1.quant"
if [ ! -f "${outcheck}" ]; then
    echo "  Sel+Support consensus-quant..."
    ./target/release/piscem-infer consensus-quant \
        -m "${manifest}" -l auto \
        --txp-selection --filter-mode support --min-ec-support 2 \
        --num-threads ${THREADS}
fi

# Sel+Adaptive
manifest="${BENCHDIR}/manifest_sel_adapt.csv"
outcheck="${BENCHDIR}/quant_sel_adapt/A_1/A_1.quant"
if [ ! -f "${outcheck}" ]; then
    echo "  Sel+Adaptive consensus-quant..."
    ./target/release/piscem-infer consensus-quant \
        -m "${manifest}" -l auto \
        --txp-selection --filter-mode support --min-ec-support 2 \
        --adaptive-ec-support \
        --num-threads ${THREADS}
fi

# Sel+Adaptive+Rescue
manifest="${BENCHDIR}/manifest_sel_adapt_rescue.csv"
outcheck="${BENCHDIR}/quant_sel_adapt_rescue/A_1/A_1.quant"
if [ ! -f "${outcheck}" ]; then
    echo "  Sel+Adaptive+Rescue consensus-quant..."
    ./target/release/piscem-infer consensus-quant \
        -m "${manifest}" -l auto \
        --txp-selection --filter-mode support --min-ec-support 2 \
        --adaptive-ec-support --condition-rescue \
        --num-threads ${THREADS}
fi

echo "=== Pipeline complete ==="
echo "Run seqc_benchmark/evaluate_titration.R to analyze results."
