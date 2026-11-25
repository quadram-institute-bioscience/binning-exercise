#!/usr/bin/env bash
# Compute coverage and run TNF-based binning for all assemblies × all samples.
#
# Pipeline
# --------
#   1. bamcountrefs (--trimmed-mean) on every BAM in backmapping/{ref}/{sample}.bam
#      → counts/{ref}/{sample}.tsv   (Contig  trimmed_mean)
#
#   2. tnf.py for every (ref, cov_sample) pair
#      → binning/{ref}_{cov_sample}_plot.png
#      → binning/{ref}_{cov_sample}_clusters.tsv
#
# Coverage metric: trimmed mean depth (5th–95th percentile positions)
#   • length-normalised  → comparable across contigs of different sizes
#   • robust             → not inflated by edge artefacts or mapping pile-ups
#   • what MetaBAT2 uses internally
#
# Usage: run from project root (binning/)

set -euo pipefail

ASSEMBLIES="assemblies"
BACKMAPPING="backmapping"
COUNTS="counts"
OUTDIR="binning"
THREADS=4
MIN_CONTIG_LEN=1500   # skip short contigs (unreliable TNF + coverage estimates)

SAMPLES=(sample1 sample2 sample3)

mkdir -p "${OUTDIR}"

# ── 1. coverage ───────────────────────────────────────────────────────────────
echo "=== Computing trimmed-mean coverage ==="

for ref in "${SAMPLES[@]}"; do
    mkdir -p "${COUNTS}/${ref}"
    for cov_sample in "${SAMPLES[@]}"; do
        bam="${BACKMAPPING}/${ref}/${cov_sample}.bam"
        tsv="${COUNTS}/${ref}/${cov_sample}.tsv"

        if [[ ! -f "${bam}" ]]; then
            echo "  [skip] ${bam} not found"
            continue
        fi

        echo "  ${cov_sample} → ${ref}  →  ${tsv}"
        bamcountrefs \
            --trimmed-mean \
            --threads "${THREADS}" \
            --mapq 20 \
            --proper-pairs \
            "${bam}" > "${tsv}"
    done
done

# ── 2. TNF + clustering ───────────────────────────────────────────────────────
echo ""
echo "=== TNF binning ==="

for ref in "${SAMPLES[@]}"; do
    fasta="${ASSEMBLIES}/${ref}.contigs.fasta"
    if [[ ! -f "${fasta}" ]]; then
        echo "  [skip] ${fasta} not found"
        continue
    fi

    for cov_sample in "${SAMPLES[@]}"; do
        tsv="${COUNTS}/${ref}/${cov_sample}.tsv"
        if [[ ! -f "${tsv}" ]]; then
            echo "  [skip] ${tsv} not found"
            continue
        fi

        echo "  assembly=${ref}  coverage=${cov_sample}"
        python scripts/tnf.py \
            "${fasta}" \
            --coverage  "${tsv}" \
            --clusters  auto \
            --min-length "${MIN_CONTIG_LEN}" \
            --output    "${OUTDIR}/${ref}_${cov_sample}"
    done
done

echo ""
echo "Done."
echo "  Plots   : ${OUTDIR}/*_plot.png"
echo "  Clusters: ${OUTDIR}/*_clusters.tsv"
