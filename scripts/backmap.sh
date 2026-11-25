#!/usr/bin/env bash
# Backmap all read datasets against every assembly with BWA-MEM
# Produces a 3×3 matrix of sorted, indexed BAM files:
#   backmapping/{assembly}/{sample}.bam
# Usage: run from project root (binning/)

set -euo pipefail

READS="reads"
ASSEMBLIES="assemblies"
OUTDIR="backmapping"
THREADS=8

SAMPLES=(sample1 sample2 sample3)

mkdir -p "${OUTDIR}"

# ── index each assembly once ─────────────────────────────────────────────────
for ref_sample in "${SAMPLES[@]}"; do
    ref="${ASSEMBLIES}/${ref_sample}.contigs.fasta"
    if [[ ! -f "${ref}.bwt" ]]; then
        echo "[index] ${ref}"
        bwa index "${ref}"
    else
        echo "[index] ${ref}  (already indexed, skipping)"
    fi
done

# ── map every sample against every assembly ───────────────────────────────────
for ref_sample in "${SAMPLES[@]}"; do
    ref="${ASSEMBLIES}/${ref_sample}.contigs.fasta"
    ref_dir="${OUTDIR}/${ref_sample}"
    mkdir -p "${ref_dir}"

    for query_sample in "${SAMPLES[@]}"; do
        r1="${READS}/${query_sample}_R1.fastq.gz"
        r2="${READS}/${query_sample}_R2.fastq.gz"
        bam="${ref_dir}/${query_sample}.bam"

        echo "[map] ${query_sample} → ${ref_sample}"

        bwa mem \
            -t "${THREADS}" \
            -R "@RG\tID:${query_sample}\tSM:${query_sample}\tPL:ILLUMINA" \
            "${ref}" "${r1}" "${r2}" \
        | samtools sort \
            -@ "${THREADS}" \
            -o "${bam}"

        samtools index "${bam}"

        reads_mapped=$(samtools flagstat "${bam}" | awk '/mapped \(/{print $1}')
        echo "  → ${bam}  (${reads_mapped} reads mapped)"
    done
done

echo ""
echo "Done. BAM files:"
find "${OUTDIR}" -name "*.bam" | sort | xargs ls -lh
