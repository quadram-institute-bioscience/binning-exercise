#!/usr/bin/env bash
# De novo metagenome assembly with MEGAHIT
# Input  : reads/sample{1,2,3}_R{1,2}.fastq.gz
# Output : assemblies/sample{1,2,3}.contigs.fasta
# Usage  : run from project root (binning/)

set -euo pipefail

# ── Apple Silicon workaround ──────────────────────────────────────────────────
# MEGAHIT v1.2.9 on ARM64 macOS crashes (SIGSEGV, exit -11) when using
# mmap-backed memory (the default --mem-flag 1). Two flags fix this:
#   --no-hw-accel  → force megahit_core base binary (no POPCNT/BMI2 paths)
#   --mem-flag 0   → use malloc instead of mmap for SdBG construction
# ─────────────────────────────────────────────────────────────────────────────

READS="reads"
OUTDIR="assemblies"
THREADS=8
MIN_CONTIG=1000   # discard short contigs (standard for binning)

mkdir -p "${OUTDIR}"

for sample in sample1 sample2 sample3; do
    r1="${READS}/${sample}_R1.fastq.gz"
    r2="${READS}/${sample}_R2.fastq.gz"
    tmpdir="${OUTDIR}/${sample}_megahit_tmp"
    final="${OUTDIR}/${sample}.contigs.fasta"

    echo "[${sample}] Assembling..."

    megahit \
        -1 "${r1}" \
        -2 "${r2}" \
        --out-dir  "${tmpdir}" \
        --out-prefix "${sample}" \
        --min-contig-len ${MIN_CONTIG} \
        --num-cpu-threads ${THREADS} \
        --no-hw-accel \
        --mem-flag 0 \
        --verbose

    cp "${tmpdir}/${sample}.contigs.fa" "${final}"
    rm -rf "${tmpdir}"

    n=$(grep -c "^>" "${final}")
    total=$(awk '/^>/{next}{sum+=length($0)}END{printf "%.1f Mb", sum/1e6}' "${final}")
    echo "  → ${final}  (${n} contigs, ${total})"
done

echo ""
echo "Done. Assemblies:"
ls -lh "${OUTDIR}"/*.contigs.fasta
