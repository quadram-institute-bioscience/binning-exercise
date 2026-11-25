#!/usr/bin/env bash
# Simulate three paired-end Illumina shotgun samples from a 5-species mock community
# Simulator : ART (art_illumina), HiSeq 2500 profile, 2×150 bp, insert ~400 bp
# Usage     : run from project root (binning/)
#
# Sample design — each number is the fold coverage for that species:
#
#   Species            | sample1 (even) | sample2 (E.coli dom.) | sample3 (Staph/Bacillus)
#   -------------------|----------------|-----------------------|-------------------------
#   E. coli K-12       |     20×        |         50×           |          5×
#   S. aureus MRSA252  |     20×        |          5×           |         40×
#   P. aeruginosa PAO1 |     20×        |         15×           |          5×
#   B. subtilis 168    |     20×        |          5×           |         50×
#   M. tuberculosis    |     20×        |         20×           |         10×

set -euo pipefail

REFDIR="ref"
OUTDIR="reads"
TMPDIR_BASE="reads/tmp"
THREADS=4
SEED=42

mkdir -p "${OUTDIR}" "${TMPDIR_BASE}"

# ── helper: simulate one genome at given coverage, append to running R1/R2 ──
simulate_genome() {
    local genome="$1"
    local coverage="$2"
    local prefix="$3"   # tmp output prefix
    local r1_out="$4"
    local r2_out="$5"

    art_illumina \
        --seqSys  HS25 \
        --in      "${REFDIR}/${genome}" \
        --paired  \
        --len     150 \
        --fcov    "${coverage}" \
        --mflen   400 \
        --sdev    80 \
        --rndSeed "${SEED}" \
        --noALN \
        --quiet \
        --out     "${prefix}"

    cat "${prefix}1.fq" >> "${r1_out}"
    cat "${prefix}2.fq" >> "${r2_out}"
    rm -f "${prefix}1.fq" "${prefix}2.fq"
}

# ── generate one sample ──────────────────────────────────────────────────────
make_sample() {
    local name="$1"   # e.g. sample1
    shift
    # remaining args: genome cov [genome cov ...]
    local r1_tmp="${TMPDIR_BASE}/${name}_R1.fastq"
    local r2_tmp="${TMPDIR_BASE}/${name}_R2.fastq"
    > "${r1_tmp}"
    > "${r2_tmp}"

    echo "[${name}] Simulating reads..."
    while [[ $# -ge 2 ]]; do
        local genome="$1"
        local cov="$2"
        shift 2
        local tag="${genome%.fa}"
        echo "  ${tag}  →  ${cov}×"
        simulate_genome "${genome}" "${cov}" \
            "${TMPDIR_BASE}/${name}_${tag}" \
            "${r1_tmp}" "${r2_tmp}"
    done

    echo "  Shuffling and compressing..."
    # Shuffle read pairs together keeping mates in sync (Python, works on macOS & Linux)
    python3 - "${r1_tmp}" "${r2_tmp}" \
              "${TMPDIR_BASE}/${name}_R1_shuf.fastq" \
              "${TMPDIR_BASE}/${name}_R2_shuf.fastq" \
              "${SEED}" <<'PYEOF'
import sys, random

def read_fastq(path):
    """Yield (name, seq, plus, qual) tuples."""
    with open(path) as fh:
        while True:
            lines = [fh.readline().rstrip() for _ in range(4)]
            if not lines[0]:
                break
            yield tuple(lines)

r1_in, r2_in, r1_out, r2_out, seed = sys.argv[1:]
pairs = list(zip(read_fastq(r1_in), read_fastq(r2_in)))
random.seed(int(seed))
random.shuffle(pairs)
with open(r1_out, 'w') as f1, open(r2_out, 'w') as f2:
    for r1, r2 in pairs:
        f1.write('\n'.join(r1) + '\n')
        f2.write('\n'.join(r2) + '\n')
PYEOF

    pigz -p "${THREADS}" -c "${TMPDIR_BASE}/${name}_R1_shuf.fastq" > "${OUTDIR}/${name}_R1.fastq.gz"
    pigz -p "${THREADS}" -c "${TMPDIR_BASE}/${name}_R2_shuf.fastq" > "${OUTDIR}/${name}_R2.fastq.gz"

    rm -f "${r1_tmp}" "${r2_tmp}" \
          "${TMPDIR_BASE}/${name}_R1_shuf.fastq" \
          "${TMPDIR_BASE}/${name}_R2_shuf.fastq"

    local n_reads
    n_reads=$(pigz -p "${THREADS}" -dc "${OUTDIR}/${name}_R1.fastq.gz" | wc -l)
    n_reads=$(( n_reads / 4 ))
    echo "  → ${OUTDIR}/${name}_R{1,2}.fastq.gz  (${n_reads} read pairs)"
}

# ── Sample 1 : even community, all species at 20× ───────────────────────────
make_sample sample1 \
    e_coli.fa        20 \
    s_aureus.fa      20 \
    p_aeruginosa.fa  20 \
    b_subtilis.fa    20 \
    m_tuberculosis.fa 20

# ── Sample 2 : E. coli dominated ────────────────────────────────────────────
make_sample sample2 \
    e_coli.fa        50 \
    s_aureus.fa       5 \
    p_aeruginosa.fa  15 \
    b_subtilis.fa     5 \
    m_tuberculosis.fa 20

# ── Sample 3 : Staph / Bacillus heavy, low E. coli ──────────────────────────
make_sample sample3 \
    e_coli.fa         5 \
    s_aureus.fa      40 \
    p_aeruginosa.fa   5 \
    b_subtilis.fa    50 \
    m_tuberculosis.fa 10

rmdir "${TMPDIR_BASE}" 2>/dev/null || true

echo ""
echo "Done. Output files:"
ls -lh "${OUTDIR}"/*.fastq.gz
