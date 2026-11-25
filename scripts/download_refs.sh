#!/usr/bin/env bash
# Download 5 reference bacterial genomes for binning tutorial
# Tool: NCBI datasets CLI (v18+)
# Usage: run from project root (binning/)

set -euo pipefail

REFDIR="ref"
mkdir -p "${REFDIR}"
cd "${REFDIR}"

# Download all 5 genomes in one call
datasets download genome accession \
    GCF_000005845.2 \
    GCF_000011505.1 \
    GCF_000006765.1 \
    GCF_000009045.1 \
    GCF_000195955.2 \
    --include genome \
    --filename genomes.zip

unzip -o genomes.zip -d genomes_raw

# Rename to human-readable names
cp genomes_raw/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna   e_coli.fa
cp genomes_raw/ncbi_dataset/data/GCF_000011505.1/GCF_000011505.1_ASM1150v1_genomic.fna  s_aureus.fa
cp genomes_raw/ncbi_dataset/data/GCF_000006765.1/GCF_000006765.1_ASM676v1_genomic.fna   p_aeruginosa.fa
cp genomes_raw/ncbi_dataset/data/GCF_000009045.1/GCF_000009045.1_ASM904v1_genomic.fna   b_subtilis.fa
cp genomes_raw/ncbi_dataset/data/GCF_000195955.2/GCF_000195955.2_ASM19595v2_genomic.fna m_tuberculosis.fa

# Write tracking manifest
cat > genomes.tsv <<'TSV'
fasta_name	species	strain	ncbi_accession
e_coli.fa	Escherichia coli	K-12 MG1655	GCF_000005845.2
s_aureus.fa	Staphylococcus aureus	MRSA252	GCF_000011505.1
p_aeruginosa.fa	Pseudomonas aeruginosa	PAO1	GCF_000006765.1
b_subtilis.fa	Bacillus subtilis	168	GCF_000009045.1
m_tuberculosis.fa	Mycobacterium tuberculosis	H37Rv	GCF_000195955.2
TSV

echo "Reference genomes written to ${REFDIR}/"
echo "Manifest: ${REFDIR}/genomes.tsv"
