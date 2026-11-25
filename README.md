# Metagenomic Binning Tutorial

A reproducible tutorial for metagenomic binning using a simulated five-species mock community.
Reads are simulated from reference genomes, assembled de novo, and contigs are binned by
combining **tetranucleotide frequency (TNF)** composition with **differential coverage** across
three samples with different species ratios.

## Mock community

| Species | Strain | GC% | Genome |
|---|---|---|---|
| *Escherichia coli* | K-12 MG1655 | 50% | GCF_000005845.2 |
| *Staphylococcus aureus* | MRSA252 | 33% | GCF_000011505.1 |
| *Pseudomonas aeruginosa* | PAO1 | 66% | GCF_000006765.1 |
| *Bacillus subtilis* | 168 | 43% | GCF_000009045.1 |
| *Mycobacterium tuberculosis* | H37Rv | 65% | GCF_000195955.2 |

## Sample design

Three paired-end Illumina samples (2×150 bp, HiSeq 2500) with different coverage ratios
to provide differential abundance signal for binning.

| Species | sample1 (even) | sample2 (E. coli dom.) | sample3 (Staph/Bacillus) |
|---|---|---|---|
| *E. coli* | 20× | 50× | 5× |
| *S. aureus* | 20× | 5× | 40× |
| *P. aeruginosa* | 20× | 15× | 5× |
| *B. subtilis* | 20× | 5× | 50× |
| *M. tuberculosis* | 20× | 20× | 10× |

## Pipeline

```
scripts/download_refs.sh      # download reference genomes → ref/
scripts/generate_shotgun.sh   # simulate reads → reads/
scripts/denovo.sh             # assemble → assemblies/
scripts/backmap.sh            # map reads to assemblies → backmapping/
scripts/binning.sh            # coverage + TNF clustering → counts/ binning/
```

Run each script from the project root with the `binning` conda environment active:

```bash
micromamba activate binning
bash scripts/download_refs.sh
bash scripts/generate_shotgun.sh
bash scripts/denovo.sh
bash scripts/backmap.sh
bash scripts/binning.sh
```

The binning step can also be run manually with custom parameters:

```bash
python scripts/tnf.py assemblies/sample1.contigs.fasta \
    --coverage counts/sample1/sample3.tsv \
    --clusters auto \
    --min-length 1500 \
    --output binning/sample1_sample3
```

## Repository layout

```
ref/genomes.tsv       species manifest (name, strain, accession)
scripts/              all pipeline scripts and the tnf.py tool
notes/results.md      benchmark results and interpretation
```

## Notes

**MEGAHIT on Apple Silicon:** MEGAHIT v1.2.9 segfaults on ARM64 macOS with its
default memory mode. `denovo.sh` includes `--no-hw-accel --mem-flag 0` to work around this.

**Coverage metric:** `bamcountrefs --trimmed-mean` (5th–95th percentile depth per position).
Raw counts are proportional to contig length and inflate long contigs from rare species;
trimmed mean depth is length-normalised and robust against coverage spikes at contig ends.

**TNF clustering:** axes (TNF PC1, log₁₀ coverage) are standardised to zero-mean unit-variance
before clustering so composition and depth contribute equally to Euclidean distance.
Without standardisation, depth (range ~3 log units) dominates over TNF PC1 (range ~0.2 units)
and clusters reduce to horizontal depth bands. See `notes/results.md` for benchmark results.
