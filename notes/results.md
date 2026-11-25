# Binning results — TNF + coverage clustering

## Mock community

| File | Species | Strain | Accession | GC% |
|---|---|---|---|---|
| e_coli.fa | *Escherichia coli* | K-12 MG1655 | GCF_000005845.2 | 50% |
| s_aureus.fa | *Staphylococcus aureus* | MRSA252 | GCF_000011505.1 | 33% |
| p_aeruginosa.fa | *Pseudomonas aeruginosa* | PAO1 | GCF_000006765.1 | 66% |
| b_subtilis.fa | *Bacillus subtilis* | 168 | GCF_000009045.1 | 43% |
| m_tuberculosis.fa | *Mycobacterium tuberculosis* | H37Rv | GCF_000195955.2 | 65% |

## Sample design (fold coverage per species)

| Species | sample1 (even) | sample2 (E. coli dom.) | sample3 (Staph/Bacillus) |
|---|---|---|---|
| *E. coli* | 20× | 50× | 5× |
| *S. aureus* | 20× | 5× | 40× |
| *P. aeruginosa* | 20× | 15× | 5× |
| *B. subtilis* | 20× | 5× | 50× |
| *M. tuberculosis* | 20× | 20× | 10× |

## Methods

- **Simulator**: ART (HiSeq 2500, 2×150 bp, insert ~400 bp)
- **Assembler**: MEGAHIT v1.2.9, min contig 1,000 bp
- **Coverage metric**: trimmed mean depth (5th–95th percentile positions, `bamcountrefs --trimmed-mean --mapq 20 --proper-pairs`)
- **TNF**: 136 canonical tetranucleotides, normalised to frequencies
- **Dimensionality reduction**: PCA, cluster on (PC1, log₁₀ trimmed-mean)
- **Axes standardised** with `StandardScaler` before clustering (zero-mean, unit-variance)
- **Clustering**: GMM with BIC-guided k selection (k = 1…10), `auto` mode

## Metrics

### Purity
For each cluster, find the dominant (most common) species among its contigs.
Purity is the fraction of all contigs that are in their cluster's dominant species:

```
purity = Σ_clusters  max_species(count in cluster)
         ─────────────────────────────────────────
                   total contigs
```

A cluster with 90 E. coli + 10 S. aureus contigs contributes 90 to the numerator.
High purity means clusters are not mixed — each bin is mostly one organism.
It does **not** penalise a species being split across multiple clusters.

### Recall
For each species, find the single cluster where it is most represented (its dominant cluster).
Recall for that species is the fraction of its contigs that land in that cluster:

```
recall(species) = contigs of species in its dominant cluster
                  ──────────────────────────────────────────
                       total contigs of that species
```

Mean recall averages this across all five species.
High recall means contigs from each species are **concentrated** in one bin.
It does **not** penalise a bin containing contigs from other species.

Together: **purity measures bin contamination; recall measures bin completeness.**
A perfect binner scores 100% on both. A run can score high purity but low recall
(many clean but incomplete bins) or high recall but low purity (complete but contaminated bins).

## Results

Sorted by purity + recall. Ground truth assigned by `bwa mem -x intractg` alignment
of assembly contigs to the five reference genomes.

| Assembly | Coverage from | k (auto) | Purity | Mean recall |
|---|---|---|---|---|
| sample3 | sample2 | 7 | **97.4%** | 77.8% |
| sample2 | sample3 | 8 | **97.4%** | 75.4% |
| sample3 | sample3 | 7 | 97.0% | **81.2%** |
| sample1 | sample3 | 8 | 90.4% | **78.3%** |
| sample2 | sample2 | 9 | 90.4% | 63.1% |
| sample3 | sample1 | 7 | 90.4% | 48.3% |
| sample2 | sample1 | 7 | 88.0% | 68.3% |
| sample1 | sample2 | 8 | 89.3% | 74.0% |
| sample1 | sample1 | 8 | 71.6% | 61.8% |
| sample1 | sample1 (k=5 fixed) | 5 | 32.5% | 27.5% |

## Key observations

**1. Standardisation matters.**
Before applying `StandardScaler`, clustering on raw (PC1, log-depth) was dominated
by depth (range ~3 log units vs ~0.2 PC1 units). `sample1/sample1` went from 28.8% → 71.6%
purity after standardising both axes to unit variance.

**2. Even coverage (sample1) is always the weakest coverage source.**
When all five species are sequenced at 20×, depth is flat across contigs.
The y-axis carries no discriminating signal and clustering relies entirely on PC1 alone.

**3. Cross-sample differential coverage rescues weak assemblies.**
`sample1/sample3` (90.4% purity) substantially outperforms `sample1/sample1` (71.6%)
using the same assembly. The assembly does not need to come from the same sample as the coverage.

**4. Differential coverage is the main driver.**
The best three runs all use a non-even coverage sample (sample2 or sample3).
This is the core principle of differential-coverage binning (MetaBAT2, CONCOCT):
species present at different abundances across samples separate along the y-axis.

**5. Recall is the remaining bottleneck.**
Purity reaches ~97% in the best runs but recall plateaus at ~75–81%.
The two high-GC species (*M. tuberculosis* ~65%, *P. aeruginosa* ~66%) overlap
in TNF PC1 and are the primary source of recall loss — their contigs land in
each other's dominant cluster.

**6. Auto k (GMM + BIC) outperforms fixed k = number of species.**
Forcing k = 5 (the true number of species) scored 32.5% purity. In a 2D
(PC1, log-depth) space, species do not necessarily each form one Gaussian,
and BIC correctly selected larger k values (7–9) that better fit the actual geometry.
