#!/usr/bin/env python3
"""
tnf.py — Tetranucleotide frequency analysis for metagenomic contigs

Modes
-----
No coverage table  →  PC1 vs PC2 scatter plot
Coverage table     →  PC1 vs log10(coverage) scatter plot

Add --clusters N (fixed k, K-Means) or --clusters auto (GMM + BIC) to
colour points by inferred cluster in either mode.

Usage
-----
    # PC1 vs PC2
    python scripts/tnf.py assembly.fa

    # PC1 vs coverage
    python scripts/tnf.py assembly.fa --coverage coverage.tsv

    # With de-novo clustering
    python scripts/tnf.py assembly.fa --coverage coverage.tsv --clusters auto
    python scripts/tnf.py assembly.fa --clusters 5 --output plot.png
"""

import argparse
import sys
from itertools import product

import numpy as np


# ── canonical 136-mer index ──────────────────────────────────────────────────

def _rc(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def canonical_tetramers() -> list[str]:
    """Return the 136 lexicographically-canonical tetranucleotides."""
    seen: set[str] = set()
    out: list[str] = []
    for bases in product("ACGT", repeat=4):
        kmer = "".join(bases)
        canon = min(kmer, _rc(kmer))
        if canon not in seen:
            seen.add(canon)
            out.append(canon)
    return sorted(out)


CANONICAL = canonical_tetramers()          # 136 tetramers, computed once
KMER_INDEX = {k: i for i, k in enumerate(CANONICAL)}


# ── TNF computation ──────────────────────────────────────────────────────────

def compute_tnf(seq: str) -> np.ndarray:
    """
    Compute the normalised canonical TNF vector for a nucleotide sequence.
    Returns a float32 array of length 136 (sums to 1.0).
    Non-ACGT characters are skipped.
    """
    seq = seq.upper()
    counts = np.zeros(136, dtype=np.float32)
    n = len(seq)
    for i in range(n - 3):
        kmer = seq[i:i + 4]
        if "N" in kmer:
            continue
        idx = KMER_INDEX.get(min(kmer, _rc(kmer)))
        if idx is not None:
            counts[idx] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return counts


# ── FASTA parser ─────────────────────────────────────────────────────────────

def load_fasta(path: str) -> dict[str, str]:
    """Return {header: sequence} for every record in a FASTA file."""
    seqs: dict[str, str] = {}
    name: str | None = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


# ── coverage table parser ────────────────────────────────────────────────────

def load_coverage(path: str) -> dict[str, float]:
    """
    Parse a two-column TSV (contig_name\\tcoverage).
    Lines starting with '#' and a header line where coverage is not numeric
    are skipped automatically.
    """
    cov: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            try:
                cov[parts[0]] = float(parts[1])
            except ValueError:
                continue          # skip header or malformed row
    return cov


# ── clustering ───────────────────────────────────────────────────────────────

def find_clusters(coords: np.ndarray, n: int | str) -> np.ndarray:
    """
    Cluster 2D plot coordinates.

    n = int    → K-Means with exactly n clusters
    n = 'auto' → Gaussian Mixture Model; k chosen by lowest BIC over 1..min(10, N//2)

    Axes are standardised (zero-mean, unit-variance) before clustering so that
    TNF and coverage contribute equally to Euclidean distance regardless of
    their original numeric ranges.

    Returns an integer label array of shape (len(coords),).
    """
    try:
        from sklearn.cluster import KMeans
        from sklearn.mixture import GaussianMixture
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        sys.exit(f"[tnf.py] Missing dependency: {exc}. "
                 "Install with: pip install scikit-learn")

    coords = StandardScaler().fit_transform(coords)

    if isinstance(n, int):
        km = KMeans(n_clusters=n, random_state=42, n_init=10)
        labels = km.fit_predict(coords)
        print(f"[tnf.py] K-Means: k={n}", file=sys.stderr)
        return labels

    # auto: BIC-guided GMM
    n_max = max(1, min(10, len(coords) // 2))
    best_bic, best_labels, best_k = np.inf, np.zeros(len(coords), dtype=int), 1
    for k in range(1, n_max + 1):
        gm = GaussianMixture(n_components=k, random_state=42, n_init=5,
                             covariance_type="full")
        gm.fit(coords)
        bic = gm.bic(coords)
        if bic < best_bic:
            best_bic = bic
            best_labels = gm.predict(coords)
            best_k = k
    print(f"[tnf.py] GMM auto: selected k={best_k}  (BIC={best_bic:.1f})",
          file=sys.stderr)
    for lbl in range(best_k):
        print(f"  cluster {lbl}: {(best_labels == lbl).sum()} contigs", file=sys.stderr)
    return best_labels


# ── plotting ─────────────────────────────────────────────────────────────────

def _cluster_colors(labels: np.ndarray):
    """Return a list of per-point colors from a qualitative palette."""
    import matplotlib.pyplot as plt
    n_clusters = labels.max() + 1
    cmap = plt.get_cmap("tab10" if n_clusters <= 10 else "tab20")
    return [cmap(int(l) % cmap.N) for l in labels]


def plot_pc1_vs_pc2(
    names: list[str],
    tnf_matrix: np.ndarray,
    output: str | None,
    labels: np.ndarray | None = None,
) -> None:
    """Scatter plot of TNF PC1 (x) vs PC2 (y), optionally coloured by cluster."""
    try:
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
    except ImportError as exc:
        sys.exit(f"[tnf.py] Missing dependency: {exc}. "
                 "Install with: pip install matplotlib scikit-learn")

    if tnf_matrix.shape[0] < 2:
        sys.exit("[tnf.py] Need at least 2 contigs to run PCA.")

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(tnf_matrix)
    var1, var2 = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(7, 6))

    if labels is not None:
        colors = _cluster_colors(labels)
        ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=18, alpha=0.70, linewidths=0)
        # legend: one entry per cluster
        import matplotlib.patches as mpatches
        cmap = plt.get_cmap("tab10" if labels.max() < 10 else "tab20")
        handles = [mpatches.Patch(color=cmap(k), label=f"Cluster {k}")
                   for k in range(labels.max() + 1)]
        ax.legend(handles=handles, fontsize=8, title="Clusters", title_fontsize=8,
                  loc="best", framealpha=0.7)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=18, alpha=0.65, linewidths=0)

    ax.set_xlabel(f"TNF PC1  ({var1:.1f}% variance)", fontsize=11)
    ax.set_ylabel(f"TNF PC2  ({var2:.1f}% variance)", fontsize=11)
    ax.set_title(f"TNF PCA  ({len(names)} contigs)", fontsize=12)
    ax.tick_params(labelsize=9)
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150)
        print(f"[tnf.py] Plot saved → {output}")
    else:
        plt.show()


def plot_pc1_vs_coverage(
    names: list[str],
    tnf_matrix: np.ndarray,
    coverages: np.ndarray,
    output: str | None,
    labels: np.ndarray | None = None,
) -> None:
    """Scatter plot of TNF PC1 (x) vs log10(coverage) (y), optionally coloured by cluster."""
    try:
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
    except ImportError as exc:
        sys.exit(f"[tnf.py] Missing dependency: {exc}. "
                 "Install with: pip install matplotlib scikit-learn")

    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(tnf_matrix)        # (n, 2)
    pc1 = coords[:, 0]
    var_explained = pca.explained_variance_ratio_[0] * 100
    log_cov = np.log10(coverages)

    fig, ax = plt.subplots(figsize=(8, 5))

    if labels is not None:
        colors = _cluster_colors(labels)
        ax.scatter(pc1, log_cov, c=colors, s=18, alpha=0.70, linewidths=0)
        import matplotlib.patches as mpatches
        cmap = plt.get_cmap("tab10" if labels.max() < 10 else "tab20")
        handles = [mpatches.Patch(color=cmap(k), label=f"Cluster {k}")
                   for k in range(labels.max() + 1)]
        ax.legend(handles=handles, fontsize=8, title="Clusters", title_fontsize=8,
                  loc="best", framealpha=0.7)
    else:
        sc = ax.scatter(pc1, log_cov, s=18, alpha=0.65, linewidths=0,
                        c=log_cov, cmap="viridis")
        cbar = fig.colorbar(sc, ax=ax, pad=0.02)
        cbar.set_label("log₁₀(coverage)", fontsize=9)

    ax.set_xlabel(f"TNF PC1  ({var_explained:.1f}% variance)", fontsize=11)
    ax.set_ylabel("log₁₀(coverage)", fontsize=11)
    ax.set_title(f"TNF composition vs coverage  ({len(names)} contigs)", fontsize=12)
    ax.tick_params(labelsize=9)
    fig.tight_layout()

    if output:
        fig.savefig(output, dpi=150)
        print(f"[tnf.py] Plot saved → {output}")
    else:
        plt.show()


# ── output helpers ───────────────────────────────────────────────────────────

def _write_clusters_tsv(names: list[str], labels: np.ndarray, path: str) -> None:
    with open(path, "w") as fh:
        fh.write("contig\tcluster\n")
        for name, lbl in zip(names, labels):
            fh.write(f"{name}\t{lbl}\n")
    print(f"[tnf.py] Clusters saved → {path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Tetranucleotide frequency (TNF) analysis for metagenomic contigs."
    )
    ap.add_argument("fasta", help="FASTA file of contig(s)")
    ap.add_argument(
        "--coverage", "-c",
        metavar="TSV",
        help="Two-column TSV: contig_name<tab>coverage. "
             "When supplied: PC1 vs log10(coverage) plot. "
             "When omitted: PC1 vs PC2 plot.",
    )
    ap.add_argument(
        "--output", "-o",
        metavar="BASENAME",
        help="Save outputs under BASENAME: writes BASENAME_plot.png and, "
             "when --clusters is set, BASENAME_clusters.tsv (contig<tab>cluster). "
             "If omitted, the plot is displayed interactively.",
    )
    ap.add_argument(
        "--min-length", "-l",
        type=int, default=0,
        metavar="INT",
        help="Skip contigs shorter than INT bp (default: 0).",
    )
    ap.add_argument(
        "--clusters", "-k",
        metavar="N|auto",
        default=None,
        help="Colour points by cluster. Pass an integer for K-Means with exactly N "
             "clusters, or 'auto' to select k automatically via GMM + BIC.",
    )
    args = ap.parse_args()

    # parse --clusters value
    clusters: int | str | None = None
    if args.clusters is not None:
        if args.clusters == "auto":
            clusters = "auto"
        else:
            try:
                clusters = int(args.clusters)
                if clusters < 1:
                    raise ValueError
            except ValueError:
                sys.exit("[tnf.py] --clusters must be a positive integer or 'auto'.")

    # Load sequences
    seqs = load_fasta(args.fasta)
    if not seqs:
        sys.exit(f"[tnf.py] No sequences found in {args.fasta}")

    if args.min_length > 0:
        before = len(seqs)
        seqs = {n: s for n, s in seqs.items() if len(s) >= args.min_length}
        dropped = before - len(seqs)
        if dropped:
            print(f"[tnf.py] Dropped {dropped} contig(s) shorter than {args.min_length} bp",
                  file=sys.stderr)

    # derive output paths from basename
    plot_path = f"{args.output}_plot.png" if args.output else None
    tsv_path  = f"{args.output}_clusters.tsv" if args.output else None

    # ── mode 1: no coverage — PCA plot PC1 vs PC2 ───────────────────────────
    if args.coverage is None:
        names = list(seqs.keys())
        print(f"[tnf.py] Computing TNF for {len(names)} contigs...", file=sys.stderr)
        tnf_matrix = np.vstack([compute_tnf(seqs[n]) for n in names])
        labels = find_clusters(tnf_matrix, clusters) if clusters is not None else None
        plot_pc1_vs_pc2(names, tnf_matrix, plot_path, labels=labels)
        if labels is not None and tsv_path:
            _write_clusters_tsv(names, labels, tsv_path)
        return

    # ── mode 2: coverage table supplied — PCA + plot ─────────────────────────
    cov_table = load_coverage(args.coverage)

    # Intersect contigs present in both FASTA and coverage table
    shared = [n for n in seqs if n in cov_table]
    missing_cov = len(seqs) - len(shared)
    if missing_cov:
        print(f"[tnf.py] {missing_cov} contig(s) have no coverage entry — skipped",
              file=sys.stderr)
    if not shared:
        sys.exit("[tnf.py] No contigs overlap between FASTA and coverage table.")

    print(f"[tnf.py] Computing TNF for {len(shared)} contigs...", file=sys.stderr)
    tnf_matrix = np.vstack([compute_tnf(seqs[n]) for n in shared])
    coverages  = np.array([cov_table[n] for n in shared], dtype=np.float64)

    if np.any(coverages <= 0):
        sys.exit("[tnf.py] Coverage values must be > 0 for log transform.")

    # cluster on (PC1, log10_cov) — the actual 2D space shown in the plot
    if clusters is not None:
        from sklearn.decomposition import PCA
        pc1 = PCA(n_components=1, random_state=42).fit_transform(tnf_matrix)
        plot_coords = np.column_stack([pc1, np.log10(coverages)])
        labels = find_clusters(plot_coords, clusters)
    else:
        labels = None

    plot_pc1_vs_coverage(shared, tnf_matrix, coverages, plot_path, labels=labels)
    if labels is not None and tsv_path:
        _write_clusters_tsv(shared, labels, tsv_path)


if __name__ == "__main__":
    main()
