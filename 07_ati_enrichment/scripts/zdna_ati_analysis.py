#!/usr/bin/env python3
"""
Z-DNA enrichment at A-to-I editing sites — O. bimaculoides v2
Data: GSE284564_Oct_editing.csv (NC_068981.1 coordinates, no liftover needed)
Z-DNA: zhunter_z-dna.bedgraph.gz

Usage: python3 zdna_ati_analysis.py
"""

import gzip
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).parent
ZDNA_FILE  = Path("/Users/dassagaripova/Project_Octopus/02_zdna_analysis/results/zhunter_z-dna.bedgraph.gz")
EDIT_FILE  = SCRIPT_DIR / "GSE284564_Oct_editing.csv"

WINDOWS   = [0, 50, 100, 200, 500]
N_PERM    = 100

MAIN_CHROMS = {f"NC_0689{80+i:02d}.1" for i in range(1, 31)} | \
              {"NC_069000.1","NC_069001.1","NC_069002.1","NC_069003.1",
               "NC_069004.1","NC_069005.1","NC_069006.1","NC_069007.1",
               "NC_069008.1","NC_069009.1","NC_069010.1"}


def load_zdna():
    """Load Z-DNA as dict of sorted numpy arrays (starts, ends) per chrom."""
    print("Loading Z-DNA predictions …")
    data = defaultdict(lambda: {"s": [], "e": []})
    chrom_sizes = {}
    with gzip.open(str(ZDNA_FILE), "rt") as f:
        for line in f:
            if line.startswith(("browser", "track", "#")):
                continue
            p = line.split()
            if len(p) < 3:
                continue
            c, s, e = p[0], int(p[1]), int(p[2])
            if s >= e:
                e = s + 1
            data[c]["s"].append(s)
            data[c]["e"].append(e)
            if chrom_sizes.get(c, 0) < e:
                chrom_sizes[c] = e

    arrays = {}
    for c, d in data.items():
        s = np.array(d["s"], dtype=np.int32)
        e = np.array(d["e"], dtype=np.int32)
        idx = np.argsort(s)
        arrays[c] = (s[idx], e[idx])

    n = sum(len(v[0]) for v in arrays.values())
    print(f"  {n:,} Z-DNA regions on {len(arrays)} chromosomes")
    return arrays, chrom_sizes


def count_overlaps_fast(chroms, positions, zdna, window=0):
    """Vectorised overlap count using searchsorted."""
    total = 0
    for c in np.unique(chroms):
        if c not in zdna:
            continue
        pos = positions[chroms == c]
        starts, ends = zdna[c]
        lo = pos - window
        hi = pos + window + 1
        # for each query [lo,hi): find zdna intervals that overlap
        # overlap iff zdna_start < hi AND zdna_end > lo
        idx = np.searchsorted(starts, hi, side="left")   # starts < hi
        # among those, check ends > lo
        for i, (l, h, ix) in enumerate(zip(lo, hi, idx)):
            if ix > 0 and np.any(ends[:ix] > l):
                total += 1
    return total


def count_overlaps(sites_df, zdna, window=0):
    chroms = sites_df["chr"].values
    positions = sites_df["pos"].values.astype(np.int32)
    return count_overlaps_fast(chroms, positions, zdna, window)


def shuffle_sites(sites_df, chrom_sizes):
    dist = Counter(sites_df["chr"])
    rows = []
    for chrom, n in dist.items():
        size = chrom_sizes.get(chrom, 0)
        if size == 0:
            continue
        pos = np.random.randint(0, size, size=n).astype(np.int32)
        rows.append(pd.DataFrame({"chr": chrom, "pos": pos}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["chr","pos"])


def load_editing_sites():
    print("Loading A-to-I editing sites …")
    df = pd.read_csv(str(EDIT_FILE))
    ati = df[(df["control_edit"] == "A->G") & (df["chr"].isin(MAIN_CHROMS))].copy()
    ati["pos"] = ati["start"].astype(int)

    cold = ati[ati["delta_editing_percent"] < 0][["chr","pos"]].reset_index(drop=True)
    warm = ati[ati["delta_editing_percent"] > 0][["chr","pos"]].reset_index(drop=True)
    all_ = ati[["chr","pos"]].reset_index(drop=True)

    print(f"  A->G (main chroms):      {len(all_):,}")
    print(f"  Cold-regulated (↑ cold): {len(cold):,}")
    print(f"  Warm-regulated (↑ warm): {len(warm):,}")
    return all_, cold, warm


def run_analysis(label, sites_df, zdna, chrom_sizes):
    if sites_df is None or len(sites_df) == 0:
        print(f"  [{label}] No sites — skipping")
        return None, None

    n = len(sites_df)
    print(f"\n[{label}]  {n:,} sites")

    obs = {w: count_overlaps(sites_df, zdna, w) for w in WINDOWS}
    for w in WINDOWS:
        print(f"  ±{w:4d}bp: {obs[w]:,}/{n:,} ({100*obs[w]/n:.2f}%)")

    print(f"  Permutations ({N_PERM}x) …")
    perm = {w: [] for w in WINDOWS}
    for i in range(N_PERM):
        if (i + 1) % 100 == 0:
            print(f"    {i+1}/{N_PERM}")
        sh = shuffle_sites(sites_df, chrom_sizes)
        for w in WINDOWS:
            perm[w].append(count_overlaps(sh, zdna, w))

    print(f"\n  {'Window':>8}  {'Obs':>7}  {'Mean_perm':>10}  {'FE':>6}  {'p-value':>9}")
    print("  " + "-"*50)
    results = []
    for w in WINDOWS:
        arr  = np.array(perm[w])
        mean = arr.mean()
        fe   = obs[w] / mean if mean > 0 else float("nan")
        pval = (arr >= obs[w]).mean()
        ps   = f"<{1/N_PERM:.4f}" if pval == 0 else f"{pval:.4f}"
        print(f"  {w:>7}bp  {obs[w]:>7,}  {mean:>10.1f}  {fe:>6.3f}  {ps:>9}")
        results.append({"label": label, "window_bp": w,
                        "n_observed": obs[w], "mean_permuted": round(mean, 1),
                        "fold_enrichment": round(fe, 4), "p_value": pval,
                        "n_sites": n, "n_permutations": N_PERM})
    return results, perm


def make_plot(all_results, all_perms):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = {"All A->G": "steelblue", "Cold↑": "royalblue", "Warm↑": "tomato"}

    ax = axes[0]
    for res in all_results:
        if res is None:
            continue
        label = res[0]["label"]
        ws = [r["window_bp"] for r in res]
        fe = [r["fold_enrichment"] for r in res]
        ax.plot(ws, fe, "o-", label=label, color=colors.get(label, "gray"), lw=2, ms=7)
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xlabel("Window (±bp)")
    ax.set_ylabel("Fold Enrichment")
    ax.set_title("Z-DNA enrichment at A-to-I editing sites")
    ax.set_xticks(WINDOWS)
    ax.legend()

    w_show = 200
    pairs = [(all_results[1], all_perms[1], "Cold↑"),
             (all_results[2], all_perms[2], "Warm↑")]
    for ax, (res, perm, label) in zip(axes[1:], pairs):
        if res is None or perm is None:
            continue
        arr = np.array(perm[w_show])
        obs_val = next(r["n_observed"] for r in res if r["window_bp"] == w_show)
        color = colors.get(label, "steelblue")
        ax.hist(arr, bins=40, color=color, alpha=0.4, edgecolor=color)
        ax.axvline(obs_val, color="red", lw=2, label=f"Observed ({obs_val:,})")
        ax.axvline(arr.mean(), color="gray", ls="--", lw=1.5,
                   label=f"Mean ({arr.mean():.0f})")
        ax.set_xlabel("Sites overlapping Z-DNA")
        ax.set_ylabel("Permutation count")
        ax.set_title(f"{label} sites  |  ±{w_show}bp")
        ax.legend(fontsize=9)

    plt.suptitle("O. bimaculoides  |  Z-DNA × A-to-I editing enrichment",
                 fontsize=13, y=1.02)
    plt.tight_layout()
    out = SCRIPT_DIR / "zdna_ati_enrichment.png"
    plt.savefig(str(out), dpi=150, bbox_inches="tight")
    print(f"\nPlot saved: {out.name}")


def main():
    print("=" * 60)
    print("Z-DNA / A-to-I editing enrichment  |  O. bimaculoides v2")
    print("=" * 60)

    zdna, chrom_sizes = load_zdna()
    all_, cold, warm = load_editing_sites()

    all_results, all_perms = [], []
    for label, sites in [("All A->G", all_), ("Cold↑", cold), ("Warm↑", warm)]:
        res, perm = run_analysis(label, sites, zdna, chrom_sizes)
        all_results.append(res)
        all_perms.append(perm)

    rows = [r for res in all_results if res for r in res]
    out_tsv = SCRIPT_DIR / "zdna_ati_results.tsv"
    pd.DataFrame(rows).to_csv(str(out_tsv), sep="\t", index=False)
    print(f"\nResults saved: {out_tsv.name}")

    make_plot(all_results, all_perms)
    print("Done.")


if __name__ == "__main__":
    main()
