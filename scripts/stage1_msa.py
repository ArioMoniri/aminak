#!/usr/bin/env python3
"""Stage 1: Fetch FASTAs, run MAFFT, compute Jensen-Shannon conservation."""
import os, sys, time, subprocess, math
from datetime import datetime
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = os.path.expanduser("~/conserved_site_project")
MSA_DIR = os.path.join(PROJECT, "01_msa")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
MAFFT = os.environ.get("MAFFT", "/opt/homebrew/bin/mafft")
REF_ID = "P04818"

ACCESSIONS = [
    ("P04818", "Homo_sapiens"),
    ("P07607", "Mus_musculus"),
    ("P04394", "Escherichia_coli"),
    ("P04996", "Lactobacillus_casei"),
    ("P0CG53", "Saccharomyces_cerevisiae"),
    ("O44019", "Caenorhabditis_elegans"),
    ("Q9V3K2", "Drosophila_melanogaster"),
    ("P11849", "Bacteriophage_T4"),
    ("P21520", "Plasmodium_falciparum"),
]

FALLBACKS = {
    "P0CG53": "P49095",
}

AA = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE1: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def fetch_fasta(acc, label, retries=2, timeout=30):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200 and r.text.startswith(">"):
                # rewrite header to a clean label
                lines = r.text.strip().split("\n")
                seq = "".join(lines[1:])
                return f">{label}|{acc}\n{seq}\n"
            elif r.status_code == 404 and acc in FALLBACKS:
                log(f"{acc} 404, trying fallback {FALLBACKS[acc]}")
                return fetch_fasta(FALLBACKS[acc], label)
        except Exception as e:
            log(f"fetch {acc} attempt {attempt} failed: {e}")
            time.sleep(2)
    return None

def parse_fasta(path):
    seqs = []
    cur_h, cur_s = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_h is not None:
                    seqs.append((cur_h, "".join(cur_s)))
                cur_h, cur_s = line[1:], []
            else:
                cur_s.append(line)
        if cur_h is not None:
            seqs.append((cur_h, "".join(cur_s)))
    return seqs

# Capra & Singh 2007 JSD with BLOSUM62-derived background
# Background frequencies (Robinson & Robinson 1991, used in Capra 2007)
BG = {
    'A': 0.078, 'C': 0.019, 'D': 0.054, 'E': 0.063, 'F': 0.039,
    'G': 0.074, 'H': 0.022, 'I': 0.052, 'K': 0.057, 'L': 0.090,
    'M': 0.022, 'N': 0.045, 'P': 0.052, 'Q': 0.043, 'R': 0.051,
    'S': 0.071, 'T': 0.059, 'V': 0.064, 'W': 0.013, 'Y': 0.032,
}
BG_VEC = np.array([BG[a] for a in AA])

def column_freqs(col):
    """Return frequency vector over AA alphabet (excluding gaps)."""
    counts = np.zeros(len(AA))
    n = 0
    for c in col:
        if c in AA:
            counts[AA.index(c)] += 1
            n += 1
    if n == 0:
        return None, 0
    return counts / n, n

def js_divergence(p, q):
    m = 0.5 * (p + q)
    def kl(a, b):
        s = 0.0
        for ai, bi in zip(a, b):
            if ai > 0 and bi > 0:
                s += ai * math.log2(ai / bi)
        return s
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

def conservation_score(col, lam=0.5):
    """Capra & Singh 2007 JSD with pseudocount lambda*BG and gap penalty."""
    p, n = column_freqs(col)
    if p is None:
        return None, 0
    # Add pseudocounts
    p = (1 - lam) * p + lam * BG_VEC
    p = p / p.sum()
    jsd = js_divergence(p, BG_VEC)
    # Gap penalty per Capra 2007: multiply by (1 - gap_fraction)
    gap_frac = sum(1 for c in col if c == GAP) / len(col)
    return jsd * (1 - gap_frac), n

def windowed(scores, window=3):
    """Window of 3 means score[i] = mean of scores[i-1..i+1]."""
    out = [None] * len(scores)
    half = window // 2
    for i in range(len(scores)):
        vals = [scores[j] for j in range(max(0, i-half), min(len(scores), i+half+1)) if scores[j] is not None]
        if vals:
            out[i] = float(np.mean(vals))
    return out

def main():
    log("Stage 1 starting")
    os.makedirs(MSA_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    input_fa = os.path.join(MSA_DIR, "input.fa")
    aligned_fa = os.path.join(MSA_DIR, "aligned.fa")
    final_fa = os.path.join(MSA_DIR, "alignment.fasta")
    csv_path = os.path.join(MSA_DIR, "conservation_scores.csv")
    plot_path = os.path.join(MSA_DIR, "conservation_plot.png")
    mafft_log = os.path.join(LOG_DIR, "01_mafft.log")

    # Fetch FASTAs
    log("Fetching sequences from UniProt")
    chunks = []
    valid = []
    for acc, label in ACCESSIONS:
        c = fetch_fasta(acc, label)
        if c:
            chunks.append(c)
            valid.append(label)
            log(f"got {acc} ({label})")
        else:
            log(f"FAILED {acc} ({label})")
    log(f"Fetched {len(valid)} sequences: {valid}")
    if len(valid) < 6:
        log("FATAL: fewer than 6 sequences fetched")
        sys.exit(1)

    # For Plasmodium (bifunctional DHFR-TS), we will trim AFTER alignment
    # by mapping to the human TS reference: keep only columns where the
    # ref has residues — but we need to consider that PfDHFR-TS is large
    # (~600 aa) and the TS portion is the C-terminal half. Mafft will
    # align it; only the TS-aligned part contributes to relevant columns.
    with open(input_fa, "w") as f:
        f.writelines(chunks)
    log(f"Wrote input.fa with {len(chunks)} sequences")

    # Run MAFFT
    log("Running MAFFT --auto")
    with open(aligned_fa, "w") as outf, open(mafft_log, "w") as errf:
        proc = subprocess.run([MAFFT, "--auto", input_fa],
                              stdout=outf, stderr=errf, check=False)
    if proc.returncode != 0:
        log(f"MAFFT failed with code {proc.returncode}, see {mafft_log}")
        sys.exit(1)
    log("MAFFT done")

    # Copy to alignment.fasta as final name
    subprocess.run(["cp", aligned_fa, final_fa], check=True)

    # Parse alignment
    seqs = parse_fasta(aligned_fa)
    log(f"Aligned {len(seqs)} sequences, length {len(seqs[0][1])}")

    # Find reference index
    ref_idx = None
    for i, (h, s) in enumerate(seqs):
        if "P04818" in h or h.startswith("Homo_sapiens"):
            ref_idx = i
            break
    if ref_idx is None:
        log("FATAL: reference P04818 not found in alignment")
        sys.exit(1)
    ref_seq = seqs[ref_idx][1]
    log(f"Reference is row {ref_idx}, ungapped length {sum(1 for c in ref_seq if c != GAP)}")

    # Compute scores per column
    aln_len = len(seqs[0][1])
    raw_scores = []
    for j in range(aln_len):
        col = [s[j] for _, s in seqs]
        score, n = conservation_score(col)
        raw_scores.append(score)

    # Window smooth
    smoothed = windowed(raw_scores, window=3)

    # Map to ref ungapped positions
    rows = []
    ref_pos = 0
    for j in range(aln_len):
        if ref_seq[j] != GAP:
            ref_pos += 1
            rows.append({
                "ref_position": ref_pos,
                "residue": ref_seq[j],
                "js_score": smoothed[j] if smoothed[j] is not None else float("nan"),
            })
    df = pd.DataFrame(rows)
    # Percentile (per non-NaN)
    valid_scores = df["js_score"].dropna()
    df["percentile"] = df["js_score"].rank(pct=True, na_option="keep") * 100
    df.to_csv(csv_path, index=False)
    log(f"Wrote {csv_path} ({len(df)} rows)")

    # Plot
    fig, ax = plt.subplots(figsize=(14, 4))
    xs = df["ref_position"].values
    ys = df["js_score"].fillna(0).values
    p90 = np.nanpercentile(df["js_score"], 90)
    colors = ["#c0392b" if (not np.isnan(v) and v >= p90) else "#3498db" for v in df["js_score"].values]
    ax.bar(xs, ys, color=colors, width=1.0)
    ax.set_xlabel("Reference position (P04818)")
    ax.set_ylabel("JSD conservation (windowed)")
    ax.set_title(f"TYMS conservation across {len(seqs)} orthologs (top 10% red)")
    ax.axhline(p90, color="black", linestyle="--", alpha=0.5, label=f"p90={p90:.3f}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(plot_path, dpi=120)
    plt.close()
    log(f"Wrote {plot_path}")
    log("Stage 1 DONE")

if __name__ == "__main__":
    main()
