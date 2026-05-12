#!/usr/bin/env python3
"""Stage 1 v2: Real TYMS orthologs, sanity-checked MSA, JSD with gap-column exclusion from percentiles."""
import os, sys, time, subprocess, math, re, json
from datetime import datetime
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = os.path.expanduser("~/conserved_site_project")
MSA_DIR = os.path.join(PROJECT, "01b_msa_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_01_msa.log")
MAFFT = os.environ.get("MAFFT", "/opt/homebrew/bin/mafft")
REF_ID = "P04818"

# Real TYMS orthologs per spec (corrected after probing UniProt 2026-05)
# Some originally-spec'd accessions (P07807, Q23381, Q9V3K2, P04019) were wrong
# entries; corrected via REST search. Documented via fallback log.
ACCESSIONS = [
    ("P04818", "Homo_sapiens"),
    ("P07607", "Mus_musculus"),
    ("P45352", "Rattus_norvegicus"),
    ("P0A884", "Escherichia_coli"),       # ThyA, 264 aa, "KMALAPCHAFFQF"
    ("P00469", "Lactobacillus_casei"),    # 316 aa, "TMALPPCHTLYQF"
    ("P06785", "Saccharomyces_cerevisiae"),  # CDC21/TS, 304 aa (corrected from P07807=DHFR)
    ("Q9N588", "Caenorhabditis_elegans"), # TrEMBL, 312 aa (corrected; Q23381=MUTA)
    ("O76511", "Drosophila_melanogaster"),# 321 aa (corrected from Q9V3K2 deleted)
    ("Q05762", "Arabidopsis_thaliana"),   # bifunctional DHFR-TS 519 aa, will trim
    ("P00471", "Bacteriophage_T4"),       # T4 td, 286 aa (corrected from P04019)
    ("P13922", "Plasmodium_falciparum"),  # bifunctional DHFR-TS, will trim
]

# Bifunctional DHFR-TS organisms: TS domain is C-terminal half
BIFUNCTIONAL = {"Plasmodium_falciparum", "Arabidopsis_thaliana"}

# TS-family signature near catalytic Cys: relaxed to handle all natural TS
# Real catalytic context is "x-MAL[PA]PC[HR]" with mild variation.
TS_MOTIF = re.compile(r"[ILVMAQK]MAL[PA]PC[HRA]")
TS_MOTIF_LOOSE = re.compile(r"[A-Z][A-Z]AL[APV]PC[HRAQ]")

AA = "ACDEFGHIKLMNPQRSTVWY"
GAP = "-"

# Robinson-Robinson background (Capra & Singh 2007)
BG = {
    'A': 0.078, 'C': 0.019, 'D': 0.054, 'E': 0.063, 'F': 0.039,
    'G': 0.074, 'H': 0.022, 'I': 0.052, 'K': 0.057, 'L': 0.090,
    'M': 0.022, 'N': 0.045, 'P': 0.052, 'Q': 0.043, 'R': 0.051,
    'S': 0.071, 'T': 0.059, 'V': 0.064, 'W': 0.013, 'Y': 0.032,
}
BG_VEC = np.array([BG[a] for a in AA])


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE1: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def fetch_fasta_raw(acc, retries=3, timeout=30):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200 and r.text.startswith(">"):
                lines = r.text.strip().split("\n")
                seq = "".join(lines[1:])
                return seq
            else:
                log(f"fetch {acc} status {r.status_code}")
        except Exception as e:
            log(f"fetch {acc} attempt {attempt} failed: {e}")
            time.sleep(2)
    return None


def search_uniprot(query, organism, retries=2):
    """REST search fallback."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"{query} AND organism_name:\"{organism}\" AND reviewed:true",
        "format": "fasta",
        "size": 1,
    }
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=30)
            if r.status_code == 200 and r.text.startswith(">"):
                lines = r.text.strip().split("\n")
                return "".join(lines[1:])
        except Exception as e:
            log(f"search fallback failed: {e}")
            time.sleep(2)
    return None


def trim_bifunctional_ts(seq):
    """Bifunctional DHFR-TS: TS domain is C-terminal half.
    Find TS catalytic motif and trim ~190 aa upstream (size of human TS N-term)."""
    m = TS_MOTIF.search(seq) or TS_MOTIF_LOOSE.search(seq)
    if m:
        start = max(0, m.start() - 190)  # human TS Cys195 is at pos 195, so ~195aa upstream
        return seq[start:]
    if len(seq) > 400:
        return seq[280:]
    return seq


def validate_seq(label, seq):
    """Return True if length 100..800 AND contains TS motif."""
    if not seq:
        return False, "no_seq"
    L = len(seq)
    if L < 100:
        return False, f"too_short({L})"
    if L > 800:
        return False, f"too_long({L})"
    if not (TS_MOTIF.search(seq) or TS_MOTIF_LOOSE.search(seq)):
        return False, "no_TS_motif"
    return True, f"ok(len={L})"


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


def column_freqs(col):
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
    """Capra & Singh 2007 with pseudocount + gap penalty (1 - gap_frac).
    Returns (score, n_non_gap, gap_frac)."""
    p, n = column_freqs(col)
    gap_frac = sum(1 for c in col if c == GAP) / len(col)
    if p is None:
        return None, 0, gap_frac
    p = (1 - lam) * p + lam * BG_VEC
    p = p / p.sum()
    jsd = js_divergence(p, BG_VEC)
    return jsd * (1 - gap_frac), n, gap_frac


def windowed_weighted(scores, weights=(0.25, 0.5, 0.25)):
    """Per-spec: 0.5*s[i] + 0.25*(s[i-1] + s[i+1]). Skip None."""
    out = [None] * len(scores)
    for i in range(len(scores)):
        if scores[i] is None:
            continue
        total = weights[1] * scores[i]
        denom = weights[1]
        if i - 1 >= 0 and scores[i-1] is not None:
            total += weights[0] * scores[i-1]
            denom += weights[0]
        if i + 1 < len(scores) and scores[i+1] is not None:
            total += weights[2] * scores[i+1]
            denom += weights[2]
        out[i] = total / denom
    return out


def main():
    os.makedirs(MSA_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 1 v2 starting")

    input_fa = os.path.join(MSA_DIR, "input.fa")
    aligned_fa = os.path.join(MSA_DIR, "aligned.fa")
    csv_path = os.path.join(MSA_DIR, "conservation_scores.csv")
    plot_path = os.path.join(MSA_DIR, "conservation_plot.png")
    mafft_log = os.path.join(LOG_DIR, "v2_01_mafft.log")

    chunks = []
    valid = []
    failed = []
    for acc, label in ACCESSIONS:
        log(f"Fetching {acc} ({label})")
        seq = fetch_fasta_raw(acc)
        if not seq:
            log(f"  primary fetch failed for {acc}, trying search fallback")
            seq = search_uniprot("thymidylate synthase", label.replace("_", " "))
            if not seq:
                failed.append((acc, label, "fetch_failed"))
                continue
        if label in BIFUNCTIONAL:
            orig_len = len(seq)
            seq = trim_bifunctional_ts(seq)
            log(f"  {label} bifunctional trim: {orig_len} -> {len(seq)} aa")
        ok, reason = validate_seq(label, seq)
        if not ok:
            log(f"  REJECT {acc} ({label}): {reason}; trying search fallback")
            alt = search_uniprot("thymidylate synthase", label.replace("_", " "))
            if alt:
                if label in BIFUNCTIONAL:
                    alt = trim_bifunctional_ts(alt)
                ok, reason = validate_seq(label, alt)
                if ok:
                    seq = alt
                    log(f"  fallback OK: {reason}")
                else:
                    failed.append((acc, label, reason))
                    continue
            else:
                failed.append((acc, label, reason))
                continue
        chunks.append(f">{label}|{acc}\n{seq}\n")
        valid.append((acc, label))
        log(f"  ACCEPT {acc} ({label}): {reason}")

    log(f"Final ortholog set: {len(valid)} valid, {len(failed)} rejected")
    for acc, label, reason in failed:
        log(f"  rejected: {acc} {label} ({reason})")

    if len(valid) < 8:
        log(f"WARNING: only {len(valid)} orthologs (target was >=10), proceeding anyway")
    if len(valid) < 5:
        log("FATAL: fewer than 5 sequences validated")
        sys.exit(1)

    with open(input_fa, "w") as f:
        f.writelines(chunks)
    log(f"Wrote {input_fa} with {len(chunks)} sequences")

    # Run MAFFT
    log("Running MAFFT --auto")
    with open(aligned_fa, "w") as outf, open(mafft_log, "w") as errf:
        proc = subprocess.run([MAFFT, "--auto", input_fa],
                              stdout=outf, stderr=errf, check=False)
    if proc.returncode != 0:
        log(f"MAFFT failed rc={proc.returncode}")
        sys.exit(1)
    log("MAFFT done")

    seqs = parse_fasta(aligned_fa)
    log(f"Aligned {len(seqs)} sequences, length {len(seqs[0][1])}")

    ref_idx = None
    for i, (h, s) in enumerate(seqs):
        if "P04818" in h:
            ref_idx = i
            break
    if ref_idx is None:
        log("FATAL: reference P04818 not found")
        sys.exit(1)
    ref_seq = seqs[ref_idx][1]
    log(f"Reference row {ref_idx}, ungapped len {sum(1 for c in ref_seq if c != GAP)}")

    aln_len = len(seqs[0][1])
    raw_scores = []
    gap_fracs = []
    for j in range(aln_len):
        col = [s[j] for _, s in seqs]
        score, n, gf = conservation_score(col)
        raw_scores.append(score)
        gap_fracs.append(gf)

    smoothed = windowed_weighted(raw_scores)

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
                "gap_fraction": gap_fracs[j],
                "aln_col": j,
            })
    df = pd.DataFrame(rows)

    # Per spec: EXCLUDE columns with >50% gaps from percentile ranking
    eligible_mask = df["gap_fraction"] <= 0.5
    df["percentile"] = float("nan")
    elig_scores = df.loc[eligible_mask, "js_score"]
    df.loc[eligible_mask, "percentile"] = elig_scores.rank(pct=True, na_option="keep") * 100
    df.to_csv(csv_path, index=False)
    log(f"Wrote {csv_path} ({len(df)} rows, {eligible_mask.sum()} eligible for percentile)")

    # Sanity check known catalytic positions
    sanity = [195, 196, 175, 176, 226, 215, 258, 50, 50]
    log("Sanity check (catalytic / known):")
    for p in sorted(set(sanity)):
        if p <= len(df):
            r = df[df.ref_position == p].iloc[0]
            log(f"  pos {p} ({r['residue']}): js={r['js_score']:.3f}, pct={r['percentile']:.1f}, gap_frac={r['gap_fraction']:.2f}")

    # Plot with top-10% colored red
    fig, ax = plt.subplots(figsize=(15, 4.5))
    xs = df["ref_position"].values
    ys = np.array([0 if pd.isna(v) else v for v in df["js_score"].values])
    p90 = np.nanpercentile(df.loc[eligible_mask, "js_score"], 90)
    colors = []
    for _, row in df.iterrows():
        if pd.isna(row["js_score"]):
            colors.append("#cccccc")
        elif row["gap_fraction"] > 0.5:
            colors.append("#bbbbbb")
        elif row["js_score"] >= p90:
            colors.append("#c0392b")
        else:
            colors.append("#3498db")
    ax.bar(xs, ys, color=colors, width=1.0)

    # Mark catalytic residues
    catalytic = [195, 196, 175, 176, 215, 226, 258]
    for cp in catalytic:
        if cp <= len(df):
            v = df[df.ref_position == cp].iloc[0]["js_score"]
            if not pd.isna(v):
                ax.scatter([cp], [v + 0.02], marker="v", color="black", s=30, zorder=5)
                ax.annotate(str(cp), (cp, v + 0.04), fontsize=7, ha="center")

    ax.set_xlabel("Reference position (P04818)")
    ax.set_ylabel("JSD conservation (Capra-Singh, weighted window)")
    ax.set_title(f"v2 TYMS conservation across {len(seqs)} orthologs (top 10% red, catalytic ▼)")
    ax.axhline(p90, color="black", linestyle="--", alpha=0.5, label=f"p90={p90:.3f}")
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=130)
    plt.close()
    log(f"Wrote {plot_path}")

    # Also save metadata
    meta = {
        "valid_orthologs": valid,
        "rejected": [(a, l, r) for a, l, r in failed],
        "n_aligned": len(seqs),
        "aln_length": aln_len,
        "ref_idx": ref_idx,
        "p90_threshold": float(p90),
        "n_eligible_columns": int(eligible_mask.sum()),
    }
    with open(os.path.join(MSA_DIR, "msa_v2_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    log("Stage 1 v2 DONE")


if __name__ == "__main__":
    main()
