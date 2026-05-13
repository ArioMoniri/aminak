#!/usr/bin/env python3
"""Task E: build a phylogenetic tree from the 10 verified TYMS orthologs."""
from __future__ import annotations
import io
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

import os
PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "12_phase7" / "05_phylogeny"
LOG = PROJECT / "logs" / "v7_task_e.log"
PIPELOG = PROJECT / "pipeline.log"
ALIGN = PROJECT / "01b_msa_v2" / "aligned.fa"

KINGDOM = {
    "Homo_sapiens": ("animal", "#d62728"),
    "Mus_musculus": ("animal", "#ff7f0e"),
    "Rattus_norvegicus": ("animal", "#bcbd22"),
    "Drosophila_melanogaster": ("animal", "#e377c2"),
    "Saccharomyces_cerevisiae": ("fungus", "#9467bd"),
    "Arabidopsis_thaliana": ("plant", "#2ca02c"),
    "Plasmodium_falciparum": ("protist", "#17becf"),
    "Escherichia_coli": ("bacterium", "#1f77b4"),
    "Lactobacillus_casei": ("bacterium", "#7f7f7f"),
    "Bacteriophage_T4": ("virus", "#8c564b"),
}


def log(msg: str) -> None:
    line = f"[V7][taskE] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as fh: fh.write(line+"\n")
    with PIPELOG.open("a") as fh: fh.write(line+"\n")


def species_of(name: str) -> str:
    return name.split("|")[0]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    log("Task E start")
    aln = AlignIO.read(str(ALIGN), "fasta")
    log(f"alignment: {len(aln)} taxa, {aln.get_alignment_length()} columns")

    # Distance matrix using identity model
    calc = DistanceCalculator("blosum62")
    dm = calc.get_distance(aln)
    constructor = DistanceTreeConstructor(calc, "nj")
    tree_nj = constructor.build_tree(aln)
    constructor_upgma = DistanceTreeConstructor(calc, "upgma")
    tree_upgma = constructor_upgma.build_tree(aln)
    # Strip the inner node labels for cleanness
    for t in (tree_nj, tree_upgma):
        for clade in t.get_nonterminals():
            clade.name = None
    nwk_path = OUT / "tymS_tree.nwk"
    Phylo.write([tree_nj], str(nwk_path), "newick")
    Phylo.write([tree_upgma], str(OUT / "tymS_tree_upgma.nwk"), "newick")
    log(f"NJ Newick -> {nwk_path}")

    # matplotlib cladogram
    # Re-shorten leaf labels (drop UniProt accession)
    for clade in tree_nj.get_terminals():
        sp = species_of(clade.name)
        clade.name = sp
    fig, ax = plt.subplots(figsize=(10, 7))
    Phylo.draw(tree_nj, do_show=False, axes=ax,
               label_func=lambda c: c.name if c.is_terminal() else "",
               label_colors=lambda n: KINGDOM.get(n, ("?", "black"))[1])
    ax.set_title("TYMS NJ tree (BLOSUM62 distance, 10 orthologs)")
    fig.tight_layout()
    png = OUT / "tymS_tree.png"
    fig.savefig(png, dpi=140)
    plt.close(fig)
    log(f"matplotlib tree -> {png}")

    # Plotly tree (manual rendering by walking the tree)
    # Compute leaf y-coords by post-order traversal, x-coords by depth
    leaf_order = tree_nj.get_terminals()
    y_of = {leaf.name: i for i, leaf in enumerate(leaf_order)}

    def assign_xy(clade, x=0.0):
        if clade.is_terminal():
            return x + (clade.branch_length or 0.0), y_of[clade.name]
        ys = []
        for child in clade.clades:
            cx, cy = assign_xy(child, x + (clade.branch_length or 0.0))
            ys.append(cy)
            child._x = cx; child._y = cy
        clade._x = x + (clade.branch_length or 0.0)
        clade._y = sum(ys) / len(ys)
        return clade._x, clade._y

    tree_nj.root._x = 0.0; tree_nj.root._y = 0.0
    assign_xy(tree_nj.root, 0.0)
    # ensure leaves have coords too
    for leaf in tree_nj.get_terminals():
        if not hasattr(leaf, "_x"):
            leaf._x = (leaf.branch_length or 0.0)
            leaf._y = y_of[leaf.name]

    edges_x, edges_y = [], []
    def walk(parent):
        for child in parent.clades:
            # vertical line at parent._x from parent._y to child._y
            edges_x.extend([parent._x, parent._x, None])
            edges_y.extend([parent._y, child._y, None])
            # horizontal line from parent._x to child._x
            edges_x.extend([parent._x, child._x, None])
            edges_y.extend([child._y, child._y, None])
            walk(child)
    walk(tree_nj.root)

    leaf_x = [c._x for c in leaf_order]
    leaf_y = [c._y for c in leaf_order]
    leaf_names = [c.name for c in leaf_order]
    leaf_colors = [KINGDOM.get(n, ("?", "black"))[1] for n in leaf_names]
    leaf_groups = [KINGDOM.get(n, ("?", "black"))[0] for n in leaf_names]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edges_x, y=edges_y, mode="lines",
                              line=dict(color="grey", width=2),
                              hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=leaf_x, y=leaf_y, mode="markers+text",
                              text=[f"{n}  [{g}]" for n, g in zip(leaf_names, leaf_groups)],
                              textposition="middle right",
                              marker=dict(color=leaf_colors, size=14, line=dict(color="black", width=1)),
                              hovertext=leaf_names, name="taxa"))
    fig.update_layout(template="plotly_white", height=600,
                      title="TYMS phylogeny (NJ, BLOSUM62)",
                      xaxis=dict(title="Branch length (substitutions/site)", showgrid=True),
                      yaxis=dict(showticklabels=False, showgrid=False),
                      margin=dict(l=20, r=300, t=60, b=40))
    html = OUT / "tymS_tree.html"
    fig.write_html(html, include_plotlyjs="cdn")
    log(f"plotly tree -> {html}")

    # Distance matrix CSV
    with (OUT / "distance_matrix.csv").open("w") as fh:
        names = [species_of(n) for n in dm.names]
        fh.write("," + ",".join(names) + "\n")
        for i, n in enumerate(names):
            row = [f"{dm[i,j]:.4f}" for j in range(len(names))]
            fh.write(n + "," + ",".join(row) + "\n")

    summary = {
        "n_taxa": len(aln),
        "alignment_length": aln.get_alignment_length(),
        "distance_model": "blosum62",
        "tree_methods": ["nj", "upgma"],
        "files": {
            "newick_nj": str(nwk_path),
            "newick_upgma": str(OUT / "tymS_tree_upgma.nwk"),
            "png": str(png),
            "html": str(html),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    log("Task E done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
