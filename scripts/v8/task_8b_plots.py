#!/usr/bin/env python3
"""Phase 8b plots: rigid Vina vs flex Vina scatter for the 8 priority mutants."""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "13_phase8" / "02_flexres"


def load_rows() -> list[dict]:
    rows = []
    with (OUT / "flexres_compare.csv").open() as fh:
        for r in csv.DictReader(fh):
            try:
                r["rigid_vina_score"] = float(r["rigid_vina_score"])
                r["flex_vina_score"] = float(r["flex_vina_score"])
            except ValueError:
                continue
            rows.append(r)
    return rows


def main() -> int:
    rows = load_rows()
    if not rows:
        print("no flex rows scorable; skip plot")
        return 0

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 6))
    xs = [r["rigid_vina_score"] for r in rows]
    ys = [r["flex_vina_score"] for r in rows]
    ax.scatter(xs, ys, c="#2ca02c", s=90, alpha=0.85,
               edgecolors="white", linewidth=0.8)
    for r in rows:
        ax.annotate(r["label"], (r["rigid_vina_score"], r["flex_vina_score"]),
                    fontsize=8, alpha=0.7, xytext=(4, 4), textcoords="offset points")
    all_vals = xs + ys
    lo, hi = min(all_vals) - 0.5, max(all_vals) + 0.5
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.4, lw=1, label="y = x (no flex effect)")
    ax.set_xlabel("Rigid Vina top affinity (kcal/mol)")
    ax.set_ylabel("Flex Vina top affinity (kcal/mol)")
    ax.set_title("Phase 8b: Flex residue (14 active-site) vs Rigid Vina")
    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    png = OUT / "flex_vs_rigid.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print(f"wrote {png}")

    try:
        import plotly.graph_objects as go
    except ImportError:
        return 0
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="markers+text",
        name="mutant",
        marker=dict(color="#2ca02c", size=12, line=dict(color="white", width=1)),
        text=[r["label"] for r in rows],
        textposition="top center",
        hovertemplate="<b>%{text}</b><br>rigid: %{x:.2f}<br>flex: %{y:.2f}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=[lo, hi], y=[lo, hi], mode="lines",
        line=dict(color="black", dash="dash", width=1),
        showlegend=False, hoverinfo="skip",
    ))
    fig.update_layout(
        title="Phase 8b: Flex vs Rigid Vina (8 priority mutants, holo)",
        xaxis_title="Rigid Vina top affinity (kcal/mol)",
        yaxis_title="Flex Vina top affinity (kcal/mol)",
        width=820, height=620, hovermode="closest",
        template="plotly_white",
    )
    html = OUT / "flex_vs_rigid.html"
    fig.write_html(str(html), include_plotlyjs="cdn")
    print(f"wrote {html}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
