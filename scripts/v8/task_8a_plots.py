#!/usr/bin/env python3
"""Phase 8a plots: Vinardo vs Vina scatter (static PNG + interactive HTML)."""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "13_phase8" / "01_alt_scoring"


def load_rows() -> list[dict]:
    rows = []
    with (OUT / "alt_scoring_results.csv").open() as fh:
        for r in csv.DictReader(fh):
            try:
                r["vina_score"] = float(r["vina_score"])
                r["vinardo_score"] = float(r["vinardo_score"])
            except ValueError:
                continue
            rows.append(r)
    return rows


def main() -> int:
    rows = load_rows()
    if not rows:
        print("no scorable rows; aborting plot")
        return 1

    apo = [r for r in rows if r["condition"] == "apo"]
    holo = [r for r in rows if r["condition"] == "holo"]

    # static PNG via matplotlib
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 6))
    if apo:
        ax.scatter([r["vina_score"] for r in apo], [r["vinardo_score"] for r in apo],
                   c="#1f77b4", label=f"apo (n={len(apo)})", s=70, alpha=0.85,
                   edgecolors="white", linewidth=0.8)
    if holo:
        ax.scatter([r["vina_score"] for r in holo], [r["vinardo_score"] for r in holo],
                   c="#d62728", label=f"holo (n={len(holo)})", s=70, alpha=0.85,
                   edgecolors="white", linewidth=0.8)
    all_vals = [r["vina_score"] for r in rows] + [r["vinardo_score"] for r in rows]
    lo, hi = min(all_vals) - 0.3, max(all_vals) + 0.3
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.4, lw=1, label="y = x")
    ax.set_xlabel("Vina score (kcal/mol)")
    ax.set_ylabel("Vinardo score (kcal/mol)")
    ax.set_title("Phase 8a: Vinardo vs Vina re-scoring on Vina top poses")
    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    png = OUT / "alt_scoring_compare.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print(f"wrote {png}")

    # interactive plotly HTML
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("plotly not available; skipping HTML")
        return 0

    fig = go.Figure()
    for cond, color in [("apo", "#1f77b4"), ("holo", "#d62728")]:
        sub = [r for r in rows if r["condition"] == cond]
        if not sub:
            continue
        fig.add_trace(go.Scatter(
            x=[r["vina_score"] for r in sub],
            y=[r["vinardo_score"] for r in sub],
            mode="markers",
            name=f"{cond} (n={len(sub)})",
            marker=dict(color=color, size=11, line=dict(color="white", width=1)),
            text=[r["label"] for r in sub],
            hovertemplate="<b>%{text}</b><br>Vina: %{x:.2f}<br>Vinardo: %{y:.2f}<extra></extra>",
        ))
    fig.add_trace(go.Scatter(
        x=[lo, hi], y=[lo, hi], mode="lines",
        line=dict(color="black", dash="dash", width=1),
        showlegend=False, hoverinfo="skip",
    ))
    fig.update_layout(
        title="Phase 8a: Vinardo vs Vina re-scoring (interactive)",
        xaxis_title="Vina score (kcal/mol)",
        yaxis_title="Vinardo score (kcal/mol)",
        width=820, height=620, hovermode="closest",
        template="plotly_white",
    )
    html = OUT / "alt_scoring_compare.html"
    fig.write_html(str(html), include_plotlyjs="cdn")
    print(f"wrote {html}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
