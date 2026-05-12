#!/usr/bin/env python3
"""
Build interactive Plotly versions of every static analysis plot.
Outputs to 11_enhanced/plotly/*.html (each ~20-40 KB, plotly.js loaded from CDN).
Hover for full per-mutant detail; click legend entries to filter; box-zoom; etc.
"""
import csv, io, pathlib, json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = pathlib.Path(".").resolve()
OUT = ROOT / "11_enhanced" / "plotly"
OUT.mkdir(parents=True, exist_ok=True)

# --- load v5 mutant data ---
CSV = ROOT / "07e_mut_docking_v5" / "mutant_results_v5.csv"
text = CSV.read_text()
lines = [l for l in text.splitlines() if not l.startswith("#")]
df = pd.read_csv(io.StringIO("\n".join(lines)))
df["delta_vina_vs_wt"] = pd.to_numeric(df["delta_vina_vs_wt"], errors="coerce")
df["top_affinity"] = pd.to_numeric(df["top_affinity"], errors="coerce")
df["rmsd_to_native"] = pd.to_numeric(df["rmsd_to_native"], errors="coerce")
df["n_modes"] = pd.to_numeric(df["n_modes"], errors="coerce")
df["mis_docked"] = df["mis_docked"].astype(str).str.lower().isin(("true","1","yes"))
df["low_confidence"] = df["low_confidence"].astype(str).str.lower().isin(("true","1","yes"))
mut_only = df[df["mutant"] != "WT"].copy()
mut_only["clean"] = ~(mut_only["mis_docked"] | mut_only["low_confidence"])

COLORS = {
    "ala_scan":"#1f77b4", "opposite":"#ff7f0e", "arg_clamp":"#2ca02c",
    "double_dyad":"#9467bd", "double_phosclamp":"#d62728",
    "double_polar_neutral":"#8c564b", "double_substrate_orient":"#e377c2",
    "double_aromatic_swap":"#7f7f7f", "control_surface":"#bcbd22", "explore_g217w":"#17becf",
}
NOISE = 0.85

PLOTLY_CONFIG = {"displaylogo": False, "responsive": True,
                 "modeBarButtonsToRemove": ["lasso2d","select2d"]}

def save(fig, name, title):
    fig.update_layout(template="plotly_white",
                      title=dict(text=title, font=dict(size=14, color="#1f3b5e")),
                      margin=dict(l=60, r=30, t=70, b=60))
    html = fig.to_html(include_plotlyjs="cdn", full_html=True, config=PLOTLY_CONFIG)
    p = OUT / f"{name}.html"
    p.write_text(html)
    print(f"  {p.name}  ({p.stat().st_size//1024} KB)")


# ========================== 1) Δ Vina by category (box) ==========================
fig = px.box(mut_only, x="category", y="delta_vina_vs_wt", color="condition",
             color_discrete_map={"apo":"#1f77b4","holo":"#ff7f0e"},
             points="all", hover_name="mutant",
             hover_data={"top_affinity":":.3f","rmsd_to_native":":.3f",
                         "n_modes":True,"mis_docked":True,"low_confidence":True})
fig.add_hline(y=NOISE, line_dash="dash", line_color="black",
              annotation_text=f"+ noise floor ({NOISE})", annotation_position="top right")
fig.add_hline(y=-NOISE, line_dash="dash", line_color="black",
              annotation_text=f"− noise floor (−{NOISE})", annotation_position="bottom right")
fig.update_layout(xaxis_tickangle=-30,
                  yaxis_title="Δ Vina score vs WT (kcal/mol)",
                  xaxis_title="Mutant category",
                  height=560)
save(fig, "delta_vina_by_category",
     "Δ Vina by mutant category (apo + holo, all mutants)")

# ========================== 2) Per-mutant Δ Vina bar (apo vs holo) ==========================
piv = (mut_only
       .pivot_table(index="mutant", columns="condition", values="delta_vina_vs_wt")
       .reset_index()
       .sort_values(by="holo", ascending=False, na_position="last"))
fig = go.Figure()
fig.add_trace(go.Bar(x=piv["mutant"], y=piv["apo"], name="apo",
                     marker_color="#1f77b4", opacity=0.85,
                     hovertemplate="<b>%{x}</b> apo<br>Δ Vina = %{y:.3f}<extra></extra>"))
fig.add_trace(go.Bar(x=piv["mutant"], y=piv["holo"], name="holo",
                     marker_color="#ff7f0e", opacity=0.85,
                     hovertemplate="<b>%{x}</b> holo<br>Δ Vina = %{y:.3f}<extra></extra>"))
fig.add_hline(y=NOISE, line_dash="dash", line_color="black",
              annotation_text=f"+ noise floor ({NOISE})", annotation_position="top right")
fig.add_hline(y=-NOISE, line_dash="dash", line_color="black",
              annotation_text=f"− noise floor (−{NOISE})", annotation_position="bottom right")
fig.update_layout(xaxis_tickangle=-45, barmode="group",
                  yaxis_title="Δ Vina score vs WT (kcal/mol)",
                  height=560, legend_title="Cofactor condition")
save(fig, "delta_vina_apo_holo",
     "Per-mutant Δ Vina — apo vs holo side-by-side (sorted by holo)")

# ========================== 3) Apo vs holo paired (scatter with y=x line) ==========================
piv2 = piv.dropna()
fig = px.scatter(piv2, x="apo", y="holo", text="mutant",
                 color_discrete_sequence=["#2b6f9c"], height=560,
                 hover_data={"mutant": True, "apo":":.3f","holo":":.3f"})
fig.update_traces(marker=dict(size=12, line=dict(width=1, color="white")),
                  textposition="top center", textfont=dict(size=10))
lo, hi = piv2[["apo","holo"]].min().min()-0.2, piv2[["apo","holo"]].max().max()+0.2
fig.add_shape(type="line", x0=lo, y0=lo, x1=hi, y1=hi,
              line=dict(dash="dash", color="#888", width=1))
fig.add_annotation(x=hi-0.1, y=hi-0.1, text="y = x", showarrow=False,
                   font=dict(size=10, color="#666"))
# noise floor band
fig.add_hrect(y0=-NOISE, y1=NOISE, fillcolor="#fafafa", opacity=0.5, line_width=0)
fig.add_vrect(x0=-NOISE, x1=NOISE, fillcolor="#fafafa", opacity=0.5, line_width=0)
fig.update_layout(xaxis_title="Δ Vina (apo) — kcal/mol",
                  yaxis_title="Δ Vina (holo) — kcal/mol")
save(fig, "delta_vina_apo_vs_holo",
     "Apo vs holo Δ Vina — concordance scatter (light band = ±0.85 noise floor)")

# ========================== 4) v3 Apo vs Holo concordance with category colour ==========================
piv3 = (mut_only.pivot_table(index=["mutant","category"], columns="condition",
                              values="delta_vina_vs_wt")
        .reset_index().dropna())
fig = px.scatter(piv3, x="apo", y="holo", color="category", text="mutant",
                 color_discrete_map=COLORS, height=560,
                 hover_data={"category":True,"apo":":.3f","holo":":.3f"})
fig.update_traces(marker=dict(size=14, line=dict(width=1, color="white")),
                  textposition="top center", textfont=dict(size=10))
lo, hi = piv3[["apo","holo"]].min().min()-0.2, piv3[["apo","holo"]].max().max()+0.2
fig.add_shape(type="line", x0=lo, y0=lo, x1=hi, y1=hi,
              line=dict(dash="dash", color="#888", width=1))
fig.add_hrect(y0=-NOISE, y1=NOISE, fillcolor="#fafafa", opacity=0.4, line_width=0)
fig.add_vrect(x0=-NOISE, x1=NOISE, fillcolor="#fafafa", opacity=0.4, line_width=0)
fig.update_layout(xaxis_title="Δ Vina (apo) — kcal/mol",
                  yaxis_title="Δ Vina (holo) — kcal/mol",
                  legend_title="Functional category")
save(fig, "apo_vs_holo_concordance",
     "Apo vs holo concordance — coloured by functional category")

# ========================== 5) v3 Mutation-effect map (Δ Vina vs pose-RMSD) ==========================
fig = make_subplots(rows=1, cols=2, subplot_titles=("Apo (well-docked + mis-docked faded)",
                                                     "Holo (well-docked + mis-docked faded)"),
                    horizontal_spacing=0.10)
for j, cond in enumerate(("apo","holo")):
    sub = mut_only[mut_only["condition"]==cond].copy()
    well = sub[~(sub["mis_docked"] | sub["low_confidence"])]
    bad  = sub[(sub["mis_docked"] | sub["low_confidence"])]
    for cat, dfc in well.groupby("category"):
        fig.add_trace(go.Scatter(x=dfc["delta_vina_vs_wt"], y=dfc["rmsd_to_native"],
                                  mode="markers+text", text=dfc["mutant"],
                                  textposition="top center", textfont=dict(size=9),
                                  marker=dict(size=12, color=COLORS.get(cat,"#888"),
                                              line=dict(color="white",width=1)),
                                  name=cat, legendgroup=cat,
                                  showlegend=(j==0),
                                  hovertemplate=f"<b>%{{text}}</b> ({cond})<br>"
                                                f"Δ Vina = %{{x:.3f}}<br>RMSD = %{{y:.3f}} Å"
                                                f"<extra>{cat}</extra>"),
                       row=1, col=j+1)
    if len(bad):
        fig.add_trace(go.Scatter(x=bad["delta_vina_vs_wt"], y=bad["rmsd_to_native"],
                                  mode="markers+text", text=bad["mutant"],
                                  textposition="top center", textfont=dict(size=9, color="#aaa"),
                                  marker=dict(size=11, color="#cccccc",
                                              symbol="circle-open"),
                                  name="mis-docked / low-conf", legendgroup="bad",
                                  showlegend=(j==0),
                                  hovertemplate="<b>%{text}</b><br>flagged<br>"
                                                "Δ Vina = %{x:.3f}<br>RMSD = %{y:.3f} Å<extra></extra>"),
                       row=1, col=j+1)
    fig.add_vline(x=0, line_color="#bbb", line_width=1, row=1, col=j+1)
    fig.add_vline(x=NOISE, line_color="black", line_dash="dash", line_width=1, row=1, col=j+1)
    fig.add_vline(x=-NOISE, line_color="black", line_dash="dash", line_width=1, row=1, col=j+1)
    fig.add_hline(y=3.0, line_color="#aaa", line_dash="dot", line_width=1, row=1, col=j+1)
    fig.update_xaxes(title_text="Δ Vina score vs WT (kcal/mol)", row=1, col=j+1)
    fig.update_yaxes(title_text="Top-pose RMSD vs crystal dUMP (Å)", row=1, col=j+1)
fig.update_layout(height=600, legend_title="Functional category")
save(fig, "mutation_effect_map",
     "Mutation-effect map — Δ Vina vs pose-RMSD (apo | holo). Dashed verticals = ±0.85 noise floor.")

# ========================== 6) Modeller — per-model quality ==========================
# Read RMSD + scores + Ramachandran stats
SC = pd.read_csv(ROOT / "10_modeller/04_modeller_run/scores.csv")
RMSD = pd.read_csv(ROOT / "10_modeller/05_comparison/rmsd_per_model.csv")
RAMA = pd.read_csv(ROOT / "10_modeller/06_validation/ramachandran_stats.csv")
mer = SC.merge(RMSD, on="model_id").merge(RAMA, on="model_id")
mer["model"] = "model " + mer["model_id"].astype(str)
fig = make_subplots(rows=2, cols=2,
                    subplot_titles=("DOPE (lower is better)", "molpdf (lower is better)",
                                    "Cα RMSD vs 1HVY (Å)", "Ramachandran % favoured"),
                    vertical_spacing=0.16, horizontal_spacing=0.12)
fig.add_trace(go.Bar(x=mer["model"], y=mer["DOPE"], marker_color="#2b6f9c",
                     hovertemplate="%{x}<br>DOPE = %{y:.1f}<extra></extra>"),
              row=1, col=1)
fig.add_trace(go.Bar(x=mer["model"], y=mer["molpdf"], marker_color="#2f8a6f",
                     hovertemplate="%{x}<br>molpdf = %{y:.1f}<extra></extra>"),
              row=1, col=2)
fig.add_trace(go.Bar(x=mer["model"], y=mer["rmsd_to_crystal"], marker_color="#caa44a",
                     hovertemplate="%{x}<br>Cα RMSD = %{y:.3f} Å<extra></extra>"),
              row=2, col=1)
fig.add_trace(go.Bar(x=mer["model"], y=mer["pct_favoured"], marker_color="#b8593c",
                     hovertemplate="%{x}<br>%favoured = %{y:.2f}<extra></extra>"),
              row=2, col=2)
fig.update_xaxes(tickangle=-30)
fig.update_layout(height=720, showlegend=False)
save(fig, "modeller_quality_overview",
     "Modeller — per-model quality overview (10 models)")

print("\nDone.")
