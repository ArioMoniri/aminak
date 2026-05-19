#!/usr/bin/env python3
"""Phase 14d evidence package for the under-explored druggable cavity 18.

Produces:
  - downloads/cavity18_apo.pdb              — apo dimer with cavity-18 residues annotated
  - downloads/cavity18_apo_pocket.pdb       — only the cavity-18 residues (B-factor = druggability)
  - downloads/cavity18_indazole_complex.pdb — apo + 1H-indazole top pose
  - downloads/cavity18_ibuprofen_complex.pdb — apo + ibuprofen top pose
  - downloads/cavity18_residues.csv         — residue × ortholog × conservation × mutation table
  - viewers/cavity18_apo.html               — 3Dmol viewer, apo, pocket coloured wheat,
                                              loop 181-197 (Anderson/Pozzi allosteric) coloured red
  - viewers/cavity18_indazole.html          — same + indazole top pose
  - viewers/cavity18_ibuprofen.html         — same + ibuprofen top pose
  - figures/cavity18_conservation.png       — JS-score bar plot, cavity residues vs whole protein
  - figures/cavity18_phylogeny_annot.png    — phylogeny tree annotated with #cavity-18 mutations
"""
from __future__ import annotations
import json, csv, re, shutil, subprocess
from pathlib import Path
from collections import defaultdict, OrderedDict

REPO = Path(__file__).resolve().parents[2]
STRAT = REPO / "14_inhibitor_design" / "04_allosteric"
OUT = STRAT / "cavity18_evidence"
DL  = OUT / "downloads"
VW  = OUT / "viewers"
FIG = OUT / "figures"
for d in (DL, VW, FIG): d.mkdir(parents=True, exist_ok=True)

APO_PDB = STRAT / "apo_for_fpocket.pdb"
INDAZOLE_POSE = STRAT / "poses" / "cav18_CID7032_pose1.pdb"
IBUPROFEN_POSE = STRAT / "poses" / "cav18_CID3672_pose1.pdb"
CONS_CSV = REPO / "01b_msa_v2" / "conservation_scores.csv"
MSA_FA = REPO / "01b_msa_v2" / "aligned.fa"
PHYLO_NWK = REPO / "12_phase7" / "05_phylogeny" / "tymS_tree.nwk"

# Cavity-18 residue set (parsed from FPocket pocket18_atm.pdb earlier)
# Format: (chain_in_PDB, resid_in_PDB)
def parse_pocket_residues():
    p = STRAT / "apo_for_fpocket_out" / "pockets" / "pocket18_atm.pdb"
    residues = set()
    for ln in p.read_text().splitlines():
        if ln.startswith("ATOM"):
            chain = ln[21].strip() or "A"
            try: rid = int(ln[22:26])
            except ValueError: continue
            residues.add((chain, rid))
    return sorted(residues)

# loop 181-197 = Anderson 2012 / Pozzi 2019 long-range allosteric communication zone
ALLOSTERIC_LOOP_RANGE = (181, 197)


def write_apo_with_chain_ids(out_pdb: Path):
    """The apo_for_fpocket.pdb has chain id 'A' for everything; split by atom count."""
    lines = APO_PDB.read_text().splitlines()
    n = sum(1 for ln in lines if ln.startswith("ATOM")) // 2
    out_lines = []
    atom_idx = 0
    for ln in lines:
        if ln.startswith("ATOM"):
            chain = "A" if atom_idx < n else "B"
            ln = ln[:21] + chain + ln[22:]
            atom_idx += 1
        if ln.strip().startswith(("ATOM","HETATM","TER","END")):
            out_lines.append(ln)
    out_pdb.write_text("\n".join(out_lines) + "\nEND\n")


def write_complex(apo_pdb: Path, ligand_pdb: Path, out_pdb: Path, lig_resname: str):
    apo_body = []
    for ln in apo_pdb.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM","TER")):
            apo_body.append(ln)
    lig_body = []
    for ln in ligand_pdb.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM")):
            # rewrite to HETATM, force chain Z, force resname
            new = "HETATM" + ln[6:17] + lig_resname.ljust(3) + " Z" + ln[22:]
            lig_body.append(new)
    out_pdb.write_text("\n".join(apo_body) + "\nTER\n" + "\n".join(lig_body) + "\nEND\n")


def write_pocket_only(apo_pdb: Path, residues: list[tuple[str,int]], out_pdb: Path, drug_score=0.994):
    """Pocket residues only, B-factor set to FPocket druggability score for shading."""
    resset = set(residues)
    out_lines = []
    for ln in apo_pdb.read_text().splitlines():
        if ln.startswith("ATOM"):
            chain = ln[21].strip() or "A"
            try: rid = int(ln[22:26])
            except ValueError: continue
            if (chain, rid) in resset:
                ln = ln[:60] + f"{drug_score*100:6.2f}" + ln[66:]
                out_lines.append(ln)
    out_pdb.write_text("\n".join(out_lines) + "\nEND\n")


def load_conservation():
    """Return dict {1-based ref position: js_score}"""
    out = {}
    with CONS_CSV.open() as f:
        rd = csv.DictReader(f)
        for r in rd:
            try:
                pos = int(r["ref_position"])
                out[pos] = float(r["js_score"])
            except (ValueError, KeyError): continue
    return out


def load_msa():
    """Return OrderedDict {ortholog_label: sequence}"""
    seqs = OrderedDict()
    cur = None
    cur_seq = []
    for ln in MSA_FA.read_text().splitlines():
        if ln.startswith(">"):
            if cur is not None: seqs[cur] = "".join(cur_seq)
            cur = ln[1:].strip(); cur_seq = []
        else:
            cur_seq.append(ln.strip())
    if cur is not None: seqs[cur] = "".join(cur_seq)
    return seqs


def aln_col_for_ref_position(human_seq: str, ref_pos: int) -> int | None:
    """Walk the aligned human sequence; return alignment column index (0-based) for the given
    1-based reference position. Returns None if ref_pos exceeds the un-gapped sequence length."""
    n = 0
    for i, ch in enumerate(human_seq):
        if ch != "-":
            n += 1
            if n == ref_pos: return i
    return None


def build_residue_table(residues: list[tuple[str,int]], cons: dict, seqs: OrderedDict):
    """For each cavity-18 residue: identity in each ortholog + JS conservation + mutation flag."""
    human_label = next(k for k in seqs if k.startswith("Homo_sapiens"))
    human_seq = seqs[human_label]
    others = [k for k in seqs if k != human_label]

    rows = []
    for chain, rid in residues:
        # Build a label like "B/L196" (1-letter code from human seq + resid)
        # Find what amino acid the residue is in the human sequence
        col = aln_col_for_ref_position(human_seq, rid) if rid <= 313 else None
        human_aa = human_seq[col] if col is not None else "?"
        # Identity in each ortholog
        ortho_aas = {}
        for k in others:
            ortho_aas[k] = seqs[k][col] if col is not None else "?"
        # Mutations = non-identical (treating - as deletion)
        n_diff = sum(1 for v in ortho_aas.values() if v not in ("-", human_aa))
        n_gap  = sum(1 for v in ortho_aas.values() if v == "-")
        # Conservation score (from Phase 1b — JS divergence; high = conserved)
        js = cons.get(rid)
        # Is this residue on the allosteric loop 181-197?
        on_allosteric_loop = ALLOSTERIC_LOOP_RANGE[0] <= rid <= ALLOSTERIC_LOOP_RANGE[1]
        rows.append({
            "chain": chain, "resid": rid, "human_aa": human_aa, "js_score": js,
            "n_orthologs_differ": n_diff, "n_orthologs_gap": n_gap,
            "on_allosteric_loop_181_197": on_allosteric_loop,
            **{f"aa_{k.split('|')[0]}": v for k, v in ortho_aas.items()}
        })
    # sort by resid then chain
    rows.sort(key=lambda r: (r["resid"], r["chain"]))
    return rows


def write_residue_csv(rows, out_csv: Path):
    if not rows: return
    fields = list(rows[0].keys())
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def plot_conservation(rows, cons_all: dict, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    # Sort cavity residues by resid
    cav_rows = sorted(rows, key=lambda r: r["resid"])
    cav_residues = sorted(set(r["resid"] for r in cav_rows))
    cav_js = [cons_all.get(r, 0.0) for r in cav_residues]
    cav_lbl = [f"{r['human_aa']}{r['resid']}" for r in [next(x for x in cav_rows if x['resid']==pos) for pos in cav_residues]]
    cav_loop = [(ALLOSTERIC_LOOP_RANGE[0] <= pos <= ALLOSTERIC_LOOP_RANGE[1]) for pos in cav_residues]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6.5),
                                    gridspec_kw={'height_ratios': [1.4, 2]})
    # Top: full-protein JS landscape, cavity residues highlighted
    all_pos = sorted(cons_all.keys())
    all_js = [cons_all[p] for p in all_pos]
    ax1.plot(all_pos, all_js, color="lightgrey", linewidth=0.8, alpha=0.7, zorder=1)
    cav_set = set(cav_residues)
    cav_x = [p for p in all_pos if p in cav_set]
    cav_y = [cons_all[p] for p in cav_x]
    ax1.scatter(cav_x, cav_y, color="orange", s=22, zorder=3, label="cavity-18 residues")
    # highlight loop 181-197
    ax1.axvspan(ALLOSTERIC_LOOP_RANGE[0], ALLOSTERIC_LOOP_RANGE[1], color="red", alpha=0.18,
                label="allosteric loop 181-197 (Anderson 2012, Pozzi 2019)")
    ax1.set_xlabel("TYMS residue (human, 1HVY numbering)")
    ax1.set_ylabel("JS conservation score\n(higher = more conserved)")
    ax1.set_title("Cavity-18 residues vs the full TYMS conservation landscape")
    ax1.legend(loc="upper right", fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Bottom: cavity residues only, bars coloured by allosteric-loop membership
    x = np.arange(len(cav_residues))
    colors = ["#c0392b" if cav_loop[i] else "#e67e22" for i in range(len(cav_residues))]
    bars = ax2.bar(x, cav_js, color=colors)
    ax2.set_xticks(x); ax2.set_xticklabels(cav_lbl, rotation=45, ha="right", fontsize=8)
    ax2.set_ylabel("JS conservation score")
    ax2.set_title("Cavity-18 residues — conservation per position\n(red = on allosteric loop 181-197, orange = other cavity walls)")
    ax2.axhline(y=np.median(all_js), ls="--", color="green", alpha=0.5,
                label=f"protein-wide median JS = {np.median(all_js):.2f}")
    ax2.axhline(y=np.percentile(all_js, 80), ls=":", color="purple", alpha=0.5,
                label=f"80th-pctile JS = {np.percentile(all_js, 80):.2f}")
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=140)
    plt.close()


def plot_phylogeny_annotated(rows, out_png: Path):
    """Read Newick, annotate each taxon with #cavity-18 mutations vs human."""
    try:
        from Bio import Phylo
    except ImportError:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not PHYLO_NWK.exists():
        return
    tree = Phylo.read(str(PHYLO_NWK), "newick")
    # mutations per taxon (vs Homo_sapiens)
    seqs = load_msa()
    human_label = next(k for k in seqs if k.startswith("Homo_sapiens"))
    human_seq = seqs[human_label]
    # collect alignment columns for cavity residues
    cav_resids = sorted(set(r["resid"] for r in rows))
    cav_cols = []
    for rid in cav_resids:
        col = aln_col_for_ref_position(human_seq, rid)
        if col is not None: cav_cols.append((rid, col))
    # taxon → list[mutations]
    muts_per_taxon = {}
    for label, seq in seqs.items():
        muts = []
        for rid, col in cav_cols:
            aa = seq[col]
            ref = human_seq[col]
            if aa not in ("-", ref):
                muts.append(f"{ref}{rid}{aa}")
        muts_per_taxon[label.split("|")[0]] = muts
    # Annotate tree leaves
    for term in tree.get_terminals():
        # canonical: term.name should match the label up to "|"
        name = term.name.split("|")[0] if term.name else ""
        n_muts = len(muts_per_taxon.get(name, []))
        term.name = f"{name} ({n_muts} mut)" if n_muts > 0 else f"{name} (=)"

    fig, ax = plt.subplots(figsize=(12, 6))
    Phylo.draw(tree, axes=ax, show_confidence=False, do_show=False)
    ax.set_title("TYMS ortholog phylogeny — #cavity-18 mutations (vs human reference)\n"
                 "'=' = identical at all cavity-18 positions; (N mut) = N substitutions")
    plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()
    # also write a mutations-per-taxon JSON
    mp = OUT / "downloads" / "cavity18_mutations_per_taxon.json"
    mp.write_text(json.dumps(muts_per_taxon, indent=2))


def build_viewer(html_path: Path, title: str, pdb_files: list[tuple[str, Path]],
                  cavity_residues: list[int], ligand_resname: str | None):
    """Build a self-contained 3Dmol.js viewer.
    pdb_files: list of (label, path) — each PDB content is embedded as a string."""
    # Inline PDB content
    pdb_blocks = []
    for label, p in pdb_files:
        txt = p.read_text()
        # escape backticks
        txt = txt.replace("\\","\\\\").replace("`","\\`")
        pdb_blocks.append(f"const pdb_{label} = `{txt}`;")
    cav_resids_js = ",".join(str(r) for r in sorted(set(cavity_residues)))
    loop_lo, loop_hi = ALLOSTERIC_LOOP_RANGE

    lig_section = ""
    show_lig_buttons = ""
    if ligand_resname:
        lig_section = f"""
  // ligand
  m.setStyle({{resn: '{ligand_resname}'}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.18}}}});
  // polar contacts (yellow dashes) ≤ 3.5 Å between ligand N/O and protein N/O
"""
        show_lig_buttons = f"""
  <button onclick="zoomLigand()">Zoom to ligand</button>
  <button onclick="hideLigand()">Hide ligand</button>
  <button onclick="showLigand()">Show ligand</button>"""

    html = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>{title}</title>
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<style>
 body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif;
        margin: 0; padding: 0 16px 12px 16px; background:#0c1116; color:#e8eef5; }}
 h1 {{ font-size: 17px; font-weight: 600; margin: 12px 0 4px; color:#9ad5ff; }}
 .meta {{ font-size: 12px; color:#aab8c8; margin: 0 0 8px; line-height:1.4; }}
 .panel {{ background:#1a2230; border:1px solid #283344; border-radius:8px;
          padding:0; overflow:hidden; }}
 .legend {{ font-size: 11.5px; color:#cad5e6; padding: 6px 12px; border-top:1px solid #283344; }}
 .legend span {{ display:inline-block; padding: 2px 8px; border-radius: 4px;
                margin: 2px 4px 2px 0; }}
 .controls {{ padding: 6px 12px; border-top:1px solid #283344;
              font-size: 12px; color:#cad5e6; }}
 .controls button {{ background:#283b53; color:#cfe6ff; border:1px solid #3a4d65;
                     border-radius: 4px; padding: 4px 10px; margin: 2px 4px 2px 0;
                     cursor: pointer; font-size: 11.5px; }}
 .controls button:hover {{ background:#3a4d65; }}
 .nav a {{ color:#9ad5ff; margin-right: 12px; font-size: 12px; text-decoration:none; }}
 .nav a:hover {{ text-decoration:underline; }}
</style>
</head><body>
<div class="nav" style="margin: 10px 0 4px;">
  <a href="../../../README.md">← Main README</a>
  <a href="../README.md">← Phase 14 README</a>
  <a href="https://github.com/ArioMoniri/aminak">GitHub</a>
</div>
<h1>{title}</h1>
<p class="meta">
TYMS dimer (apo, Phase-6c-hardened receptor). <strong style="color:#e8b250">Wheat sticks</strong>: cavity-18 residues
(FPocket druggability 0.994). <strong style="color:#e85050">Red sticks</strong>: residues that are *also* on the
published allosteric communication loop 181–197 (Anderson 2012, Pozzi 2019). Click & drag to rotate, scroll to zoom.
</p>

<div class="panel">
 <div id="gl" style="height:560px; width:100%;"></div>
 <div class="legend">
  <strong>Legend:</strong>
  <span style="background:#cccccc; color:#000">grey cartoon: protein backbone</span>
  <span style="background:#e8b250; color:#000">wheat sticks: cavity-18 residues</span>
  <span style="background:#e85050; color:#fff">red sticks: allosteric loop 181-197 ∩ cavity 18</span>
  {("<span style='background:#22cccc; color:#000'>cyan: ligand</span>" if ligand_resname else "")}
 </div>
 <div class="controls">
  <strong>Toggle:</strong>
  <button onclick="hideSurface()">Hide surface</button>
  <button onclick="showSurface()">Show surface</button>
  <button onclick="zoomPocket()">Zoom to cavity 18</button>
  <button onclick="zoomLoop()">Zoom to loop 181-197</button>
  <button onclick="spinOn()">Spin</button>
  <button onclick="spinOff()">Stop spin</button>{show_lig_buttons}
 </div>
</div>

<script>
{chr(10).join(pdb_blocks)}
const CAV_RESIDS = [{cav_resids_js}];
const LOOP_LO = {loop_lo}, LOOP_HI = {loop_hi};
const LIG_RESN = {f"'{ligand_resname}'" if ligand_resname else "null"};

function _withV(fn) {{ if (!window._v) return; try {{ fn(window._v); }} catch(e) {{ console.error(e); }} }}
function hideSurface() {{ _withV(v => {{ v.removeAllSurfaces(); v.render(); }}); }}
function showSurface() {{ _withV(v => {{
  // Whole protein surface in light grey
  v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.30, color:'lightgrey'}}, {{polymer:true, not:{{resi: CAV_RESIDS}}}}).then(()=>{{
    // Cavity-18 surface patch in wheat
    v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.85, color:'#e8b250'}}, {{resi: CAV_RESIDS, not:{{resi: [`${{LOOP_LO}}-${{LOOP_HI}}`]}}}}).then(()=>{{
      // Allosteric-loop subset in red
      v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.95, color:'#e85050'}}, {{resi: [`${{LOOP_LO}}-${{LOOP_HI}}`], and:{{resi: CAV_RESIDS}}}}).then(()=> v.render());
    }});
  }});
}}); }}
function zoomPocket() {{ _withV(v => {{ v.zoomTo({{resi: CAV_RESIDS}}); v.zoom(0.85); v.render(); }}); }}
function zoomLoop()   {{ _withV(v => {{ v.zoomTo({{resi: [`${{LOOP_LO}}-${{LOOP_HI}}`]}}); v.zoom(0.85); v.render(); }}); }}
function zoomLigand() {{ if (!LIG_RESN) return; _withV(v => {{ v.zoomTo({{resn: LIG_RESN}}); v.zoom(0.8); v.render(); }}); }}
function hideLigand() {{ if (!LIG_RESN) return; _withV(v => {{ v.setStyle({{resn: LIG_RESN}}, {{}}); v.render(); }}); }}
function showLigand() {{ if (!LIG_RESN) return; _withV(v => {{
  v.setStyle({{resn: LIG_RESN}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.18}}}}); v.render();
}}); }}
function spinOn()  {{ _withV(v => v.spin(true)); }}
function spinOff() {{ _withV(v => v.spin(false)); }}

document.addEventListener("DOMContentLoaded", function () {{
  const v = $3Dmol.createViewer("gl", {{backgroundColor: "#0c1116"}});
  v.addModel(pdb_apo, "pdb");
  // Cartoon: chain A grey, chain B lightgrey
  v.setStyle({{chain: 'A'}}, {{cartoon: {{color: '#888888', opacity: 0.7}}}});
  v.setStyle({{chain: 'B'}}, {{cartoon: {{color: '#bbbbbb', opacity: 0.7}}}});
  // Cavity-18 residues: wheat sticks (heavy atoms)
  v.addStyle({{resi: CAV_RESIDS, not: {{atom: 'H'}}}}, {{stick: {{colorscheme: 'wheatCarbon', radius: 0.16}}}});
  // Allosteric-loop overlap: red sticks (override)
  const loopResids = [];
  for (let r = LOOP_LO; r <= LOOP_HI; r++) if (CAV_RESIDS.indexOf(r) !== -1) loopResids.push(r);
  v.addStyle({{resi: loopResids, not: {{atom: 'H'}}}}, {{stick: {{colorscheme: 'redCarbon', radius: 0.20}}}});
  // Ligand model (separate addModel for the complex viewers)
  if (LIG_RESN) {{
    v.addModel(pdb_lig, "pdb");
    v.setStyle({{resn: LIG_RESN}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.18}}}});
    // polar-contact dashes
    v.addLine({{start: {{resn: LIG_RESN}}, end: {{resi: CAV_RESIDS, atom: ['N','O','OD1','OD2','OE1','OE2','OG','OG1','OH','NE','NH1','NH2','ND1','NE2','NZ']}}, dashed: true, color: 'yellow', dashLength: 0.20, gapLength: 0.20}});
  }}
  v.zoomTo({{resi: CAV_RESIDS}});
  v.zoom(0.85);
  // initial: show cavity surface tinted
  v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.30, color:'lightgrey'}}, {{polymer:true, not:{{resi: CAV_RESIDS}}}}).then(()=>{{
    v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.85, color:'#e8b250'}}, {{resi: CAV_RESIDS}}).then(()=>{{
      v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.95, color:'#e85050'}}, {{resi: loopResids}}).then(()=>{{
        v.render(); window._v = v;
      }});
    }});
  }});
}});
</script>
</body></html>
"""
    html_path.write_text(html)


def main():
    print("=== Phase 14d cavity-18 evidence package ===")
    residues = parse_pocket_residues()
    print(f"  cavity-18 residues parsed: {len(residues)} from FPocket pocket18_atm.pdb")
    cav_resids_int = sorted(set(r for _, r in residues))

    # 1. APO + chains
    apo_chained = DL / "cavity18_apo.pdb"
    write_apo_with_chain_ids(apo_chained)
    print(f"  → {apo_chained}")
    # 2. Pocket-only with druggability in B-factor
    pocket_only = DL / "cavity18_apo_pocket.pdb"
    write_pocket_only(apo_chained, residues, pocket_only, drug_score=0.994)
    print(f"  → {pocket_only}")
    # 3. Complexes
    if INDAZOLE_POSE.exists():
        cx1 = DL / "cavity18_indazole_complex.pdb"
        write_complex(apo_chained, INDAZOLE_POSE, cx1, "IND")
        print(f"  → {cx1}")
    if IBUPROFEN_POSE.exists():
        cx2 = DL / "cavity18_ibuprofen_complex.pdb"
        write_complex(apo_chained, IBUPROFEN_POSE, cx2, "IBU")
        print(f"  → {cx2}")
    # 4. Residue + conservation + mutation table
    cons = load_conservation()
    seqs = load_msa()
    print(f"  conservation: {len(cons)} positions; MSA: {len(seqs)} orthologs")
    table_rows = build_residue_table(residues, cons, seqs)
    res_csv = DL / "cavity18_residues.csv"
    write_residue_csv(table_rows, res_csv)
    print(f"  → {res_csv} ({len(table_rows)} rows)")
    # 5. Conservation plot
    cons_png = FIG / "cavity18_conservation.png"
    plot_conservation(table_rows, cons, cons_png)
    print(f"  → {cons_png}")
    # 6. Phylogeny annotated
    phylo_png = FIG / "cavity18_phylogeny_annot.png"
    try:
        plot_phylogeny_annotated(table_rows, phylo_png)
        print(f"  → {phylo_png}")
    except Exception as e:
        print(f"  ! phylogeny annotation failed: {e}")
    # 7. Viewers — three of them
    cav_loop_resids = [r for r in cav_resids_int if ALLOSTERIC_LOOP_RANGE[0] <= r <= ALLOSTERIC_LOOP_RANGE[1]]
    print(f"  cavity-18 ∩ loop 181-197 residues: {cav_loop_resids}")
    # apo viewer
    build_viewer(VW / "cavity18_apo.html",
                 "Cavity 18 — apo (no ligand)",
                 [("apo", apo_chained)], cav_resids_int, None)
    # indazole viewer
    if INDAZOLE_POSE.exists():
        cx1 = DL / "cavity18_indazole_complex.pdb"
        build_viewer(VW / "cavity18_indazole.html",
                     "Cavity 18 + 1H-indazole (−7.52 kcal/mol)",
                     [("apo", apo_chained), ("lig", INDAZOLE_POSE.with_suffix(".pdb"))],
                     cav_resids_int, "UNL")
    # ibuprofen viewer
    if IBUPROFEN_POSE.exists():
        build_viewer(VW / "cavity18_ibuprofen.html",
                     "Cavity 18 + ibuprofen (−7.28 kcal/mol)",
                     [("apo", apo_chained), ("lig", IBUPROFEN_POSE.with_suffix(".pdb"))],
                     cav_resids_int, "UNL")
    print(f"  → {VW}/")
    print("Done.")

if __name__ == "__main__":
    main()
