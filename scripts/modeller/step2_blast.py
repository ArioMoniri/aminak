"""Step 2 — Remote BLAST chain-A against PDB; pick 3 templates with 30-95% identity.

Then download each chosen template's PDB from RCSB and extract the relevant chain.
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
from pathlib import Path

from Bio import PDB, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

sys.path.insert(0, str(Path(__file__).parent))
from _common import STEP1_DIR, STEP2_DIR, TEMPLATES_DIR, setup_logger  # noqa: E402

LOG = setup_logger("step2_blast")

MAX_HITS = 250
MIN_IDENT = 30.0
MAX_IDENT = 95.0
MIN_COVERAGE = 0.80
MAX_RES = 2.5  # Angstroms

# 1HVY itself — exclude (source structure, would be 100% identity)
EXCLUDE_PDB = {"1HVY"}


def parse_pdb_resolution(pdb_path: Path) -> float | None:
    """Parse REMARK 2 RESOLUTION from a PDB file."""
    try:
        with open(pdb_path) as fh:
            for line in fh:
                if line.startswith("REMARK   2 RESOLUTION."):
                    m = re.search(r"([\d\.]+)\s*ANGSTROM", line)
                    if m:
                        try:
                            return float(m.group(1))
                        except ValueError:
                            return None
                if line.startswith("ATOM"):
                    return None  # passed header
    except Exception:
        return None
    return None


def download_pdb(pdb_id: str, dest: Path) -> bool:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    LOG.info("Downloading %s", url)
    try:
        urllib.request.urlretrieve(url, str(dest))
        return dest.exists() and dest.stat().st_size > 1000
    except Exception as e:
        LOG.warning("Download failed for %s: %s", pdb_id, e)
        return False


def extract_chain(src_pdb: Path, chain_id: str, dst_pdb: Path) -> bool:
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(src_pdb.stem, str(src_pdb))
    except Exception as e:
        LOG.warning("PDB parse failed for %s: %s", src_pdb, e)
        return False

    class ChainSelect(PDB.Select):
        def accept_chain(self, chain):
            return chain.id == chain_id

        def accept_residue(self, residue):
            return residue.id[0].strip() == ""

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(str(dst_pdb), ChainSelect())
    return dst_pdb.exists() and dst_pdb.stat().st_size > 500


def run_blast(query_seq: str, xml_out: Path) -> None:
    LOG.info("Running NCBIWWW.qblast on pdb DB (this may take 1-3 minutes)...")
    handle = NCBIWWW.qblast(
        program="blastp",
        database="pdb",
        sequence=query_seq,
        hitlist_size=MAX_HITS,
        expect=10.0,
        word_size=3,
    )
    xml = handle.read()
    xml_out.write_text(xml)
    LOG.info("BLAST result XML written: %s (%d bytes)", xml_out, xml_out.stat().st_size)


def parse_hits(xml_path: Path, query_len: int) -> list[dict]:
    """Return list of hit dicts: pdb_id, chain, identity_pct, coverage, hit_def, hsp_score."""
    hits: list[dict] = []
    with open(xml_path) as fh:
        records = list(NCBIXML.parse(fh))
    if not records:
        return hits
    rec = records[0]
    for aln in rec.alignments:
        # PDB hit_def looks like "pdb|2RD8|A Chain A, ..." for older format,
        # or "Chain A, Thymidylate Synthase [pdb|2RD8|A]" — extract PDB ID and chain.
        hit_def = aln.hit_def
        # Try new fasta-style header in hit_id, e.g. "pdb|2RD8|A"
        m = re.search(r"pdb\|([0-9A-Za-z]{4})\|([A-Za-z0-9])", aln.hit_id)
        if not m:
            m = re.search(r"pdb\|([0-9A-Za-z]{4})\|([A-Za-z0-9])", hit_def)
        if not m:
            # try "1HVY_A" style
            m = re.search(r"\b([0-9][A-Za-z0-9]{3})_([A-Za-z0-9])\b", hit_def)
        if not m:
            continue
        pdb_id = m.group(1).upper()
        chain = m.group(2).upper()
        if pdb_id in EXCLUDE_PDB:
            continue
        # take best HSP
        if not aln.hsps:
            continue
        hsp = aln.hsps[0]
        identity_pct = 100.0 * hsp.identities / hsp.align_length
        coverage = hsp.align_length / query_len
        hits.append({
            "pdb_id": pdb_id,
            "chain": chain,
            "identity_pct": round(identity_pct, 2),
            "coverage": round(coverage, 3),
            "align_length": hsp.align_length,
            "hsp_score": hsp.score,
            "hsp_evalue": hsp.expect,
            "hit_def": hit_def[:160],
        })
    return hits


def main() -> int:
    STEP2_DIR.mkdir(parents=True, exist_ok=True)
    TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)

    fasta = STEP1_DIR / "1hvy_chainA.fasta"
    if not fasta.exists():
        LOG.error("Missing %s", fasta)
        return 2
    rec = next(SeqIO.parse(str(fasta), "fasta"))
    query_seq = str(rec.seq)
    LOG.info("Query length: %d", len(query_seq))

    xml_path = STEP2_DIR / "blast.xml"
    if not xml_path.exists() or xml_path.stat().st_size < 200:
        run_blast(query_seq, xml_path)
    else:
        LOG.info("Reusing cached BLAST XML: %s", xml_path)

    hits = parse_hits(xml_path, len(query_seq))
    LOG.info("Parsed %d BLAST hits", len(hits))

    # Write tsv of all hits
    tsv = STEP2_DIR / "blast.tsv"
    with open(tsv, "w") as fh:
        fh.write("pdb_id\tchain\tidentity_pct\tcoverage\talign_length\tscore\tevalue\thit_def\n")
        for h in hits:
            fh.write("\t".join(str(h[k]) for k in [
                "pdb_id", "chain", "identity_pct", "coverage",
                "align_length", "hsp_score", "hsp_evalue", "hit_def"]) + "\n")
    LOG.info("Wrote %s", tsv)

    # Filter by identity, coverage; resolution will be checked after download
    passed: list[dict] = []
    seen_pdb = set()
    for h in hits:
        if h["pdb_id"] in seen_pdb:
            continue
        if not (MIN_IDENT <= h["identity_pct"] <= MAX_IDENT):
            continue
        if h["coverage"] < MIN_COVERAGE:
            continue
        seen_pdb.add(h["pdb_id"])
        passed.append(h)
    LOG.info("After identity/coverage filter: %d candidates", len(passed))

    # Download candidates and check resolution; pick top 3 by score
    selected: list[dict] = []
    for h in passed:
        if len(selected) >= 3:
            break
        pdb_path = TEMPLATES_DIR / f"{h['pdb_id']}.pdb"
        if not pdb_path.exists() or pdb_path.stat().st_size < 1000:
            ok = download_pdb(h["pdb_id"], pdb_path)
            if not ok:
                LOG.warning("Skip %s: download failed", h["pdb_id"])
                continue
        res = parse_pdb_resolution(pdb_path)
        h["resolution"] = res
        if res is not None and res > MAX_RES:
            LOG.info("Skip %s: resolution %.2f > %.2f", h["pdb_id"], res, MAX_RES)
            continue
        # Extract chain
        chain_pdb = TEMPLATES_DIR / f"{h['pdb_id']}_{h['chain']}.pdb"
        if not extract_chain(pdb_path, h["chain"], chain_pdb):
            LOG.warning("Skip %s: chain %s extraction failed", h["pdb_id"], h["chain"])
            continue
        h["chain_pdb"] = str(chain_pdb)
        selected.append(h)
        LOG.info("Selected %s_%s identity=%.2f%% cov=%.2f res=%s",
                 h["pdb_id"], h["chain"], h["identity_pct"], h["coverage"], res)
        time.sleep(0.5)  # be polite to RCSB

    if len(selected) < 3:
        LOG.warning("Only %d templates selected (wanted 3)", len(selected))
        if len(selected) == 0:
            LOG.error("No usable templates — abort")
            return 4

    # Write selected_templates.json
    out = {
        "criteria": {
            "min_identity_pct": MIN_IDENT,
            "max_identity_pct": MAX_IDENT,
            "min_coverage": MIN_COVERAGE,
            "max_resolution_A": MAX_RES,
            "excluded_pdbs": sorted(EXCLUDE_PDB),
        },
        "note": (
            "Templates were chosen with %identity 30-95% deliberately, for "
            "educational purposes; in production one would use the highest-quality "
            "template available, including 100% identity matches (e.g. 1HVY itself)."
        ),
        "selected": selected,
        "n_total_hits": len(hits),
        "n_passed_filters": len(passed),
    }
    sel_path = STEP2_DIR / "selected_templates.json"
    sel_path.write_text(json.dumps(out, indent=2))
    LOG.info("Wrote %s", sel_path)
    LOG.info("STEP 2 OK: %d templates selected", len(selected))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
