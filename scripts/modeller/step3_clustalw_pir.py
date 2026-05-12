"""Step 3 — ClustalW MSA of target + templates, then convert to Modeller PIR (.ali)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from Bio import AlignIO, PDB, SeqIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    CLUSTALW_BIN, STEP1_DIR, STEP2_DIR, STEP3_DIR, TEMPLATES_DIR, setup_logger,
)  # noqa: E402

LOG = setup_logger("step3_clustalw_pir")


def chain_seq_and_range(pdb_path: Path, chain_id: str) -> tuple[str, int, int]:
    """Return (seq, first_resnum, last_resnum) for the chain."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    chain = None
    for model in structure:
        for ch in model:
            if ch.id == chain_id:
                chain = ch
                break
        if chain is not None:
            break
    if chain is None:
        raise ValueError(f"Chain {chain_id} not found in {pdb_path}")
    ppb = PPBuilder()
    seq_parts = []
    res_nums = []
    for pp in ppb.build_peptides(chain):
        seq_parts.append(str(pp.get_sequence()))
        res_nums.extend([r.id[1] for r in pp])
    if not seq_parts:
        raise ValueError(f"No peptides built from {pdb_path}:{chain_id}")
    seq = "".join(seq_parts)
    return seq, min(res_nums), max(res_nums)


def main() -> int:
    STEP3_DIR.mkdir(parents=True, exist_ok=True)
    sel_path = STEP2_DIR / "selected_templates.json"
    sel = json.loads(sel_path.read_text())
    templates = sel["selected"]
    LOG.info("Templates: %s", [(t["pdb_id"], t["chain"]) for t in templates])

    # Target sequence
    target_fa = STEP1_DIR / "1hvy_chainA.fasta"
    target_rec = next(SeqIO.parse(str(target_fa), "fasta"))
    target_seq = str(target_rec.seq)
    # Target residue range — use 1HVY chain A actual residue numbering
    target_pdb = STEP1_DIR / "1hvy_chainA.pdb"
    _, t_first, t_last = chain_seq_and_range(target_pdb, "A")
    LOG.info("Target residue range %d-%d (len=%d)", t_first, t_last, len(target_seq))

    # Build combined input FASTA
    records: list[SeqRecord] = []
    records.append(SeqRecord(Seq(target_seq), id="target", description="P04818 TYMS chain A"))
    template_meta: list[dict] = []
    for t in templates:
        pdb_id = t["pdb_id"]
        chain = t["chain"]
        chain_pdb = TEMPLATES_DIR / f"{pdb_id}_{chain}.pdb"
        seq, first, last = chain_seq_and_range(chain_pdb, chain)
        rec_id = f"{pdb_id}"
        records.append(SeqRecord(Seq(seq), id=rec_id, description=f"{pdb_id} chain {chain}"))
        template_meta.append({
            "pdb_id": pdb_id,
            "chain": chain,
            "first": first,
            "last": last,
            "len": len(seq),
        })
        LOG.info("Template %s_%s: range %d-%d, len=%d", pdb_id, chain, first, last, len(seq))

    input_fa = STEP3_DIR / "input.fasta"
    SeqIO.write(records, str(input_fa), "fasta")
    LOG.info("Wrote %s (%d records)", input_fa, len(records))

    # Run ClustalW
    aln_out = STEP3_DIR / "aligned.aln"
    if aln_out.exists():
        aln_out.unlink()
    cmd = [
        CLUSTALW_BIN,
        f"-INFILE={input_fa}",
        "-ALIGN",
        f"-OUTFILE={aln_out}",
        "-OUTPUT=CLUSTAL",
    ]
    LOG.info("Running ClustalW: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(STEP3_DIR))
    if proc.returncode != 0:
        LOG.error("ClustalW failed: %s", proc.stderr)
        return 3
    LOG.info("ClustalW done; written %s (%d bytes)", aln_out, aln_out.stat().st_size)

    aln = AlignIO.read(str(aln_out), "clustal")
    # Map record id -> aligned sequence (with gaps)
    aligned: dict[str, str] = {r.id: str(r.seq) for r in aln}
    LOG.info("Aligned columns: %d", aln.get_alignment_length())

    # Build PIR
    pir_lines: list[str] = []
    # Target
    target_aligned = aligned.get("target")
    if target_aligned is None:
        LOG.error("target sequence missing from alignment")
        return 4
    pir_lines.append(">P1;target")
    pir_lines.append(f"sequence:target:{t_first}:A:{t_last}:A:::-1.00:-1.00")
    # PIR requires '*' terminator and lines wrapped at 75 chars but Modeller is lenient
    pir_lines.append(target_aligned + "*")

    for tm in template_meta:
        pdb_id = tm["pdb_id"]
        chain = tm["chain"]
        first = tm["first"]
        last = tm["last"]
        seq = aligned.get(pdb_id)
        if seq is None:
            LOG.error("template %s missing", pdb_id)
            return 5
        pir_lines.append(f">P1;{pdb_id}_{chain}")
        pir_lines.append(
            f"structureX:{pdb_id}_{chain}:{first}:{chain}:{last}:{chain}:::-1.00:-1.00"
        )
        pir_lines.append(seq + "*")

    pir_text = "\n".join(pir_lines) + "\n"
    ali_path = STEP3_DIR / "alignment.ali"
    pir_path = STEP3_DIR / "alignment.pir"
    ali_path.write_text(pir_text)
    pir_path.write_text(pir_text)
    LOG.info("Wrote %s and %s", ali_path, pir_path)

    # Persist knowns list for Step 4
    knowns = [f"{tm['pdb_id']}_{tm['chain']}" for tm in template_meta]
    (STEP3_DIR / "knowns.json").write_text(json.dumps({
        "knowns": knowns,
        "sequence_code": "target",
        "target_first": t_first,
        "target_last": t_last,
    }, indent=2))
    LOG.info("knowns: %s", knowns)
    LOG.info("STEP 3 OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
