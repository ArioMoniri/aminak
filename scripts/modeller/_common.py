"""Common utilities for Phase 6 (Modeller homology modelling).

Logging, paths, and small helpers used across step scripts.
"""
from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

PROJECT_DIR = Path(os.environ.get("PROJECT_DIR", str(Path.home() / "conserved_site_project")))
PHASE_DIR = PROJECT_DIR / "10_modeller"
LOGS_DIR = PROJECT_DIR / "logs"
PIPELINE_LOG = PROJECT_DIR / "pipeline.log"

# Step folders
STEP1_DIR = PHASE_DIR / "01_clean_pdb"
STEP2_DIR = PHASE_DIR / "02_blast"
TEMPLATES_DIR = STEP2_DIR / "templates"
STEP3_DIR = PHASE_DIR / "03_alignment"
STEP4_DIR = PHASE_DIR / "04_modeller_run"
MODELS_DIR = STEP4_DIR / "models"
STEP5_DIR = PHASE_DIR / "05_comparison"
STEP6_DIR = PHASE_DIR / "06_validation"
STEP7_DIR = PHASE_DIR / "07_viewers"
VIEWERS_DIR = PROJECT_DIR / "viewers"

PYMOL_BIN = "/opt/homebrew/bin/pymol"
CLUSTALW_BIN = "/opt/homebrew/bin/clustalw"
BLASTP_BIN = "/opt/homebrew/bin/blastp"


def setup_logger(step_name: str) -> logging.Logger:
    """Configure a logger that writes to logs/modeller_<step>.log AND pipeline.log."""
    LOGS_DIR.mkdir(exist_ok=True, parents=True)
    logger = logging.getLogger(f"modeller.{step_name}")
    logger.setLevel(logging.INFO)
    # avoid duplicate handlers on re-import
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    # per-step log
    fh_step = logging.FileHandler(LOGS_DIR / f"modeller_{step_name}.log", mode="a")
    fh_step.setFormatter(fmt)
    logger.addHandler(fh_step)
    # pipeline.log with [MODELLER] prefix
    fmt_pipeline = logging.Formatter(
        f"%(asctime)s [MODELLER:{step_name}] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh_pipeline = logging.FileHandler(PIPELINE_LOG, mode="a")
    fh_pipeline.setFormatter(fmt_pipeline)
    logger.addHandler(fh_pipeline)
    # stdout
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    return logger


def png_ok(path: Path, min_size: int = 50_000) -> bool:
    """Return True if PNG file decodes and is large enough."""
    try:
        from PIL import Image
    except Exception:
        return path.exists() and path.stat().st_size >= min_size
    if not path.exists():
        return False
    if path.stat().st_size < min_size:
        return False
    try:
        with Image.open(path) as im:
            im.verify()
        return True
    except Exception:
        return False
