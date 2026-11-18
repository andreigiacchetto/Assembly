#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Automated Bacterial Genome Assembly & Annotation Pipeline (Offline / Conda)
#                    -------------------------------------------------
#                    |                                               |
#                    |         ANDREI GIACCHETTO FELICE              |
#                    |                                               |
#                    |                PhD Student                    |
#                    |  Laboratory of Immunology and Bioinformatics  |
#                    |    Federal University of Triangulo Mineiro    |
#                    -------------------------------------------------
#
#                               VERSION: 2.0
#                               LAST UPDATE: 18/11/2026
#
#                    -------------------------------------------------

import argparse
import glob
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# =============================================================================
# Configuration & Tool Registry
# =============================================================================

ALL_STEPS = [
    "fastqc_raw",
    "trimmomatic",
    "fastqc_trimmed",
    "unicycler",
    "ragtag",
    "quast",
    "prokka",
    "multiqc",
]

IGNORED_DIRS = {
    "fastqc_raw",
    "fastqc_trimmed",
    "trimmomatic_result",
    "unicycler_result",
    "ragtag_result",
    "quast_result",
    "unicycler_quast_result",
    "prokka_result",
    "unicycler_prokka_result",
    "multiqc_report",
}


@dataclass
class ToolSpec:
    name: str
    executables: List[str]
    conda_pkg: str
    description: str
    fallback_paths: List[str] = field(default_factory=list)
    optional: bool = False


TOOLS: Dict[str, ToolSpec] = {
    "fastqc": ToolSpec(
        name="FastQC",
        executables=["fastqc"],
        conda_pkg="fastqc",
        description="Quality control for raw and trimmed FASTQ reads",
        fallback_paths=["/home/aluno/Programas/FastQC/fastqc"],
    ),
    "trimmomatic": ToolSpec(
        name="Trimmomatic",
        executables=["trimmomatic"],
        conda_pkg="trimmomatic",
        description="Illumina adapter removal and Phred quality trimming",
    ),
    "unicycler": ToolSpec(
        name="Unicycler",
        executables=["unicycler"],
        conda_pkg="unicycler",
        description="De novo bacterial genome assembler",
    ),
    "ragtag": ToolSpec(
        name="RagTag",
        executables=["ragtag.py", "ragtag"],
        conda_pkg="ragtag",
        description="Reference-guided contig scaffolding (optional)",
        optional=True,
    ),
    "quast": ToolSpec(
        name="QUAST",
        executables=["quast.py", "quast"],
        conda_pkg="quast",
        description="Genome assembly quality evaluation",
    ),
    "prokka": ToolSpec(
        name="Prokka",
        executables=["prokka"],
        conda_pkg="prokka",
        description="Rapid prokaryotic genome annotation",
    ),
    "multiqc": ToolSpec(
        name="MultiQC",
        executables=["multiqc"],
        conda_pkg="multiqc",
        description="Aggregate QC/assembly/annotation reports into one HTML (optional)",
        optional=True,
    ),
}


@dataclass
class SamplePair:
    sample_id: str
    r1: Path
    r2: Path
    output_root: Path
    already_trimmed: bool = False


@dataclass
class SampleResult:
    sample_id: str
    directory: Path
    status: str  # "OK", "FAILED", "SKIPPED"
    duration_sec: float
    details: str = ""


# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(log_file: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger("run_assembly")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    console_level = logging.DEBUG if verbose else logging.INFO
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(console_level)
    ch.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(ch)

    if log_file is not None:
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(fh)

    return logger


# =============================================================================
# Dependency & Adapter Discovery (Offline-First)
# =============================================================================

def find_executable(spec: ToolSpec) -> Optional[str]:
    for exe in spec.executables:
        found = shutil.which(exe)
        if found:
            return found
    for fb in spec.fallback_paths:
        if os.path.isfile(fb) and os.access(fb, os.X_OK):
            return fb
    return None


def find_trimmomatic_adapters(custom_path: Optional[str] = None) -> Optional[Path]:
    """Locate Illumina adapter FASTA file offline from custom path or Conda/system dirs."""
    if custom_path:
        p = Path(custom_path).expanduser().resolve()
        return p if p.is_file() else None

    candidate_names = ["TruSeq3-PE-2.fa", "TruSeq3-PE.fa", "NexteraPE-PE.fa"]
    search_dirs: List[Path] = []

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        for pattern in [
            Path(conda_prefix) / "share" / "trimmomatic*" / "adapters",
            Path(conda_prefix) / "share" / "trimmomatic" / "adapters",
        ]:
            for match in glob.glob(str(pattern)):
                search_dirs.append(Path(match))

    trimmomatic_bin = shutil.which("trimmomatic")
    if trimmomatic_bin:
        bin_dir = Path(trimmomatic_bin).resolve().parent
        for rel in ["../share/trimmomatic/adapters", "../share/trimmomatic*/adapters"]:
            for match in glob.glob(str(bin_dir / rel)):
                search_dirs.append(Path(match))

    search_dirs.extend([
        Path("/usr/share/trimmomatic"),
        Path("/usr/local/share/trimmomatic/adapters"),
    ])

    for d in search_dirs:
        if not d.is_dir():
            continue
        for name in candidate_names:
            candidate = d / name
            if candidate.is_file():
                return candidate
    return None


def get_conda_installer() -> Optional[str]:
    for mgr in ("mamba", "micromamba", "conda"):
        if shutil.which(mgr):
            return mgr
    return None


def check_dependencies(
    active_steps: List[str],
    has_reference: bool,
    auto_install: bool,
    logger: logging.Logger,
    check_only: bool = False,
) -> Tuple[Dict[str, str], bool]:
    """
    Check required and optional tools.
    Returns (resolved_paths, all_required_ok).
    """
    step_to_tool = {
        "fastqc_raw": "fastqc",
        "fastqc_trimmed": "fastqc",
        "trimmomatic": "trimmomatic",
        "unicycler": "unicycler",
        "ragtag": "ragtag",
        "quast": "quast",
        "prokka": "prokka",
        "multiqc": "multiqc",
    }

    needed_tools = set()
    for step in active_steps:
        if step == "ragtag" and not has_reference:
            continue
        tool_key = step_to_tool.get(step)
        if tool_key:
            needed_tools.add(tool_key)

    if check_only:
        needed_tools = set(TOOLS.keys())

    resolved: Dict[str, str] = {}
    missing_required: List[ToolSpec] = []
    missing_optional: List[ToolSpec] = []

    logger.info("\n=== Dependency Check ===")
    for key, spec in TOOLS.items():
        if key not in needed_tools:
            continue
        path = find_executable(spec)
        if path:
            resolved[key] = path
            logger.info(f"  [OK]      {spec.name:<12} -> {path}")
        else:
            is_opt = spec.optional and not (key == "ragtag" and has_reference and "ragtag" in active_steps)
            tag = "[OPT-MISS]" if is_opt else "[MISSING] "
            logger.info(f"  {tag} {spec.name:<12} -> not found in PATH")
            if is_opt:
                missing_optional.append(spec)
            else:
                missing_required.append(spec)

    adapters = find_trimmomatic_adapters()
    if "trimmomatic" in needed_tools:
        if adapters:
            logger.info(f"  [OK]      {'Adapters':<12} -> {adapters}")
        else:
            logger.info(f"  [WARN]    {'Adapters':<12} -> TruSeq3-PE.fa not auto-detected (use --adapters)")

    if missing_optional:
        opt_pkgs = " ".join(s.conda_pkg for s in missing_optional)
        logger.info(f"\n[INFO] Optional tools not installed ({', '.join(s.name for s in missing_optional)}).")
        logger.info(f"       To install: conda install -c bioconda -c conda-forge {opt_pkgs}")

    if not missing_required:
        return resolved, True

    req_pkgs = " ".join(s.conda_pkg for s in missing_required)
    installer = get_conda_installer()
    install_cmd = (
        f"{installer or 'conda'} install -y -c bioconda -c conda-forge {req_pkgs}"
    )

    logger.error("\n[ERROR] Missing required tools:")
    for spec in missing_required:
        logger.error(f"  - {spec.name} ({spec.description}) -> conda pkg: {spec.conda_pkg}")

    logger.error("\n[HINT] Run the following command to install them:")
    logger.error(f"  {install_cmd}")

    if auto_install and installer:
        logger.info(f"\n[INFO] --auto-install enabled. Running: {install_cmd}")
        try:
            subprocess.run(shlex.split(install_cmd), check=True)
            # Re-verify after installation
            still_missing = []
            for spec in missing_required:
                key = next(k for k, v in TOOLS.items() if v is spec)
                path = find_executable(spec)
                if path:
                    resolved[key] = path
                else:
                    still_missing.append(spec)
            if not still_missing:
                logger.info("[OK] All missing required tools were installed successfully.")
                return resolved, True
        except subprocess.CalledProcessError as e:
            logger.error(f"[ERROR] Auto-install failed: {e}")
    elif auto_install and not installer:
        logger.error("[ERROR] --auto-install requested, but neither conda, mamba, nor micromamba was found.")

    return resolved, False


# =============================================================================
# Sample Discovery
# =============================================================================

RAW_R1_PATTERNS = [
    re.compile(r"^(.+?)(?:_L\d{3})?_R1(?:_001)?\.(fastq|fq)\.gz$"),
    re.compile(r"^(.+?)_1\.(fastq|fq)\.gz$"),
]

TRIMMED_R1_PATTERN = re.compile(r"^(.+?)_1P\.(fastq|fq)\.gz$")


def find_sample_pairs_in_dir(directory: Path) -> List[SamplePair]:
    """Detect paired-end FASTQ files inside a single directory."""
    if not directory.is_dir():
        return []

    files = sorted(
        p for p in directory.iterdir()
        if p.is_file() and (p.name.endswith(".fastq.gz") or p.name.endswith(".fq.gz"))
    )
    if not files:
        return []

    file_names = {p.name: p for p in files}
    pairs: List[SamplePair] = []
    used_files = set()

    # 1. Check raw paired-end patterns first
    for p in files:
        if p.name in used_files:
            continue
        if any(tag in p.name for tag in ("_1P.", "_2P.", "_1U.", "_2U.")):
            continue

        for idx, pat in enumerate(RAW_R1_PATTERNS):
            m = pat.match(p.name)
            if not m:
                continue
            sample_id = m.group(1)
            if idx == 0:
                r2_name = re.sub(r"_R1(_001)?\.(fastq|fq)\.gz$", r"_R2\1.\2.gz", p.name)
            else:
                r2_name = re.sub(r"_1\.(fastq|fq)\.gz$", r"_2.\1.gz", p.name)

            if r2_name in file_names:
                r2_path = file_names[r2_name]
                used_files.add(p.name)
                used_files.add(r2_name)
                pairs.append(
                    SamplePair(
                        sample_id=sample_id,
                        r1=p,
                        r2=r2_path,
                        output_root=directory,
                        already_trimmed=False,
                    )
                )
            break

    # 2. Fallback: if no raw pairs found, check already-trimmed *_1P / *_2P files
    if not pairs:
        for p in files:
            m = TRIMMED_R1_PATTERN.match(p.name)
            if not m:
                continue
            sample_id = m.group(1)
            r2_name = re.sub(r"_1P\.(fastq|fq)\.gz$", r"_2P.\1.gz", p.name)
            if r2_name in file_names:
                pairs.append(
                    SamplePair(
                        sample_id=sample_id,
                        r1=p,
                        r2=file_names[r2_name],
                        output_root=directory,
                        already_trimmed=True,
                    )
                )

    # If multiple samples live in the same directory, isolate each sample's outputs
    if len(pairs) > 1:
        for pair in pairs:
            pair.output_root = directory / f"{pair.sample_id}_results"

    return pairs


def discover_samples(input_dir: Path, recursive: bool = True) -> List[SamplePair]:
    """Discover all paired-end samples under input_dir."""
    all_pairs: List[SamplePair] = []

    if not recursive:
        return find_sample_pairs_in_dir(input_dir)

    for root, dirs, _ in os.walk(input_dir):
        # Prune tool output directories so we never re-scan trimmed/assembled outputs
        dirs[:] = sorted(
            d for d in dirs
            if d not in IGNORED_DIRS and not d.endswith("_results") and not d.startswith(".")
        )
        root_path = Path(root)
        all_pairs.extend(find_sample_pairs_in_dir(root_path))

    return all_pairs


# =============================================================================
# Command Execution Helper
# =============================================================================

def exec_cmd(
    cmd: List[str],
    cwd: Path,
    dry_run: bool,
    logger: logging.Logger,
) -> bool:
    cmd_str = " ".join(shlex.quote(str(c)) for c in cmd)
    logger.info(f"  [CMD] {cmd_str}")
    if dry_run:
        return True
    try:
        subprocess.run([str(c) for c in cmd], cwd=str(cwd), check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"  [ERROR] Command failed with exit code {e.returncode}: {cmd_str}")
        return False


# =============================================================================
# Sample Pipeline
# =============================================================================

def process_sample(
    pair: SamplePair,
    active_steps: List[str],
    tools: Dict[str, str],
    args: argparse.Namespace,
    adapters_path: Optional[Path],
    reference_path: Optional[Path],
    logger: logging.Logger,
) -> SampleResult:
    start_time = time.time()
    out_root = pair.output_root
    if not args.dry_run:
        out_root.mkdir(parents=True, exist_ok=True)

    logger.info(f"\n============================================================")
    logger.info(f"[SAMPLE] {pair.sample_id}")
    logger.info(f"  R1:     {pair.r1}")
    logger.info(f"  R2:     {pair.r2}")
    logger.info(f"  Output: {out_root}")
    logger.info(f"============================================================")

    threads = str(args.threads)
    force = args.force

    # Paths for intermediate/final outputs
    fastqc_raw_dir = out_root / "fastqc_raw"
    trim_dir = out_root / "trimmomatic_result"
    fastqc_trim_dir = out_root / "fastqc_trimmed"
    unicycler_dir = out_root / "unicycler_result"
    ragtag_dir = out_root / "ragtag_result"
    quast_dir = out_root / "quast_result"
    prokka_dir = out_root / "prokka_result"

    p1 = trim_dir / f"{pair.sample_id}_1P.fastq.gz"
    u1 = trim_dir / f"{pair.sample_id}_1U.fastq.gz"
    p2 = trim_dir / f"{pair.sample_id}_2P.fastq.gz"
    u2 = trim_dir / f"{pair.sample_id}_2U.fastq.gz"

    if pair.already_trimmed:
        p1, p2 = pair.r1, pair.r2

    # 1. FastQC on raw reads
    if "fastqc_raw" in active_steps and not pair.already_trimmed:
        r1_html = fastqc_raw_dir / (pair.r1.name.replace(".fastq.gz", "_fastqc.html").replace(".fq.gz", "_fastqc.html"))
        if r1_html.exists() and not force:
            logger.info("  [SKIP] FastQC (raw) already completed.")
        else:
            logger.info("  [STEP] Running FastQC on raw reads...")
            if not args.dry_run:
                fastqc_raw_dir.mkdir(parents=True, exist_ok=True)
            ok = exec_cmd(
                [tools.get("fastqc", "fastqc"), "-t", threads, "-o", str(fastqc_raw_dir), str(pair.r1), str(pair.r2)],
                cwd=out_root,
                dry_run=args.dry_run,
                logger=logger,
            )
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "fastqc_raw")

    # 2. Trimmomatic (Phred >= 20 + Illumina adapter clipping)
    if "trimmomatic" in active_steps and not pair.already_trimmed:
        if p1.exists() and p2.exists() and not force:
            logger.info("  [SKIP] Trimmomatic outputs already exist.")
        else:
            logger.info(f"  [STEP] Running Trimmomatic (Phred >= {args.phred}, MinLen >= {args.minlen})...")
            if not args.dry_run:
                trim_dir.mkdir(parents=True, exist_ok=True)
            trim_cmd = [
                tools.get("trimmomatic", "trimmomatic"),
                "PE",
                "-threads", threads,
                "-phred33",
                str(pair.r1),
                str(pair.r2),
                str(p1),
                str(u1),
                str(p2),
                str(u2),
            ]
            if adapters_path:
                trim_cmd.append(f"ILLUMINACLIP:{adapters_path}:2:30:10:2:True")
            else:
                logger.warning("  [WARN] No adapter FASTA found; running quality trimming without ILLUMINACLIP.")
            trim_cmd.extend([
                f"LEADING:{args.phred}",
                f"TRAILING:{args.phred}",
                f"SLIDINGWINDOW:4:{args.phred}",
                f"MINLEN:{args.minlen}",
            ])
            ok = exec_cmd(trim_cmd, cwd=out_root, dry_run=args.dry_run, logger=logger)
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "trimmomatic")
    elif not pair.already_trimmed and not (p1.exists() and p2.exists()):
        # Trimmomatic was skipped by user; use raw reads for downstream steps
        p1, p2 = pair.r1, pair.r2

    # 3. FastQC on trimmed reads
    if "fastqc_trimmed" in active_steps:
        if not args.dry_run and not (p1.exists() and p2.exists()):
            logger.error("  [ERROR] Paired reads for FastQC (trimmed) not found.")
            return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "missing trimmed reads")
        p1_html = fastqc_trim_dir / (p1.name.replace(".fastq.gz", "_fastqc.html").replace(".fq.gz", "_fastqc.html"))
        if p1_html.exists() and not force:
            logger.info("  [SKIP] FastQC (trimmed) already completed.")
        else:
            logger.info("  [STEP] Running FastQC on trimmed reads...")
            if not args.dry_run:
                fastqc_trim_dir.mkdir(parents=True, exist_ok=True)
            ok = exec_cmd(
                [tools.get("fastqc", "fastqc"), "-t", threads, "-o", str(fastqc_trim_dir), str(p1), str(p2)],
                cwd=out_root,
                dry_run=args.dry_run,
                logger=logger,
            )
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "fastqc_trimmed")

    # 4. Unicycler (De novo assembly)
    assembly_fasta = unicycler_dir / "assembly.fasta"
    if "unicycler" in active_steps:
        if assembly_fasta.exists() and not force:
            logger.info("  [SKIP] Unicycler assembly already exists.")
        else:
            logger.info("  [STEP] Running Unicycler...")
            ok = exec_cmd(
                [tools.get("unicycler", "unicycler"), "-1", str(p1), "-2", str(p2), "-o", str(unicycler_dir), "-t", threads],
                cwd=out_root,
                dry_run=args.dry_run,
                logger=logger,
            )
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "unicycler")

    # 5. RagTag (Optional reference-guided scaffolding)
    ragtag_fasta = ragtag_dir / "ragtag.scaffold.fasta"
    if "ragtag" in active_steps and reference_path is not None:
        if not args.dry_run and not assembly_fasta.exists():
            logger.error("  [ERROR] Unicycler assembly.fasta not found for RagTag.")
            return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "missing assembly.fasta")
        if ragtag_fasta.exists() and not force:
            logger.info("  [SKIP] RagTag scaffold already exists.")
        else:
            logger.info("  [STEP] Running RagTag scaffolding...")
            ok = exec_cmd(
                [
                    tools.get("ragtag", "ragtag.py"),
                    "scaffold",
                    str(reference_path),
                    str(assembly_fasta),
                    "-o", str(ragtag_dir),
                    "-t", threads,
                ],
                cwd=out_root,
                dry_run=args.dry_run,
                logger=logger,
            )
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "ragtag")

    # Choose target FASTA for QUAST & Prokka: RagTag scaffold if available, else Unicycler assembly
    target_fasta = ragtag_fasta if (reference_path is not None and "ragtag" in active_steps and (args.dry_run or ragtag_fasta.exists())) else assembly_fasta

    # 6. QUAST (Assembly QC)
    quast_report = quast_dir / "report.txt"
    if "quast" in active_steps:
        if not args.dry_run and not target_fasta.exists():
            logger.error(f"  [ERROR] Target FASTA not found for QUAST: {target_fasta}")
            return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "missing target fasta")
        if quast_report.exists() and not force:
            logger.info("  [SKIP] QUAST report already exists.")
        else:
            logger.info(f"  [STEP] Running QUAST on {target_fasta.name}...")
            quast_cmd = [tools.get("quast", "quast.py"), "-o", str(quast_dir), "-t", threads]
            if reference_path is not None:
                quast_cmd.extend(["-r", str(reference_path)])
            quast_cmd.append(str(target_fasta))
            ok = exec_cmd(quast_cmd, cwd=out_root, dry_run=args.dry_run, logger=logger)
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "quast")

    # 7. Prokka (Genome Annotation)
    prokka_gbk = prokka_dir / f"{pair.sample_id}.gbk"
    if "prokka" in active_steps:
        if not args.dry_run and not target_fasta.exists():
            logger.error(f"  [ERROR] Target FASTA not found for Prokka: {target_fasta}")
            return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "missing target fasta")
        if prokka_gbk.exists() and not force:
            logger.info("  [SKIP] Prokka annotation already exists.")
        else:
            logger.info(f"  [STEP] Running Prokka on {target_fasta.name}...")
            prokka_cmd = [
                tools.get("prokka", "prokka"),
                "--outdir", str(prokka_dir),
                "--prefix", pair.sample_id,
                "--cpus", threads,
                "--addgenes",
                "--force",
            ]
            if args.prokka_rnammer:
                prokka_cmd.append("--rnammer")
            prokka_cmd.append(str(target_fasta))
            ok = exec_cmd(prokka_cmd, cwd=out_root, dry_run=args.dry_run, logger=logger)
            if not ok:
                return SampleResult(pair.sample_id, out_root, "FAILED", time.time() - start_time, "prokka")

    return SampleResult(pair.sample_id, out_root, "OK", time.time() - start_time)


# =============================================================================
# CLI & Main Entrypoint
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    default_threads = min(8, os.cpu_count() or 4)
    parser = argparse.ArgumentParser(
        prog="run_assembly.py",
        description="Automated offline bacterial genome assembly & annotation pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Check installed tools & Conda packages
  python run_assembly.py --check-deps

  # Full automated pipeline (FastQC -> Trimmomatic Phred 20 -> Unicycler -> QUAST -> Prokka -> MultiQC)
  python run_assembly.py -i ./samples -t 8

  # Include RagTag reference-guided scaffolding
  python run_assembly.py -i ./samples -r reference.fasta -t 8

  # Dry-run to preview commands without executing
  python run_assembly.py -i ./samples --dry-run

  # Run only specific steps
  python run_assembly.py -i ./samples --steps trimmomatic,unicycler
        """,
    )

    parser.add_argument("-i", "--input", type=str, default=".", help="Input directory containing FASTQ files or sample subdirectories (default: .)")
    parser.add_argument("-r", "--reference", type=str, default=None, help="Reference genome FASTA (enables RagTag scaffolding & QUAST reference comparison)")
    parser.add_argument("-t", "--threads", type=int, default=default_threads, help=f"Number of CPU threads (default: {default_threads})")
    parser.add_argument("--no-recursive", action="store_true", help="Only scan the top-level input directory")

    # Trimmomatic tuning
    parser.add_argument("--phred", type=int, default=20, help="Minimum Phred quality score for Trimmomatic LEADING/TRAILING/SLIDINGWINDOW (default: 20)")
    parser.add_argument("--minlen", type=int, default=30, help="Minimum read length after trimming (default: 30)")
    parser.add_argument("--adapters", type=str, default=None, help="Custom path to Illumina adapters FASTA (e.g. TruSeq3-PE.fa)")

    # Step control
    parser.add_argument("--steps", type=str, default=None, help=f"Comma-separated steps to run. Available: {','.join(ALL_STEPS)}")
    parser.add_argument("--skip", type=str, default=None, help="Comma-separated steps to skip")
    parser.add_argument("--only-fastqc", action="store_true", help="Run only FastQC steps")
    parser.add_argument("--only-trimmomatic", action="store_true", help="Run only Trimmomatic")
    parser.add_argument("--only-unicycler", action="store_true", help="Run only Unicycler")
    parser.add_argument("--only-quast", action="store_true", help="Run only QUAST")
    parser.add_argument("--only-prokka", action="store_true", help="Run only Prokka")
    parser.add_argument("--resume", action="store_true", help="Resume pipeline, skipping already-completed steps (default behavior unless --force)")
    parser.add_argument("--force", action="store_true", help="Re-run steps even if output files already exist")

    # Prokka / Environment / Utility
    parser.add_argument("--prokka-rnammer", action="store_true", help="Pass --rnammer to Prokka (requires RNAmmer installed)")
    parser.add_argument("--check-deps", action="store_true", help="Check if required bioinformatics tools are installed and exit")
    parser.add_argument("--auto-install", action="store_true", help="Attempt to install missing tools via mamba/conda (requires internet)")
    parser.add_argument("--dry-run", action="store_true", help="Print commands that would be executed without running them")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")

    return parser


def resolve_active_steps(args: argparse.Namespace) -> List[str]:
    if args.steps:
        requested = [s.strip() for s in args.steps.split(",") if s.strip()]
        invalid = [s for s in requested if s not in ALL_STEPS]
        if invalid:
            sys.exit(f"[ERROR] Unknown step(s): {', '.join(invalid)}. Valid: {', '.join(ALL_STEPS)}")
        steps = requested
    else:
        only_map = {
            "only_fastqc": ["fastqc_raw", "fastqc_trimmed"],
            "only_trimmomatic": ["trimmomatic"],
            "only_unicycler": ["unicycler"],
            "only_quast": ["quast"],
            "only_prokka": ["prokka"],
        }
        selected = []
        for attr, mapped in only_map.items():
            if getattr(args, attr):
                selected.extend(mapped)
        steps = selected if selected else list(ALL_STEPS)

    if args.skip:
        to_skip = {s.strip() for s in args.skip.split(",") if s.strip()}
        steps = [s for s in steps if s not in to_skip]

    return steps


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    input_dir = Path(args.input).expanduser().resolve()
    log_path = None if (args.check_deps or args.dry_run) else input_dir / f"run_assembly_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    if log_path is not None and not input_dir.exists():
        log_path = None

    logger = setup_logging(log_file=log_path, verbose=args.verbose)
    active_steps = resolve_active_steps(args)

    reference_path: Optional[Path] = None
    if args.reference:
        reference_path = Path(args.reference).expanduser().resolve()
        if not reference_path.is_file():
            logger.error(f"[ERROR] Reference genome file not found: {reference_path}")
            sys.exit(1)

    resolved_tools, deps_ok = check_dependencies(
        active_steps=active_steps,
        has_reference=reference_path is not None,
        auto_install=args.auto_install,
        logger=logger,
        check_only=args.check_deps,
    )

    if args.check_deps:
        sys.exit(0 if deps_ok else 1)

    if not deps_ok and not args.dry_run:
        logger.error("\n[ABORT] Install the missing tools above (or pass --dry-run to preview commands).")
        sys.exit(1)

    if not input_dir.is_dir():
        logger.error(f"[ERROR] Input directory does not exist: {input_dir}")
        sys.exit(1)

    adapters_path = find_trimmomatic_adapters(args.adapters)
    if args.adapters and not adapters_path:
        logger.error(f"[ERROR] Custom adapters file not found: {args.adapters}")
        sys.exit(1)

    logger.info(f"\n[INFO] Scanning for paired-end FASTQ samples in: {input_dir}")
    samples = discover_samples(input_dir, recursive=not args.no_recursive)

    if not samples:
        logger.warning("[WARN] No paired-end FASTQ samples (*_R1/*_R2 or *_1/*_2) found.")
        sys.exit(0)

    logger.info(f"[INFO] Found {len(samples)} sample(s): {', '.join(s.sample_id for s in samples)}")
    logger.info(f"[INFO] Active steps: {', '.join(active_steps)}")

    results: List[SampleResult] = []
    for pair in samples:
        res = process_sample(
            pair=pair,
            active_steps=active_steps,
            tools=resolved_tools,
            args=args,
            adapters_path=adapters_path,
            reference_path=reference_path,
            logger=logger,
        )
        results.append(res)

    # 8. MultiQC summary across all samples
    if "multiqc" in active_steps and ("multiqc" in resolved_tools or args.dry_run):
        multiqc_out = input_dir / "multiqc_report"
        logger.info("\n[STEP] Running MultiQC across all samples...")
        exec_cmd(
            [resolved_tools.get("multiqc", "multiqc"), str(input_dir), "-o", str(multiqc_out), "--force"],
            cwd=input_dir,
            dry_run=args.dry_run,
            logger=logger,
        )

    # Final Summary
    logger.info("\n" + "=" * 60)
    logger.info("PIPELINE SUMMARY")
    logger.info("=" * 60)
    failed = 0
    for r in results:
        if r.status == "FAILED":
            failed += 1
        detail_str = f" ({r.details})" if r.details else ""
        logger.info(f"  {r.sample_id:<25} {r.status:<8} {r.duration_sec:>6.1f}s{detail_str}")

    if log_path is not None:
        logger.info(f"\n[INFO] Full log saved to: {log_path}")

    if failed > 0:
        logger.error(f"\n[DONE] Finished with {failed} failed sample(s).")
        sys.exit(1)
    else:
        logger.info("\n[DONE] All samples completed successfully!")


if __name__ == "__main__":
    main()
