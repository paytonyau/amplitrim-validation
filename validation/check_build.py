#!/usr/bin/env python
"""One place for everything around the validation: check the build and the machine, see what is left to test, fetch the data that is missing, run it all.

  python validation\\check_build.py --overnight                 THE ONE COMMAND: find your folders, check the machine, fetch and prepare V4 and the inter-laboratory mock, run every step, write MORNING_REPORT.md
  python validation\\check_build.py --todo                      the whole board: what is done, what is ready, the decision, the manuscript and release items, what to do next
  python validation\\check_build.py --discover --write          find the `ready` folders (by their number of samples), input files and figaro.py, and save them
  python validation\\check_build.py --decide                    the rule for the default of --allow-unmergeable, applied to the saved coverage results
  python validation\\check_build.py --table                     the test-stage table, one for 16S and one for ITS, from the saved results (also written to TEST_TABLE.md)
  python validation\\check_build.py --init                      write amplitrim_paths.json (your folders; edit it once)
  python validation\\check_build.py                             check: the build id (and, with --manifest, which files differ), the tools, the data
  python validation\\check_build.py --plan                      every test and download that is left, what each needs, and the next command
  python validation\\check_build.py --list-runs mockB          the runs of that project by kind of sample name (so that --select can pick the right ones)
  python validation\\check_build.py --download mockB            fetch the reads of a new dataset from ENA and describe them (mockB, v4v5, v3v4); --select REGEX, --min-reads N, --dry-run
  python validation\\check_build.py --prepare mockB             remove the primers (run_dataset.py) and register the reads; --samples A,B --mixed-orientation --min-reads N
  python validation\\check_build.py --make-reference mockB=ZymoBIOMICS.STD.refseq.v2.zip   cut full-length 16S sequences between the dataset's primers, collapse copies, register them
  python validation\\check_build.py --set-reference mockB=zymo_16S.fasta      register expected sequences that are already amplicon-length
  python validation\\check_build.py --run --quick               every ready step on small slices (about 40 minutes): check that it all works
  python validation\\check_build.py --run                       the real run (a few hours; steps already finished are not repeated)

--manifest FILE compares the files of this copy with the list of another build (the lines written by the earlier check_build.py). --run passes
--only, --skip, --quick and --fresh on to overnight.py. Paths live in amplitrim_paths.json next to the package (--config chooses another).
"""

import argparse
import contextlib
import csv
import hashlib
import importlib.util
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))

DEFAULT_CONFIG = {
    "data": "C:\\data",
    "out": "C:\\data\\overnight",
    "mockA_high": "C:\\data\\results3\\mockA_high\\ready",
    "mockA_all": "C:\\data\\results2\\mockA_ee26\\ready",
    "runinfo": "C:\\data\\mockA\\runs_info.tsv",
    "asv_table": "C:\\data\\mockA\\Mock_staggered-A_ASVs.csv",
    "reference": "%USERPROFILE%\\Downloads\\mockA_reference.fasta",
    "its": "C:\\data\\results4\\cryobank_its3\\ready",
    "v1v2": "",
    "figaro": "",
    "extras": {},
}

# datasets that can be fetched; the expected sequences of a mock community are not downloadable by script and are registered with --set-reference
CATALOG = {
    "mockB": {"short": "mockB (independent mock, V1-V3)", "title": "an independent mock community (Zymo D6321, V1-V3)", "project": "PRJEB67827", "primer_f": "AGAGTTTGATYMTGGCTCAG", "primer_r": "ATTACCGCGGCTGCTGG",
              "reference": "Zymo's D6321 reference sequences (D6321.refseq.zip): python validation\\check_build.py --make-reference mockB=D6321.refseq.zip", "gives": "accuracy against a reference that does not come from a DADA2 run",
              "select": r"(53[6-9]|54[4-7])-BiomPsy", "mixed": True, "min_reads": 1000, "max_runs": 8, "auto": True,
              "note": "PRJEB67827 holds three mock communities and eight extraction protocols (94 runs); the mock B of Hulpusch, Rauer et al. 2023 is the three-species Zymo spike-in D6321 (Truepera radiovictrix, Imtechella halotolerans, Allobacillus halotolerans) in the samples 536-539 (1.1e5 cells) and 544-547 (5.55e3 cells); 552-555 are negative controls, left out"},
    "v4": {"short": "v4 (PRJNA574030, V4)", "title": "V4, 36 soil samples (PRJNA574030): used in development, so run again with the final build", "project": "PRJNA574030", "primer_f": "GTGYCAGCMGCCGCGGTAA", "primer_r": "GGACTACNVGGGTWTCTAAT",
           "reference": None, "auto": True, "max_runs": 36, "gives": "a V4 re-test of the final build on data that informed the development (not a held-out test)"},
    "v4v5": {"short": "v4v5 (PRJNA971391, V4-V5)", "title": "V4-V5 reads (PRJNA971391)", "project": "PRJNA971391", "primer_f": None, "primer_r": None, "reference": None,
             "gives": "a second region; exactness and calibration; coverage levels"},
    "mcm_v46": {"short": "mcm_v46 (PRJEB34919, V4-V6, real tight overlap)", "title": "inter-laboratory mock communities MCM2a and MCM2b, V4-V6 (about 564 bp amplicon, MiSeq 2x300: a real tight-overlap case; 6 samples)", "project": "PRJEB34919",
                "primer_f": "GTGCCAGCAGYYGCGGTAATAC", "primer_r": "CACRACACGAGCTGACGACA",
                "instrument": "MiSeq", "mixed": True, "min_reads": 5000, "max_runs": 12, "auto": True,
                "reference": "the expected 16S sequences of the 15 bacteria of Table 1 of O'Sullivan et al. 2021 (a FASTA built from their genomes; not downloadable by script)",
                "gives": "tight overlap on real reads with an independent truth (digital PCR); only 6 samples, so no paired intervals and no per-sample rank; the shotgun runs of the project are left out by the MiSeq filter"},
    "mcm_v12": {"short": "mcm_v12 (PRJEB34919, V1-V2)", "title": "inter-laboratory mock communities MCM2a and MCM2b, V1-V2 (about 337 bp amplicon, MiSeq 2x300; 6 samples)", "project": "PRJEB34919",
                "primer_f": "GYTCAGRDYKAACGCTGGCGG", "primer_r": "ACTGCTGCCTCCCGTAGGAGT",
                "instrument": "MiSeq", "mixed": True, "min_reads": 5000, "max_runs": 12, "auto": True, "raw_from": "mcm_v46",
                "reference": "the expected 16S sequences of the 15 bacteria of Table 1 of O'Sullivan et al. 2021 (a FASTA built from their genomes; not downloadable by script)",
                "gives": "a V1-V2 truth with an independent reference; only 6 samples; it shares the downloaded reads of mcm_v46 (the runs of the other region have no primers of this one and are left out)"},
    "sputum_v1v2": {"dropped": True, "short": "sputum_v1v2 (PRJNA915241, V1-V2)", "title": "sputum, region V1-V2 (PRJNA915241, QIAseq 16S/ITS panel, MiSeq 2x300; a held-out study: use --list-runs, then --select the runs of this region)",
                  "project": "PRJNA915241", "primer_f": None, "primer_r": None, "reference": None,
                  "gives": "a held-out test of V1-V2 on clinical samples (33 sputum samples per region); the primers are not given in the paper, so --prepare needs --primer-f and --primer-r read off the reads by the describer"},
    "sputum_v3v4": {"dropped": True, "short": "sputum_v3v4 (PRJNA915241, V3-V4)", "title": "sputum, region V3-V4 (PRJNA915241, QIAseq 16S/ITS panel, MiSeq 2x300; a held-out study: use --list-runs, then --select the runs of this region)",
                  "project": "PRJNA915241", "primer_f": None, "primer_r": None, "reference": None,
                  "gives": "a held-out test of V3-V4 on clinical samples (33 sputum samples per region); the primers are not given in the paper, so --prepare needs --primer-f and --primer-r read off the reads by the describer"},
    "sputum_v5v7": {"dropped": True, "short": "sputum_v5v7 (PRJNA915241, V5-V7)", "title": "sputum, region V5-V7 (PRJNA915241, QIAseq 16S/ITS panel, MiSeq 2x300; a held-out study: use --list-runs, then --select the runs of this region)",
                  "project": "PRJNA915241", "primer_f": None, "primer_r": None, "reference": None,
                  "gives": "a held-out test of V5-V7 on clinical samples (33 sputum samples per region); the primers are not given in the paper, so --prepare needs --primer-f and --primer-r read off the reads by the describer"},
    "sputum_v7v9": {"dropped": True, "short": "sputum_v7v9 (PRJNA915241, V7-V9)", "title": "sputum, region V7-V9 (PRJNA915241, QIAseq 16S/ITS panel, MiSeq 2x300; a held-out study: use --list-runs, then --select the runs of this region)",
                  "project": "PRJNA915241", "primer_f": None, "primer_r": None, "reference": None,
                  "gives": "a held-out test of V7-V9 on clinical samples (33 sputum samples per region); the primers are not given in the paper, so --prepare needs --primer-f and --primer-r read off the reads by the describer"},
    "v3v4": {"short": "v3v4 (V3-V4, tight overlap)", "title": "V3-V4 2x250: tight overlap on real reads (the candidate runs of PRJEB58189)", "project": None, "primer_f": None, "primer_r": None, "reference": None,
             "gives": "the tight-overlap case on real data; AmpliTrim against FIGARO"},
}
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
       "&fields=run_accession,sample_alias,library_layout,instrument_model,read_count&format=tsv")

PACKAGES = [("numpy", "numpy", True), ("pandas", "pandas", True), ("rapidfuzz", "rapidfuzz", True), ("Bio", "biopython", True), ("matplotlib", "matplotlib", False)]
TOOLS = [("vsearch", "vsearch", "bioconda"), ("cutadapt", "cutadapt", "bioconda")]


# ------------------------------------------------------------------------------------------------------------- config
def config_path(args) -> Path:
    return Path(args.config) if args.config else ROOT / "amplitrim_paths.json"


def load_config(path: Path) -> dict:
    cfg = dict(DEFAULT_CONFIG)
    if path.exists():
        cfg.update(json.loads(path.read_text(encoding="utf-8")))
    return cfg


def save_config(path: Path, cfg: dict) -> None:
    path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")


def p(cfg: dict, key: str):
    v = cfg.get(key)
    return Path(os.path.expandvars(v)) if v else None


def conda_env() -> str:
    return os.environ.get("CONDA_DEFAULT_ENV") or Path(sys.prefix).name


# ------------------------------------------------------------------------------------------------------------- the build
def tracked_files() -> list:
    return sorted(list((ROOT / "amplitrim").glob("*.py")) + list(HERE.glob("*.py")) + list(HERE.glob("*.R")), key=lambda p: (p.parent.name, p.name))


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()[:10]


def compare_manifest(manifest: Path) -> dict:
    """{'same': n, 'differ': [...], 'missing': [...], 'extra': [...]} between this copy and a manifest of lines '<hash>  <relative path>'."""
    wanted = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2 and len(parts[0]) == 10 and not line.startswith("build "):
            wanted[parts[1].strip()] = parts[0]
    have = {f.relative_to(ROOT).as_posix(): file_hash(f) for f in tracked_files()}
    return {"same": sum(1 for k, v in wanted.items() if have.get(k) == v), "differ": sorted(k for k, v in wanted.items() if k in have and have[k] != v),
            "missing": sorted(k for k in wanted if k not in have), "extra": sorted(k for k in have if k not in wanted)}


SHIPPED_MANIFEST = HERE / "build_manifest.txt"          # written when the zip is made: the hashes of the files it contained


def file_lines() -> list:
    return [f"{file_hash(f)}  {f.relative_to(ROOT).as_posix()}" for f in tracked_files()]


def check_build(args) -> list:
    from run_dataset import build_id
    lines = [f"build {build_id()} ({len(tracked_files())} files hashed: amplitrim/*.py, validation/*.py, validation/*.R)"]
    manifest = Path(args.manifest) if args.manifest else (SHIPPED_MANIFEST if SHIPPED_MANIFEST.is_file() else None)
    if args.manifest and not manifest.is_file():
        lines.append(f"manifest file not found: {manifest}")
        lines.append("  give its full path (for example C:\\Users\\Admin\\Downloads\\build_manifest_xxxxxxxxxx.txt); without it, --files prints this copy's hashes to compare by hand")
        manifest = SHIPPED_MANIFEST if SHIPPED_MANIFEST.is_file() else None
    if manifest is not None:
        first = manifest.read_text(encoding="utf-8").splitlines()[:1]
        shipped = first[0].split(":")[0] if first and first[0].startswith("build ") else "?"
        c = compare_manifest(manifest)
        lines.append(f"against {manifest.name} ({shipped}): {c['same']} files identical")
        for key, text in (("differ", "differs from that build (edited or an older version)"), ("missing", "is in that build but not here"),
                          ("extra", "is here but not in that build (it changes the build id)")):
            for f in c[key]:
                lines.append(f"  {f}  {text}")
        if not (c["differ"] or c["missing"] or c["extra"]):
            lines.append("  every file matches: this copy is exactly that build")
    else:
        lines.append("no manifest to compare with (the zip's validation/build_manifest.txt is missing); --files prints the hashes of this copy")
    return lines


# ------------------------------------------------------------------------------------------------------------- the machine
def tool_version(cmd: list):
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if res.returncode != 0:
            return None
        lines = [x.strip() for x in (res.stdout + "\n" + res.stderr).splitlines() if x.strip()]
        lines = [x for x in lines if not re.search(r"doi|https?://", x, re.I)]                              # a citation is not a version
        named = next((x for x in lines if re.search(r"\bv?\d+\.\d+", x) and x.lower().startswith(Path(cmd[0]).stem.lower())), None)
        return named or next((x for x in lines if re.search(r"\bv?\d+\.\d+", x)), lines[0] if lines else None)
    except (OSError, subprocess.SubprocessError):
        return None


def check_machine() -> list:
    """[(item, ok, detail, how to fix)]"""
    env = conda_env()
    out = [("python", True, sys.version.split()[0] + f" (environment {env})", "")]
    for module, conda_name, needed in PACKAGES:
        found = importlib.util.find_spec(module) is not None
        out.append((module, found or not needed, "installed" if found else ("MISSING" if needed else "missing (only the figures need it)"),
                    "" if found else f"conda install -n {env} -c conda-forge {conda_name}"))
    for exe, conda_name, channel in TOOLS:
        v = tool_version([exe, "--version"]) if shutil.which(exe) else None
        out.append((exe, v is not None, v or "MISSING", "" if v else f"conda install -n {env} -c conda-forge -c {channel} {conda_name}"))
    import rtools
    rscript = rtools.find_rscript()
    if rscript is None:
        out.append(("R", False, "MISSING", f"conda install -n {env} -c conda-forge r-base   (or install R from cran.r-project.org and put it on PATH)"))
    else:
        out.append(("R", True, tool_version([str(rscript), "--version"]) or "found", ""))
        v = tool_version([str(rscript), "-e", 'cat(as.character(packageVersion("dada2")))'])
        out.append(("R + dada2", v is not None, f"dada2 {v}" if v else "R found, dada2 MISSING",
                    "" if v else f"conda install -n {env} -c conda-forge -c bioconda bioconductor-dada2"))
    return out


def check_data(cfg: dict) -> list:
    from amplitrim.subsample import find_pairs
    out = []
    for key in ("mockA_high", "mockA_all", "its", "v1v2"):
        path = p(cfg, key)
        if path is None:
            out.append((key, None, "not set in amplitrim_paths.json"))
        elif not path.is_dir():
            out.append((key, False, f"{path} does not exist"))
        else:
            try:
                out.append((key, True, f"{path}: {len(find_pairs(path))} samples"))
            except Exception as err:                                             # a folder without paired reads
                out.append((key, False, f"{path}: {str(err)[:80]}"))
    for key in ("runinfo", "asv_table", "reference", "figaro"):
        path = p(cfg, key)
        out.append((key, None if path is None else path.is_file(), "not set" if path is None else (str(path) if path.is_file() else f"{path} does not exist")))
    for name, e in cfg.get("extras", {}).items():
        if not e.get("ready"):
            out.append((f"extra {name}", False, f"reference registered, reads not prepared (--prepare {name})"))
            continue
        ready = Path(os.path.expandvars(e["ready"]))
        out.append((f"extra {name}", ready.is_dir(), f"{ready}" + (f" with reference {e['reference']}" if e.get("reference") else " (no reference)")))
    return out


# ------------------------------------------------------------------------------------------------------------- downloads
def sample_kind(alias: str) -> str:
    """A sample name with its numbers replaced, so that 103-KitA-1-1-8_S103_L001 and 99-KitA-1-1-4_S99_L001 fall in one group."""
    return re.sub(r"_S#(_L#)?$", "", re.sub(r"\d+", "#", alias))


def fetch_run_table(project: str, dest: Path, opener=urllib.request.urlopen, select: str = None, min_reads: int = 0, instrument: str = None) -> list:
    """The paired runs of an ENA/SRA project: writes runs.tsv (the whole table) and runs.txt (one accession per line, only the runs whose sample name matches `select` and that
    have at least `min_reads` reads) in dest; returns those accessions."""
    dest.mkdir(parents=True, exist_ok=True)
    with opener(ENA.format(acc=project), timeout=120) as r:
        text = r.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows or "run_accession" not in rows[0]:
        raise SystemExit(f"ENA returned no run table for {project}")
    paired = [x for x in rows if (x.get("library_layout") or "PAIRED").upper() == "PAIRED"]
    if not paired:
        raise SystemExit(f"{project} has no paired runs ({len(rows)} runs found)")
    (dest / "runs.tsv").write_text(text, encoding="utf-8")
    if instrument:
        paired = [x for x in paired if instrument.lower() in (x.get("instrument_model") or "").lower()]
        if not paired:
            raise SystemExit(f"{project} has no paired run on a {instrument} instrument (its runs.tsv lists the instruments)")
    chosen = [x for x in paired if (not select or re.search(select, x.get("sample_alias") or "", re.I)) and int(float(x.get("read_count") or 0)) >= min_reads]
    if not chosen:
        raise SystemExit(f"no run of {project} matches --select {select!r} with at least {min_reads} reads ({len(paired)} paired runs; --list-runs shows the sample names)")
    accessions = [x["run_accession"] for x in chosen]
    (dest / "runs.txt").write_text("\n".join(accessions) + "\n", encoding="utf-8")
    return accessions


def list_runs(cfg: dict, name: str, opener=urllib.request.urlopen) -> list:
    """The runs of a dataset's project grouped by kind of sample name (what --select can choose), with read counts and what is already downloaded."""
    entry = CATALOG[name]
    d = dataset_dirs(cfg, name)
    if not entry["project"]:
        raise SystemExit(f"{name} has no project; its runs come from --runs-file")
    if not (d["base"] / "runs.tsv").exists():
        fetch_run_table(entry["project"], d["base"], opener)
    rows = [r for r in csv.DictReader((d["base"] / "runs.tsv").read_text(encoding="utf-8").splitlines(), delimiter="\t") if (r.get("library_layout") or "PAIRED").upper() == "PAIRED"]
    groups = {}
    for r in rows:
        groups.setdefault(sample_kind(r.get("sample_alias") or "?"), []).append(r)
    lines = [f"{entry['project']}: {len(rows)} paired runs in {len(groups)} kinds of sample name; choose with --select REGEX (a match on the sample name, any case) and --min-reads N", ""]
    lines.append(f"{'kind of sample name':40s} {'runs':>5s} {'reads min / median / max':>28s} {'downloaded':>10s}  example")
    for kind, rs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        reads = sorted(int(float(r.get("read_count") or 0)) for r in rs)
        have = sum((d["raw"] / f"{r['run_accession']}_R1.fastq.gz").exists() for r in rs)
        lines.append(f"{kind[:40]:40s} {len(rs):5d} {f'{reads[0]:,} / {reads[len(reads) // 2]:,} / {reads[-1]:,}':>28s} {have:10d}  {rs[0].get('sample_alias')} ({rs[0]['run_accession']})")
    return lines


def dataset_dirs(cfg: dict, name: str) -> dict:
    base = Path(os.path.expandvars(cfg["data"])) / "new" / name
    raw_base = Path(os.path.expandvars(cfg["data"])) / "new" / CATALOG.get(name, {}).get("raw_from", name)          # a dataset can share the downloaded reads of another
    return {"base": base, "raw": raw_base / "raw", "prepared": base / "prepared", "ready": base / "prepared" / name / "ready"}


def downloaded_aliases(d: dict) -> list:
    """One line per downloaded run with its sample name and read count (from runs.tsv), so that the kinds of sample in a download can be told apart."""
    table = d["base"] / "runs.tsv"
    if not table.exists():
        return []
    rows = {r["run_accession"]: r for r in csv.DictReader(table.read_text(encoding="utf-8").splitlines(), delimiter="\t")}
    have = sorted(f.name.split("_R1")[0] for f in d["raw"].glob("*_R1.fastq.gz")) if d["raw"].is_dir() else []
    return [f"  {run}  {rows[run].get('sample_alias', '?'):36s} {int(float(rows[run].get('read_count') or 0)):>9,} read pairs" for run in have if run in rows]


def download(cfg: dict, name: str, args, run=subprocess.run, opener=urllib.request.urlopen) -> int:
    entry = CATALOG[name]
    d = dataset_dirs(cfg, name)
    if entry.get("raw_from") and d["raw"].is_dir() and any(d["raw"].glob("*_R1.fastq.gz")) and not args.runs_file:
        print(f"{name} uses the reads that were downloaded for {entry['raw_from']} ({d['raw']}): nothing to download")
        return 0
    if args.runs_file:
        runs = Path(args.runs_file)
        if not runs.is_file():
            raise SystemExit(f"--runs-file {runs} does not exist")
    elif entry["project"]:
        select = getattr(args, "select", None) or entry.get("select")                      # the catalog knows which samples of a shared project belong to the dataset
        args.select = select
        accessions = fetch_run_table(entry["project"], d["base"], opener, select, getattr(args, "min_reads", 0) or 0, entry.get("instrument"))
        runs = d["base"] / "runs.txt"
        chosen = f" match --select {args.select!r}" if getattr(args, "select", None) else " listed"
        print(f"{entry['project']}: {len(accessions)} paired runs{chosen}; downloading {min(args.max_runs, len(accessions))} (spread evenly; --max-runs changes it)")
    else:
        raise SystemExit(f"{name} has no project: give the runs with --runs-file (one run accession per line, for example the candidate list of PRJEB58189)")
    d["raw"].mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(HERE / "download_ena.py"), str(runs), "--out", str(d["raw"]), "--max-runs", str(args.max_runs)]
    if args.dry_run:
        cmd.append("--dry-run")
    if run(cmd).returncode != 0:
        return 1
    if not args.dry_run:
        aliases = downloaded_aliases(d)
        if aliases:
            print("\nThe sample names of the downloaded runs (--list-runs shows all kinds; --select chooses):\n" + "\n".join(aliases))
        print("\nWhat the reads look like (the describer lists the primers it finds):")
        run([sys.executable, str(HERE / "describe_reads.py"), str(d["raw"]), "--lengths", "--positions", "15"])
        print(f"\nNext: python validation\\check_build.py --prepare {name}" + ("" if entry["primer_f"] else " --primer-f <forward primer> --primer-r <reverse primer>"))
        if entry["reference"]:
            print(f"Also needed for the accuracy tests: {entry['reference']}; register it with --set-reference {name}=<file.fasta>")
    return 0


def prepare(cfg: dict, name: str, args, cfg_path: Path, run=subprocess.run) -> int:
    entry = CATALOG[name]
    d = dataset_dirs(cfg, name)
    if not d["raw"].is_dir() or not any(d["raw"].glob("*_R1.fastq.gz")):
        raise SystemExit(f"no reads in {d['raw']}: run --download {name} first")
    pf, pr = args.primer_f or entry["primer_f"], args.primer_r or entry["primer_r"]
    no_primers = getattr(args, "no_primers", False)
    mixed = getattr(args, "mixed_orientation", False) or bool(entry.get("mixed"))                  # a library made by adapter ligation is sequenced from either end
    min_reads = (getattr(args, "min_reads", 0) or 0) or entry.get("min_reads", 0)                   # samples without this dataset's primers (the other region's runs) fall below it
    if not (pf and pr) and not no_primers:
        raise SystemExit(f"the primers of {name} are not known: give --primer-f and --primer-r (the describer's list above shows what the reads contain), or --no-primers if the reads are already free of primers")
    cmd = [sys.executable, str(HERE / "run_dataset.py"), name, str(d["raw"]), "--out", str(d["prepared"]), "--error-rate", "0.15", "--pairs", "100000"]
    if not no_primers:
        cmd += ["--primer-f", pf, "--primer-r", pr]
    if getattr(args, "samples", None):
        cmd += ["--only", args.samples]
    if mixed:
        cmd.append("--mixed-orientation")
    if min_reads > 0:
        cmd += ["--min-reads", str(min_reads)]
    if run(cmd).returncode != 0:
        return 1
    if not d["ready"].is_dir():
        print(f"run_dataset.py finished but {d['ready']} was not found")
        return 1
    extras = cfg.setdefault("extras", {})
    extras.setdefault(name, {})["ready"] = str(d["ready"])
    extras[name].setdefault("reference", "")
    save_config(cfg_path, cfg)
    print(f"\nregistered {name}: {d['ready']}" + ("" if extras[name]["reference"] else "  (no reference yet)"))
    return 0


def set_reference(cfg: dict, spec: str, cfg_path: Path) -> int:
    name, _, path = spec.partition("=")
    if name not in CATALOG:
        raise SystemExit(f"{name} is not one of {', '.join(sorted(CATALOG))}")
    if not Path(os.path.expandvars(path)).is_file():
        raise SystemExit(f"{path} does not exist")
    cfg.setdefault("extras", {}).setdefault(name, {"ready": "", "reference": ""})["reference"] = path
    save_config(cfg_path, cfg)
    print(f"{name}: reference set to {path}" + ("" if cfg["extras"][name].get("ready") else f"  (the reads are not prepared yet: --prepare {name})"))
    return 0


def find_reference_candidates(cfg: dict, name: str, limit: int = 8) -> list:
    """Zip files, folders and FASTA files in the usual places whose name looks like the one asked for (a ZymoBIOMICS refseq, D6300, D6321, or the same stem)."""
    stem = re.sub(r"\.(zip|fa|fasta|fna)$", "", name, flags=re.I).lower()
    pattern = re.compile(r"refseq|d63\d\d|zymo" + (("|" + re.escape(stem)) if stem else ""), re.I)
    home = Path(os.path.expanduser("~"))
    found = []
    for root in (home / "Downloads", home / "Desktop", home / "Documents", Path(os.path.expandvars(cfg.get("data", "")))):
        if not root.is_dir():
            continue
        base = len(root.parts)
        for dirpath, dirnames, files in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and len(Path(dirpath).parts) - base < 3]
            for n in dirnames + files:
                if pattern.search(n) and (n.lower().endswith((".zip", ".fa", ".fasta", ".fna")) or n in dirnames):
                    found.append(Path(dirpath) / n)
            if len(found) >= limit:
                return found[:limit]
    return found


def compare_reference_for(cfg: dict, spec: str, args) -> int:
    """Cut a published full-length reference between the primers and compare it with the expected sequences already registered for the dataset, without changing anything:
    which strains the two share, and which only one of them has."""
    import make_reference as mr
    name, _, source = spec.partition("=")
    data = Path(os.path.expandvars(cfg["data"]))
    if name == "mockA":
        pf, pr, base, ref = (getattr(args, "primer_f", None) or "AGAGTTTGATYMTGGCTCAG", getattr(args, "primer_r", None) or "ATTACCGCGGCTGCTGG", data / "mockA", p(cfg, "reference"))
    elif name in CATALOG:
        pf, pr = getattr(args, "primer_f", None) or CATALOG[name]["primer_f"], getattr(args, "primer_r", None) or CATALOG[name]["primer_r"]
        base = dataset_dirs(cfg, name)["base"]
        ref = (lambda r: Path(os.path.expandvars(r)) if r else None)((cfg.get("extras", {}).get(name) or {}).get("reference"))
    else:
        raise SystemExit("use --compare-reference NAME=INPUT with NAME mockA or one of " + ", ".join(sorted(CATALOG)))
    if not (pf and pr):
        raise SystemExit(f"the primers of {name} are not known: give --primer-f and --primer-r")
    if ref is None or not Path(ref).is_file():
        raise SystemExit(f"no expected sequences are registered for {name} (--set-reference {name}=FILE, or the reference of amplitrim_paths.json for mockA)")
    src = Path(os.path.expandvars(source))
    if not src.exists():
        hits = find_reference_candidates(cfg, src.name)
        raise SystemExit(f"{src} does not exist. (A path such as C:\\path\\... in an example is a placeholder: give the real location of the file.)"
                         + ("\nFiles or folders with a similar name that were found:\n  " + "\n  ".join(str(h) for h in hits) if hits else ""))
    fasta, lines = mr.build(mr.read_inputs([src]), pf, pr)
    print(f"primers used for {name}: {pf} / {pr}\n" + "\n".join(lines))
    if not fasta:
        raise SystemExit("no sequence has sites for both primers: check the primers (--primer-f, --primer-r) and the reference")
    out = base / "reference_check" / f"{src.stem}_amplicons.fasta"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(fasta, encoding="utf-8")
    labels = {}
    map_file = Path(ref).with_name(Path(ref).stem + "_map.tsv")
    if map_file.is_file():
        with open(map_file, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("strain_id") and row.get("species"):
                    labels[row["strain_id"]] = row["species"]
    print("\n=== comparison with " + str(ref))
    print("\n".join(mr.compare_references(out, Path(ref), labels=labels)))
    print(f"written: {out}")
    return 0


def make_reference_for(cfg: dict, spec: str, args, cfg_path: Path) -> int:
    """Cut a full-length reference (a FASTA, a folder or a zip such as ZymoBIOMICS.STD.refseq) between the primers of the dataset, collapse identical copies, and register the result."""
    import make_reference as mr
    name, _, source = spec.partition("=")
    if name not in CATALOG or not source:
        raise SystemExit("use --make-reference NAME=INPUT, for example --make-reference mockB=ZymoBIOMICS.STD.refseq.v2.zip (NAME is one of " + ", ".join(sorted(CATALOG)) + ")")
    entry = CATALOG[name]
    pf, pr = getattr(args, "primer_f", None) or entry["primer_f"], getattr(args, "primer_r", None) or entry["primer_r"]
    if not (pf and pr):
        raise SystemExit(f"the primers of {name} are not known: give --primer-f and --primer-r")
    src = Path(os.path.expandvars(source))
    if not src.exists():
        hits = find_reference_candidates(cfg, src.name)
        raise SystemExit(f"{src} does not exist. (A path such as C:\\path\\... in an example is a placeholder: give the real location of the file.)"
                         + ("\nFiles or folders with a similar name that were found:\n  " + "\n  ".join(str(h) for h in hits) if hits else
                            "\nNothing with a similar name was found in Downloads, Desktop, Documents or the data folder."))
    fasta, lines = mr.build(mr.read_inputs([src]), pf, pr)
    print("\n".join(lines))
    if not fasta:
        raise SystemExit("no sequence has sites for both primers: check the primers (--primer-f, --primer-r) and the reference")
    out = dataset_dirs(cfg, name)["base"] / "reference.fasta"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(fasta, encoding="utf-8")
    print(f"written: {out}")
    return set_reference(cfg, f"{name}={out}", cfg_path)


# ------------------------------------------------------------------------------------------------------------- finding what is already on this machine
SKIP_DIRS = {"miniconda3", "anaconda3", "mambaforge", "AppData", "node_modules", ".git", "site-packages", "Windows", "Program Files", "Program Files (x86)", "$Recycle.Bin", "__pycache__"}
SIZES = {6: "mockA_high", 23: "mockA_all", 33: "v1v2", 36: "v4", 37: "its"}          # the number of samples of each dataset in `ready` folders
FILES = {"runinfo": "runs_info.tsv", "asv_table": "Mock_staggered-A_ASVs.csv", "reference": "mockA_reference.fasta"}


def find_ready_folders(root: Path, depth: int = 6) -> list:
    """[(folder, samples, modified)] for every folder called `ready` under root that holds paired reads."""
    from amplitrim.subsample import find_pairs
    out = []
    root = Path(root)
    if not root.is_dir():
        return out
    base = len(root.parts)
    for dirpath, dirnames, _files in os.walk(root):
        here = Path(dirpath)
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        if len(here.parts) - base >= depth:
            dirnames[:] = []
        if here.name == "ready":
            try:
                n = len(find_pairs(here))
            except Exception:
                n = 0
            if n:
                out.append((here, n, here.stat().st_mtime))
            dirnames[:] = []
    return out


def find_file(name: str, roots: list, depth: int = 5):
    for root in roots:
        root = Path(root)
        if not root.is_dir():
            continue
        base = len(root.parts)
        for dirpath, dirnames, files in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
            if len(Path(dirpath).parts) - base >= depth:
                dirnames[:] = []
            if name in files:
                return Path(dirpath) / name
    return None


def discover(cfg: dict) -> dict:
    """{config key: (path, how it was recognised)} for what is missing from the config and can be found: the reads by their number of samples, the files by
    name, FIGARO by figaro.py."""
    data = Path(os.path.expandvars(cfg["data"]))
    home = Path(os.path.expanduser("~"))
    found = {}
    by_size = {}
    for folder, n, mtime in find_ready_folders(data):
        if n in SIZES:
            by_size.setdefault(SIZES[n], []).append((mtime, folder, n))
    for key, items in by_size.items():
        cur = p(cfg, key) if key != "v4" else (Path(os.path.expandvars(cfg.get("extras", {}).get("v4", {}).get("ready") or "")) if cfg.get("extras", {}).get("v4", {}).get("ready") else None)
        if cur is not None and cur.is_dir():
            continue
        items.sort(reverse=True)                                                  # the newest folder with that many samples
        preferred = [i for i in items if re.search(r"42528|v1v2|v1-v2", str(i[1]), re.I)] if key == "v1v2" else []
        mtime, folder, n = (preferred or items)[0]
        found[key] = (folder, f"a `ready` folder with {n} samples" + (f" ({len(items) - 1} other with that many)" if len(items) > 1 else ""))
    for key, name in FILES.items():
        cur = p(cfg, key)
        if cur is None or not cur.is_file():
            hit = find_file(name, [data, home / "Downloads", home])
            if hit:
                found[key] = (hit, f"found by name ({name})")
    cur = p(cfg, "figaro")
    if cur is None or not cur.is_file():
        hit = find_file("figaro.py", [Path("C:/figaro"), home / "figaro", home / "Downloads", data, home], depth=5)
        if hit:
            found["figaro"] = (hit, "found by name (figaro.py)")
    return found


def apply_discovery(cfg: dict, found: dict) -> dict:
    for key, (path, _how) in found.items():
        if key == "v4":                                                              # a dataset that is tested like a new one
            cfg.setdefault("extras", {}).setdefault("v4", {"ready": "", "reference": ""})["ready"] = str(path)
        else:
            cfg[key] = str(path)
    return cfg


# ------------------------------------------------------------------------------------------------------------- what the results say
def step_states(out: Path) -> dict:
    res = {}
    for f in sorted(Path(out).glob("*/result.json")):
        try:
            r = json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            continue
        r["_when"] = time.strftime("%Y-%m-%d %H:%M", time.localtime(f.stat().st_mtime))
        res[f.parent.name] = r
    return res


def collect_coverage_data(out: Path) -> list:
    data = []
    r = step_states(out).get("coverage_levels")
    if r and r.get("status") == "ok":
        data += r.get("data", [])
    return data


def decide_coverage(datasets: list, level: float = 0.002, max_median_cost: float = 2.0, max_cost: float = 5.0) -> dict:
    """The rule for changing the default of --allow-unmergeable from 1% to 0.2%, stated before the evidence was in:
    on at least two datasets with expected sequences, covering the longest 0.2% instead of 1% recovers at least one more strain on at least one of them and fewer
    variants or strains on none; and across all datasets it costs a median of at most 2 points of merged pairs and at most 5 on any."""
    entries = []
    for ds in datasets:
        rows = {round(r["level"], 4): r for r in ds["rows"]}
        a, b = rows.get(0.01), rows.get(level)
        if not a or not b:
            continue
        e = {"name": ds["name"], "same": (a["tf"], a["tr"]) == (b["tf"], b["tr"]), "from": f"{a['tf']}/{a['tr']}", "to": f"{b['tf']}/{b['tr']}", "cost": 100 * (a["merged"] - b["merged"])}
        if ds.get("has_truth") and a.get("strains_found") is not None:
            e.update({"d_strains": b["strains_found"] - a["strains_found"], "d_variants": b["variants"] - a["variants"], "d_outside": b["outside"] - a["outside"]})
        entries.append(e)
    truth = [e for e in entries if "d_strains" in e]
    lines = [f"  {e['name']:12s} {e['from']} -> {e['to']}: merged pairs {-e['cost']:+.1f} points" + (f"; variants {e['d_variants']:+d}, strains {e['d_strains']:+d}, ASVs outside the reference {e['d_outside']:+d}" if "d_strains" in e else "; no expected sequences, cost only")
             for e in entries]
    if len(entries) < 2 or len(truth) < 2:
        need = []
        if len(truth) < 2:
            need.append("a second dataset with expected sequences (an independent mock community: --download mockB, --prepare, --set-reference)")
        if len(entries) - len(truth) < 1:
            need.append("a dataset without expected sequences for the cost (V1-V2: set v1v2 in amplitrim_paths.json)")
        return {"verdict": "UNDETERMINED", "lines": lines, "why": "needs " + " and ".join(need), "entries": entries}
    c_benefit = any(e["d_strains"] >= 1 for e in truth) and all(e["d_strains"] >= 0 and e["d_variants"] >= 0 for e in truth)
    c_median = statistics.median(e["cost"] for e in entries) <= max_median_cost
    c_max = max(e["cost"] for e in entries) <= max_cost
    why = [f"benefit on the datasets with expected sequences: {'yes' if c_benefit else 'NO'}", f"median cost {statistics.median(e['cost'] for e in entries):.1f} points (limit {max_median_cost}): {'yes' if c_median else 'NO'}",
           f"largest cost {max(e['cost'] for e in entries):.1f} points (limit {max_cost}): {'yes' if c_max else 'NO'}"]
    return {"verdict": "ADOPT 0.2% as the default" if (c_benefit and c_median and c_max) else "KEEP 1%", "lines": lines, "why": "; ".join(why), "entries": entries}


MANUSCRIPT = [
    "Title block: affiliation(s) and the target journal.",
    "Table 1: the platform and read length of the V4 dataset (PRJNA574030).",
    "Figures 3 and 5 and Table 9: their values were transcribed from the comparison outputs; regenerate them from the saved compare_methods.tsv files (V1-V2, ITS at 5,000 pairs and at full depth, and the battery folders for the other datasets).",
    "ITS: reconcile the cost of a 200/200 truncation in the final-build battery (+16.1, +14.3 to +17.9 for AmpliTrim minus 200/200) with the earlier run (+23.5, +21.8 to +25.2, Table 9 and criterion 4) before either number is cited.",
    "Availability: repository URL, licence, the final build identifier, operating systems.",
    "Declarations: competing interests, the AI-tools statement (adapt to the journal's policy), funding, author contributions, acknowledgements, ethics.",
    "References: check all of them against the sources (10 is from the paper itself); add the titles of 2, 5, 7 and 8; add the reference for PRJEB42528 (9).",
    "After each batch of results: update Results, Table 10, the Discussion and the Abstract (send MORNING_REPORT.md and OVERNIGHT_SUMMARY.md).",
    "A last full run on the final build (python validation\\check_build.py --overnight, after --fresh on the steps), so that the paper cites one build id; say that earlier ids depended on the file order.",
]
RELEASE = [
    "Run the CI workflow (.github/workflows) once on GitHub and fix what it finds (it has never run).",
    "Push the repository to GitHub (needs your credentials) and put its URL in the manuscript and in CITATION.cff.",
    "Choose a licence and add the LICENSE file (the mock A reference is derived from a GPL-3.0 repository and is not bundled).",
    "Archive the release for a DOI (for example Zenodo from the GitHub release) and cite it.",
    "Test the install on a second machine or system (pip install -e . and python validation\\check_build.py).",
]


def todo(cfg: dict, args) -> list:
    import overnight
    out = Path(os.path.expandvars(cfg["out"]))
    states = step_states(out)
    extras = [f"{n}={os.path.expandvars(e['ready'])}" + (f"|{os.path.expandvars(e['reference'])}" if e.get("reference") else "") for n, e in cfg.get("extras", {}).items() if e.get("ready")]
    ctx = SimpleNamespace(args=SimpleNamespace(extra=extras), p=lambda name: p(cfg, name))
    lines = ["== 1. Tests (overnight.py); results in " + str(out)]
    ready = []
    for sid, tier, minutes, _q, needs, tools, _fn in overnight.STEPS:
        st = states.get(sid)
        if st and st.get("status") == "ok":
            notmet = sum(1 for v in st.get("verdicts", []) if v["outcome"] == "NOT MET")
            stale = " (saved before the decision rule existed: it reruns)" if sid == "coverage_levels" and "data" not in st else ""
            lines.append(f"  DONE    {sid:16s} {st['_when']}" + (f"  {notmet} criterion not met" if notmet else "") + stale)
            if stale:
                ready.append(sid)
        elif st and st.get("status") == "failed":
            lines.append(f"  FAILED  {sid:16s} {st.get('error', '')[:100]}" + ("   (round 4: --round4)" if tier == 4 else ""))
            if tier <= 3:
                ready.append(sid)
        else:
            problems = overnight.readiness(ctx, needs, tools)
            lines.append(f"  {'READY  ' if not problems else 'BLOCKED'} {sid:16s} " + ("about %d min" % minutes if not problems else "; ".join(problems)) + ("   (round 4: --round4)" if tier == 4 else ""))
            if not problems and tier <= 3:                                   # the steps of round 4 are run by --round4, not by --overnight
                ready.append(sid)
    lines += ["", "== 2. Decisions"]
    d = decide_coverage(collect_coverage_data(out))
    lines += [f"  Default of --allow-unmergeable (1% -> 0.2%): {d['verdict']}"] + d["lines"] + [f"  {d['why']}"]
    lines += ["", "== 3. New data"]
    for name, e in CATALOG.items():
        if e.get("dropped"):
            continue
        reg = cfg.get("extras", {}).get(name, {})
        ref = ", reference ok" if reg.get("reference") else (", reference MISSING" if e["reference"] else "")
        lines.append(f"  {name:6s} {'registered' if reg.get('ready') else 'not prepared'}{ref}   ({e['title']})")
    lines += ["", "== 4. Manuscript (things only you can complete)"] + [f"  [ ] {x}" for x in MANUSCRIPT]
    lines += ["", "== 5. Release"] + [f"  [ ] {x}" for x in RELEASE]
    nxt, found = [], discover(cfg)
    if found:
        nxt.append("python validation\\check_build.py --discover --write      (found: " + ", ".join(sorted(found)) + ")")
    if ready:
        nxt.append("python validation\\check_build.py --overnight              (runs " + ", ".join(ready) + ")")
    if not cfg.get("extras", {}).get("mockB", {}).get("ready"):
        nxt.append("python validation\\check_build.py --list-runs mockB         (then --download mockB --select ..., --prepare mockB ..., --make-reference mockB=D6321.refseq.zip)")
    nxt.append("python validation\\check_build.py --round4 --quick   (then without --quick): the open questions of round 4 in phases, with breaks that go on by themselves (criteria: TEST_REGISTER.md, section 9)")
    nxt.append("then send OVERNIGHT_SUMMARY.md, and work through sections 4 and 5")
    lines += ["", "== Next, in this order"] + [f"  {i}. {x}" for i, x in enumerate(nxt, 1)]
    return lines


AUTO_ORDER = ["v4", "mockB", "mcm_v46", "mcm_v12"]                  # the datasets the overnight run fetches and prepares without help (mcm_v12 uses the reads of mcm_v46)


class Tee:
    """What the overnight run prints, to the screen and to a file."""
    def __init__(self, stream, path: Path):
        self.stream, self.f = stream, open(path, "a", encoding="utf-8")

    def write(self, s):
        self.stream.write(s)
        self.f.write(s)
        self.f.flush()

    def flush(self):
        self.stream.flush()

    def close(self):
        self.f.close()


@contextlib.contextmanager
def keep_awake():
    """Windows goes to sleep after a while and stops the run; this asks it not to for as long as the run lasts (nothing happens elsewhere)."""
    windll = None
    if os.name == "nt":
        try:
            import ctypes
            windll = ctypes.windll.kernel32
            windll.SetThreadExecutionState(0x80000000 | 0x00000001)          # ES_CONTINUOUS | ES_SYSTEM_REQUIRED
        except Exception:
            windll = None
    try:
        yield
    finally:
        if windll is not None:
            windll.SetThreadExecutionState(0x80000000)


def auto_dataset(cfg: dict, name: str, cfg_path: Path, download_fn, prepare_fn) -> str:
    """Download (if the reads are not there), trim and register one dataset; the outcome as a line of the report. Never raises: a dataset that cannot be fetched does not stop the run."""
    entry = CATALOG[name]
    registered = (cfg.get("extras", {}).get(name) or {}).get("ready")
    if registered and Path(os.path.expandvars(registered)).is_dir():
        return "already registered"
    dl = SimpleNamespace(runs_file=None, select=entry.get("select"), min_reads=0, max_runs=entry.get("max_runs", 12), dry_run=False)
    pr = SimpleNamespace(primer_f=None, primer_r=None, no_primers=False, samples=None, mixed_orientation=False, min_reads=entry.get("min_reads", 0))
    try:
        d = dataset_dirs(cfg, name)
        rc = download_fn(cfg, name, dl)                                      # resumable: files already verified are skipped, the missing runs are tried again
        have = d["raw"].is_dir() and any(d["raw"].glob("*_R1.fastq.gz"))
        if rc and not have:
            return f"download failed (exit {rc}); run --download {name} to see why"
        rc_prepare = prepare_fn(cfg, name, pr, cfg_path)
        if rc_prepare:
            return f"preparation failed (exit {rc_prepare}); run --prepare {name} to see why"
        return "downloaded, trimmed and registered" + (f" (the download was not complete, so some runs are missing: run --download {name} again, then --overnight)" if rc else "")
    except SystemExit as err:
        return f"stopped: {err}"
    except Exception as err:
        return f"failed: {type(err).__name__}: {err}"


def morning_report(cfg: dict, args, phases: list, code: int, started: str) -> Path:
    """One file to read in the morning: what ran, the step verdicts, the test tables, the decision, what is left and what needs you."""
    import overnight
    out = Path(os.path.expandvars(cfg["out"]) + ("_quick" if args.quick else ""))
    states = step_states(out)
    L = ["# AmpliTrim overnight run: morning report", "", f"Started {started}, finished {time.strftime('%Y-%m-%d %H:%M')}; overnight.py exit code {code} (0 = no step failed).", "",
         "## What the run did", "", "| phase | outcome |", "| --- | --- |"] + [f"| {p} | {o} |" for p, o in phases]
    L += ["", "## Steps", "", "| step | status | minutes | criteria not met | datasets |", "| --- | --- | --- | --- | --- |"]
    for sid, *_rest in overnight.STEPS:
        st = states.get(sid)
        if st:
            notmet = sum(1 for v in st.get("verdicts", []) if v["outcome"] == "NOT MET")
            L.append(f"| {sid} | {st.get('status')} | {st.get('seconds', 0) / 60:.1f} | {notmet} | {', '.join(st.get('datasets') or []) or '-'} |")
    failed = [sid for sid, st in states.items() if st.get("status") == "failed"]
    if failed:
        L += ["", "Steps that failed: " + ", ".join(failed) + ". The reason of each is in OVERNIGHT_SUMMARY.md."]
    L += ["", "## Test tables", ""] + table_lines(cfg)
    d = decide_coverage(collect_coverage_data(out))
    L += ["", "## Decision: the default of --allow-unmergeable (1% or 0.2%)", "", f"**{d['verdict']}**", ""] + [f"- {x.strip()}" for x in d["lines"]] + ["", d["why"]]
    L += ["", "## What still needs you", "",
          "- mockB (PRJEB67827): its eight mock samples are fetched by the run; if it has no expected sequences yet, `--make-reference mockB=D6321.refseq.zip` (Zymo's three species) registers them.",
          "- mcm_v46 and mcm_v12: the expected sequences of the 15 bacteria of Table 1 of the paper as a FASTA (`--set-reference mcm_v46=FILE`), for the truth columns.",
          "- V4-V5 (PRJNA971391): the primers, then `--download v4v5`, `--prepare v4v5`.",
          "- The manuscript and release items on the board below.", "", "## The board", "", "```"] + todo(cfg, args) + ["```", "",
          f"Every step's own report: {out / 'OVERNIGHT_SUMMARY.md'}. The log of this run: {out / 'overnight_all.log'} and {out / 'overnight.log'}.", ""]
    path = out / "MORNING_REPORT.md"
    path.write_text("\n".join(L), encoding="utf-8")
    return path


def overnight_all(cfg: dict, args, cfg_path: Path, download_fn=None, prepare_fn=None, run_fn=None, machine_fn=None, discover_fn=None) -> int:
    """Everything that can run without you, in one command: find the folders, check the machine, fetch and prepare the datasets that need no choice (V4 and the two regions
    of the inter-laboratory mock), run every step (a step that covered other datasets is repeated, the rest is kept), and write the morning report. The machine is kept awake."""
    download_fn, prepare_fn, run_fn = download_fn or download, prepare_fn or prepare, run_fn or run_overnight
    machine_fn, discover_fn = machine_fn or check_machine, discover_fn or discover
    started = time.strftime("%Y-%m-%d %H:%M")
    out = Path(os.path.expandvars(cfg["out"]) + ("_quick" if args.quick else ""))
    out.mkdir(parents=True, exist_ok=True)
    phases, code = [], 1
    old, tee = sys.stdout, Tee(sys.stdout, out / "overnight_all.log")
    sys.stdout = tee
    try:
        with keep_awake():
            found = discover_fn(cfg)
            if found:
                apply_discovery(cfg, found)
                save_config(cfg_path, cfg)
            phases.append(("find the folders and FIGARO", ", ".join(f"{k}" for k in sorted(found)) if found else "nothing new"))
            bad = [m for m in machine_fn() if not m[1]]
            if bad:
                phases.append(("check the machine", "NOT READY: " + "; ".join(f"{m[0]} ({m[3]})" for m in bad)))
                print("the machine is not ready:\n" + "\n".join(f"  {m[0]}: {m[2]}  fix: {m[3]}" for m in bad))
            else:
                phases.append(("check the machine", "ok"))
                left_out = {n.strip() for n in (getattr(args, "without", None) or "").split(",") if n.strip()}
                for name in AUTO_ORDER:
                    if name in left_out:
                        phases.append((f"dataset {name}", "left out (--without)"))
                        continue
                    print(f"\n=== {name}: {CATALOG[name]['title']}")
                    phases.append((f"dataset {name}", auto_dataset(cfg, name, cfg_path, download_fn, prepare_fn)))
                    cfg = load_config(cfg_path)                             # prepare registered the dataset in the file
                code = run_fn(cfg, args)
                phases.append(("overnight steps", f"exit code {code}"))
    finally:
        sys.stdout = old
        tee.close()
    path = morning_report(cfg, args, phases, code, started)
    print(f"\nmorning report: {path}")
    return code


def _pt(x, sign: bool = True) -> str:
    return "-" if x is None else (f"{x:+.1f}" if sign else f"{x:.1f}")


def _ratio(x) -> str:
    return "-" if x is None else f"{x:.2f}"


def _ci(d) -> str:
    return "-" if not d else f"{d['mean']:+.1f} ({d['lo']:+.1f} to {d['hi']:+.1f})"


def _across(r: dict) -> str:
    c = r.get("calibration")
    return "-" if not c else f"{c['offset_mean']:+.1f} ({c['offset_lo']:+.1f} to {c['offset_hi']:+.1f})"


def _rank(r: dict) -> str:
    c = r.get("calibration")
    return "-" if not c or c.get("rank") is None else f"{c['rank']:+.2f} ({c['n']})"


def _coverage_cell(cov) -> str:
    if not cov:
        return "-"
    rows = {round(r["level"], 4): r for r in cov["rows"]}
    a, b = rows.get(0.01), rows.get(0.002)
    if not a or not b:
        return "-"
    text = f"{100 * (b['merged'] - a['merged']):+.1f}"
    if cov.get("has_truth") and a.get("strains_found") is not None:
        text += f", strains {b['strains_found'] - a['strains_found']:+d}"
    return text


def render_tables(rows: list, coverage: dict, rhos: dict, pending: list, landscape: dict = None, stability: dict = None) -> list:
    """Markdown: one table for the 16S datasets and one for ITS, each row a dataset with the numbers of the standard battery."""
    out = []
    s16 = [r for r in rows if r["marker"] == "16s"]
    out += ["### 16S rRNA datasets", "",
            "| Dataset | Region | Samples | AmpliTrim setting | Filter: AmpliTrim - DADA2 (points) | Yield estimate - merged (points) | Merged / estimate | Across maxEE settings: offset mean (range) | "
            "Across maxEE settings: rank correlation (settings) | AmpliTrim - no truncation (points, 95% interval) | "
            "AmpliTrim - FIGARO at its maxEE (points, 95% interval) | AmpliTrim default - FIGARO as shipped (points, 95% interval) | Coverage 1% -> 0.2%: merged (points) | "
            "Landscape: points below the best DADA2 setting | Stability: truncation spread (nt) | Per-sample rank |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in s16:
        rho = rhos.get(r["name"])
        rho_text = "-" if rho is None else "%+.2f" % rho
        out.append(f"| {r['name']} | {r['region']} | {r['samples']} | {r['setting']} (maxEE {r['maxEE']}) | {_pt(r['filter_diff'])} | {_pt(r['yield_offset'])} | {_ratio(r['ratio'])} | {_across(r)} | {_rank(r)} | "
                   f"{_ci(r.get('vs_other'))} | {_ci(r.get('vs_figaro'))} | {_ci(r.get('vs_figaro_end'))} | {_coverage_cell(coverage.get(r['name']))} | "
                   f"{_pt((landscape or {}).get(r['name']), sign=False)} | {'-' if (stability or {}).get(r['name']) is None else '%d' % (stability or {})[r['name']]} | {rho_text} |")
    for name in pending:
        out.append(f"| {name} | | | not yet run | | | | | | | | | | | | |")
    its = [r for r in rows if r["marker"] == "its"]
    out += ["", "### ITS datasets", "", "| Dataset | Region | Samples | maxEE | Filter: AmpliTrim - DADA2 (points) | Yield estimate - merged (points) | Merged / estimate | AmpliTrim (maxEE 2) - 200/200 truncation (points, 95% interval) |",
            "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in its:
        for i, x in enumerate(r.get("its_rows", [])):
            first = i == 0
            out.append(f"| {r['name'] if first else ''} | {r['region'] if first else ''} | {r['samples'] if first else ''} | {x['maxEE']} | {_pt(x['filter_diff'])} | {_pt(x['yield_offset'])} | {_ratio(x.get('ratio'))} | "
                       f"{_ci(r.get('vs_other')) if first else ''} |")
    if not its:
        out.append("| not yet run | | | | | | | |")
    truth = [(r, x) for r in rows for x in r.get("truth", [])]
    if truth:                                                                            # the same methods, scored against the expected sequences
        out += ["", "### Datasets with expected sequences (the same methods, scored against the expected sequences)", "",
                "| Dataset | Method | Truncation | maxEE | Pairs merged | Variants recovered | Strains recovered | Reads in identical ASVs | ASVs outside the reference |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for r, x in truth:
            merged = "-" if x["merged"] is None else f"{100 * x['merged']:.1f}%"
            ident = "-" if x["identical"] is None else f"{100 * x['identical']:.2f}%"
            out.append(f"| {r['name']} | {x['method']} | {x['setting']} | {x['maxEE']} | {merged} | {'-' if x['found'] is None else int(x['found'])} | {x['strains'] or '-'} | {ident} | {'-' if x['extra'] is None else int(x['extra'])} |")
    return out


def table_lines(cfg: dict) -> list:
    out = Path(os.path.expandvars(cfg["out"]))
    states = step_states(out)
    rows = []
    for name in ("battery", "extra_datasets"):
        r = states.get(name)
        if r and r.get("status") == "ok":
            rows += r.get("data", [])
    coverage = {d["name"]: d for d in collect_coverage_data(out)}
    rhos = {d["name"]: d["rho"] for d in (states.get("per_sample_rank", {}) or {}).get("data", [])} if states.get("per_sample_rank", {}).get("status") == "ok" else {}
    have = {r["name"] for r in rows}
    pending = [e["short"] for n, e in CATALOG.items() if n not in have and not e.get("dropped")]
    def data_of(step: str) -> list:
        r = states.get(step) or {}
        return r.get("data", []) if r.get("status") == "ok" else []
    landscape = {d["name"]: d["regret"] for d in data_of("landscape")}
    stability = {d["name"]: d["truncation_spread"] for d in data_of("stability")}
    return render_tables(rows, coverage, rhos, pending, landscape, stability)


# ------------------------------------------------------------------------------------------------------------- plan and run
def overnight_args(cfg: dict, args) -> list:
    out = ["--out", os.path.expandvars(cfg["out"]) + ("_quick" if getattr(args, "quick", False) else "")]
    for flag, key in (("--mockA-high", "mockA_high"), ("--mockA-all", "mockA_all"), ("--runinfo", "runinfo"), ("--asv-table", "asv_table"), ("--reference", "reference"),
                      ("--its", "its"), ("--v1v2", "v1v2"), ("--figaro-script", "figaro")):
        if cfg.get(key):
            out += [flag, os.path.expandvars(cfg[key])]
    for name, e in cfg.get("extras", {}).items():
        if not e.get("ready"):
            continue                                                       # a reference registered before the reads were prepared
        out += ["--extra", f"{name}={os.path.expandvars(e['ready'])}" + (f"|{os.path.expandvars(e['reference'])}" if e.get("reference") else "")]
    for flag, attr in (("--only", "only"), ("--skip", "skip")):
        if getattr(args, attr, None):
            out += [flag, getattr(args, attr)]
    for flag, attr in (("--quick", "quick"), ("--fresh", "fresh")):
        if getattr(args, attr, False):
            out.append(flag)
    return out


def plan(cfg: dict, args) -> list:
    import overnight
    lines = ["A. Tests (overnight.py): ready, done, or what they are missing"]
    import io
    import contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        overnight.main(overnight_args(cfg, args) + ["--list"])
    lines += ["  " + ln for ln in buf.getvalue().splitlines() if ln.strip()]
    lines += ["", "B. Data that can be fetched (check_build.py --download NAME, then --prepare NAME)"]
    for name, e in CATALOG.items():
        d = dataset_dirs(cfg, name)
        state = "registered for the tests" if cfg.get("extras", {}).get(name, {}).get("ready") else ("downloaded, not prepared" if d["raw"].is_dir() and any(d["raw"].glob("*_R1.fastq.gz")) else "not downloaded")
        ref = "" if not e["reference"] else (" | reference: " + (cfg["extras"][name]["reference"] if cfg.get("extras", {}).get(name, {}).get("reference") else f"MISSING ({e['reference']})"))
        lines.append(f"  {name:6s} {e['title']}: {state}{ref}\n         gives: {e['gives']}")
    lines += ["", "C. By hand",
              "  - FIGARO: download it (github.com/Zymo-Research/figaro), put the path of figaro.py in amplitrim_paths.json (key figaro), and rerun: it unlocks figaro_mockA, figaro_v1v2 and the FIGARO rows of the others",
              "  - the PRJEB42528 (V1-V2) folder: put its `ready` folder in amplitrim_paths.json (key v1v2): it unlocks figaro_v1v2, v1v2_fulldepth and the V1-V2 rows of coverage_levels, landscape, per_sample_rank",
              "  - the original authors' ITS parameters (agmicrobiomebase): compare them with AmpliTrim's table and with the chimera result of its_chimeras",
              "  - after coverage_levels: decide the default of --allow-unmergeable from its table (see the summary)",
              "  - release: run the CI workflow, push to GitHub, choose a licence, archive for a DOI, check the references of the manuscript"]
    lines += ["", "Next: python validation\\check_build.py --run --quick   (then without --quick)"]
    return lines


def run_overnight(cfg: dict, args) -> int:
    cmd = [sys.executable, str(HERE / "overnight.py")] + overnight_args(cfg, args)
    print("$ " + " ".join(cmd) + "\n", flush=True)
    return subprocess.call(cmd)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", help="paths file (default: amplitrim_paths.json next to the package)")
    ap.add_argument("--init", action="store_true", help="write the paths file with the default folders")
    ap.add_argument("--manifest", help="compare the files of this copy with a manifest of another build (default: validation/build_manifest.txt, shipped in the zip)")
    ap.add_argument("--files", action="store_true", help="print the hash of every file the build id covers (to send for comparison) and stop")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--todo", action="store_true", help="the whole board: tests done and ready, decisions, new data, the manuscript and release items, and what to do next")
    ap.add_argument("--table", action="store_true", help="the test-stage table, one per region (16S, ITS), from the saved results; also written to TEST_TABLE.md in the output folder")
    ap.add_argument("--summary", action="store_true", help="rewrite OVERNIGHT_SUMMARY.md from the saved results of every step (nothing is run)")
    ap.add_argument("--decide", action="store_true", help="apply the rule for the default of --allow-unmergeable to the saved coverage results")
    ap.add_argument("--discover", action="store_true", help="find the `ready` folders (by their number of samples), the input files and figaro.py that the config lacks")
    ap.add_argument("--write", action="store_true", help="with --discover: save what was found in the config")
    ap.add_argument("--overnight", action="store_true", help="EVERYTHING in one run: find the folders, check the machine, fetch and prepare V4 and the inter-laboratory mock, run every step, write MORNING_REPORT.md; keeps Windows awake")
    ap.add_argument("--finish", action="store_true", help="the same as --overnight")
    ap.add_argument("--round4", action="store_true", help="round 4 in phases with breaks: the open questions of the validation (see validation/round4.py and TEST_REGISTER.md, section 9); goes on by itself when nothing is pressed")
    ap.add_argument("--break-seconds", type=int, default=90, help="with --round4: the break between phases (default 90 s; 0 = none); Enter continues, S skips the next phase, Q stops, P holds")
    ap.add_argument("--phases", metavar="A,B,C", help="with --round4: only these phases (A quick analyses, B simulations, C new datasets)")
    ap.add_argument("--without", metavar="NAMES", help="with --overnight: datasets of the run to leave out, comma separated (for example --without mockB); one that is already registered still takes part in the steps")
    ap.add_argument("--download", choices=sorted(CATALOG), metavar="NAME", help=f"fetch a dataset ({', '.join(sorted(CATALOG))})")
    ap.add_argument("--prepare", choices=sorted(CATALOG), metavar="NAME", help="remove the primers of a downloaded dataset and register it")
    ap.add_argument("--make-reference", metavar="NAME=INPUT", help="build the expected sequences of a mock community from full-length 16S sequences (a FASTA, a folder or a zip such as ZymoBIOMICS.STD.refseq) by cutting them between the primers of the dataset, and register them")
    ap.add_argument("--compare-reference", metavar="NAME=INPUT", help="cut a published full-length reference (a zip, folder or FASTA) between the primers of the dataset and compare it with the expected sequences already registered (NAME is mockA or a catalog dataset); changes nothing")
    ap.add_argument("--set-reference", metavar="NAME=FASTA")
    ap.add_argument("--runs-file", help="with --download: a file of run accessions, one per line (needed for v3v4)")
    ap.add_argument("--list-runs", choices=sorted(CATALOG), metavar="NAME", help="list the runs of a dataset's project by kind of sample name (what --select can choose)")
    ap.add_argument("--select", metavar="REGEX", help="with --download: only the runs whose sample name matches (see --list-runs)")
    ap.add_argument("--min-reads", type=int, default=0, help="with --download or --prepare: leave out runs or samples with fewer read pairs")
    ap.add_argument("--samples", metavar="NAMES", help="with --prepare: only these samples (run accessions, comma-separated)")
    ap.add_argument("--no-primers", action="store_true", help="with --prepare: the reads are already free of primers and adapters (used as they are)")
    ap.add_argument("--mixed-orientation", action="store_true", help="with --prepare: some pairs were sequenced the other way round (the describer says so)")
    ap.add_argument("--max-runs", type=int, default=12, help="with --download: how many runs, spread evenly (default 12)")
    ap.add_argument("--dry-run", action="store_true", help="with --download: list the runs and sizes only")
    ap.add_argument("--primer-f"), ap.add_argument("--primer-r")
    ap.add_argument("--run", action="store_true", help="run overnight.py with the paths of the config")
    ap.add_argument("--quick", action="store_true"), ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--only"), ap.add_argument("--skip")
    args = ap.parse_args(argv)
    path = config_path(args)
    if args.init:
        if path.exists():
            print(f"{path} exists; edit it, or delete it to start again")
        else:
            save_config(path, DEFAULT_CONFIG)
            print(f"written {path}: edit the paths (and v1v2, figaro) once")
        return 0
    if args.files:
        print("\n".join(file_lines()))
        return 0
    cfg = load_config(path)
    if args.round4:
        import round4
        return round4.run_round4(cfg, args)
    if args.discover:
        found = discover(cfg)
        print("\n".join(f"  {k:11s} {v[0]}   ({v[1]})" for k, v in found.items()) if found else "  nothing new found (every path in the config exists, or nothing matches)")
        if found and args.write:
            save_config(path, apply_discovery(cfg, found))
            print(f"saved to {path}")
        elif found:
            print("add --write to save these to the config")
        return 0
    if args.finish or args.overnight:
        return overnight_all(cfg, args, path)
    if args.summary:
        return subprocess.call([sys.executable, str(HERE / "overnight.py"), "--out", os.path.expandvars(cfg["out"]) + ("_quick" if args.quick else ""), "--summary"])
    if args.todo:
        print("\n".join(todo(cfg, args)))
        return 0
    if args.table:
        lines = table_lines(cfg)
        target = Path(os.path.expandvars(cfg["out"])) / "TEST_TABLE.md"
        if target.parent.is_dir():
            target.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("\n".join(lines) + (f"\n\nwritten to {target}" if target.parent.is_dir() else ""))
        return 0
    if args.decide:
        d = decide_coverage(collect_coverage_data(Path(os.path.expandvars(cfg["out"]))))
        print("\n".join(["The rule (fixed before the evidence): change the default of --allow-unmergeable from 1% to 0.2% if, on at least two datasets with expected sequences,",
                          "it recovers at least one more strain on one and fewer variants or strains on none, and across all datasets it costs a median of at most 2 points of merged",
                          "pairs and at most 5 on any.", "", f"Decision: {d['verdict']}"] + d["lines"] + [d["why"]]))
        return 0
    if args.list_runs:
        print("\n".join(list_runs(cfg, args.list_runs)))
        return 0
    if args.download:
        return download(cfg, args.download, args)
    if args.prepare:
        return prepare(cfg, args.prepare, args, path)
    if args.compare_reference:
        return compare_reference_for(cfg, args.compare_reference, args)
    if args.make_reference:
        return make_reference_for(cfg, args.make_reference, args, path)
    if args.set_reference:
        return set_reference(cfg, args.set_reference, path)
    if args.run:
        return run_overnight(cfg, args)
    if args.plan:
        print("\n".join(plan(cfg, args)))
        return 0
    print("\n".join(check_build(args)))
    print("\nMachine:")
    bad = 0
    for item, ok, detail, fix in check_machine():
        bad += not ok
        print(f"  {'ok     ' if ok else 'MISSING'} {item:11s} {detail}" + (f"\n          fix: {fix}" if fix else ""))
    print("\nData (from " + (str(path) if path.exists() else "the defaults; run --init to write amplitrim_paths.json") + "):")
    for key, ok, detail in check_data(cfg):
        print(f"  {'ok     ' if ok else ('not set' if ok is None else 'MISSING')} {key:12s} {detail}")
    print("\nNext: python validation\\check_build.py --plan")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
