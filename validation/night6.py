#!/usr/bin/env python
"""Night 6: an independent real dataset, PRJNA933120 (Buetas et al. 2024, BMC Genomics 25:310): the same nine people sequenced as 16S V3-V4 on an Illumina
MiSeq (26 samples) and as full-length 16S on a PacBio Sequel II (the 26 twin runs). The Illumina half tests the RELEASED package on data that was never used
to design it; the PacBio half is described (no AmpliTrim mode exists for it yet). The criteria are in TEST_REGISTER.md, section 25, and the margins below.

  python validation/night6.py --check                                   the machine, the tools, the release package and ENA, then stop
  python validation/night6.py --out C:\\data\\runs\\night6_pilot --package C:\\path\\to\\amplitrim --limit 2      a pilot: two pairs of runs, about half an hour
  python validation/night6.py --out C:\\data\\runs\\night6_PRJNA933120 --package C:\\path\\to\\amplitrim         the whole run, to leave overnight

Steps, in order (--list shows what each needs and whether it is done):
  download    the 26 Illumina pairs and the 26 PacBio runs from ENA (about 1.7 GB), every file checked against ENA's MD5 sum; interrupted transfers resume
  describe    the Illumina reads: lengths, quality, primers, orientation (describe_reads.py)
  prepare     remove the primers (cutadapt), then thin to --pairs-per-sample pairs per sample (5,000, as in the nights)
  amplitrim   the RELEASED package (--package) on the thinned reads: its recommendation, the cut it would choose at the default maxEE 2 and at 2,5, and its REPORT.md
  dada2       DADA2 on the same reads at those three settings (compare_methods.py): the share of pairs passing the filter, merged, and after chimera removal
  longread    the PacBio reads: primers, orientation, lengths, expected errors, retention at each maxEE, quality at the ends against the middle
  summary     the pre-registered checks (P46 to P53), written to NIGHT6_SUMMARY.md

A step that finished is not repeated when you start the command again (--fresh repeats everything, --redo STEP repeats one step); a step that fails is recorded and
the run goes on with the steps that do not need it. A step counts as finished only when its results exist: the dada2 step fails, and says why, if DADA2 gave no result. Needs: Python 3.9+, numpy, cutadapt, vsearch on the PATH, R with dada2, and access to www.ebi.ac.uk and ftp.sra.ebi.ac.uk.

The margins that decide "met" and "not met" are in MARGINS below. Edit them BEFORE the run, not after seeing results.
"""
import argparse
import csv
import datetime
import http.client
import importlib
import json
import shutil
import subprocess
import sys
import time
import traceback
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import download_ena  # noqa: E402
import longread_describe  # noqa: E402

PAIRS_TABLE = HERE / "PRJNA933120_pairs.csv"
ILLUMINA_PRIMERS = ("CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC")        # 341F and 805R (V3-V4), as in the paper
PACBIO_PRIMERS = ("AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT")          # 27F and 1492R (full length), as in the paper
DEFAULT_PAIRS_PER_SAMPLE = 5000                                           # what DADA2 and AmpliTrim got in the nights
STEPS = ["download", "describe", "prepare", "amplitrim", "dada2", "longread", "summary"]
NEEDED_GB = 12                                                             # raw 1.7 GB + prepared reads + DADA2's working files, with room

MARGINS = {
    "filter_estimate": 0.02,           # P46: |AmpliTrim's estimate - DADA2| of the share of pairs passing the filter, at the released and at the default setting
    "rule_gain": 0.02,                 # P47: if the rule relaxes maxEE, the share merged must be at least this much higher than at the default setting
    "both_primers": 0.90,              # P50: share of PacBio reads with both primers
    "median_length": (1430, 1480),     # P51: the median length between the primers of the PacBio reads
    "within_10pct": 0.95,              # P51: share of those reads within 10% of the median
    "distinct_q": 40,                  # P52: distinct quality values, in every PacBio file
    "max_q": 60,                       # P52: the highest quality value, in every PacBio file
    "end_vs_middle": 3.0,              # P53: the first and last 50 bases within this many Q of the middle half
}


class Skip(Exception):
    """A step that cannot run because what it needs is missing."""


class Ctx:
    def __init__(self, args):
        self.args = args
        self.out = Path(args.out)
        self.raw_il, self.raw_pb = self.out / "raw_illumina", self.out / "raw_pacbio"
        self.ready, self.thin = self.out / "illumina_ready", self.out / "illumina_thin"
        self.amp_dir, self.dada2_dir, self.long_dir = self.out / "amplitrim", self.out / "dada2", self.out / "longread"
        self.state_file, self.log_file = self.out / "night6_state.json", self.out / "night6.log"
        self.state = {"steps": {}}
        self.pairs = []
        self.package = None


LOG = None


def say(message="") -> None:
    line = f"{datetime.datetime.now():%H:%M:%S}  {message}" if message else ""
    print(line, flush=True)
    if LOG is not None:
        with open(LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def read_pairs(path: Path, limit: int = 0) -> list:
    """The volunteer / body site / PacBio run / Illumina run table, or `limit` rows of it spread over the list."""
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return download_ena.evenly_spaced(rows, limit) if limit else rows


# ----------------------------------------------------------------------------- downloads: ENA, single-end (the paired ones use download_ena.download)

def ena_listing(accession: str, opener=None, attempts: int = 4, wait: float = 3.0) -> dict:
    """ENA's file report of one run (its FASTQ URLs, MD5 sums and sizes)."""
    opener = opener or urllib.request.urlopen
    url = download_ena.API.format(acc=accession)
    problem = ""
    for attempt in range(1, attempts + 1):
        try:
            with opener(url, timeout=60) as response:
                text = response.read().decode("utf-8", "replace")
            rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
            if rows and rows[0].get("fastq_ftp"):
                return rows[0]
            problem = f"ENA lists no FASTQ file for {accession}"
            break
        except (OSError, http.client.HTTPException) as err:
            problem = f"could not ask ENA about {accession} ({err})"
            if attempt < attempts:
                time.sleep(wait * attempt)
    raise SystemExit(problem)


def single_file(accession: str, opener=None) -> tuple:
    """(url, md5, bytes) of the one FASTQ file of a single-end run."""
    row = ena_listing(accession, opener, wait=download_ena.PAUSE * 3.0)
    urls = row["fastq_ftp"].split(";")
    sums, sizes = (row.get("fastq_md5") or "").split(";"), (row.get("fastq_bytes") or "").split(";")
    exact = [i for i, u in enumerate(urls) if u.endswith(f"/{accession}.fastq.gz")]
    i = exact[0] if exact else (0 if len(urls) == 1 else None)
    if i is None:
        raise SystemExit(f"{accession}: expected one FASTQ file (a single-end run), ENA lists {urls}")
    return urls[i], (sums[i] if i < len(sums) else ""), (int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else 0)


def download_single(accession: str, out_dir: Path, retries: int = 8, log=say, opener=None, retriever=None) -> str:
    """Download one single-end run as <run>.fastq.gz, verified against ENA's MD5 sum; returns the 'md5  filename' line. Resumes interrupted transfers."""
    out_dir.mkdir(parents=True, exist_ok=True)
    url, md5, size = single_file(accession, opener)
    target = out_dir / f"{accession}.fastq.gz"
    if target.exists() and (not md5 or download_ena.md5_of(target) == md5):
        log(f"{accession}: already downloaded, checksum ok")
        return f"{md5}  {target.name}"
    part = out_dir / f"{accession}.part"
    interrupted = False
    for attempt in range(1, retries + 1):
        try:
            if retriever is None:
                download_ena.fetch_resumable("https://" + url, part, size)
            else:
                retriever("https://" + url, part)
        except download_ena.ServerAnswer as err:
            part.unlink(missing_ok=True)
            raise SystemExit(f"{accession}: ENA answered with a short message instead of the {size:,} byte file ({err}); the run may be missing or restricted at ENA")
        except (OSError, http.client.HTTPException) as err:
            interrupted = True
            log(f"{accession}: interrupted ({err}); resuming (attempt {attempt} of {retries})")
            continue
        interrupted = False
        if not md5 or download_ena.md5_of(part) == md5:
            part.replace(target)
            log(f"{accession}: downloaded, checksum ok")
            return f"{md5}  {target.name}"
        log(f"{accession}: checksum mismatch, starting again (attempt {attempt} of {retries})")
        part.unlink(missing_ok=True)
    if not interrupted:
        part.unlink(missing_ok=True)
    raise SystemExit(f"{accession}: no verified download after {retries} attempts" + ("; the partial file is kept, so running the same command again resumes it" if interrupted else ""))


# ----------------------------------------------------------------------------- the released package

def resolve_package(given: Path):
    """The folder that contains the `amplitrim` package, from what was given: that folder, the package folder itself, or the folder above the repository."""
    given = Path(given)
    if (given / "amplitrim" / "__init__.py").is_file():
        return given
    if (given / "__init__.py").is_file() and given.name == "amplitrim":
        return given.parent
    if (given / "amplitrim" / "amplitrim" / "__init__.py").is_file():
        return given / "amplitrim"
    return None


def load_release(package: Path):
    """Make the released package the one that `import amplitrim` finds, and return its pipeline module."""
    root = resolve_package(package) if package else None
    if root is None:
        raise SystemExit(f"The folder given with --package ({package}) has no 'amplitrim' package in it: give the folder of the clone of the software repository.")
    sys.path.insert(0, str(Path(root).resolve()))
    for m in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:
        del sys.modules[m]
    pl = importlib.import_module("amplitrim.pipeline")
    if not hasattr(pl, "WARN_ABOVE"):
        raise SystemExit(f"That package ({Path(pl.__file__).parent}) is not the release without a cap (it has no WARN_ABOVE): give --package the folder of the clone of the software repository, up to date.")
    return pl


def ee_text(f: float, r: float) -> str:
    return "%g" % f if f == r else "%g,%g" % (f, r)


def run_amplitrim(pl, reads: Path, work: Path, report_dir: Path, say=say) -> dict:
    """The released package on `reads`: its own recommendation (the maxEE rule), and the cut it chooses at maxEE 2 and at 2,5."""
    from amplitrim.config import ee_pair
    from amplitrim.report import write_outputs
    out = {"package": str(Path(pl.__file__).parent), "version": getattr(importlib.import_module("amplitrim"), "__version__", "?")}
    for name, ee in (("released", None), ("default", 2.0), ("g25", (2.0, 5.0))):
        say(f"AmpliTrim, {name} ...")
        (work / name).mkdir(parents=True, exist_ok=True)
        res = pl.run(pl.Settings(input_dir=reads, work_dir=work / name, **({} if ee is None else {"max_ee": ee})))
        rec = res.recommendation
        if rec is None:
            out[name] = {"error": "no recommendation"}
            continue
        ef, er = ee_pair(res.cfg.max_ee)
        out[name] = {"trim_f": int(rec.trim_f), "trim_r": int(rec.trim_r), "max_ee_f": float(ef), "max_ee_r": float(er), "max_ee": ee_text(float(ef), float(er)),
                     "short_lost": getattr(res, "short_lost", None), "floor": getattr(res, "floor", None)}
        if name == "released":
            write_outputs(res, report_dir)
        say(f"  {name}: truncLenF {rec.trim_f}, truncLenR {rec.trim_r}, maxEE {out[name]['max_ee']}")
    return out


# ----------------------------------------------------------------------------- DADA2 through compare_methods.py

def dada2_settings(amp: dict) -> list:
    """[(name, truncF, truncR, maxEE text)] of the settings DADA2 runs, with identical settings run once (released kept under its own name last)."""
    seen, out = {}, []
    for name in ("default", "g25", "released"):
        s = amp.get(name)
        if not s or "trim_f" not in s:
            continue
        key = (s["trim_f"], s["trim_r"], s["max_ee"])
        if key in seen:
            continue
        seen[key] = name
        out.append((name, *key))
    return out


def compare_command(thin: Path, out: Path, settings: list, figaro: Path = None, min_overlap: int = 12) -> list:
    cmd = [sys.executable, str(HERE / "compare_methods.py"), str(thin), "--out", str(out), "--amplitrim-max-ee", "2", "--min-overlap", str(min_overlap), "--baseline", settings[0][0]]
    for name, f, r, ee in settings:
        cmd += ["--setting", f"{name}={f},{r},{ee}"]
    if figaro:
        cmd += ["--figaro-script", str(figaro)]
    return cmd


# ----------------------------------------------------------------------------- the pre-registered checks (TEST_REGISTER.md, section 25)

def fraction(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v / 100.0 if v > 1.5 else v


def read_table(path: Path) -> dict:
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}


def illumina_quality_values(path: Path) -> list:
    """The number of distinct quality values of each sample, from the `Qvals` column of describe_illumina.txt; [] when the file or the column is missing."""
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    head = next((i for i, line in enumerate(lines) if line.startswith("sample") and "Qvals" in line), None)
    if head is None:
        return []
    column = lines[head].split().index("Qvals")                 # the first four columns (sample, reads, len, Qvals) never contain a space
    out = []
    for line in lines[head + 1:]:
        if not line.strip():
            break
        parts = line.split()
        if len(parts) > column and parts[column].isdigit():
            out.append(int(parts[column]))
    return out


def simplified_qualities(values: list) -> bool:
    """True when every sample has two or fewer distinct quality values: the per-base quality is missing, as in the SRA 'lite' format, and DADA2 cannot learn an error model."""
    return bool(values) and max(values) <= 2


def dada2_problems(rows: dict, names: list) -> list:
    """['name: why'] for every setting that was to be run by DADA2 and has no measured result: no row at all, an error that compare_methods.py recorded, or no share of pairs passing."""
    out = []
    for n in names:
        r = rows.get(n)
        if r is None:
            out.append(f"{n}: not in the table")
        elif r.get("error"):
            out.append(f"{n}: {r['error']}")
        elif fraction(r.get("pass")) is None or fraction(r.get("merged")) is None:
            out.append(f"{n}: DADA2 gave no share of pairs passing or merged")
    return out


def evaluate_illumina(rows: dict, amp: dict, margins=MARGINS) -> list:
    """[(id, verdict, text)] for P46 to P49 from compare_methods.tsv (rows by name) and the released package's results."""
    out = []

    def row(name):
        """The DADA2 row of a setting: under its own name, or under the name of the identical setting that was run instead."""
        if name in rows:
            return rows[name]
        s = amp.get(name)
        if not s or "trim_f" not in s:
            return None
        key = (s["trim_f"], s["trim_r"], s["max_ee"])
        return next((rows[n] for n, f, r, e in dada2_settings(amp) if (f, r, e) == key and n in rows), None)

    gaps, missing = [], []
    for name in ("released", "default"):
        r = row(name)
        est, got = (fraction(r.get("est_pass")), fraction(r.get("pass"))) if r else (None, None)
        if est is None or got is None:
            missing.append(name)
        else:
            gaps.append((name, est - got))
    if missing:
        out.append(("P46", "not decided", f"no estimate and DADA2 measurement for: {', '.join(missing)}"))
    else:
        ok = all(abs(g) <= margins["filter_estimate"] for _, g in gaps)
        out.append(("P46", "MET" if ok else "NOT MET", "estimate minus DADA2 of the share of pairs passing the filter: " + "; ".join(f"{n} {100 * g:+.1f} points" for n, g in gaps) + f" (margin {100 * margins['filter_estimate']:.0f} points)"))
    rel, dflt = amp.get("released"), amp.get("default")
    if not rel or not dflt or "trim_f" not in rel or "trim_f" not in dflt:
        out.append(("P47", "not decided", "the released package or the default setting gave no recommendation"))
    elif (rel["trim_f"], rel["trim_r"], rel["max_ee"]) == (dflt["trim_f"], dflt["trim_r"], dflt["max_ee"]):
        out.append(("P47", "not applicable", f"the rule kept the default setting ({rel['trim_f']}/{rel['trim_r']}, maxEE {rel['max_ee']})"))
    else:
        a, b = row("released"), row("default")
        ma, mb = (fraction(a.get("merged")) if a else None), (fraction(b.get("merged")) if b else None)
        if ma is None or mb is None:
            out.append(("P47", "not decided", "no merged share for the released or the default setting"))
        else:
            ok = ma - mb >= margins["rule_gain"]
            out.append(("P47", "MET" if ok else "NOT MET", f"merged {100 * ma:.1f}% at the released setting ({rel['trim_f']}/{rel['trim_r']}, maxEE {rel['max_ee']}) against {100 * mb:.1f}% at the default "
                                                           f"({dflt['trim_f']}/{dflt['trim_r']}, maxEE {dflt['max_ee']}): {100 * (ma - mb):+.1f} points (needs at least {100 * margins['rule_gain']:+.0f})"))
    shares = []
    for name, label in (("default", "default (maxEE 2)"), ("g25", "2,5"), ("released", "released")):
        r = row(name)
        v = fraction(r.get("merged")) if r else None
        shares.append(f"{label} {'-' if v is None else '%.1f%%' % (100 * v)}")
    out.append(("P48", "descriptive", "share of pairs merged by DADA2: " + ", ".join(shares)))
    dev = next((r for k, r in rows.items() if k.startswith("AmpliTrim")), None)
    if dev and dflt and "trim_f" in dflt:
        same = (int(float(dev["truncF"])), int(float(dev["truncR"]))) == (dflt["trim_f"], dflt["trim_r"])
        out.append(("P49", "descriptive", f"the development build chose {dev['truncF']}/{dev['truncR']} at maxEE 2, the released package {dflt['trim_f']}/{dflt['trim_r']}: " + ("the same cut" if same else "DIFFERENT cuts")))
    else:
        out.append(("P49", "not decided", "no development-build row in the table"))
    return out


def evaluate_pacbio(rows: list, margins=MARGINS) -> list:
    """[(id, verdict, text)] for P50 to P53 from the per-file rows of longread_describe (the pooled summary)."""
    if not rows:
        return [(i, "not decided", "no PacBio description") for i in ("P50", "P51", "P52", "P53")]
    p = longread_describe.pooled(rows)
    out = []
    ok = p["share_both"] is not None and p["share_both"] >= margins["both_primers"]
    out.append(("P50", "as expected" if ok else "NOT as expected", f"{100 * (p['share_both'] or 0):.1f}% of {p['reads']:,} reads have both primers (expected at least {100 * margins['both_primers']:.0f}%)"))
    lo, hi = margins["median_length"]
    med, within = p["len_median_of_medians"], p["share_within_10pct"]
    ok = med is not None and lo <= med <= hi and within is not None and within >= margins["within_10pct"]
    out.append(("P51", "as expected" if ok else "NOT as expected", f"median length between the primers {0 if med is None else med:.0f} nt (expected {lo} to {hi}); {100 * (within or 0):.1f}% within 10% of the median (expected at least {100 * margins['within_10pct']:.0f}%)"))
    ok = p["distinct_q_min"] is not None and p["distinct_q_min"] >= margins["distinct_q"] and p["max_q_min"] is not None and p["max_q_min"] >= margins["max_q"]
    out.append(("P52", "as expected" if ok else "NOT as expected", f"fewest distinct quality values in a file {p['distinct_q_min']} (expected at least {margins['distinct_q']}); lowest maximum Q{p['max_q_min']} (expected at least Q{margins['max_q']})"
                                                             + ("" if ok else "; the SRA copy may have simplified qualities, so expected errors are not reliable")))
    if p["q_mid"] is None:
        out.append(("P53", "not decided", "no quality at the ends and the middle"))
    else:
        d = max(abs(p["q_start"] - p["q_mid"]), abs(p["q_end"] - p["q_mid"]))
        out.append(("P53", "as expected" if d <= margins["end_vs_middle"] else "NOT as expected", f"first 50 bases Q{p['q_start']:.1f}, middle half Q{p['q_mid']:.1f}, last 50 bases Q{p['q_end']:.1f} (expected within {margins['end_vs_middle']:g})"))
    return out


# ----------------------------------------------------------------------------- the steps

def need(path: Path, what: str, pattern: str = "*.fastq.gz"):
    if not path.is_dir() or not any(path.glob(pattern)):
        raise Skip(f"{what} is missing: {path}")


def step_download(c: Ctx) -> str:
    failures, done = [], 0
    il = [r["illumina_run"] for r in c.pairs]
    pb = [r["pacbio_run"] for r in c.pairs]
    for acc in il:
        try:
            lines = download_ena.download(acc, c.raw_il, retries=c.args.retries, log=say)
            with open(c.raw_il / "fastq_checksums.md5", "a", encoding="utf-8") as fh:
                fh.write("\n".join(lines) + "\n")
            done += 1
        except SystemExit as err:
            failures.append(f"{acc}: {err}")
            say(f"FAILED {acc}: {err}")
    for acc in pb:
        try:
            line = download_single(acc, c.raw_pb, retries=c.args.retries)
            with open(c.raw_pb / "fastq_checksums.md5", "a", encoding="utf-8") as fh:
                fh.write(line + "\n")
            done += 1
        except SystemExit as err:
            failures.append(f"{acc}: {err}")
            say(f"FAILED {acc}: {err}")
    with open(c.out / "samples.tsv", "w", encoding="utf-8") as fh:
        fh.write("volunteer\tsite\tillumina_run\tpacbio_run\n")
        for r in c.pairs:
            fh.write(f"{r['volunteer']}\t{r['site']}\t{r['illumina_run']}\t{r['pacbio_run']}\n")
    if failures:
        raise SystemExit(f"{len(failures)} of {len(il) + len(pb)} downloads failed (start the command again: finished files are kept and interrupted ones resume):\n  " + "\n  ".join(failures[:10]))
    return f"{done} runs downloaded and checked"


def step_describe(c: Ctx) -> str:
    need(c.raw_il, "the Illumina reads", "*_R1.fastq.gz")
    cmd = [sys.executable, str(HERE / "describe_reads.py"), str(c.raw_il), "--primer", ILLUMINA_PRIMERS[0], "--primer", ILLUMINA_PRIMERS[1], "--lengths", "--csv", str(c.out / "describe_illumina.csv")]
    res = subprocess.run(cmd, capture_output=True, text=True)
    (c.out / "describe_illumina.txt").write_text(res.stdout + res.stderr, encoding="utf-8")
    if res.returncode != 0:
        raise SystemExit(f"describe_reads.py failed: {(res.stderr or res.stdout)[-400:]}")
    return "describe_illumina.txt written"


def step_prepare(c: Ctx) -> str:
    need(c.raw_il, "the Illumina reads", "*_R1.fastq.gz")
    have = len(list(c.raw_il.glob("*_R1.fastq.gz")))
    if have < len(c.pairs) and not c.args.allow_partial:
        raise Skip(f"only {have} of {len(c.pairs)} Illumina runs are downloaded (run the download step again, or give --allow-partial)")
    load_release(c.package)
    from real_data import thin_reads, trim_reads
    if c.ready.is_dir():
        shutil.rmtree(c.ready)
    if c.thin.is_dir():
        shutil.rmtree(c.thin)
    say(f"removing the primers {ILLUMINA_PRIMERS[0]} and {ILLUMINA_PRIMERS[1]} (cutadapt) ...")
    trim_reads(c.raw_il, c.ready, primer_f=ILLUMINA_PRIMERS[0], primer_r=ILLUMINA_PRIMERS[1], mixed=c.args.mixed_orientation)
    say(f"thinning to at most {c.args.pairs_per_sample} pairs per sample ...")
    thin_reads(c.ready, c.thin, c.args.pairs_per_sample)
    return f"{len(list(c.thin.glob('*_R1.fastq.gz')))} samples ready in {c.thin}"


def step_amplitrim(c: Ctx) -> str:
    need(c.thin, "the thinned Illumina reads", "*_R1.fastq.gz")
    pl = load_release(c.package)
    c.amp_dir.mkdir(parents=True, exist_ok=True)
    result = run_amplitrim(pl, c.thin, c.amp_dir / "work", c.amp_dir / "report")
    (c.amp_dir / "amplitrim_results.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
    shutil.rmtree(c.amp_dir / "work", ignore_errors=True)
    r = result.get("released", {})
    return f"released package {result['version']}: " + (f"truncLenF {r['trim_f']}, truncLenR {r['trim_r']}, maxEE {r['max_ee']}" if "trim_f" in r else "no recommendation")


def step_dada2(c: Ctx) -> str:
    need(c.thin, "the thinned Illumina reads", "*_R1.fastq.gz")
    q = illumina_quality_values(c.out / "describe_illumina.txt")
    if simplified_qualities(q):
        raise Skip(f"the Illumina reads have at most {max(q)} distinct quality values in each of {len(q)} samples (a simplified-quality copy, as in the SRA 'lite' format): "
                   "DADA2 cannot learn an error model from them. Get files with the original quality scores.")
    f = c.amp_dir / "amplitrim_results.json"
    if not f.is_file():
        raise Skip("the amplitrim step has not finished: DADA2 needs the settings it chose")
    settings = dada2_settings(json.loads(f.read_text(encoding="utf-8")))
    if not settings:
        raise Skip("the released package gave no setting to run")
    c.dada2_dir.mkdir(parents=True, exist_ok=True)
    cmd = compare_command(c.thin, c.dada2_dir, settings, c.args.figaro_script)
    say("running: " + " ".join(f'"{x}"' if " " in x else x for x in cmd))
    with open(c.out / "dada2_run.log", "w", encoding="utf-8") as fh:
        res = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=str(ROOT))
    if res.returncode != 0 or not (c.dada2_dir / "compare_methods.tsv").is_file():
        tail = (c.out / "dada2_run.log").read_text(encoding="utf-8", errors="replace")[-600:]
        raise SystemExit(f"compare_methods.py failed (exit {res.returncode}); the end of {c.out / 'dada2_run.log'}:\n{tail}")
    problems = dada2_problems(read_table(c.dada2_dir / "compare_methods.tsv"), [n for n, *_ in settings])
    if problems:
        raise SystemExit(f"DADA2 gave no result for {len(problems)} of {len(settings)} setting(s); compare_methods.py itself ended without an error:\n  " + "\n  ".join(p[:400] for p in problems[:4])
                         + f"\nThe whole table is in {c.dada2_dir / 'compare_methods.txt'} and the output of R in {c.out / 'dada2_run.log'}. Fix the cause, then repeat just this step with --redo dada2.")
    return f"DADA2 ran at {len(settings)} setting(s): " + ", ".join(f"{n} {f}/{r} maxEE {e}" for n, f, r, e in settings)


def step_longread(c: Ctx) -> str:
    need(c.raw_pb, "the PacBio reads")
    rows = longread_describe.describe_folder(c.raw_pb, c.long_dir, PACBIO_PRIMERS[0], PACBIO_PRIMERS[1], log=say)
    p = longread_describe.pooled(rows)
    return f"{p['files']} files, {p['reads']:,} reads, {100 * (p['share_both'] or 0):.1f}% with both primers"


def read_longread_rows(folder: Path) -> list:
    path = folder / "longread_describe.tsv"
    if not path.is_file():
        return []
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            row = {}
            for k, v in r.items():
                if v == "" or v is None:
                    continue
                row[k] = v if k == "sample" else (int(v) if v.lstrip("-").isdigit() else float(v))
            rows.append(row)
    return rows


def step_summary(c: Ctx) -> str:
    c.out.mkdir(parents=True, exist_ok=True)
    L = [f"# Night 6: PRJNA933120 (Buetas et al. 2024)", "", f"Run on {datetime.date.today()}; folder {c.out}. Criteria: TEST_REGISTER.md, section 25. Margins: {json.dumps({k: v for k, v in MARGINS.items()})}", ""]
    amp_file = c.amp_dir / "amplitrim_results.json"
    amp = json.loads(amp_file.read_text(encoding="utf-8")) if amp_file.is_file() else None
    L += ["## Illumina, 16S V3-V4: the released package against DADA2", ""]
    if amp:
        L.append(f"The released package (version {amp['version']}, {amp['package']}) on the reads, thinned to at most {c.args.pairs_per_sample} pairs per sample:")
        L += [f"- {n}: " + (f"truncLenF {s['trim_f']}, truncLenR {s['trim_r']}, maxEE {s['max_ee']}" if "trim_f" in s else s.get("error", "-")) for n, s in ((k, amp[k]) for k in ("released", "default", "g25") if k in amp)]
        L.append("")
    q = illumina_quality_values(c.out / "describe_illumina.txt")
    if q:
        L.append(f"Distinct quality values per sample in the Illumina reads: {min(q)} to {max(q)} ({len(q)} samples).")
        L.append("")
    table = c.dada2_dir / "compare_methods.tsv"
    if simplified_qualities(q):
        L += [f"- **P46 to P49 cannot be decided on this copy of the reads.** Every sample has at most {max(q)} distinct quality values, so the per-base quality is missing (as in the SRA 'lite' format) "
              "and DADA2 cannot learn an error model from it. AmpliTrim's recommendation above was computed from the same scores, so it is no more informative than they are, and P47 'does not apply' "
              "carries no evidence about the rule. The criteria of section 25 stay as they are for a copy of these reads, or another dataset, with the original quality scores."]
    elif amp and table.is_file():
        rows = read_table(table)
        problems = dada2_problems(rows, [n for n, *_ in dada2_settings(amp)])
        if problems:
            L += [f"- **DADA2 gave no result for {len(problems)} setting(s)**, so P46 and P48 cannot be decided:"] + [f"  - {p[:300]}" for p in problems[:4]] + [""]
        verdicts = evaluate_illumina(rows, amp)
        L += [f"- **{i}** {v}: {t}" for i, v, t in verdicts]
        L += ["", "The whole DADA2 table is in `dada2/compare_methods.txt`."]
    else:
        L.append("Not decided: " + ("the AmpliTrim step has not finished." if not amp else "DADA2 has not finished (dada2/compare_methods.tsv is missing)."))
    L += ["", "## PacBio, full-length 16S: description only", ""]
    rows = read_longread_rows(c.long_dir)
    if rows:
        verdicts = evaluate_pacbio(rows)
        L += [f"- **{i}** {v}: {t}" for i, v, t in verdicts]
        L += ["", "The table of every file, with what is unusual, is in `longread/LONGREAD_DESCRIBE.md`."]
    else:
        L.append("Not decided: the longread step has not finished.")
    L += ["", "## Steps", ""]
    for name in STEPS[:-1]:
        s = c.state["steps"].get(name)
        L.append(f"- {name}: " + (f"{s['status']} in {s['seconds'] // 60} min" + (f" ({s['detail'].splitlines()[0]})" if s.get("detail") else "") if s else "not run"))
    (c.out / "NIGHT6_SUMMARY.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    say("\n".join(L))
    return "NIGHT6_SUMMARY.md written"


STEP_FUNCTIONS = {"download": step_download, "describe": step_describe, "prepare": step_prepare, "amplitrim": step_amplitrim, "dada2": step_dada2, "longread": step_longread, "summary": step_summary}


def save_state(c: Ctx) -> None:
    c.out.mkdir(parents=True, exist_ok=True)
    c.state_file.write_text(json.dumps(c.state, indent=1), encoding="utf-8")


def run_step(c: Ctx, name: str, force: bool = False) -> None:
    done = c.state["steps"].get(name, {})
    if done.get("status") == "done" and not c.args.fresh and not force and name != "summary":
        say(f"{name}: finished earlier, skipped (--fresh repeats it)")
        return
    say(f"=== {name} ===")
    t0 = time.time()
    try:
        detail, status = STEP_FUNCTIONS[name](c), "done"
    except Skip as err:
        detail, status = str(err), "skipped"
    except SystemExit as err:
        detail, status = str(err), "failed"
    except Exception:                                                  # a bug or a tool that crashed: record it and go on with the other steps
        detail, status = traceback.format_exc(limit=6), "failed"
    c.state["steps"][name] = {"status": status, "seconds": int(time.time() - t0), "detail": (detail or "")[:3000], "finished": datetime.datetime.now().isoformat(timespec="seconds")}
    save_state(c)
    say(f"{name}: {status}" + (f" - {detail.splitlines()[0]}" if detail else ""))
    if status == "failed" and detail:
        say(detail)


# ----------------------------------------------------------------------------- --check and --list

def tool_version(cmd: list, timeout: int = 60):
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as err:
        return None, str(err)
    text = (res.stdout or res.stderr).strip().splitlines()
    return (res.returncode == 0), (text[0] if text else "")


def check_environment(c: Ctx, opener=None) -> list:
    """[(what, ok, detail)] for everything the run needs."""
    out = []
    out.append(("Python 3.9 or later", sys.version_info >= (3, 9), sys.version.split()[0]))
    try:
        import numpy
        out.append(("numpy", True, numpy.__version__))
    except ImportError:
        out.append(("numpy", False, "pip install numpy"))
    ok, v = tool_version([sys.executable, "-m", "cutadapt", "--version"])
    out.append(("cutadapt", bool(ok), v if ok else "pip install cutadapt"))
    ok, v = tool_version(["vsearch", "--version"])
    out.append(("vsearch on the PATH", bool(ok), v if ok else "vsearch is not on the PATH (see the README, Install)"))
    try:
        from rtools import NOT_FOUND, find_rscript
        rs = find_rscript()
    except Exception as err:                                           # noqa: BLE001
        rs, NOT_FOUND = None, None
        out.append(("Rscript", False, f"could not look for it ({err})"))
    else:
        if rs in (None, NOT_FOUND):
            out.append(("Rscript", False, "R was not found"))
        else:
            ok, v = tool_version([str(rs), "-e", "cat(as.character(packageVersion('dada2')))"])
            out.append(("R with dada2", bool(ok), f"dada2 {v}" if ok else f"dada2 is not installed for {rs}"))
    try:
        load_release(c.package)
        import amplitrim
        out.append(("released package (--package)", True, f"version {getattr(amplitrim, '__version__', '?')} at {Path(sys.modules['amplitrim'].__file__).parent}"))
    except SystemExit as err:
        out.append(("released package (--package)", False, str(err)))
    try:
        pairs = read_pairs(Path(c.args.pairs_table))
        out.append(("table of the 26 pairs", len(pairs) == 26, f"{len(pairs)} rows in {c.args.pairs_table}"))
    except OSError as err:
        out.append(("table of the 26 pairs", False, str(err)))
    probe = Path(c.args.out)
    while not probe.exists() and probe != probe.parent:
        probe = probe.parent
    free = shutil.disk_usage(probe).free / 1e9
    out.append((f"free disk space at {probe} (needs about {NEEDED_GB} GB)", free >= NEEDED_GB, f"{free:.0f} GB"))
    try:
        ena_listing(c.pairs[0]["illumina_run"] if c.pairs else "SRR23380961", opener, attempts=1)
        out.append(("ENA reachable (www.ebi.ac.uk)", True, "answered"))
    except SystemExit as err:
        out.append(("ENA reachable (www.ebi.ac.uk)", False, str(err)))
    return out


def list_steps(c: Ctx) -> None:
    ready = {"download": "internet (ENA)", "describe": "the Illumina download", "prepare": "the Illumina download, cutadapt", "amplitrim": "the prepared reads, --package",
             "dada2": "the AmpliTrim step, vsearch, R with dada2", "longread": "the PacBio download", "summary": "whatever has finished"}
    for name in STEPS:
        s = c.state["steps"].get(name)
        say(f"  {name:10s} needs: {ready[name]:46s} " + (f"[{s['status']}, {s['seconds'] // 60} min]" if s else "[not run]"))


# ----------------------------------------------------------------------------- main

def main(argv=None) -> int:
    global LOG
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("night6_PRJNA933120"), help="the folder for everything this run makes (use a different one for a pilot)")
    ap.add_argument("--package", type=Path, default=ROOT.parent / "amplitrim", help="the folder of the clone of the software repository, the package under test (default: next to this archive)")
    ap.add_argument("--pairs-table", type=Path, default=PAIRS_TABLE, help="the 26 volunteer/site pairs of runs (default: PRJNA933120_pairs.csv next to this script)")
    ap.add_argument("--limit", type=int, default=0, metavar="N", help="a pilot: only N pairs of runs, spread over the list")
    ap.add_argument("--pairs-per-sample", type=int, default=DEFAULT_PAIRS_PER_SAMPLE, help=f"pairs per sample given to AmpliTrim and DADA2 (default {DEFAULT_PAIRS_PER_SAMPLE}, as in the nights)")
    ap.add_argument("--only", help="comma-separated steps to run")
    ap.add_argument("--skip", help="comma-separated steps to leave out")
    ap.add_argument("--redo", help="comma-separated steps to repeat even though they finished; the others are left as they are, and the summary is written again")
    ap.add_argument("--fresh", action="store_true", help="repeat steps that already finished")
    ap.add_argument("--check", action="store_true", help="check the machine, the tools, the release package and ENA, then stop")
    ap.add_argument("--list", action="store_true", help="show the steps and which have finished, then stop")
    ap.add_argument("--figaro-script", type=Path, help="path to figaro.py: also run DADA2 at FIGARO's recommendation (optional)")
    ap.add_argument("--mixed-orientation", action="store_true", help="if describe_illumina.txt shows pairs in both orientations")
    ap.add_argument("--allow-partial", action="store_true", help="go on with fewer Illumina runs than the table lists")
    ap.add_argument("--retries", type=int, default=8, help="attempts per file (default 8)")
    a = ap.parse_args(argv)
    c = Ctx(a)
    c.package = a.package
    try:
        c.pairs = read_pairs(Path(a.pairs_table), a.limit)
    except OSError as err:
        print(f"Cannot read the table of pairs: {err}")
        return 2
    if a.check:
        results = check_environment(c)
        for what, ok, detail in results:
            print(f"  [{'ok' if ok else 'NOT READY'}] {what}: {detail}")
        bad = [w for w, ok, _ in results if not ok]
        print("\nEverything is ready." if not bad else f"\n{len(bad)} thing(s) to fix before the run: " + "; ".join(bad))
        return 0 if not bad else 1
    c.out.mkdir(parents=True, exist_ok=True)
    LOG = c.log_file
    if c.state_file.is_file() and not a.fresh:
        c.state = json.loads(c.state_file.read_text(encoding="utf-8"))
    if a.fresh:
        c.state = {"steps": {}}
    if a.list:
        list_steps(c)
        return 0
    wanted = [s.strip() for s in a.only.split(",")] if a.only else list(STEPS)
    skipped = [s.strip() for s in a.skip.split(",")] if a.skip else []
    redo = [s.strip() for s in a.redo.split(",")] if a.redo else []
    if redo and not a.only:
        wanted = redo + ["summary"]
    unknown = [s for s in wanted + skipped + redo if s not in STEPS]
    if unknown:
        print(f"Unknown step(s): {', '.join(unknown)}. The steps are: {', '.join(STEPS)}")
        return 2
    say(f"Night 6, {len(c.pairs)} pairs of runs" + (f" (a pilot of {a.limit})" if a.limit else "") + f", folder {c.out}")
    t0 = time.time()
    uses = {"download": ["describe", "prepare", "longread"], "prepare": ["amplitrim", "dada2"], "amplitrim": ["dada2"]}
    for name in redo:
        stale = [x for x in uses.get(name, []) if x not in redo]
        if stale:
            say(f"Note: {', '.join(stale)} use what {name} makes and are not repeated; add them to --redo if {name} changed.")
    for name in STEPS:
        if name in wanted and name not in skipped:
            run_step(c, name, force=name in redo)
    say(f"Finished in {(time.time() - t0) / 60:.0f} minutes. Read {c.out / 'NIGHT6_SUMMARY.md'} first.")
    failed = [n for n, s in c.state["steps"].items() if s["status"] == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
