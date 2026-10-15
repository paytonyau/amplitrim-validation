"""Night 7: the long-read checks on the PacBio Kinnex-16S data of the ATCC MSA-1003 mock community, in one run.

Steps (--list shows which have finished, --redo STEP repeats one, --check tests the machine and the files first):
  genomes         downloads the 20 strains' genomes from NCBI (resumable, each checked against NCBI's MD5 sum)
  amplicons       genome_amplicons.py: the expected 16S variants and the expected share of reads of each strain
  filter_monomer  test L4 on three Monomer files: AmpliTrim's estimate of what a maxEE filter keeps against DADA2's own filter, on the same reads
  filter_kinnex   test L4 on three Kinnex files (the first N reads of each)
  summary         NIGHT7_SUMMARY.md: what ran, the gate, the warnings, and what to send back

    python night7.py --check
    python night7.py --out C:\\data\\runs\\night7_longread

It needs the files of the earlier steps in --data: sq2_monomer_fastq (the Monomer files, at least the three named below) and sq2_kinnex_fastq (the three Kinnex files).
Needs Python 3.9+, numpy, cutadapt (pip install cutadapt) and R with DADA2.
"""

import argparse
import datetime
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import download_ena  # noqa: E402
import genome_amplicons as ga  # noqa: E402
import longread_filter_check as lf  # noqa: E402
from rtools import NOT_FOUND, find_rscript  # noqa: E402

PRIMERS = ("AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT")                  # 27F and 1492R
STEPS = ["genomes", "amplicons", "filter_monomer", "filter_kinnex", "summary"]
GENOME_TABLE = HERE / "msa1003_genomes.tsv"
DEFAULT_DATA = Path(r"C:\data\longread\pacbio_kinnex16s")
MONOMER_PAIRS = ("Fwd_07--Kinnex16S_Rev_30", "Fwd_09--Kinnex16S_Rev_29", "Fwd_04--Kinnex16S_Rev_22")     # the smallest, the median and the largest Monomer file
NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
NEEDED_GB = 2
GRID = lf.DEFAULT_GRID


class Skip(Exception):
    """A step that cannot run because what it needs is missing."""


class Ctx:
    def __init__(self, args):
        self.args = args
        self.out, self.data = Path(args.out), Path(args.data)
        self.genomes = self.data / "ref" / "genomes"
        self.monomer, self.kinnex = self.data / "sq2_monomer_fastq", self.data / "sq2_kinnex_fastq"
        self.amplicons = self.out / "genome_amplicons"
        self.monomer_in = self.out / "l4_monomer_input"
        self.monomer_out, self.kinnex_out = self.out / "filter_check_monomer", self.out / "filter_check_kinnex"
        self.state_file, self.log_file = self.out / "night7_state.json", self.out / "night7.log"
        self.state = {"steps": {}}
        self.table = Path(args.table)


LOG = None


def say(message="") -> None:
    line = f"{datetime.datetime.now():%H:%M:%S}  {message}" if message else ""
    print(line, flush=True)
    if LOG is not None:
        with open(LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


# ----------------------------------------------------------------------------- the genomes

def read_genome_rows(table: Path) -> list:
    return [r for r in ga.read_table(table).values()]


def genome_file_name(row: dict) -> str:
    return f"{ga.label_of(row)}_{row['accession'].strip()}.fna.gz"


def genome_folder_url(accession: str, assembly: str) -> str:
    """The NCBI folder of an assembly: GCF_000005845.2 with ASM584v2 is genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/."""
    prefix, rest = accession.split("_")
    digits = rest.split(".")[0]
    return f"{NCBI}{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{accession}_{assembly}/"


def gzip_ok(path: Path) -> bool:
    """True when the whole file reads as gzip (a truncated or half-written download does not)."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 22):
                pass
        return True
    except (OSError, EOFError):
        return False


def expected_md5(folder_url: str, file_name: str, opener=None):
    """NCBI's MD5 sum of one file of an assembly folder, or None when the folder's md5checksums.txt cannot be read or does not list it."""
    opener = opener or urllib.request.urlopen
    try:
        with opener(folder_url + "md5checksums.txt", timeout=60) as response:
            text = response.read().decode("utf-8", "replace")
    except OSError:
        return None
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[1].lstrip("./") == file_name:
            return fields[0]
    return None


def download_genome(row: dict, folder: Path, retries: int = 6, log=say, opener=None) -> str:
    folder.mkdir(parents=True, exist_ok=True)
    acc, asm = row["accession"].strip(), row["assembly_name"].strip()
    remote = f"{acc}_{asm}_genomic.fna.gz"
    base = genome_folder_url(acc, asm)
    target = folder / genome_file_name(row)
    if target.exists() and gzip_ok(target):
        log(f"{target.name}: already here, reads as gzip")
        return "present"
    md5 = expected_md5(base, remote, opener)
    if md5 is None:
        log(f"{target.name}: NCBI's checksum could not be read, so only the gzip check will be made")
    part = folder / (target.name + ".part")
    for attempt in range(1, retries + 1):
        try:
            download_ena.fetch_resumable(base + remote, part, 0, opener)
        except download_ena.ServerAnswer as err:
            part.unlink(missing_ok=True)
            raise SystemExit(f"{target.name}: NCBI answered with a short message instead of the file ({err}); check {base + remote}")
        except OSError as err:
            log(f"{target.name}: interrupted ({err}); resuming (attempt {attempt} of {retries})")
            continue
        if (md5 is None or download_ena.md5_of(part) == md5) and gzip_ok(part):
            part.replace(target)
            log(f"{target.name}: downloaded, " + ("checksum ok" if md5 else "gzip ok"))
            return "downloaded"
        log(f"{target.name}: checksum or gzip check failed, starting again (attempt {attempt} of {retries})")
        part.unlink(missing_ok=True)
    raise SystemExit(f"{target.name}: no verified download after {retries} attempts ({base + remote})")


def step_genomes(c: Ctx) -> str:
    rows = read_genome_rows(c.table)
    if not rows:
        raise SystemExit(f"No rows with an accession in {c.table}")
    results = [download_genome(r, c.genomes, c.args.retries, say) for r in rows]
    return f"{len(rows)} genomes in {c.genomes} ({results.count('downloaded')} downloaded, {results.count('present')} already there)"


# ----------------------------------------------------------------------------- the other steps

def step_amplicons(c: Ctx) -> str:
    if not c.genomes.is_dir() or not ga.find_genomes(c.genomes):
        raise Skip(f"no genome files in {c.genomes}; run the genomes step first")
    result = ga.run(c.genomes, c.table, PRIMERS[0], PRIMERS[1], c.amplicons, log=say)
    found = [k for k, v in result["strains"].items() if v["variants"]]
    return (f"{len(found)} of {len(result['strains'])} genomes give an amplicon, {sum(len(result['strains'][k]['variants']) for k in found)} expected variants; "
            f"{len(result['warnings'])} thing(s) to check in {c.amplicons / 'EXPECTED_AMPLICONS.md'}")


def link_or_copy(src: Path, dst: Path) -> None:
    """A hard link where the drive allows it (no second copy of the data; organise_data.py counts the two names as one file), else a copy."""
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def need_files(folder: Path, what: str) -> list:
    files = lf.find_files(folder) if folder.is_dir() else []
    if not files:
        raise Skip(f"no FASTQ files in {folder} ({what})")
    return files


def verdict_text(out: Path) -> str:
    met, worst = lf.verdict(read_check_rows(out), lf.TOLERANCE)
    return "gate L4 " + ("MET" if met else "NOT MET") + ("" if worst is None else f" (largest difference {100 * abs(worst['difference']):.2f} points)")


def read_check_rows(out: Path) -> list:
    """[(sample, [row, ...])] from a filter_check.tsv written by longread_filter_check.py."""
    table = Path(out) / "filter_check.tsv"
    if not table.is_file():
        return []
    by = {}
    for r in lf.read_table(table):
        def num(x):
            return None if x in ("", "None") else float(x)
        by.setdefault(r["sample"], []).append({"maxEE": float(r["maxEE"]), "estimate": num(r["estimate"]), "dada2": num(r["dada2"]), "difference": num(r["difference"]),
                                               "within": r["within_tolerance"] == "True"})
    return list(by.items())


def step_filter_monomer(c: Ctx) -> str:
    files = need_files(c.monomer, "the Monomer files")
    chosen = [p for p in files if any(pair in p.name for pair in MONOMER_PAIRS)]
    if len(chosen) < len(MONOMER_PAIRS):
        raise Skip(f"only {len(chosen)} of the {len(MONOMER_PAIRS)} Monomer files for the check are in {c.monomer} ({', '.join(MONOMER_PAIRS)})")
    c.monomer_in.mkdir(parents=True, exist_ok=True)
    try:
        for p in chosen:
            link_or_copy(p, c.monomer_in / p.name)
        lf.check_tools()
        lf.check_folder(c.monomer_in, c.monomer_out, PRIMERS[0], PRIMERS[1], GRID, c.args.reads_monomer, 0.15, c.args.tolerance, False, say)
    finally:                                                          # the three names are removed afterwards: the data stay where they were, and no copy is left for a duplicate search to find
        shutil.rmtree(c.monomer_in, ignore_errors=True)
    return verdict_text(c.monomer_out)


def step_filter_kinnex(c: Ctx) -> str:
    need_files(c.kinnex, "the Kinnex files")
    lf.check_tools()
    lf.check_folder(c.kinnex, c.kinnex_out, PRIMERS[0], PRIMERS[1], GRID, c.args.reads, 0.15, c.args.tolerance, False, say)
    return verdict_text(c.kinnex_out)


# ----------------------------------------------------------------------------- the summary

def tool_version(cmd: list, timeout: int = 60):
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired) as err:
        return None, str(err)
    text = (res.stdout or res.stderr).strip().splitlines()
    return (res.returncode == 0), (text[0] if text else "")


def versions() -> dict:
    out = {"Python": sys.version.split()[0]}
    try:
        import numpy
        out["numpy"] = numpy.__version__
    except ImportError:
        out["numpy"] = "not installed"
    ok, text = tool_version([sys.executable, "-m", "cutadapt", "--version"])
    out["cutadapt"] = text if ok else "not available"
    rscript = find_rscript()
    if rscript:
        ok, text = tool_version([rscript, "-e", "cat(R.version.string, '| DADA2', as.character(packageVersion('dada2')))"])
        out["R and DADA2"] = text if ok else "DADA2 is not installed in this R"
    else:
        out["R and DADA2"] = "Rscript not found"
    return out


def md_lines(path: Path, start: str, stop: str = None) -> list:
    """The lines of a markdown file from the line that starts with `start` up to (not including) the next line that starts with `stop`."""
    if not Path(path).is_file():
        return []
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    begin = next((i for i, l in enumerate(lines) if l.startswith(start)), None)
    if begin is None:
        return []
    end = next((i for i in range(begin + 1, len(lines)) if stop and lines[i].startswith(stop)), len(lines))
    return lines[begin:end]


def step_summary(c: Ctx) -> str:
    v = versions()
    lines = ["# Night 7: long-read checks on the PacBio Kinnex-16S data of ATCC MSA-1003", "",
             f"Run on {datetime.date.today()}; folder {c.out}. Gates and predictions: FUTURE.md of the archive (written before the run). "
             "Tools: " + "; ".join(f"{k} {x}" for k, x in v.items()) + ".", "", "## Steps", ""]
    for name in STEPS:
        if name == "summary":
            continue
        s = c.state["steps"].get(name)
        lines.append(f"- {name}: " + (f"{s['status']} in {max(1, round(s['seconds'] / 60))} min" + (f" ({s['detail'].splitlines()[0]})" if s["detail"] else "") if s else "not run"))
    for title, out in (("Monomer files (Sequel II, three barcode pairs)", c.monomer_out), ("Kinnex files (Sequel II, three barcode pairs, the first reads of each)", c.kinnex_out)):
        lines += ["", f"## L4, filter agreement: {title}", ""]
        gate = md_lines(out / "FILTER_CHECK.md", "**Gate")
        if not gate:
            lines.append("No result (the step did not finish; see the steps above).")
            continue
        lines += [gate[0], ""] + md_lines(out / "FILTER_CHECK.md", "| sample | maxEE", "## ")[:40] + [""] + md_lines(out / "FILTER_CHECK.md", "## Reads that survive")
    lines += ["", "## Expected amplicons of the 20 strains", ""]
    exp = c.amplicons / "EXPECTED_AMPLICONS.md"
    head = md_lines(exp, "Primers")
    if not head:
        lines.append("No result (the step did not finish; see the steps above).")
    else:
        lines += md_lines(exp, "Primers", "## ") + md_lines(exp, "## To check", "## What make_reference")
    lines += ["", "## Send back", "", f"- this file ({c.out / 'NIGHT7_SUMMARY.md'})", f"- {c.monomer_out / 'FILTER_CHECK.md'} and {c.kinnex_out / 'FILTER_CHECK.md'}",
              f"- {exp}", f"- {c.log_file} if a step failed"]
    (c.out / "NIGHT7_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return "NIGHT7_SUMMARY.md written"


STEP_FUNCTIONS = {"genomes": step_genomes, "amplicons": step_amplicons, "filter_monomer": step_filter_monomer, "filter_kinnex": step_filter_kinnex, "summary": step_summary}


def save_state(c: Ctx) -> None:
    c.out.mkdir(parents=True, exist_ok=True)
    c.state_file.write_text(json.dumps(c.state, indent=1), encoding="utf-8")


def run_step(c: Ctx, name: str, force: bool = False) -> None:
    done = c.state["steps"].get(name, {})
    if done.get("status") == "done" and not c.args.fresh and not force and name != "summary":
        say(f"{name}: finished earlier, skipped (--redo {name} repeats it)")
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

def check_environment(c: Ctx, opener=None) -> list:
    """[(what, ok, detail)] for everything the run needs."""
    out = [("Python 3.9 or later", sys.version_info >= (3, 9), sys.version.split()[0])]
    try:
        import numpy
        out.append(("numpy", True, numpy.__version__))
    except ImportError:
        out.append(("numpy", False, "pip install numpy"))
    ok, text = tool_version([sys.executable, "-m", "cutadapt", "--version"])
    out.append(("cutadapt", bool(ok), ("version " + text) if ok else "pip install cutadapt"))
    rscript = find_rscript()
    out.append(("Rscript", bool(rscript), rscript or NOT_FOUND.splitlines()[0]))
    if rscript:
        ok, text = tool_version([rscript, "-e", "cat(as.character(packageVersion('dada2')))"])
        out.append(("DADA2 in R", bool(ok), ("version " + text) if ok else "install it in R: BiocManager::install('dada2')"))
    try:
        free = shutil.disk_usage(c.out if c.out.exists() else (c.out.parent if c.out.parent.exists() else Path("."))).free / 1e9
        out.append((f"{NEEDED_GB} GB free", free >= NEEDED_GB, f"{free:.0f} GB free"))
    except OSError as err:
        out.append((f"{NEEDED_GB} GB free", False, str(err)))
    out.append(("strain table", c.table.is_file(), str(c.table)))
    monomer = lf.find_files(c.monomer) if c.monomer.is_dir() else []
    found = [p for p in monomer if any(pair in p.name for pair in MONOMER_PAIRS)]
    out.append(("Monomer files for L4", len(found) == len(MONOMER_PAIRS), f"{len(found)} of {len(MONOMER_PAIRS)} in {c.monomer}"))
    kinnex = lf.find_files(c.kinnex) if c.kinnex.is_dir() else []
    out.append(("Kinnex files for L4", bool(kinnex), f"{len(kinnex)} in {c.kinnex}"))
    have = len(ga.find_genomes(c.genomes)) if c.genomes.is_dir() else 0
    wanted = len(read_genome_rows(c.table)) if c.table.is_file() else 0
    if wanted and have >= wanted:
        out.append(("NCBI reachable", True, "not needed, all genomes are here"))
    else:
        try:
            with (opener or urllib.request.urlopen)(NCBI, timeout=30):
                out.append(("NCBI reachable", True, NCBI))
        except OSError as err:
            out.append(("NCBI reachable", False, f"{err} (needed for the genomes step)"))
    return out


def list_steps(c: Ctx) -> None:
    for name in STEPS:
        s = c.state["steps"].get(name)
        print(f"  {name:15s} " + (f"{s['status']} ({s['finished']})" if s else "not run"))


def main(argv=None) -> int:
    global LOG
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("night7_longread"), help="the folder for everything this run makes")
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA, help=f"the folder of the downloaded data (default {DEFAULT_DATA})")
    ap.add_argument("--table", type=Path, default=GENOME_TABLE, help="the table of the 20 strains and their genomes (default: msa1003_genomes.tsv next to this script)")
    ap.add_argument("--reads", type=int, default=50000, help="reads per Kinnex file for test L4 (default 50000; 0 = all)")
    ap.add_argument("--reads-monomer", type=int, default=0, help="reads per Monomer file for test L4 (default: all)")
    ap.add_argument("--tolerance", type=float, default=lf.TOLERANCE, help="the gate of test L4 as a share (default 0.02 = 2 points)")
    ap.add_argument("--only", help="comma-separated steps to run")
    ap.add_argument("--skip", help="comma-separated steps to leave out")
    ap.add_argument("--redo", help="comma-separated steps to repeat even though they finished; the summary is written again")
    ap.add_argument("--fresh", action="store_true", help="repeat every step")
    ap.add_argument("--check", action="store_true", help="check the machine, the tools and the files, then stop")
    ap.add_argument("--list", action="store_true", help="show the steps and which have finished, then stop")
    ap.add_argument("--retries", type=int, default=6, help="attempts per genome (default 6)")
    a = ap.parse_args(argv)
    c = Ctx(a)
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
    say(f"Night 7, folder {c.out}")
    t0 = time.time()
    if "genomes" in redo and "amplicons" not in redo:
        say("Note: amplicons uses the genomes and is not repeated; add it to --redo if the genomes changed.")
    for name in STEPS:
        if name in wanted and name not in skipped:
            run_step(c, name, force=name in redo)
    say(f"Finished in {(time.time() - t0) / 60:.0f} minutes. Read {c.out / 'NIGHT7_SUMMARY.md'} first.")
    failed = [n for n, s in c.state["steps"].items() if s["status"] == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
