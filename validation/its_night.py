#!/usr/bin/env python
"""ITS overnight run (TEST_REGISTER.md, the section of the ITS overnight run; criteria P83 to P93): one resumable night for every dataset and every arm.

  python validation\\its_night.py --out C:\\data\\runs\\its_night --config C:\\amplitrim_020\\amplitrim_paths.json
  python validation\\its_night.py --out C:\\data\\runs\\its_night_smoke --config C:\\amplitrim_020\\amplitrim_paths.json --smoke      (300 pairs per sample, every dataset)
  python validation\\its_night.py --out ... --dry-run          (lists the datasets and steps; runs nothing)

Run it from the validation archive's folder, or with PYTHONPATH set to its validation folder: it uses the archive's own helpers (real_data.thin_reads, night10's dataset
list and depth rules, its_estimates' AmpliTrim call, rtools). It needs R with DADA2 (found as run 25 found it), vsearch on PATH (AmpliTrim's), and AmpliTrim 0.2.0.

ITS datasets (folders under --base; --dataset LABEL=FOLDER changes one): the ATCC standard (atcc_its1, atcc_its2), the piglet isolates (arfken_its1, arfken_its2) and mixes
(arfken_mix_its1, arfken_mix_its2) of Arfken et al. 2023, the second ITS study (its2) and the first (its). 16S controls: the 16S datasets of amplitrim_paths.json (--config),
as the closing run read them (night10.datasets), at AmpliTrim's cut; --no-16s leaves them out.

Every dataset is thinned with real_data.thin_reads to --pairs read pairs per sample (default 5,000; every pair of mock A, as in the closing run). For an ITS dataset the steps are:
  thin, sw15 (a sliding-window trim of the read ends: window 4, mean quality 15, cut before the first window below it), mfirst (vsearch merges the raw pairs:
  --fastq_minovlen 12 --fastq_allowmergestagger, otherwise its default allowance), arms (its_night.R: std, ee23, ee24, ee25, ee33, ee55, tq10, tq20, sw15, fwd, mo20, mo30,
  mfirst), amplitrim (AmpliTrim 0.2.0 --its on the thinned reads), estimates (vsearch on the reads that passed std's filter, allowing 10, 5, 2 and 0 differences).
For a 16S control: thin, amplitrim (its recommended cut), arms (std only, at that cut, maxEE 2,2), estimates. When AmpliTrim cannot analyse a control's thinned reads,
the cut that the final run recommended for the whole dataset is used (--final-run), and amplitrim\\CUT_FROM_FINAL_RUN.txt says so.
A step that finished is not run again (night_status.json); a step that fails is recorded and the rest of that dataset is skipped, while the other datasets go on.
Next, after the night: python validation\\its_night_report.py --out <the same folder>.
"""
import argparse
import gzip
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)]

ITS_FOLDERS = {
    "atcc_its1": "20_mock_ITS1_arfken_atcc_PRJNA693350/raw", "atcc_its2": "21_mock_ITS2_arfken_atcc_PRJNA693350/raw",
    "arfken_its1": "16_mock_ITS1_arfken_isolates_PRJNA693350/raw", "arfken_its2": "17_mock_ITS2_arfken_isolates_PRJNA693350/raw",
    "arfken_mix_its1": "18_mock_ITS1_arfken_mixes_PRJNA693350/raw", "arfken_mix_its2": "19_mock_ITS2_arfken_mixes_PRJNA693350/raw",
    "its2": "15_real_ITS_fullITS_PRJEB108994/prepared/its_soil2/ready", "its": "09_real_ITS_cryobank_PRJEB58189/prepared/all/ready",
}
COLS = ["arm", "mode", "reads", "truncF", "truncR", "maxEEF", "maxEER", "truncQ", "minOverlap"]
ITS_ARMS = [  # arm, mode, reads (a folder name under the dataset, or the arm reused), truncF, truncR, maxEEF, maxEER, truncQ, minOverlap
    ("std", "paired", "reads", 0, 0, 2, 2, 2, 12),
    ("ee23", "paired", "reads", 0, 0, 2, 3, 2, 12), ("ee24", "paired", "reads", 0, 0, 2, 4, 2, 12), ("ee25", "paired", "reads", 0, 0, 2, 5, 2, 12),
    ("ee33", "paired", "reads", 0, 0, 3, 3, 2, 12), ("ee55", "paired", "reads", 0, 0, 5, 5, 2, 12),
    ("tq10", "paired", "reads", 0, 0, 2, 2, 10, 12), ("tq20", "paired", "reads", 0, 0, 2, 2, 20, 12),
    ("sw15", "paired", "reads_sw15", 0, 0, 2, 2, 2, 12),
    ("fwd", "fwd", "std", 0, 0, 2, 2, 2, 12), ("mo20", "overlap", "std", 0, 0, 2, 2, 2, 20), ("mo30", "overlap", "std", 0, 0, 2, 2, 2, 30),
    ("mfirst", "single", "reads_mfirst", 0, 0, 2, 2, 2, 12),
]
DIFFS = (10, 5, 2, 0)
WINDOW, WINDOW_Q = 4, 15
# --plan rules (the ITS rules night, register section 35): truncQ 3 to 10, the trim at four strengths, and maxEE computed from the reads (the arm eer, added at run time)
RULES_ARMS = [("std", "paired", "reads", 0, 0, 2, 2, 2, 12)] + [(f"tq{q}", "paired", "reads", 0, 0, 2, 2, q, 12) for q in (3, 4, 5, 6, 8, 10)] \
    + [(f"sw{q}", "paired", f"reads_sw{q}", 0, 0, 2, 2, 2, 12) for q in (10, 15, 20, 25)]
SW_STRENGTHS = (10, 15, 20, 25)
# --plan bakker (the third mock community, register section 37): the 8 ITS datasets and Bakker's 3, the advised trim, its weaker form, the computed maxEE,
# the trim with maxEE computed on the trimmed reads (sweer), and truncQ = 10 as the known-unsafe reference
BAKKER_FOLDERS = {"bakker_pcr_its1": "22_mock_ITS1_bakker_pcr_PRJNA377530/prepared/its1/ready",
                  "bakker_pol_its1": "23_mock_ITSmix_bakker_polymerase_PRJNA377530/prepared/its1/ready",
                  "bakker_pol_its2": "23_mock_ITSmix_bakker_polymerase_PRJNA377530/prepared/its2/ready"}
BAKKER_ARMS = [("std", "paired", "reads", 0, 0, 2, 2, 2, 12), ("sw10", "paired", "reads_sw10", 0, 0, 2, 2, 2, 12),
               ("sw15", "paired", "reads_sw15", 0, 0, 2, 2, 2, 12), ("tq10", "paired", "reads", 0, 0, 2, 2, 10, 12)]
BAKKER_STRENGTHS = (10, 15)
RULE_Q, TRUNC_Q, MIN_LEN = 0.75, 2, 50


# ----------------------------------------------------------------------------------------------------------- pure helpers (unit-tested)
def sliding_cut(qual: str, window: int = WINDOW, q: float = WINDOW_Q, offset: int = 33) -> int:
    """The length to keep: everything before the first window of `window` bases whose mean quality is below q (the whole read when there is none)."""
    s = [ord(c) - offset for c in qual]
    for i in range(0, max(0, len(s) - window + 1)):
        if sum(s[i:i + window]) / window < q:
            return i
    return len(s)


def fastq_records(path: Path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, p, qq = fh.readline().rstrip("\n"), fh.readline().rstrip("\n"), fh.readline().rstrip("\n")
            yield h.rstrip("\n"), s, p, qq


def trim_folder(src: Path, dst: Path, q: float = WINDOW_Q) -> int:
    """The sliding-window trim (sw15 by default; sw10, sw20, sw25 with q): every read of every FASTQ in src cut by sliding_cut, written to dst under the same name
    (both mates, so the pairs stay in order). Returns the files written."""
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in sorted(src.glob("*.fastq*")):
        out = dst / f.name
        part = out.with_name(out.name + ".part")
        with gzip.open(part, "wt", newline="\n") if f.name.endswith(".gz") else open(part, "w", newline="\n") as o:
            for h, s, p, qq in fastq_records(f):
                k = sliding_cut(qq, q=q)
                o.write(f"{h}\n{s[:k]}\n{p}\n{qq[:k]}\n")
        part.replace(out)
        n += 1
    return n


def read_ee(qual: str, trunc_q: int = TRUNC_Q, offset: int = 33) -> tuple:
    """(length kept, expected errors) of a read after DADA2's truncQ cut: the read stops before its first base of quality truncQ or below; EE = sum of 10^(-Q/10)."""
    s = [ord(c) - offset for c in qual]
    k = next((i for i, v in enumerate(s) if v <= trunc_q), len(s))
    return k, sum(10 ** (-v / 10) for v in s[:k])


def ee_rule(folder: Path, min_len: int = MIN_LEN) -> dict:
    """AmpliTrim's maxEE rule (pipeline._rule_max_ee: numpy's 75th percentile of the expected errors per direction, int + 1, at least 2), applied to ITS reads
    without truncation: the expected errors of each read after the truncQ = 2 cut, over the reads that keep at least min_len bases."""
    import numpy as np
    ef, er = [], []
    for _, r1, r2 in pairs_in(folder):
        for f, acc in ((r1, ef), (r2, er)):
            for _, _, _, q in fastq_records(f):
                k, e = read_ee(q)
                if k >= min_len:
                    acc.append(e)
    if not ef or not er:
        raise RuntimeError(f"no read keeps {min_len} bases in {folder}")
    pct = (float(np.percentile(ef, 100 * RULE_Q)), float(np.percentile(er, 100 * RULE_Q)))
    raw = (int(pct[0]) + 1, int(pct[1]) + 1)
    return {"percentile": list(pct), "raw": list(raw), "max_ee": [max(2, v) for v in raw], "reads": [len(ef), len(er)]}


def pairs_in(folder: Path) -> list:
    """(sample, R1, R2) of a folder of <sample>_R1 / <sample>_R2 FASTQs."""
    out = []
    for r1 in sorted(folder.glob("*_R1*.f*q*")):
        r2 = r1.with_name(r1.name.replace("_R1", "_R2", 1))
        if r2.is_file():
            out.append((r1.name.split("_R1")[0], r1, r2))
    return out


def plain(path: Path, folder: Path) -> Path:
    """An uncompressed copy for vsearch (whose Windows builds may not read gzip)."""
    if not str(path).endswith(".gz"):
        return path
    folder.mkdir(parents=True, exist_ok=True)
    out = folder / path.name[:-3]
    with gzip.open(path, "rb") as i, open(out, "wb") as o:
        shutil.copyfileobj(i, o)
    return out


def vsearch_merge_cmd(r1: Path, r2: Path, out: Path, max_diffs=None) -> list:
    cmd = ["vsearch", "--fastq_mergepairs", str(r1), "--reverse", str(r2), "--fastq_minovlen", "12", "--fastq_allowmergestagger", "--fastqout", str(out), "--quiet"]
    if max_diffs is not None:
        cmd += ["--fastq_maxdiffs", str(max_diffs)]
    return cmd


def count_fastq(path: Path) -> int:
    if not Path(path).is_file():
        return 0
    with open(path, "rb") as fh:
        return sum(1 for _ in fh) // 4


def arms_rows(dataset_dir: Path, its: bool, cut=(0, 0)) -> list:
    """The rows of its_night.R's arms table: every ITS arm, or for a 16S control the default arm at AmpliTrim's cut."""
    def row(arm, mode, reads, tf, tr, ef, er, tq, mo):
        where = str(dataset_dir / reads) if mode in ("paired", "single") else reads
        return [arm, mode, where, tf, tr, ef, er, tq, mo]
    if its:
        return [row(*a) for a in ITS_ARMS]
    return [row("std", "paired", "reads", cut[0], cut[1], 2, 2, 2, 12)]


# ----------------------------------------------------------------------------------------------------------- the night
class Night:
    def __init__(self, out: Path, run=subprocess.run, rscript=None, thin=None, log=print):
        self.out, self.run, self.rscript, self.thin, self.log = Path(out), run, rscript, thin, log
        self.status_path = self.out / "night_status.json"
        self.status = json.loads(self.status_path.read_text(encoding="utf-8")) if self.status_path.is_file() else {}

    def save(self):
        self.out.mkdir(parents=True, exist_ok=True)
        self.status_path.write_text(json.dumps(self.status, indent=1), encoding="utf-8")

    def call(self, cmd, cwd=None):
        p = self.run([str(c) for c in cmd], capture_output=True, text=True, cwd=str(cwd) if cwd else None)
        if p.returncode != 0:
            raise RuntimeError(f"{Path(str(cmd[0])).name} stopped: {((p.stderr or '') + (p.stdout or '')).strip()[-600:]}")
        return p

    def dataset(self, label: str, steps: list):
        """Run the steps of one dataset in order; a step done before is skipped; after a failure the rest of the dataset is skipped."""
        for name, fn in steps:
            key = f"{label}/{name}"
            if self.status.get(key, {}).get("state") == "done":
                self.log(f"  {key}: done before")
                continue
            t = time.time()
            try:
                fn()
            except Exception as e:                                     # recorded; the other datasets go on
                self.status[key] = {"state": "failed", "error": str(e)[-800:], "seconds": round(time.time() - t, 1)}
                self.save()
                self.log(f"  {key}: FAILED: {str(e)[-300:]}")
                for later, _ in steps[[n for n, _ in steps].index(name) + 1:]:
                    self.status[f"{label}/{later}"] = {"state": "skipped", "error": f"{name} failed"}
                self.save()
                return False
            self.status[key] = {"state": "done", "seconds": round(time.time() - t, 1)}
            self.save()
            self.log(f"  {key}: done in {time.time() - t:.0f} s")
        return True

    # --- the steps
    def step_thin(self, src: Path, d: Path, pairs: int):
        reads = d / "reads"
        if pairs == 0:                                                 # every pair (mock A): a copy of the folder's pairs under the usual names
            reads.mkdir(parents=True, exist_ok=True)
            for s, r1, r2 in pairs_in(src):
                for r, m in ((r1, 1), (r2, 2)):
                    shutil.copyfile(r, reads / f"{s}_R{m}.fastq{'.gz' if r.name.endswith('.gz') else ''}")
        else:
            self.thin(src, reads, pairs)
        if not pairs_in(reads):
            raise RuntimeError(f"no read pairs in {src}")

    def step_sw15(self, d: Path):
        trim_folder(d / "reads", d / "reads_sw15")

    def step_trims(self, d: Path, strengths=SW_STRENGTHS):
        for q in strengths:
            trim_folder(d / "reads", d / f"reads_sw{q}", q)

    def step_eer(self, d: Path):
        (d / "eer.json").write_text(json.dumps(ee_rule(d / "reads"), indent=1), encoding="utf-8")

    def step_eer_sw15(self, d: Path):
        (d / "eer_sw15.json").write_text(json.dumps(ee_rule(d / "reads_sw15"), indent=1), encoding="utf-8")

    def step_mfirst(self, d: Path):
        out = d / "reads_mfirst"
        out.mkdir(parents=True, exist_ok=True)
        tmp = d / "tmp_plain"
        for s, r1, r2 in pairs_in(d / "reads"):
            self.call(vsearch_merge_cmd(plain(r1, tmp), plain(r2, tmp), out / f"{s}_R1.fastq"))
        shutil.rmtree(tmp, ignore_errors=True)

    def step_arms(self, d: Path, rows: list):
        arms = d / "arms.tsv"
        arms.write_text("\t".join(COLS) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows), encoding="utf-8")
        if not self.rscript:
            raise RuntimeError("Rscript was not found (R with DADA2 is needed)")
        self.call([self.rscript, HERE / "its_night.R", arms, d / "dada2", "TRUE"])

    def step_amplitrim(self, d: Path, its: bool, call_text: str, fallback: Path = None):
        """AmpliTrim on the thinned reads. For a 16S control that AmpliTrim cannot analyse (too few pairs overlap at full read length in the thinned reads), the cut
        that the final run recommended for the whole dataset is used instead, when `fallback` (its summary.json) is there; the file CUT_FROM_FINAL_RUN.txt says so."""
        rep = d / "amplitrim"
        try:
            self.call([sys.executable, "-c", call_text, "--input-dir", d / "reads", "--work-dir", rep / "work", "--out-dir", rep] + (["--its"] if its else []), cwd=self.out)
        except RuntimeError as e:
            if its or fallback is None or not Path(fallback).is_file():
                raise
            rep.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(fallback, rep / "summary.json")
            (rep / "CUT_FROM_FINAL_RUN.txt").write_text(f"AmpliTrim could not analyse the thinned reads ({str(e)[-300:]}); the cut is the final run's for the whole dataset: {fallback}\n",
                                                        encoding="utf-8")
            self.log(f"    AmpliTrim could not analyse the thinned reads; the final run's cut is used ({fallback})")
        if not (rep / "summary.json").is_file():
            raise RuntimeError("AmpliTrim wrote no summary.json")

    def step_estimates(self, d: Path):
        filt = d / "dada2" / "std" / "filtered"
        tmp = d / "tmp_plain"
        lines = ["sample\tfiltered_pairs\t" + "\t".join(f"merged_d{k}" for k in DIFFS)]
        for s, r1, r2 in pairs_in(filt):
            p1, p2 = plain(r1, tmp), plain(r2, tmp)
            n = count_fastq(p1)
            got = []
            for k in DIFFS:
                merged = tmp / f"{s}_d{k}.fastq"
                self.call(vsearch_merge_cmd(p1, p2, merged, max_diffs=k))
                got.append(count_fastq(merged))
            lines.append(f"{s}\t{n}\t" + "\t".join(str(g) for g in got))
        shutil.rmtree(tmp, ignore_errors=True)
        (d / "estimates.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")


def cut_of(d: Path) -> tuple:
    r = json.loads((d / "amplitrim" / "summary.json").read_text(encoding="utf-8"))["recommended"]
    return int(r["truncLenF"]), int(r["truncLenR"])


def rules_rows(d: Path) -> list:
    """The arms of the rules night for one dataset: RULES_ARMS and the arm eer at the maxEE that ee_rule computed (eer.json)."""
    ef, er = json.loads((d / "eer.json").read_text(encoding="utf-8"))["max_ee"]
    rows = [[arm, mode, str(d / reads), tf, tr, a, b, tq, mo] for arm, mode, reads, tf, tr, a, b, tq, mo in RULES_ARMS]
    return rows + [["eer", "paired", str(d / "reads"), 0, 0, ef, er, 2, 12]]


def bakker_rows(d: Path) -> list:
    """The arms of the Bakker night for one dataset: BAKKER_ARMS, eer (maxEE computed on the reads) and sweer (the trim, then maxEE computed on the trimmed reads)."""
    rows = [[arm, mode, str(d / reads), tf, tr, a, b, tq, mo] for arm, mode, reads, tf, tr, a, b, tq, mo in BAKKER_ARMS]
    e = json.loads((d / "eer.json").read_text(encoding="utf-8"))["max_ee"]
    w = json.loads((d / "eer_sw15.json").read_text(encoding="utf-8"))["max_ee"]
    return rows + [["eer", "paired", str(d / "reads"), 0, 0, e[0], e[1], 2, 12], ["sweer", "paired", str(d / "reads_sw15"), 0, 0, w[0], w[1], 2, 12]]


def plan(a, sixteen) -> list:
    """(label, source folder, ITS or not, pairs per sample)."""
    its = {k: Path(a.base) / v for k, v in ITS_FOLDERS.items()}
    if getattr(a, "plan", "night") == "bakker":
        its.update({k: Path(a.base) / v for k, v in BAKKER_FOLDERS.items()})
    for item in a.dataset:
        k, sep, v = item.partition("=")
        if not sep or k not in its:
            raise SystemExit(f"--dataset wants LABEL=FOLDER with LABEL one of {sorted(its)}")
        its[k] = Path(v)
    for item in getattr(a, "extra", []):
        k, sep, v = item.partition("=")
        if not sep or not k or k in its:
            raise SystemExit(f"--extra wants a new LABEL=FOLDER, not {item!r}")
        its[k] = Path(v)
    out = [(k, its[k], True, a.pairs) for k in its]
    if not a.no_16s:
        for ds in sixteen:
            out.append((f"16s_{ds['name']}", Path(ds["reads"]), False, 0 if ds["pairs"] == 0 else a.pairs))
    only = {x for x in a.only.split(",") if x}
    return [p for p in out if not only or p[0] in only]


def main(argv=None, run=subprocess.run, thin=None, rscript=None, sixteen=None, amplitrim_call=None, installed=None, log=print) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--base", type=Path, default=Path(r"C:\data"))
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json, for the 16S controls")
    ap.add_argument("--dataset", action="append", default=[], metavar="LABEL=FOLDER")
    ap.add_argument("--extra", action="append", default=[], metavar="LABEL=FOLDER", help="one more ITS dataset (for example a third mock community)")
    ap.add_argument("--plan", choices=["night", "rules", "bakker"], default="night",
                    help="night: the ITS overnight run (section 33); rules: the ITS rules night (section 35); bakker: the third mock community (section 37)")
    ap.add_argument("--pairs", type=int, default=5000)
    ap.add_argument("--only", default="", help="labels to run, e.g. atcc_its1,its2 (default: all)")
    ap.add_argument("--no-16s", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--smoke", action="store_true", help="the smoke run: every dataset, mock A included, thinned to 300 read pairs per sample")
    ap.add_argument("--expect-version", default="0.2.0")
    ap.add_argument("--final-run", type=Path, help="the final run's command-line reports (results/24_night10_closing_2026-11-06/cli of the archive; found by itself when the archive's "
                                                     "helpers are used): the fallback cut of a 16S control that AmpliTrim cannot analyse")
    a = ap.parse_args(argv)
    if a.plan in ("rules", "bakker"):
        a.no_16s = True                                                 # no question of the rules night concerns 16S
    if sixteen is None:
        sixteen = []
        if not a.no_16s:
            import night10                                             # the archive's: the 16S datasets as the closing run read them, and its depth rules
            for ds in night10.datasets(night10.load_config(a.config)):
                if ds["marker"] == "16s":
                    sixteen.append({"name": ds["name"], "reads": ds["reads"], "pairs": night10.PAIRS.get(ds["name"], night10.DEFAULT_PAIRS)})
            if not sixteen:
                print("No 16S dataset was found in --config; give it, or --no-16s.")
                return 2
    if a.smoke:
        a.pairs = 300
        sixteen = [dict(ds, pairs=300) for ds in sixteen]
    todo = plan(a, sixteen)
    for label, src, its, pairs in todo:
        log(f"{label:22s} {'ITS' if its else '16S'}  {'every pair' if pairs == 0 else f'{pairs:,} pairs per sample'}  {src}" + ("" if Path(src).is_dir() else "   <- NOT FOUND"))
    if a.dry_run:
        log(f"\n{len(todo)} datasets. Nothing was run.")
        return 0
    missing = [str(src) for _, src, _, _ in todo if not Path(src).is_dir()]
    if missing:
        log("Not found: " + ", ".join(missing))
        return 2
    a.out.mkdir(parents=True, exist_ok=True)
    if a.plan == "night" and (installed is None or amplitrim_call is None):
        import its_estimates                                            # the archive's: AmpliTrim started from --out, so the archive's early copy cannot shadow it
        installed = installed or its_estimates.installed
        amplitrim_call = amplitrim_call or its_estimates.CALL
    version = installed(run, a.out) if a.plan == "night" else ("not used", "this night does not run AmpliTrim")
    if a.plan == "night" and a.expect_version != "any" and version[0] != a.expect_version:
        log(f"The AmpliTrim in this Python is version {version[0]} ({version[1]}), not {a.expect_version}. Nothing was run.")
        return 2
    if thin is None:
        from real_data import thin_reads as thin                        # the archive's: the same thinning as runs 24 and 25
    if rscript is None:
        import rtools
        rscript = rtools.find_rscript()
    final = a.final_run
    if final is None:
        try:
            import night10
            final = Path(night10.__file__).resolve().parents[1] / "results" / "24_night10_closing_2026-11-06" / "cli"
        except ImportError:
            final = None
    night = Night(a.out, run=run, rscript=rscript, thin=thin, log=log)
    log(f"AmpliTrim {version[0]} ({version[1]}); Rscript {rscript}\n")
    for label, src, its, pairs in todo:
        d = a.out / label
        log(label)
        if its and a.plan == "bakker":
            steps = [("thin", lambda d=d, s=src, p=pairs: night.step_thin(s, d, p)), ("trims", lambda d=d: night.step_trims(d, BAKKER_STRENGTHS)),
                     ("eer", lambda d=d: night.step_eer(d)), ("eer_sw15", lambda d=d: night.step_eer_sw15(d)), ("arms", lambda d=d: night.step_arms(d, bakker_rows(d)))]
        elif its and a.plan == "rules":
            steps = [("thin", lambda d=d, s=src, p=pairs: night.step_thin(s, d, p)), ("trims", lambda d=d: night.step_trims(d)), ("eer", lambda d=d: night.step_eer(d)),
                     ("arms", lambda d=d: night.step_arms(d, rules_rows(d)))]
        elif its:
            steps = [("thin", lambda d=d, s=src, p=pairs: night.step_thin(s, d, p)), ("sw15", lambda d=d: night.step_sw15(d)),
                     ("mfirst", lambda d=d: night.step_mfirst(d)), ("arms", lambda d=d: night.step_arms(d, arms_rows(d, True))),
                     ("amplitrim", lambda d=d: night.step_amplitrim(d, True, amplitrim_call)), ("estimates", lambda d=d: night.step_estimates(d))]
        else:
            fb = Path(final) / label[len("16s_"):] / "report" / "summary.json" if final else None
            steps = [("thin", lambda d=d, s=src, p=pairs: night.step_thin(s, d, p)), ("amplitrim", lambda d=d, fb=fb: night.step_amplitrim(d, False, amplitrim_call, fb)),
                     ("arms", lambda d=d: night.step_arms(d, arms_rows(d, False, cut_of(d)))), ("estimates", lambda d=d: night.step_estimates(d))]
        night.dataset(label, steps)
    bad = {k: v for k, v in night.status.items() if v.get("state") != "done"}
    lines = ["# ITS overnight run: status", "", "| Step | State | Seconds | Note |", "|---|---|---|---|"]
    lines += [f"| {k} | {v.get('state')} | {v.get('seconds', '')} | {v.get('error', '')[:200]} |" for k, v in sorted(night.status.items())]
    (a.out / "NIGHT_STATUS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"\n{len(night.status) - len(bad)} steps done, {len(bad)} not. written: {a.out / 'NIGHT_STATUS.md'}")
    log((f"Next: python validation\\its_bakker_report.py --out {a.out}" if a.plan == "bakker" else f"Next: python validation\\its_rules_report.py --out {a.out}" if a.plan == "rules" else f"Next: python validation\\its_night_report.py --out {a.out}") if not bad else "Start the same command again to retry what failed; finished steps are kept.")
    return 0 if not bad else 2


if __name__ == "__main__":
    sys.exit(main())
