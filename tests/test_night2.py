"""validation/night2.py with stand-ins for compare_methods, AmpliTrim's own run and the quality table. The user's own config is never read: a temporary one is given."""
import contextlib, gzip, io, json, shutil, sys
from pathlib import Path
R = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(R / "validation"), str(R)]
import night2 as n2
from amplitrim.trimmer import ReadPairs

def write_fq(path, reads):
    with gzip.open(path, "wt") as fh:
        for i, (seq, q) in enumerate(reads): fh.write(f"@p{i}\n{seq}\n+\n{q}\n")

def world(tmp):
    shutil.rmtree(tmp, ignore_errors=True); tmp.mkdir(parents=True)
    data = tmp / "data"
    q37, q10, q15 = chr(33 + 37), chr(33 + 10), chr(33 + 15)
    write_fq(tmp / "a_R1.fq.gz", [("A" * 150, q37 * 100 + q10 * 50)] * 200); write_fq(tmp / "a_R2.fq.gz", [("C" * 150, q37 * 60 + q15 * 90)] * 200)
    pairs = ReadPairs.load(tmp / "a_R1.fq.gz", tmp / "a_R2.fq.gz")
    cfg = {"data": str(data), "figaro": "C:/figaro/figaro.py", "mockA_high": str(data / "mockA_high/ready"), "reference": str(data / "mockA_ref.fasta"), "v1v2": str(data / "v1v2/ready"), "its": "",
           "extras": {"mockB": {"ready": str(data / "mockB/ready"), "reference": str(data / "mockB/ref.fasta")}, "d6300_V4-V5": {"ready": str(data / "d6300_V4-V5/ready"), "reference": str(data / "ref.fasta")},
                      "swabs": {"ready": str(data / "swabs/ready"), "reference": ""}, "broken": {"ready": str(data / "broken/ready"), "reference": ""}}}
    for d in ("mockA_high", "v1v2", "mockB", "d6300_V4-V5", "swabs", "broken"): (data / d / "ready").mkdir(parents=True)
    (tmp / "cfg.json").write_text(json.dumps(cfg))
    return data, pairs

INFO = {   # what AmpliTrim's own run would say for each stand-in dataset
    "mockB": {"picks": {0.5: {"ee": (2, 4), "reached": True, "yield": 0.5}, 0.8: {"ee": (3, 5), "reached": True, "yield": 0.8}}, "trim": (100, 150), "longest": 515, "typical": 512, "table": []},
    "d6300_V4-V5": {"picks": {0.5: {"ee": (2, 2), "reached": True, "yield": 0.64}, 0.8: {"ee": (2, 5), "reached": True, "yield": 0.8}}, "trim": (100, 150), "longest": 550, "typical": 374, "table": []},
    "swabs": {"picks": {0.5: {"ee": (2, 2), "reached": True, "yield": 0.85}, 0.8: {"ee": (2, 2), "reached": True, "yield": 0.85}}, "trim": (100, 150), "longest": 430, "typical": 428, "table": []},
    "v1v2": {"picks": {0.5: {"ee": (2, 2), "reached": True, "yield": 0.9}, 0.8: {"ee": (2, 2), "reached": True, "yield": 0.9}}, "trim": (100, 150), "longest": 346, "typical": 345, "table": []},
    "mockA_high": {"picks": {0.5: {"ee": (2, 2), "reached": True, "yield": 0.8}, 0.8: {"ee": (2, 2), "reached": True, "yield": 0.8}}, "trim": (100, 150), "longest": 508, "typical": 507, "table": []},
}
COLS = ["name", "truncF", "truncR", "maxEE", "own", "own_what", "est_yield", "est_pass", "est_strict", "strains", "pass", "merged", "nonchim", "asvs", "top15", "n01", "tail", "found", "identical", "extra", "error"]

def tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows: fh.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")

def row(name, tf, tr, ee, merged, found="", strains="", identical="", err=None):
    d = {"name": name, "truncF": tf, "truncR": tr, "maxEE": ee, "pass": "" if merged == "" else merged + 0.02, "merged": merged, "nonchim": "" if merged == "" else merged - 0.01, "found": found, "strains": strains, "identical": identical}
    if err: d["error"] = err
    return d

TABLES = {
    "mockB": lambda: [row("AmpliTrim (maxEE 2)", 279, 248, "2", 0.20, 6, "3/3", 0.687), row("AmpliTrim (maxEE 2,4)", 278, 249, "2,4", 0.45, 6, "3/3", 0.70), row("AmpliTrim (maxEE 3,5)", 277, 250, "3,5", 0.60, 6, "3/3", 0.72),
                      row("AmpliTrim (maxEE 1,4)", 279, 248, "1,4", 0.30, 6, "3/3", 0.69), row("FIGARO (amplicon 515)", 279, 248, "6,9", 0.743, 7, "3/3", 0.7828),
                      row("AmpliTrim (FIGARO's maxEE 6,9)", 274, 253, "6,9", 0.767, 6, "3/3", 0.7526), row("full", 0, 0, "2", 0.02, 4, "3/3", 0.54), row("q30", 120, 90, "2", 0.01, 2, "2/3", 0.4)],
    "d6300_V4-V5": lambda: [row("AmpliTrim (maxEE 2)", 282, 280, "2", 0.572), row("AmpliTrim (maxEE 2,5)", 282, 280, "2,5", 0.70), row("AmpliTrim (maxEE 1,4)", 282, 280, "1,4", 0.55),
                            row("AmpliTrim (maxEE 2, leave out 1.0%)", 200, 200, "2", 0.69), row("AmpliTrim (maxEE 2, leave out 5.0%)", 190, 190, "2", 0.70), row("full", 0, 0, "2", 0.639), row("q30", 100, 60, "2", 0.0)],
    "swabs": lambda: [row("AmpliTrim (maxEE 2)", 224, 218, "2", 0.77), row("AmpliTrim (maxEE 1,4)", 224, 218, "1,4", 0.71), row("FIGARO (amplicon 430)", 224, 218, "1", 0.711), row("full", 0, 0, "2", "", err="DADA2 could not run: kmer"), row("q30", 120, 90, "2", 0.5)],
    "v1v2": lambda: [row("AmpliTrim (maxEE 2)", 201, 157, "2", 0.842), row("AmpliTrim (maxEE 1,4)", 201, 157, "1,4", 0.80), row("FIGARO (amplicon 346)", 225, 133, "2,1", 0.82), row("AmpliTrim (FIGARO's maxEE 2,1)", 211, 147, "2,1", 0.823), row("full", 0, 0, "2", 0.35), row("q30", 120, 90, "2", 0.2)],
    "mockA_high": lambda: [row("AmpliTrim (maxEE 2)", 258, 262, "2", 0.742, 30, "11/13", 0.9), row("AmpliTrim (maxEE 1,4)", 258, 262, "1,4", 0.70, 28, "10/13", 0.9), row("FIGARO (amplicon 508)", 273, 247, "3,4", 0.758, 32, "11/13", 0.9),
                           row("AmpliTrim (FIGARO's maxEE 3,4)", 272, 248, "3,4", 0.758, 31, "11/13", 0.9), row("full", 0, 0, "2", 0.70, 29, "10/13", 0.9), row("q30", 120, 90, "2", 0.3, 5, "3/13", 0.5)],
}

def make(tmp, calls, dry=False, setup=True, **kw):
    data, pairs = world(tmp) if setup else (tmp / "data", None)
    if pairs is None:
        pairs = ReadPairs.load(tmp / "a_R1.fq.gz", tmp / "a_R2.fq.gz")
    def explore(reads, work, marker):
        calls["explore"].append(Path(reads).parent.name)
        if Path(reads).parent.name == "broken": raise RuntimeError("AmpliTrim could not merge these reads")
        return INFO[Path(reads).parent.name]
    def run(cmd):
        cmd = [str(c) for c in cmd]; calls["cmd"].append(cmd)
        name = Path(cmd[2]).parent.name; out = Path(cmd[cmd.index("--out") + 1])
        tsv(out / "compare_methods.tsv", TABLES[name]())
        return 0, "ok"
    # shares: merged fraction per sample, 12 samples around each arm's pooled value (a fixed small spread), read from the 'merged' column of the table
    def shares(track):
        key = Path(track).name; comp = Path(track).parent
        rows = list(__import__("csv").DictReader(open(comp / "compare_methods.tsv"), delimiter="\t"))
        from compare_methods import method_key
        r = next(x for x in rows if "track_" + method_key(x["name"]) + ".tsv" == key)
        return {f"s{i}": float(r["merged"]) + 0.001 * (i - 6) for i in range(12)}
    class P:    # the quality table: looser maxEE keeps noisier reads
        def __init__(s, base): s.b = base
        def error_load(s, tf, tr, ee): return (0.1 * ee[0], 0.1 * ee[1])
    real = pairs
    class Both(P):
        _phred1, _phred2, len1, len2 = real._phred1, real._phred2, real.len1, real.len2
    fake_stable = lambda ds, info, out: {"picks": {"AmpliTrim": [{"trim": [100, 150], "ee": "2"}] * 3, "FIGARO": [{"trim": [100, 150], "ee": "1"}, {"trim": [100, 150], "ee": "1"}, {"trim": [120, 150], "ee": "1"}]},
                                            "stable": {"AmpliTrim": True, "FIGARO": False}}
    n = n2.Night2(R, dry, run=run, explore=explore, pairs=lambda reads: Both(None), shares=shares, stable=fake_stable, disk_free=lambda p: kw.get("disk", 50e9), now=lambda f: "2026-10-10", fresh=kw.get("fresh", False), cfg_path=tmp / "cfg.json")
    return n, data

def go(n):
    with contextlib.redirect_stdout(io.StringIO()): return n.go()

def test_night2_end_to_end(tmp_path):
    tmp = tmp_path
    calls = {"explore": [], "cmd": []}
    n, data = make(tmp, calls); assert go(n) == 0
    assert n.logf is None, "the log file must be closed when the run ends: Windows cannot delete a file that is still open"
    out = data / "runs/night2_2026-10-10"; rep = (out / "NIGHT2_REPORT.md").read_text(); prot = (out / "NIGHT2_PROTOCOL.md").read_text()
    # the protocol was written first and states what was fixed
    for need in ("before any analysis", "Predictions, written before the run", "P2.", "(5,5)", "-m 12", "FIGARO failure is reported as 'no result'", "stability 10 nt"): assert need in prot, need
    # a dataset that fails is noted and the others finish
    assert "broken: failed: RuntimeError" in rep and "Problems: 1" in rep and not (out / "broken/night2_result.json").exists()
    res = {d: json.loads((out / d / "night2_result.json").read_text()) for d in ("mockB", "d6300_V4-V5", "swabs", "v1v2", "mockA_high")}
    # the commands
    cmd = {Path(c[2]).parent.name: c for c in calls["cmd"]}
    mb = cmd["mockB"]; assert mb[mb.index("--amplitrim-max-ee") + 1: mb.index("--amplitrim-max-ee") + 5] == ["2", "2,4", "3,5", "1,4"]           # A, B, B2, B4: four runs
    i = mb.index("--figaro-amplicon"); assert mb[i + 1] == "515" and mb[i + 2].startswith("--") and "--allow-unmergeable" not in mb and "--reference" in mb      # one amplicon length, no tail, a mock with truth
    d6 = cmd["d6300_V4-V5"]; i = d6.index("--figaro-amplicon"); assert d6[i + 1: i + 3] == ["550", "374"] and d6[d6.index("--allow-unmergeable") + 1: d6.index("--allow-unmergeable") + 3] == ["0.01", "0.05"]   # a tail: FIGARO at both lengths, coverage rows
    assert "q30=100,60,2" in mb and "full=0,0,2" in mb and cmd["mockA_high"][cmd["mockA_high"].index("--pairs") + 1] == "0" and cmd["swabs"][cmd["swabs"].index("--pairs") + 1] == "5000"
    # mock B: AmpliTrim's pick merges far fewer pairs than FIGARO's, with cleaner reads, and is less accurate
    b = res["mockB"]["verdicts"]; assert res["mockB"]["against"] == "C"
    assert b["B"]["retention"] == "fewer" and b["B"]["noise"] == "cleaner" and b["B"]["verdict"] == "trade-off: fewer pairs, cleaner" and b["B"]["accuracy"] == "less accurate", b["B"]
    assert b["D"]["retention"] == "more" and b["D"]["noise"] == "about the same" and b["D"]["verdict"] == "better"        # D has FIGARO's maxEE: 0.767 against 0.743 is +2.4, beyond the 2-point margin, at the same noise
    assert abs(b["D"]["diff"]["mean"] - 2.4) < 1e-6 and b["D"]["diff"]["n"] == 12 and b["D"]["diff"]["lo"] is not None                               # 12 samples: an interval
    # swabs: FIGARO's stricter maxEE; A merges more with noisier reads; a failed no-truncation arm is simply absent
    s = res["swabs"]; assert s["verdicts"]["A"]["verdict"] == "trade-off: more pairs, noisier" and "E" not in s["arms"] and s["errors"]
    # v1v2: the control D matches; mock A: found/strains compared
    assert res["v1v2"]["verdicts"]["D"]["verdict"] == "matches"
    ma_ = res["mockA_high"]["verdicts"]; assert ma_["A"]["accuracy"] == "less accurate" and ma_["D"]["accuracy"] == "as accurate"      # 30 against 32 variants is beyond the margin of 1; 31 against 32 is not
    # FIGARO gave no result on the tail dataset: arms are classified against A, the failure is not a criterion
    t = res["d6300_V4-V5"]; assert t["against"] == "A" and set(t["figaro_failed"]) == {"FIGARO (amplicon 550)", "FIGARO (amplicon 374)"} and "FIGARO gave no result at" in rep and "NOT MET" not in rep
    assert [c["name"] for c in t["coverage"]][0].endswith("(default 0.2%)") and len(t["coverage"]) == 3 and t["tail"]
    # stability is reported per tool; the scorecard separates mock and real
    assert "| FIGARO | 100/150; 1 ; 100/150; 1 ; 120/150; 1 | no |" in rep and "| AmpliTrim | 100/150; 2 ; 100/150; 2 ; 100/150; 2 | yes |" in rep
    assert "**mock communities (2 datasets)**" in rep and "**real communities (2 datasets)**" in rep and "trade-off: fewer pairs, cleaner" in rep
    print("full run: arms, commands, verdicts, FIGARO no-result, coverage, stability, scorecard OK")
    # restart: finished datasets are kept (AmpliTrim is not run again), the protocol is not rewritten, and --fresh repeats them
    calls2 = {"explore": [], "cmd": []}
    again, _ = make(tmp, calls2, setup=False); assert go(again) == 0
    assert calls2["explore"] == ["broken"] and not calls2["cmd"] and (out / "NIGHT2_PROTOCOL.md").read_text() == prot, (calls2["explore"], len(calls2["cmd"]))     # only the one that failed is tried again
    calls3 = {"explore": [], "cmd": []}
    fresh, _ = make(tmp, calls3, setup=False, fresh=True); assert go(fresh) == 0 and len(calls3["cmd"]) == 5
    print("restart: finished datasets kept, the failed one retried, --fresh repeats all")
    # a dry run runs nothing and writes nothing
    shutil.rmtree(data / "runs"); calls4 = {"explore": [], "cmd": []}
    dry, _ = make(tmp, calls4, dry=True, setup=False); assert go(dry) == 0
    assert not (data / "runs").exists() and not calls4["explore"] and not calls4["cmd"] or all("compare_methods.py" in " ".join(c) for c in calls4["cmd"])
    print("dry run: nothing written, AmpliTrim not run")
    # the disk guard stops before anything is created
    low, _ = make(tmp, {"explore": [], "cmd": []}, setup=False, disk=4e9); assert go(low) == 2 and not (data / "runs").exists()
    print("disk guard: stops below 6 GB and creates nothing")
