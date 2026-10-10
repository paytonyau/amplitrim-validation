"""validation/night3.py (arms G, R and F2) with stand-ins for compare_methods, AmpliTrim's own run and the quality table, and F2's rule checked by hand."""
import contextlib, gzip, io, json, shutil, sys
from pathlib import Path
R = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(Path(__file__).resolve().parent), str(R / "validation"), str(R)]
import maxee_arms as ma
import night3 as n3
import test_night2 as t2
from amplitrim.trimmer import ReadPairs


def pairs_of(tmp, q_f, q_r, lengths=None, n=200):
    """Reads whose quality is given by strings: q_f for the forward and q_r for the reverse read (one string per read length)."""
    tmp.mkdir(parents=True, exist_ok=True)
    f1, f2 = tmp / "p_R1.fq.gz", tmp / "p_R2.fq.gz"
    lengths = lengths or [(len(q_f), len(q_r))] * n
    with gzip.open(f1, "wt") as a, gzip.open(f2, "wt") as b:
        for i, (lf, lr) in enumerate(lengths):
            a.write(f"@p{i}\n{'A' * lf}\n+\n{q_f[:lf]}\n")
            b.write(f"@p{i}\n{'C' * lr}\n+\n{q_r[:lr]}\n")
    return ReadPairs.load(f1, f2)


def q(score, n):
    return chr(33 + score) * n


def test_f2_follows_the_four_steps_of_the_register(tmp_path):
    P = pairs_of(tmp_path / "a", q(37, 100) + q(10, 50), q(37, 60) + q(15, 90))                       # median quality falls below 30 at 100 (forward) and 60 (reverse)
    assert ma.q30_cut(P) == (100, 60)
    # 1. the amplicon fits in the Q30 cut: nothing is added; the pairs pass at (2,2), so maxEE stays 2
    a = ma.f2_plan(P, 100)
    assert a["cut"] == (100, 60) and a["reaches"] and a["need"] == 112 - 0 and a["ee"] == "2" and a["pass22"] == 1.0 and a["caps"] == (150, 150)
    # 3. too short for the amplicon: the missing 52 nt are split in halves (26 + 26), and the longer reads then fail at (2,2), so maxEE is relaxed to (2,5)
    b = ma.f2_plan(P, 200)
    assert b["cut"] == (126, 86) and b["reaches"] and b["pass22"] < ma.F2_RELAX_BELOW and b["ee"] == "2,5"
    # the forward read can only take 50 more, so the reverse takes the rest of an odd split: 102 missing -> forward 50 (its cap), reverse 52
    c = ma.f2_plan(P, 250)
    assert c["cut"] == (150, 112) and c["cut"][0] + c["cut"][1] == 262 and c["reaches"]
    # an amplicon the reads cannot reach: the cut stays at the caps and is flagged
    d = ma.f2_plan(P, 515)
    assert d["cut"] == (150, 150) and not d["reaches"] and d["need"] == 527
    # an odd number of missing nucleotides: the forward read takes the odd one
    e = ma.f2_plan(P, 153)                                         # need 165, Q30 total 160, 5 missing
    assert e["cut"] == (103, 62)


def test_f2_caps_a_cut_at_the_length_that_90_percent_of_the_reads_reach(tmp_path):
    lengths = [(150, 150)] * 160 + [(100, 120)] * 40                  # 20% of the pairs are shorter, so the 10th percentile is the short length
    P = pairs_of(tmp_path / "b", q(37, 150), q(37, 150), lengths=lengths)
    plan = ma.f2_plan(P, 200)
    assert plan["caps"] == (100, 120) and plan["q30"] == (150, 150)
    assert plan["cut"] == (100, 120)                                  # Q30 says 150, the reads do not all reach it; 220 >= 212, so nothing is added


# ------------------------------------------------------------------------------------------------ the run, with stand-ins
def A(merged, **kw): return t2.row("AmpliTrim (maxEE 2)", kw.get("tf", 250), kw.get("tr", 200), "2", merged, kw.get("found", ""), kw.get("strains", ""), kw.get("identical", ""))
def G(merged, **kw): return t2.row("AmpliTrim (maxEE 2,5)", kw.get("tf", 248), kw.get("tr", 202), "2,5", merged, kw.get("found", ""), kw.get("strains", ""), kw.get("identical", ""))
def F2(merged, **kw): return t2.row("f2", 150, 150, "2,5", merged, kw.get("found", ""), kw.get("strains", ""), kw.get("identical", ""))
def FIG(amp, tf, tr, ee, merged, **kw): return t2.row(f"FIGARO (amplicon {amp})", tf, tr, ee, merged, kw.get("found", ""), kw.get("strains", ""), kw.get("identical", ""))


TABLES3 = {
    "mockB": lambda: [A(0.20, found=6, strains="3/3", identical=0.687), G(0.50, found=6, strains="3/3", identical=0.73), F2(0.30, found=5, strains="3/3", identical=0.60), FIG(515, 279, 248, "6,9", 0.743, found=7, strains="3/3", identical=0.7828)],
    "d6300_V4-V5": lambda: [A(0.572, found=9, strains="8/9", identical=0.988), G(0.70, found=9, strains="8/9", identical=0.988), F2(0.62, found=9, strains="8/9", identical=0.988)],     # FIGARO failed at 550
    "swabs": lambda: [A(0.77), G(0.79), t2.row("f2", 150, 150, "2", "", err="DADA2 could not run: all reads removed"), FIG(430, 224, 218, "1", 0.711)],
    "v1v2": lambda: [A(0.842), G(0.86), F2(0.70), FIG(346, 225, 133, "2,1", 0.82)],
    "mockA_high": lambda: [A(0.742, found=30, strains="11/13", identical=0.9), G(0.78, found=31, strains="11/13", identical=0.9), F2(0.60, found=20, strains="8/13", identical=0.8), FIG(508, 273, 247, "3,4", 0.758, found=32, strains="11/13", identical=0.9)],
}


class FakePairs:
    """A quality table whose pass rate and noise are simple: the real pairs give lengths and qualities; looser maxEE keeps noisier reads."""
    def __init__(s, real): s.r = real; s._phred1, s._phred2, s.len1, s.len2 = real._phred1, real._phred2, real.len1, real.len2
    def pass_table(s, ee=2.0): return s.r.pass_table(ee)
    def error_load(s, tf, tr, ee): return (0.1 * ee[0], 0.1 * ee[1])


def make3(tmp, calls, dry=False, setup=True, **kw):
    data, real = t2.world(tmp) if setup else (tmp / "data", None)
    if real is None:
        real = ReadPairs.load(tmp / "a_R1.fq.gz", tmp / "a_R2.fq.gz")

    def explore(reads, work, marker):
        calls["explore"].append(Path(reads).parent.name)
        if Path(reads).parent.name == "broken": raise RuntimeError("AmpliTrim could not merge these reads")
        return t2.INFO[Path(reads).parent.name]

    def run(cmd):
        cmd = [str(c) for c in cmd]; calls["cmd"].append(cmd)
        name = Path(cmd[2]).parent.name; out = Path(cmd[cmd.index("--out") + 1])
        t2.tsv(out / "compare_methods.tsv", TABLES3[name]())
        return 0, "ok"

    def shares(track):
        key = Path(track).name; comp = Path(track).parent
        rows = list(__import__("csv").DictReader(open(comp / "compare_methods.tsv"), delimiter="\t"))
        from compare_methods import method_key
        r = next(x for x in rows if "track_" + method_key(x["name"]) + ".tsv" == key)
        return {f"s{i}": float(r["merged"]) + 0.001 * (i - 6) for i in range(12)}
    n = n3.Night3(R, dry, False, kw.get("fresh", False), run=run, explore=explore, pairs=lambda reads: FakePairs(real), shares=shares, disk_free=lambda p: kw.get("disk", 50e9),
                  now=lambda f: "2026-10-10", cfg_path=tmp / "cfg.json")
    return n, data


def go3(n):
    with contextlib.redirect_stdout(io.StringIO()): return n.go()


def test_night3_end_to_end(tmp_path):
    calls = {"explore": [], "cmd": []}
    n, data = make3(tmp_path, calls); assert go3(n) == 0
    assert n.logf is None, "the log file must be closed when the run ends: Windows cannot delete a file that is still open"
    out = data / "runs/night3_2026-10-10"; rep = (out / "NIGHT3_REPORT.md").read_text(); prot = (out / "NIGHT3_PROTOCOL.md").read_text()
    assert (out / "night3.log").is_file() and not (out / "NIGHT2_PROTOCOL.md").exists()
    for need in ("before any analysis", "section 14", "P6.", "P7.", "P9.", "P12.", "F2, exactly", "-m 12", "if P7 and P9 hold"): assert need in prot, need
    assert "broken: failed: RuntimeError" in rep and "Problems: 1" in rep and not (out / "broken/night3_result.json").exists()
    res = {d: json.loads((out / d / "night3_result.json").read_text()) for d in ("mockB", "d6300_V4-V5", "swabs", "v1v2", "mockA_high")}
    # the commands: A and G through AmpliTrim, F2 as a fixed setting, FIGARO at one amplicon length, nothing from night 2's extras
    cmd = {Path(c[2]).parent.name: c for c in calls["cmd"]}
    mb = cmd["mockB"]; i = mb.index("--amplitrim-max-ee"); assert mb[i + 1: i + 3] == ["2", "2,5"]
    assert "--allow-unmergeable" not in mb and not any(x.startswith(("q30=", "full=")) for x in mb) and mb[mb.index("--figaro-amplicon") + 1] == "515" and mb[mb.index("--figaro-amplicon") + 2].startswith("--")
    f2 = [x for x in mb if x.startswith("f2=")]; assert f2 == ["f2=150,150,2,5"]                      # the stand-in reads are 150 nt and cannot reach a 515 nt amplicon: the cut stays at the caps, and almost nothing passes at (2,2)
    assert cmd["mockA_high"][cmd["mockA_high"].index("--pairs") + 1] == "0" and cmd["swabs"][cmd["swabs"].index("--pairs") + 1] == "5000"
    # F2's plan is recorded and flagged where the reads cannot reach the amplicon
    assert res["mockB"]["f2"]["reaches"] is False and res["mockB"]["f2"]["cut"] == [150, 150] and "DOES NOT REACH THE AMPLICON" in rep
    # mock B: A and G both merge fewer pairs than FIGARO, with cleaner reads; G gains 30 points over A and is more accurate
    b = res["mockB"]; assert b["against"] == "C"
    assert b["verdicts"]["A"]["verdict"] == "trade-off: fewer pairs, cleaner" and b["verdicts"]["G"]["verdict"] == "trade-off: fewer pairs, cleaner"
    assert abs(b["vs_A"]["G"]["diff"]["mean"] - 30.0) < 1e-6 and b["vs_A"]["G"]["noise"] == "noisier" and b["vs_A"]["G"]["accuracy"] == "more accurate" and b["vs_A"]["G"]["diff"]["lo"] is not None
    assert abs(b["verdicts"]["G"]["diff"]["mean"] - (50.0 - 74.3)) < 1e-6
    # rule R: (2,5) where A passes under 70% (mock B 22%, V4-V5 59%), A elsewhere
    assert b["rule"]["uses"] == "G" and res["d6300_V4-V5"]["rule"]["uses"] == "G" and res["v1v2"]["rule"]["uses"] == "A" and res["swabs"]["rule"]["uses"] == "A" and res["mockA_high"]["rule"]["uses"] == "A"
    assert b["rule"]["verdict"] == "trade-off: fewer pairs, cleaner" and res["v1v2"]["rule"]["verdict"] == (res["v1v2"]["verdicts"]["A"]["verdict"])
    # FIGARO gave no result at 550: arms are classified against A, and the failure is not a criterion
    t = res["d6300_V4-V5"]; assert t["against"] == "A" and t["figaro_failed"] == ["FIGARO (amplicon 550)"] and "NOT MET" not in rep and "FIGARO gave no result at: FIGARO (amplicon 550)" in rep
    assert "A" not in t["verdicts"] and abs(t["vs_A"]["G"]["diff"]["mean"] - 12.8) < 1e-6
    # an arm DADA2 could not run is absent, not a zero
    s = res["swabs"]; assert "F2" not in s["arms"] and any("all reads removed" in e for e in s["errors"])
    # mock A: G is as accurate as the default (one more variant, inside the margin of 1)
    assert res["mockA_high"]["vs_A"]["G"]["accuracy"] == "as accurate" and res["mockA_high"]["vs_A"]["F2"]["accuracy"] == "less accurate"
    # the report: scorecard with the rule, and the gain table, mock and real apart
    assert "**mock communities (2 datasets)**" in rep and "**real communities (2 datasets)**" in rep and "| R (uses) |" in rep and "## Gain over the default" in rep
    print("full run: arms, commands, F2 plan, verdicts, rule R, FIGARO no-result, gain table OK")
    # restart: finished datasets kept, the failed one retried, --fresh repeats all
    calls2 = {"explore": [], "cmd": []}
    again, _ = make3(tmp_path, calls2, setup=False); assert go3(again) == 0
    assert calls2["explore"] == ["broken"] and not calls2["cmd"] and (out / "NIGHT3_PROTOCOL.md").read_text() == prot
    calls3 = {"explore": [], "cmd": []}
    fresh, _ = make3(tmp_path, calls3, setup=False, fresh=True); assert go3(fresh) == 0 and len(calls3["cmd"]) == 5
    print("restart: finished datasets kept, the failed one retried, --fresh repeats all")
    # a dry run runs nothing and writes nothing
    shutil.rmtree(data / "runs"); calls4 = {"explore": [], "cmd": []}
    dry, _ = make3(tmp_path, calls4, dry=True, setup=False); assert go3(dry) == 0
    assert not calls4["cmd"] and not (data / "runs").exists(); print("dry run: nothing written, AmpliTrim not run")
    # the disk guard stops before anything is created
    low, _ = make3(tmp_path, {"explore": [], "cmd": []}, setup=False, disk=5e9); assert go3(low) == 2 and not (data / "runs").exists(); print("disk guard: stops below 6 GB and creates nothing")


def test_night3_leaves_night2_s_names_alone(tmp_path):
    assert n3.Night3.TAG == "night3" and n3.Night3.HEAD.startswith("night3.py")
    import night2
    assert night2.Night2.TAG == "night2" and night2.Night2.HEAD == "night2.py: the maxEE head-to-head"
