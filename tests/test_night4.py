"""validation/night4.py (arm H, the percentile rule) with stand-ins for compare_methods, AmpliTrim's own run and the quality table."""
import contextlib, io, json, sys
from pathlib import Path
R = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(Path(__file__).resolve().parent), str(R / "validation"), str(R)]
import night4 as n4
import test_night2 as t2
import test_night3 as t3
from amplitrim.trimmer import ReadPairs


class FakePairs(t3.FakePairs):
    """The night 3 stand-in with the number of pairs, which the rule needs."""
    def __init__(s, real): super().__init__(real); s.n = real.n


def skewed_pairs(tmp):
    """156 clean pairs and 44 whose reverse read is noisy: the 75th percentile of the reverse expected errors is clean, the 80th is not."""
    q37, q15 = chr(33 + 37), chr(33 + 15)
    rows1 = [("A" * 150, q37 * 150)] * 200
    rows2 = [("C" * 150, q37 * 150)] * 156 + [("C" * 150, q37 * 60 + q15 * 90)] * 44
    t2.write_fq(tmp / "s_R1.fq.gz", rows1); t2.write_fq(tmp / "s_R2.fq.gz", rows2)
    return ReadPairs.load(tmp / "s_R1.fq.gz", tmp / "s_R2.fq.gz")


def noisy_pairs(tmp):
    """Every reverse read carries 15 expected errors at 150 nt: the uncapped rule goes to 15 on it, the capped one to 5."""
    t2.write_fq(tmp / "n_R1.fq.gz", [("A" * 150, chr(33 + 37) * 150)] * 200); t2.write_fq(tmp / "n_R2.fq.gz", [("C" * 150, chr(33 + 10) * 150)] * 200)
    return ReadPairs.load(tmp / "n_R1.fq.gz", tmp / "n_R2.fq.gz")


A = t3.A; G = t3.G; FIG = t3.FIG
def H(ee, merged, **kw): return t2.row(f"AmpliTrim (maxEE {ee})", kw.get("tf", 248), kw.get("tr", 202), ee, merged, kw.get("found", ""), kw.get("strains", ""), kw.get("identical", ""))

TABLES4 = {      # mock B (very noisy reverse reads): H is (2,5), the cap, which is G's setting, and H80 is (2,15); V4-V5: both are (2,3); mock A: H is (2,2), so it is A, and H80 is (2,3)
    "mockB": lambda: [A(0.20, found=6, strains="3/3", identical=0.687), G(0.50, found=6, strains="3/3", identical=0.73), H("2,15", 0.70, found=6, strains="3/3", identical=0.76), FIG(515, 279, 248, "6,9", 0.743, found=7, strains="3/3", identical=0.7828)],
    "d6300_V4-V5": lambda: [A(0.572, found=9, strains="8/9", identical=0.988), G(0.70, found=9, strains="8/9", identical=0.988), H("2,3", 0.69, found=9, strains="8/9", identical=0.988)],        # FIGARO failed
    "mockA_high": lambda: [A(0.742, found=30, strains="11/13", identical=0.9), G(0.78, found=31, strains="11/13", identical=0.9), H("2,3", 0.77, found=31, strains="11/13", identical=0.9), FIG(508, 273, 247, "3,4", 0.78, found=31, strains="11/13", identical=0.9)],
}
INFO4 = {"mockB": dict(t2.INFO["mockB"], trim=(100, 150)), "d6300_V4-V5": dict(t2.INFO["d6300_V4-V5"], trim=(100, 150)), "mockA_high": dict(t2.INFO["mockA_high"], trim=(100, 150)),
         "v1v2": dict(t2.INFO["v1v2"], trim=(100, 60)), "swabs": dict(t2.INFO["swabs"], trim=(100, 60))}     # a clean cut: the rule gives (2,2)


def make4(tmp, calls, dry=False, setup=True, **kw):
    data, real = t2.world(tmp) if setup else (tmp / "data", None)
    if real is None:
        real = ReadPairs.load(tmp / "a_R1.fq.gz", tmp / "a_R2.fq.gz")
    skew, noisy = skewed_pairs(tmp), noisy_pairs(tmp)

    def explore(reads, work, marker):
        name = Path(reads).parent.name
        calls["explore"].append(name)
        if name == "broken": raise RuntimeError("AmpliTrim could not merge these reads")
        return INFO4[name]

    def run(cmd):
        cmd = [str(c) for c in cmd]; calls["cmd"].append(cmd)
        name = Path(cmd[2]).parent.name; out = Path(cmd[cmd.index("--out") + 1])
        t2.tsv(out / "compare_methods.tsv", TABLES4[name]())
        return 0, "ok"

    def shares(track):
        key = Path(track).name; comp = Path(track).parent
        rows = list(__import__("csv").DictReader(open(comp / "compare_methods.tsv"), delimiter="\t"))
        from compare_methods import method_key
        r = next(x for x in rows if "track_" + method_key(x["name"]) + ".tsv" == key)
        return {f"s{i}": float(r["merged"]) + 0.001 * (i - 6) for i in range(12)}
    n = n4.Night4(R, dry, False, kw.get("fresh", False), run=run, explore=explore, pairs=lambda reads: FakePairs({"mockA_high": skew, "mockB": noisy}.get(Path(reads).parent.name, real)), shares=shares,
                  disk_free=lambda p: 50e9, now=lambda f: "2026-10-10", cfg_path=tmp / "cfg.json")
    n.only = kw.get("only")
    return n, data


def go4(n):
    with contextlib.redirect_stdout(io.StringIO()): return n.go()


def test_night4_end_to_end(tmp_path):
    calls = {"explore": [], "cmd": []}
    n, data = make4(tmp_path, calls); assert go4(n) == 0
    assert n.logf is None, "the log file must be closed when the run ends: Windows cannot delete a file that is still open"
    out = data / "runs/night4_2026-10-10"; rep = (out / "NIGHT4_REPORT.md").read_text(); prot = (out / "NIGHT4_PROTOCOL.md").read_text()
    for need in ("before any analysis", "section 18", "P19.", "P20.", "P21.", "P22 (exploratory).", "H80", "never above 5", "The rule, exactly", "'mixed' fails", "expected to FAIL this on mock A", "Rule R is the rule of the release", "what the cap costs"): assert need in prot, need
    assert "H85" not in prot
    res = {d: json.loads((out / d / "night4_result.json").read_text()) for d in ("mockB", "d6300_V4-V5", "mockA_high", "v1v2", "swabs")}
    # the rule's maxEE is computed from the reads at AmpliTrim's cut, written into the plan, and given to compare_methods; arms are deduplicated
    cmd = {Path(c[2]).parent.name: c for c in calls["cmd"]}
    assert set(cmd) == {"mockB", "d6300_V4-V5", "mockA_high"}, "clean cuts and ITS run nothing"
    mb = cmd["mockB"]; i = mb.index("--amplitrim-max-ee"); assert mb[i + 1: i + 4] == ["2", "2,5", "2,15"] and mb[i + 4].startswith("--")          # H is (2,5) (the cap: G's setting, one arm), H80 is (2,15)
    ma_ = cmd["mockA_high"]; j = ma_.index("--amplitrim-max-ee"); assert ma_[j + 1: j + 4] == ["2", "2,5", "2,3"]                                    # q 0.75 gives (2,2) (arm A), q 0.80 gives (2,3)
    v4 = cmd["d6300_V4-V5"]; k = v4.index("--amplitrim-max-ee"); assert v4[k + 1: k + 4] == ["2", "2,5", "2,3"] and v4[k + 4].startswith("--")          # both rules give (2,3): one arm
    assert ma_[ma_.index("--figaro-amplicon") + 1] == "508" and cmd["d6300_V4-V5"][cmd["d6300_V4-V5"].index("--figaro-amplicon") + 1] == "550"     # FIGARO is asked at each longest amplicon (it fails at 550)
    assert all("--reference" in cmd[d] for d in cmd) and cmd["mockA_high"][cmd["mockA_high"].index("--pairs") + 1] == "0" and cmd["mockB"][cmd["mockB"].index("--pairs") + 1] == "5000"
    assert res["mockB"]["plan"]["H"] == "2,5" and res["mockB"]["plan"]["H80"] == "2,15" and res["mockB"]["plan"]["cut"] == [100, 150]           # the cap: 16 becomes 5
    assert res["mockA_high"]["plan"]["H"] == "2" and res["mockA_high"]["plan"]["H80"] == "2,3"
    assert res["mockA_high"]["arms"]["H"]["maxEE"] == "2" and res["mockA_high"]["arms"]["H80"]["maxEE"] == "2,3"                                   # H is the default's row on mock A
    assert res["d6300_V4-V5"]["plan"]["H"] == "2,3" and res["d6300_V4-V5"]["plan"]["H80"] == "2,3"
    # clean cuts: the rule gives (2,2) for both q, so the arm is A by construction and nothing is run
    for d in ("v1v2", "swabs"): assert "skipped" in res[d] and "identical to A by construction" in res[d]["skipped"]
    assert "swabs: the rule gives (2,2)" in rep and "## Not run" in rep
    # the differences from G: mock B 0 (H is G's setting), V4-V5 -1 (passes), mock A -3.8 (a miss beyond the margin of 2: reported, not hidden)
    assert abs(res["mockB"]["vs_G"]["H"]["diff"]["mean"]) < 1e-6 and abs(res["d6300_V4-V5"]["vs_G"]["H"]["diff"]["mean"] + 1.0) < 1e-6 and abs(res["mockA_high"]["vs_G"]["H"]["diff"]["mean"] + 3.8) < 1e-6
    assert abs(res["mockB"]["cap_cost"]["mean"] - 20.0) < 1e-6 and "What the cap costs (H80 minus H, merged pairs): +20.0 points" in rep       # uncapped (2,15) merges 20 points more than the capped (2,5)
    assert res["mockB"]["adopt"] == {"merged": True, "accuracy": True} and res["d6300_V4-V5"]["adopt"] == {"merged": True, "accuracy": True}
    assert res["mockA_high"]["adopt"]["merged"] is False and res["mockA_high"]["adopt"]["accuracy"] is True          # one variant fewer is inside the margin of 1
    # accuracy against G: mock B as accurate (the same setting), V4-V5 as accurate; FIGARO gave no result at 550, so V4-V5 is classified against A
    assert res["mockB"]["vs_G"]["H"]["accuracy"] == "as accurate" and res["d6300_V4-V5"]["vs_G"]["H"]["accuracy"] == "as accurate" and res["mockB"]["vs_G"]["H80"]["accuracy"] == "more accurate"
    assert res["d6300_V4-V5"]["against"] == "A" and res["d6300_V4-V5"]["figaro_failed"] == ["FIGARO (amplicon 550)"] and "FIGARO gave no result at: FIGARO (amplicon 550)" in rep
    # the predictions, decided by the data: P19 fails on mock A alone, P20 holds, P21 cannot be decided with one of its two datasets, P22 holds on mock B
    pred = {p: line for p in ("P19", "P20", "P21", "P22") for line in rep.splitlines() if line.startswith(f"- **{p}**")}
    assert "NOT MET" in pred["P19"] and "mockA_high -3.8" in pred["P19"] and "mockB +0.0" in pred["P19"]
    assert "MET" in pred["P20"] and "NOT MET" not in pred["P20"] and "not decided" in pred["P21"] and "MET" in pred["P22"] and "NOT MET" not in pred["P22"] and "mockB +20.0" in pred["P22"]
    assert "## The rule's maxEE and the gain over rule R's setting" in rep and "| mockB | 2,5 / 2,15 |" in rep and "| mockA_high | 2 / 2,3 |" in rep
    assert "broken: failed: RuntimeError" in rep and "Problems: 1" in rep and not (out / "broken/night4_result.json").exists()
    print("full run: the rule's maxEE from the reads, arms and commands, skips, adoption check, predictions OK")
    # restart: finished datasets (including the skipped ones) are kept, the failed one is retried, --fresh repeats all
    calls2 = {"explore": [], "cmd": []}
    again, _ = make4(tmp_path, calls2, setup=False); assert go4(again) == 0
    assert calls2["explore"] == ["broken"] and not calls2["cmd"] and (out / "NIGHT4_PROTOCOL.md").read_text() == prot
    calls3 = {"explore": [], "cmd": []}
    fresh, _ = make4(tmp_path, calls3, setup=False, fresh=True); assert go4(fresh) == 0 and len(calls3["cmd"]) == 3


def test_night4_dry_run_prints_commands_and_writes_nothing(tmp_path):
    calls = {"explore": [], "cmd": []}
    n, data = make4(tmp_path, calls, dry=True); assert go4(n) == 0
    assert calls["explore"] == [] and not (data / "runs").exists()
    assert len(calls["cmd"]) == 0                       # dry runs print through run_fn, which does not call the executor
    n2, _ = make4(tmp_path, {"explore": [], "cmd": []}, dry=True, setup=False, only=["mockB"])
    assert [d["name"] for d in n2.datasets()] == ["mockB"]


def test_night4_only_runs_the_datasets_it_is_given(tmp_path):
    calls = {"explore": [], "cmd": []}
    n, data = make4(tmp_path, calls, only=["mockB"]); assert go4(n) == 0
    assert calls["explore"] == ["mockB"] and [Path(c[2]).parent.name for c in calls["cmd"]] == ["mockB"]


def test_the_rule_in_the_script_is_the_one_the_register_fixed():
    assert (n4.Q_RULE, n4.CAP, n4.Q_ALT) == (0.75, 5, 0.80) and n4.PASSING_ACCURACY == ("as accurate", "more accurate")
    import maxee_rules as mr
    assert mr.QGRID == (0.75, 0.80, 0.85, 0.90) and n4.Q_RULE in mr.QGRID and n4.Q_ALT in mr.QGRID
