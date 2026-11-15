"""validation/its_night.py: the helpers, and the whole night with stand-ins for R, vsearch, AmpliTrim and the thinning (no real program is run)."""
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent


def load():
    spec = importlib.util.spec_from_file_location("its_night", HERE.parent / "validation" / "its_night.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


N = load()


def fq(path, records, gz=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f"@r{i}\n{s}\n+\n{q}\n" for i, (s, q) in enumerate(records))
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


# ------------------------------------------------------------------ helpers
def test_the_sliding_window_cuts_before_the_first_window_below_15():
    q = "I" * 10 + "#" * 4 + "I" * 6                                        # I = 40, # = 2
    assert N.sliding_cut(q) == 9          # windows at 7 and 8 have means 30.5 and 21; the window at 9 holds 40, 2, 2, 2: mean 11.5, the first below 15
    assert N.sliding_cut("I" * 20) == 20 and N.sliding_cut("####") == 0 and N.sliding_cut("II") == 2


def test_trim_folder_cuts_every_read_and_keeps_the_pairs_in_order(tmp_path):
    fq(tmp_path / "in" / "s_R1.fastq.gz", [("ACGTACGTAC", "IIIIII####"), ("ACGT", "IIII")])
    fq(tmp_path / "in" / "s_R2.fastq.gz", [("TTTTTTTTTT", "IIIIIIIIII"), ("GGGG", "####")])
    assert N.trim_folder(tmp_path / "in", tmp_path / "out") == 2
    r1 = list(N.fastq_records(tmp_path / "out" / "s_R1.fastq.gz"))
    r2 = list(N.fastq_records(tmp_path / "out" / "s_R2.fastq.gz"))
    assert [(r[1], r[3]) for r in r1] == [("ACGTA", "IIIII"), ("ACGT", "IIII")] and [r[1] for r in r2] == ["TTTTTTTTTT", ""]
    assert [r[0] for r in r1] == ["@r0", "@r1"] == [r[0] for r in r2]


def test_the_arms_table_runs_std_first_and_the_reusing_arms_after_it(tmp_path):
    rows = N.arms_rows(tmp_path / "x", True)
    arms = [r[0] for r in rows]
    assert arms[0] == "std" and len(arms) == 13 and arms.index("fwd") > 0 and arms.index("mo30") > 0
    by = {r[0]: r for r in rows}
    assert by["std"][2] == str(tmp_path / "x" / "reads") and by["sw15"][2] == str(tmp_path / "x" / "reads_sw15")
    assert by["mfirst"][1:3] == ["single", str(tmp_path / "x" / "reads_mfirst")] and by["fwd"][1:3] == ["fwd", "std"]
    assert by["ee25"][5:7] == [2, 5] and by["tq20"][7] == 20 and by["mo20"][8] == 20
    assert N.arms_rows(tmp_path / "y", False, (240, 160)) == [["std", "paired", str(tmp_path / "y" / "reads"), 240, 160, 2, 2, 2, 12]]


def test_the_vsearch_command_allows_stagger_and_sets_the_differences_only_when_asked():
    c = N.vsearch_merge_cmd(Path("a"), Path("b"), Path("o"))
    assert "--fastq_allowmergestagger" in c and c[c.index("--fastq_minovlen") + 1] == "12" and "--fastq_maxdiffs" not in c
    assert N.vsearch_merge_cmd(Path("a"), Path("b"), Path("o"), 0)[-2:] == ["--fastq_maxdiffs", "0"]


# ------------------------------------------------------------------ the whole night with stand-ins
class World:
    """Stand-ins: R writes DONE for every arm and copies the reads as std's filtered reads; vsearch merges fewer pairs with fewer differences; AmpliTrim writes a
    summary with a cut. fail: (dataset label, program) pairs that stop with an error."""

    MERGED = {None: 8, 10: 8, 5: 7, 2: 5, 0: 2}

    def __init__(self, fail=()):
        self.calls, self.fail = [], set(fail)

    def __call__(self, cmd, capture_output=True, text=True, cwd=None):
        self.calls.append((cmd, cwd))
        prog = Path(cmd[0]).name
        if any(lab in " ".join(cmd) and p == prog for lab, p in self.fail):
            return SimpleNamespace(returncode=1, stdout="", stderr=f"{prog} broke")
        if prog == "vsearch":
            k = int(cmd[cmd.index("--fastq_maxdiffs") + 1]) if "--fastq_maxdiffs" in cmd else None
            fq(Path(cmd[cmd.index("--fastqout") + 1]), [("ACGT", "IIII")] * self.MERGED[k], gz=False)
        elif prog == "Rscript-fake":
            arms, out = Path(cmd[2]), Path(cmd[3])
            reads = arms.parent / "reads"
            for line in arms.read_text().splitlines()[1:]:
                arm = line.split("\t")[0]
                (out / arm).mkdir(parents=True, exist_ok=True)
                (out / arm / "DONE").write_text("")
                if arm == "std":
                    (out / arm / "filtered").mkdir(exist_ok=True)
                    for f in reads.glob("*.fastq.gz"):
                        (out / arm / "filtered" / f.name).write_bytes(f.read_bytes())
        elif "-c" in cmd:                                                  # AmpliTrim
            rep = Path(cmd[cmd.index("--out-dir") + 1])
            rep.mkdir(parents=True, exist_ok=True)
            (rep / "summary.json").write_text(json.dumps({"recommended": {"truncLenF": 0 if "--its" in cmd else 240, "truncLenR": 0 if "--its" in cmd else 160}}))
        return SimpleNamespace(returncode=0, stdout="", stderr="")


def fake_thin(src, dst, pairs):
    for f in sorted(Path(src).glob("*.fastq.gz")):
        recs = list(N.fastq_records(f))[:pairs]
        fq(Path(dst) / f.name, [(r[1], r[3]) for r in recs])


def data(base, labels):
    for lab in labels:
        for s in ("A", "B"):
            for m in (1, 2):
                fq(base / lab / f"{s}_R{m}.fastq.gz", [("ACGTACGTAC", "IIIIIIIIII")] * 6)


def go(tmp_path, world, extra=(), sixteen=None, version="0.2.0"):
    base = tmp_path / "data"
    its = {"atcc_its1": base / "a1", "atcc_its2": base / "a2"}
    data(base, ["a1", "a2", "m16"])
    args = ["--out", str(tmp_path / "night"), "--base", str(base), "--only", "atcc_its1,atcc_its2,16s_mockA_high", "--pairs", "4"]
    for k, v in its.items():
        args += ["--dataset", f"{k}={v}"]
    logged = []
    rc = N.main(args + list(extra), run=world, thin=fake_thin, rscript="Rscript-fake",
                sixteen=sixteen if sixteen is not None else [{"name": "mockA_high", "reads": base / "m16", "pairs": 0}],
                amplitrim_call="CALL", installed=lambda run, cwd: (version, "/site/amplitrim"), log=logged.append)
    return rc, logged


def test_the_whole_night_runs_every_step_and_leaves_its_files(tmp_path):
    world = World()
    rc, logged = go(tmp_path, world)
    night = tmp_path / "night"
    assert rc == 0, "\n".join(logged)
    status = json.loads((night / "night_status.json").read_text())
    assert all(v["state"] == "done" for v in status.values()) and len(status) == 6 + 6 + 4
    a1 = night / "atcc_its1"
    assert len(list((a1 / "reads").glob("*.fastq.gz"))) == 4 and len(list(N.fastq_records(a1 / "reads" / "A_R1.fastq.gz"))) == 4      # thinned to 4 pairs
    assert (a1 / "reads_sw15" / "A_R2.fastq.gz").is_file() and (a1 / "reads_mfirst" / "A_R1.fastq").is_file()
    est = (a1 / "estimates.tsv").read_text().splitlines()
    assert est[0] == "sample\tfiltered_pairs\tmerged_d10\tmerged_d5\tmerged_d2\tmerged_d0" and est[1] == "A\t4\t8\t7\t5\t2"
    rows16 = (night / "16s_mockA_high" / "arms.tsv").read_text().splitlines()
    assert len(rows16) == 2 and rows16[1].split("\t")[3:5] == ["240", "160"]                     # the 16S control at AmpliTrim's cut
    assert len(list((night / "16s_mockA_high" / "reads").glob("*.fastq.gz"))) == 4 and len(list(N.fastq_records(night / "16s_mockA_high" / "reads" / "A_R1.fastq.gz"))) == 6
    amp = [(c, cwd) for c, cwd in world.calls if "-c" in c]
    assert amp and all(cwd == str(night) for _, cwd in amp)                                     # started from --out, never from the archive
    assert sum("--its" in c for c, _ in amp) == 2 and len(amp) == 3
    assert "Next: python validation\\its_night_report.py" in logged[-1] and (night / "NIGHT_STATUS.md").is_file()


def test_a_second_start_repeats_nothing(tmp_path):
    world = World()
    assert go(tmp_path, world)[0] == 0
    n = len(world.calls)
    rc, logged = go(tmp_path, world)
    assert rc == 0 and len(world.calls) == n and sum("done before" in x for x in logged) == 16


def test_a_failure_stops_only_its_dataset_and_a_second_start_finishes_the_rest(tmp_path):
    world = World(fail={("atcc_its2", "Rscript-fake")})
    rc, logged = go(tmp_path, world)
    status = json.loads((tmp_path / "night" / "night_status.json").read_text())
    assert rc == 2 and status["atcc_its2/arms"]["state"] == "failed" and "broke" in status["atcc_its2/arms"]["error"]
    assert status["atcc_its2/amplitrim"]["state"] == status["atcc_its2/estimates"]["state"] == "skipped"
    assert status["atcc_its1/estimates"]["state"] == "done" and status["16s_mockA_high/estimates"]["state"] == "done"
    world.fail.clear()
    n = len(world.calls)
    assert go(tmp_path, world)[0] == 0
    redone = [c for c, _ in world.calls[n:]]
    assert redone and all("atcc_its2" in " ".join(c) for c in redone)                          # only what was not done


def test_another_amplitrim_version_stops_the_night_before_anything_runs(tmp_path):
    world = World()
    rc, logged = go(tmp_path, world, version="0.1.0")
    assert rc == 2 and world.calls == [] and "not 0.2.0" in logged[-1]


def test_dry_run_lists_and_runs_nothing(tmp_path):
    world = World()
    rc, logged = go(tmp_path, world, extra=["--dry-run"])
    assert rc == 0 and world.calls == [] and "3 datasets. Nothing was run." in logged[-1]
    assert any("16S" in x and "every pair" in x for x in logged)


def test_a_missing_folder_is_named(tmp_path):
    world = World()
    rc, logged = go(tmp_path, world, extra=["--dataset", f"atcc_its2={tmp_path / 'nowhere'}"])
    assert rc == 2 and "nowhere" in logged[-1] and world.calls == []
    rc, logged = go(tmp_path, world, extra=["--dataset", f"atcc_its2={tmp_path / 'nowhere'}", "--dry-run"])
    assert rc == 0 and any("nowhere" in x and "NOT FOUND" in x for x in logged)                  # the dry run marks it


def test_a_wrong_dataset_label_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="LABEL=FOLDER"):
        go(tmp_path, World(), extra=["--dataset", "nonsense=C:/x"])



def test_the_smoke_run_thins_every_dataset_mock_a_included(tmp_path):
    rc, logged = go(tmp_path, World(), extra=["--smoke", "--dry-run"])
    assert rc == 0 and sum("300 pairs per sample" in l for l in logged) == 3 and not any("every pair" in l for l in logged)


def test_a_16s_control_that_amplitrim_cannot_analyse_takes_the_final_runs_cut(tmp_path):
    final = tmp_path / "final_cli"
    (final / "mockA_high" / "report").mkdir(parents=True)
    (final / "mockA_high" / "report" / "summary.json").write_text(json.dumps({"recommended": {"truncLenF": 270, "truncLenR": 220}}))
    py = Path(sys.executable).name
    world = World(fail={("16s_mockA_high", py)})
    rc, logged = go(tmp_path, world, extra=["--final-run", str(final)])
    d = tmp_path / "night" / "16s_mockA_high"
    assert rc == 0 and (d / "amplitrim" / "CUT_FROM_FINAL_RUN.txt").is_file()
    assert (d / "arms.tsv").read_text().splitlines()[1].split("\t")[3:5] == ["270", "220"]
    assert any("the final run's cut is used" in l for l in logged)
    world2 = World(fail={("16s_mockA_high", py)})                                              # no fallback file: the step fails as before
    rc, _ = go(tmp_path / "b", world2, extra=["--final-run", str(tmp_path / "nowhere")])
    status = json.loads((tmp_path / "b" / "night" / "night_status.json").read_text())
    assert rc == 2 and status["16s_mockA_high/amplitrim"]["state"] == "failed"


# ------------------------------------------------------------------ the rules night (--plan rules, register section 35)
def test_the_trim_strength_cuts_more_the_higher_it_is(tmp_path):
    q = "I" * 20 + "5" * 8 + "+" * 8 + "I" * 4                                     # 40, then 20, then 10, then 40
    fq(tmp_path / "in" / "s_R1.fastq.gz", [("A" * 40, q)])
    fq(tmp_path / "in" / "s_R2.fastq.gz", [("A" * 40, q)])
    kept = {}
    for strength in N.SW_STRENGTHS:
        N.trim_folder(tmp_path / "in", tmp_path / f"sw{strength}", strength)
        kept[strength] = len(list(N.fastq_records(tmp_path / f"sw{strength}" / "s_R1.fastq.gz"))[0][1])
    assert kept[10] >= kept[15] >= kept[20] >= kept[25] and kept[10] > kept[25]
    assert (kept[10], kept[15], kept[20], kept[25]) == (40, 27, 25, 20)                 # 25: the window at base 20 averages 20; 20: at 25, 17.5; 15: at 27, 12.5; 10: none falls below 10


def test_the_computed_maxee_is_amplitrims_rule_on_reads_cut_at_their_first_q2():
    import numpy as np
    assert N.read_ee("I" * 10 + "#" + "I" * 10) == (10, 10 * 10 ** -4)              # cut before the first Q2 ('#'); EE of the 10 bases kept
    assert N.read_ee("+" * 3)[0] == 3 and abs(N.read_ee("+" * 3)[1] - 0.3) < 1e-9     # Q10: 0.1 each


def test_ee_rule_on_a_folder(tmp_path):
    import numpy as np
    quals_f = ["I" * 60, "5" * 60, "+" * 60, "+" * 30 + "I" * 30]                     # EE 0.006, 0.6, 6, 3.003
    quals_r = ["+" * 60, "+" * 60, "+" * 60, "+" * 40]                               # the last keeps only 40 bases: left out (minLen 50)
    fq(tmp_path / "r" / "a_R1.fastq.gz", [("A" * len(q), q) for q in quals_f])
    fq(tmp_path / "r" / "a_R2.fastq.gz", [("A" * len(q), q) for q in quals_r])
    res = N.ee_rule(tmp_path / "r")
    ef = [N.read_ee(q)[1] for q in quals_f]
    expected_f = max(2, int(float(np.percentile(ef, 75))) + 1)
    assert res["reads"] == [4, 3] and res["max_ee"][0] == expected_f and res["max_ee"][1] == 7     # reverse: 6, 6, 6 -> P75 6.0 -> 7


class RulesWorld(World):
    pass


def test_the_rules_night_runs_its_steps_and_writes_the_arms_with_the_computed_maxee(tmp_path):
    world = RulesWorld()
    base = tmp_path / "data"
    for lab in ("a1", "a2"):                                                          # reads of 60 bases: the computed maxEE needs reads of at least 50
        for smp in ("A", "B"):
            for m in (1, 2):
                fq(base / lab / f"{smp}_R{m}.fastq.gz", [("ACGT" * 15, "I" * 50 + "+" * 10)] * 6)
    logged = []
    rc = N.main(["--out", str(tmp_path / "night"), "--base", str(base), "--plan", "rules", "--only", "atcc_its1,mock3", "--pairs", "4",
                 "--dataset", f"atcc_its1={base / 'a1'}", "--extra", f"mock3={base / 'a2'}"],
                run=world, thin=fake_thin, rscript="Rscript-fake", sixteen=[], log=logged.append)
    assert rc == 0, "\n".join(logged)
    d = tmp_path / "night" / "atcc_its1"
    rows = [l.split("\t") for l in (d / "arms.tsv").read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["std", "tq3", "tq4", "tq5", "tq6", "tq8", "tq10", "sw10", "sw15", "sw20", "sw25", "eer"]
    eer = json.loads((d / "eer.json").read_text())
    assert rows[-1][5:7] == [str(eer["max_ee"][0]), str(eer["max_ee"][1])] and rows[3][7] == "5"
    assert all((d / f"reads_sw{q}").is_dir() for q in N.SW_STRENGTHS)
    status = json.loads((tmp_path / "night" / "night_status.json").read_text())
    assert sorted(k.split("/")[1] for k in status if k.startswith("atcc_its1/")) == ["arms", "eer", "thin", "trims"] and "mock3/arms" in status
    assert not any("-c" in c for c, _ in world.calls)                                # no AmpliTrim in the rules night
    assert "its_rules_report.py" in logged[-1]


# ------------------------------------------------------------------ the Bakker night (--plan bakker, register section 37)
def test_the_bakker_night_runs_eleven_datasets_six_arms_and_computes_maxee_on_the_trimmed_reads_for_sweer(tmp_path):
    world = World()
    base = tmp_path / "data"
    for lab in ("a1", "b3"):
        for smp in ("A", "B"):
            for m in (1, 2):
                fq(base / lab / f"{smp}_R{m}.fastq.gz", [("ACGT" * 27 + "AC", "I" * 60 + "5" * 10 + "$" * 40)] * 6)   # a Q3 tail: about 20 expected errors untrimmed, 0.1 trimmed
    logged = []
    rc = N.main(["--out", str(tmp_path / "night"), "--base", str(base), "--plan", "bakker", "--pairs", "4", "--dry-run"], run=world, thin=fake_thin, rscript="Rscript-fake",
                sixteen=[], log=logged.append)
    assert rc == 0 and "11 datasets. Nothing was run." in logged[-1] and any(l.startswith("bakker_pol_its2") for l in logged)
    import argparse
    a = argparse.Namespace(base=base, dataset=[f"atcc_its1={base / 'a1'}"], extra=[], no_16s=True, only="", pairs=4, plan="bakker")
    labels = [p[0] for p in N.plan(a, [])]
    assert labels[-3:] == ["bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2"] and len(labels) == 11
    night = tmp_path / "night2"
    rc = N.main(["--out", str(night), "--base", str(base), "--plan", "bakker", "--only", "atcc_its1", "--pairs", "4", "--dataset", f"atcc_its1={base / 'a1'}"],
                run=world, thin=fake_thin, rscript="Rscript-fake", sixteen=[], log=logged.append)
    assert rc == 0, "\n".join(logged)
    d = night / "atcc_its1"
    rows = [l.split("\t") for l in (d / "arms.tsv").read_text().splitlines()[1:]]
    assert [r[0] for r in rows] == ["std", "sw10", "sw15", "tq10", "eer", "sweer"]
    sweer, eer = rows[-1], rows[-2]
    assert sweer[2].endswith("reads_sw15") and eer[2].endswith("reads")
    assert eer[5:7] == ["21", "21"] and sweer[5:7] == ["2", "2"]                          # the computed maxEE follows the reads it is computed on
    assert sweer[5:7] == [str(x) for x in json.loads((d / "eer_sw15.json").read_text())["max_ee"]]
    assert eer[5:7] == [str(x) for x in json.loads((d / "eer.json").read_text())["max_ee"]]
    assert (d / "reads_sw10").is_dir() and (d / "reads_sw15").is_dir() and not (d / "reads_sw20").exists()
    status = json.loads((night / "night_status.json").read_text())
    assert sorted(k.split("/")[1] for k in status) == ["arms", "eer", "eer_sw15", "thin", "trims"]
    assert "its_bakker_report.py" in logged[-1]
