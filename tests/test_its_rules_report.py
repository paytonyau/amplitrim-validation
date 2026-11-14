"""validation/its_rules_report.py: P94 to P104, safety on every truth, the merging mechanism, the consistency gate and the advice, on a synthetic rules night
whose answers are known by construction."""
import importlib.util
import random
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("its_rules_report", HERE.parent / "validation" / "its_rules_report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = load()
random.seed(9)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
A_, B_, X_, K_, C_ = rnd(250), rnd(260), rnd(255), rnd(220), rnd(230)
JUNK, JUNK2 = rnd(250), rnd(240)


def edit(s, n):
    s = list(s)
    for i in range(n):
        s[10 + 7 * i] = "A" if s[10 + 7 * i] != "A" else "C"
    return "".join(s)


BASE = {A_: 50, edit(B_, 1): 30, edit(X_, 1): 10, JUNK: 10}                          # true exact, true near, unjudged, spurious
LOST_B = {A_: 50, edit(X_, 1): 10}                                                     # B lost, no spurious read
FLOWS = {  # input, filtered, denoised, merged, nonchim, asvs (per sample) on the ATCC standard
    "std": (1000, 800, 750, 700, 600, 4), "tq3": (1000, 850, 800, 790, 700, 4), "tq4": (1000, 900, 850, 840, 800, 3), "tq5": (1000, 880, 820, 810, 750, 4),
    "tq6": (1000, 820, 760, 712, 620, 4), "tq8": (1000, 900, 850, 840, 800, 4), "tq10": (1000, 900, 850, 840, 800, 4), "sw10": (1000, 830, 780, 700, 650, 4),
    "sw15": (1000, 900, 870, 860, 800, 4), "sw20": (1000, 920, 900, 890, 850, 4), "sw25": (1000, 950, 920, 910, 900, 3), "eer": (1000, 850, 800, 780, 720, 4)}
ATCC_SEQS = {arm: (LOST_B if arm in ("tq4", "sw25") else BASE) for arm in FLOWS}
STUDY_ASV = {"std": 10, "tq3": 11, "tq6": 20, "eer": 12}                              # median ASVs on the two studies; others 11


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows))


def arm(out, lab, arm_, seqs, flow):
    d = out / lab / "dada2" / arm_
    write(d / "samples.tsv", ["sample", "input", "filtered", "denoised", "merged", "nonchim", "asvs"], [(s,) + tuple(flow) for s in seqs])
    write(d / "seqs.tsv", ["sample", "sequence", "abundance"], [(s, q, a) for s, v in seqs.items() for q, a in v.items()])


def refs(folder):
    folder.mkdir(parents=True, exist_ok=True)
    for m in ("ITS1", "ITS2"):
        (folder / f"references_{m}.fasta").write_text(f">A|sp|G1|strain_genome|copy1x1|250nt\n{A_}\n>B|sp|G2|strain_record|copy1x1|260nt\n{B_}\n>X|sp|G3|species_genome|copy1x1|255nt\n{X_}\n")
    return folder


def build(tmp):
    out = tmp / "rules"
    for lab in ("atcc_its1", "atcc_its2"):
        for a, f in FLOWS.items():
            arm(out, lab, a, {f"ATCC-{lab[-4:]}a": dict(ATCC_SEQS[a]), f"ATCC-{lab[-4:]}b": dict(ATCC_SEQS[a])}, f)
    for lab in ("arfken_its1", "arfken_its2"):
        arm(out, lab, "std", {"K-ITSa": {K_: 100}, "C-ITSa": {C_: 100}}, (1000, 900, 850, 840, 800, 1))
    for lab in ("arfken_mix_its1", "arfken_mix_its2"):
        for a in FLOWS:
            seqs = {K_: 60, edit(C_, 2): 30, JUNK: 10}
            if a == "tq8":
                seqs[JUNK2] = 5                                                        # safe on the ATCC standard, not in the mixes
            arm(out, lab, a, {"MixB-ITSa": seqs}, (1000, 900, 850, 840, 800, 3))
    for lab in ("its2", "its"):
        for a in FLOWS:
            arm(out, lab, a, {"s1": {rnd(300): 1}}, (1000, 800, 700, 600, 500, STUDY_ASV.get(a, 11)))
    prev = tmp / "run26"
    for lab_dir in out.iterdir():
        for a in ("std", "sw15"):
            src = lab_dir / "dada2" / a / "samples.tsv"
            if src.is_file():
                (prev / lab_dir.name / "dada2" / a).mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, prev / lab_dir.name / "dada2" / a / "samples.tsv")
    return out, prev


def run(tmp, extra=()):
    out, prev = build(tmp)
    rc = R.main(["--out", str(out), "--references", str(refs(tmp / "refs")), "--previous", str(prev)] + list(extra))
    return out, prev, rc


def verdict(out, cid):
    for l in (out / "ITS_RULES_SUMMARY.md").read_text().splitlines():
        if l.startswith(f"- **{cid}** "):
            return l
    raise AssertionError(cid)


def test_the_verdicts_and_the_advice(tmp_path):
    out, _, rc = run(tmp_path)
    assert rc == 0
    met = lambda c: verdict(out, c).startswith(f"- **{c}** MET")
    assert met("P94 tq3") and met("P95 tq3") and met("P96 tq3") and met("P97 tq3")
    assert not met("P94 tq4") and "recall -50.0" in verdict(out, "P94 tq4")
    assert met("P94 tq5") and met("P95 tq5")
    assert met("P94 tq6") and not met("P95 tq6") and not met("P96 tq6") and not met("P97 tq6")        # safe, not useful, merging not better, ASVs doubled
    assert not met("P94 tq8") and "arfken_mix_its1 spurious +4.29" in verdict(out, "P94 tq8")          # caught only by the second truth
    assert met("P98 sw10") and met("P99 sw10") and "+5.0" in verdict(out, "P99 sw10")                  # exactly 5 points: met
    assert met("P98 sw20") and not met("P98 sw25") and met("P101")
    assert met("P102") and met("P103") and met("P104")
    s = (out / "ITS_RULES_SUMMARY.md").read_text()
    assert "`truncQ` option: `tq5`" in s and "The trim strength: `sw20`" in s and "Advised for poor read ends: `sw20`" in s and "may be used for ITS" in s
    assert "reproduces run 26 on 8 datasets exactly" in s


def test_a_changed_count_in_the_default_arm_stops_every_verdict(tmp_path):
    out, prev = build(tmp_path)
    p = prev / "its" / "dada2" / "std" / "samples.tsv"
    p.write_text(p.read_text().replace("\t1000\t800\t700", "\t1000\t801\t700"))
    rc = R.main(["--out", str(out), "--references", str(refs(tmp_path / "refs")), "--previous", str(prev)])
    s = (out / "ITS_RULES_SUMMARY.md").read_text()
    assert rc == 2 and "FAILED, so no criterion is decided" in s and "its: 1 sample(s) differ" in s
    assert "- **P94 tq3** not decided (the consistency check failed)" in s and "truncQ stays 2" in s and "ITS keeps maxEE 2,2" in s


def test_sw15_must_reproduce_run_26(tmp_path):
    out, prev = build(tmp_path)
    p = prev / "atcc_its2" / "dada2" / "sw15" / "samples.tsv"
    p.write_text(p.read_text().replace("\t870\t", "\t871\t"))
    R.main(["--out", str(out), "--references", str(refs(tmp_path / "refs")), "--previous", str(prev)])
    assert verdict(out, "P101").startswith("- **P101** NOT MET") and "atcc_its2: 2 sample(s) differ" in verdict(out, "P101")


def test_a_further_mock_given_with_truth_must_be_safe_too(tmp_path):
    out, prev = build(tmp_path)
    for a, f in FLOWS.items():                                                         # mock3: the computed maxEE loses member B there
        arm(out, "mock3", a, {"M3a": dict(LOST_B if a == "eer" else BASE)}, f)
    rc = R.main(["--out", str(out), "--references", str(refs(tmp_path / "refs")), "--previous", str(prev), "--truth", f"mock3={refs(tmp_path / 'refs3')}"])
    assert rc == 0 and verdict(out, "P102").startswith("- **P102** NOT MET") and "mock3 spurious" in verdict(out, "P102")
    assert "ITS keeps maxEE 2,2" in (out / "ITS_RULES_SUMMARY.md").read_text()


def test_the_report_holds_the_merge_share_and_the_computed_maxee(tmp_path):
    out, prev = build(tmp_path)
    import json
    (out / "atcc_its2" / "eer.json").write_text(json.dumps({"percentile": [2.4, 6.7], "raw": [3, 7], "max_ee": [3, 7], "reads": [9000, 8500]}))
    R.main(["--out", str(out), "--references", str(refs(tmp_path / "refs")), "--previous", str(prev)])
    rep = (out / "ITS_RULES_REPORT.md").read_text()
    assert "| atcc_its2 | 2.40 / 6.70 | 3,7 | 9,000 / 8,500 |" in rep and "| atcc_its2 | tq3 | 2,000 | 85.00% | 79.00% | 70.00% | 98.75% | 4 |" in rep
    assert (out / "its_rules_metrics.tsv").read_text().splitlines()[0].startswith("dataset\tarm\tinput\tpass\tmerged\tyield\tmerge_share")
