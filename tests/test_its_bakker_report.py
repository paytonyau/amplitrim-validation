"""validation/its_bakker_report.py: P105 to P114 on a synthetic Bakker night whose answers are known by construction, and one variant per decision."""
import importlib.util
import random
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("its_bakker_report", HERE.parent / "validation" / "its_bakker_report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B = load()
random.seed(21)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
A_, X_, K_, C_, M1, M2, M3 = (rnd(250) for _ in range(7))
JUNK = rnd(250)
ARMS = ["std", "sw10", "sw15", "tq10", "eer", "sweer"]
ATCC_SEQ = {A_: 80, JUNK: 10, X_: 10}                       # A true (strain), X unjudged (species), JUNK spurious
BAKKER_SEQ = {M1: 50, M2: 30, M3: 20}                       # all three members of the (synthetic) Bakker truth
SAMPLES = {"bakker_pcr_its1": ["FMockE1", "FMockE_HC1"], "bakker_pol_its1": ["MockE_GoTaq_1", "MockE_Phusion_1"], "bakker_pol_its2": ["MockE_GoTaq_1", "MockE_PreAmp_1"]}


def flows(lab, arm):
    """input, filtered, denoised, merged, nonchim, asvs per sample; the yield is nonchim / 1000."""
    low = lab == "atcc_its2"
    base = {"std": (270 if low else 740, 250 if low else 700)}
    y = {"atcc_its2": {"std": 200, "sw10": 300, "sw15": 400, "tq10": 450, "eer": 210, "sweer": 430},
         "atcc_its1": {"std": 600, "sw10": 620, "sw15": 650, "tq10": 700, "eer": 600, "sweer": 650},
         "bakker_pcr_its1": {"std": 800, "sw10": 800, "sw15": 790, "tq10": 820, "eer": 800, "sweer": 790},
         "bakker_pol_its1": {"std": 800, "sw10": 805, "sw15": 810, "tq10": 830, "eer": 795, "sweer": 810},
         "bakker_pol_its2": {"std": 780, "sw10": 790, "sw15": 800, "tq10": 820, "eer": 780, "sweer": 795}}.get(lab, {}).get(arm, 600)
    filt = 270 if (low and arm in ("std", "eer")) else 900 if lab != "atcc_its2" else 850
    if lab in B.BAKKER or lab.startswith("arfken") or lab == "atcc_its1":
        filt = 900
    return (1000, filt, y + 20, y + 10, y, {"std": 10, "sw15": 11, "eer": 10, "sweer": 12}.get(arm, 10))


def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows))


def arm(out, lab, a, seqs):
    d = out / lab / "dada2" / a
    write(d / "samples.tsv", ["sample", "input", "filtered", "denoised", "merged", "nonchim", "asvs"], [(s,) + flows(lab, a) for s in seqs])
    write(d / "seqs.tsv", ["sample", "sequence", "abundance"], [(s, q, n) for s, v in seqs.items() for q, n in v.items()])


def refs(folder, pairs):
    folder.mkdir(parents=True, exist_ok=True)
    for m in ("ITS1", "ITS2"):
        (folder / f"references_{m}.fasta").write_text("".join(f">{n}|sp|G|{k}|copy1x1|250nt\n{s}\n" for n, k, s in pairs))
    return folder


def build(tmp, change=None):
    """change: {(lab, arm): seqs} to override."""
    out = tmp / "night"
    change = change or {}
    for lab in ("atcc_its1", "atcc_its2"):
        for a in ARMS:
            arm(out, lab, a, {f"ATCC-{lab[-4:]}{r}": dict(change.get((lab, a), ATCC_SEQ)) for r in "ab"})
    for lab in ("arfken_its1", "arfken_its2"):
        arm(out, lab, "std", {"K-ITSa": {K_: 100}, "C-ITSa": {C_: 100}})
    for lab in B.MIXES:
        for a in ARMS:
            arm(out, lab, a, {"MixB-ITSa": dict(change.get((lab, a), {K_: 60, C_: 30, JUNK: 10}))})
    for lab in B.BAKKER:
        for a in ARMS:
            seqs = change.get((lab, a), {M1: 50, M3: 20, JUNK: 30} if a == "tq10" and lab == "bakker_pcr_its1" else BAKKER_SEQ)
            arm(out, lab, a, {s: dict(seqs) for s in SAMPLES[lab]})
    for lab in B.STUDIES:
        for a in ARMS:
            arm(out, lab, a, {"s1": {rnd(300): 1}})
    prev26, prev27 = tmp / "run26", tmp / "run27"
    for lab_dir in out.iterdir():
        for a, prev in (("std", prev26), ("sw15", prev26), ("eer", prev27)):
            src = lab_dir / "dada2" / a / "samples.tsv"
            if src.is_file() and lab_dir.name in B.OLD:
                (prev / lab_dir.name / "dada2" / a).mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, prev / lab_dir.name / "dada2" / a / "samples.tsv")
    r = refs(tmp / "refs", [("A", "strain_genome", A_), ("X", "species_genome", X_)])
    br = tmp / "brefs"
    br.mkdir(parents=True, exist_ok=True)
    for lab in B.BAKKER:                                                                # one reference file per Bakker dataset
        (br / f"references_{lab}.fasta").write_text("".join(f">{n}|sp|G|strain_record|copy1x1|250nt\n{q}\n" for n, q in (("M1", M1), ("M2", M2), ("M3", M3))))
    return out, r, br, prev26, prev27


def run(tmp, change=None, tweak=None):
    out, r, br, p26, p27 = build(tmp, change)
    if tweak:
        tweak(p26, p27)
    rc = B.main(["--out", str(out), "--references", str(r), "--bakker-references", str(br), "--previous", str(p26), "--previous-rules", str(p27)])
    return out, rc, (out / "ITS_BAKKER_SUMMARY.md").read_text()


def v(summary, cid):
    return next(l for l in summary.splitlines() if l.startswith(f"- **{cid}** "))


met = lambda s, c: v(s, c).startswith(f"- **{c}** MET")


def test_everything_holds_in_the_base_night(tmp_path):
    out, rc, s = run(tmp_path)
    assert rc == 0, s
    for c in ("P105", "P106", "P107", "P108", "P109", "P110", "P111", "P112 sw15", "P112 eer", "P112 sweer", "P113", "P114"):
        assert met(s, c), v(s, c)
    assert "atcc_its2 +20.0" in v(s, "P106") and "bakker_pcr_its1 safe, yield -1.0" in v(s, "P107")
    assert "holds on a third mock community, with other read lengths" in s and "0.3 may combine the trim" in s and "the 50% trigger is not needed for safety" in s


def test_a_member_lost_only_in_bakkers_its2_makes_the_trim_unsafe(tmp_path):
    out, rc, s = run(tmp_path, {("bakker_pol_its2", "sw15"): {M1: 70, M3: 30}})
    assert not met(s, "P105") and "bakker_pol_its2 spurious +0.00, recall -33.3" in v(s, "P105") and "restricted to what it was safe on" in s


def test_the_trigger_stays_when_the_trim_costs_yield_on_well_passing_data(tmp_path):
    import its_bakker_report as BB
    out, r, br, p26, p27 = build(tmp_path)
    p = out / "bakker_pcr_its1" / "dada2" / "sw15" / "samples.tsv"
    p.write_text(p.read_text().replace("\t790\t", "\t770\t"))                          # -3 points
    rc = B.main(["--out", str(out), "--references", str(r), "--bakker-references", str(br), "--previous", str(p26), "--previous-rules", str(p27)])
    s = (out / "ITS_BAKKER_SUMMARY.md").read_text()
    assert not met(s, "P107") and "bakker_pcr_its1 safe, yield -3.0" in v(s, "P107") and "the 50% trigger stays" in s


def test_a_combination_no_better_than_the_trim_stays_separate(tmp_path):
    out, r, br, p26, p27 = build(tmp_path)
    p = out / "atcc_its2" / "dada2" / "sweer" / "samples.tsv"
    p.write_text(p.read_text().replace("\t430\t", "\t405\t"))                          # +0.5 over the trim, not +2
    B.main(["--out", str(out), "--references", str(r), "--bakker-references", str(br), "--previous", str(p26), "--previous-rules", str(p27)])
    s = (out / "ITS_BAKKER_SUMMARY.md").read_text()
    assert met(s, "P110") and not met(s, "P111") and "stay separate" in s


def test_reproduction_and_the_gate(tmp_path):
    def eer_differs(p26, p27):
        q = p27 / "atcc_its1" / "dada2" / "eer" / "samples.tsv"
        q.write_text(q.read_text().replace("\t600\t", "\t601\t"))
    out, rc, s = run(tmp_path / "a", tweak=eer_differs)
    assert not met(s, "P113") and "atcc_its1: 2 sample(s) differ" in v(s, "P113")

    def std_differs(p26, p27):
        q = p26 / "its" / "dada2" / "std" / "samples.tsv"
        q.write_text(q.read_text().replace("\t600\t", "\t599\t"))
    out, rc, s = run(tmp_path / "b", tweak=std_differs)
    assert rc == 2 and "FAILED, so no criterion is decided" in s and "- **P105** not decided (the consistency check failed)" in s


def test_bakkers_samples_are_grouped_by_condition(tmp_path):
    assert [B.condition(x) for x in ("FMockE1", "FMockA_HC2", "FMockB_LE3", "MockE_GoTaq_1", "MockA_Phusion_2", "MockB_PreAmp_3")] == \
        ["standard", "high cycles", "long extension", "GoTaq", "Phusion", "PreAmp"]
    out, rc, s = run(tmp_path)
    rep = (out / "ITS_BAKKER_REPORT.md").read_text()
    assert "| bakker_pcr_its1 | tq10 | high cycles | 30.00% | 66.67% | 2.0 | 1 |" in rep and "| bakker_pol_its2 | sw15 | PreAmp | 0.00% | 100.00% | 3.0 | 1 |" in rep
    assert "| bakker_pcr_its1 | sw15 | 2,000 | 90.00% | 80.00% | 79.00% | 98.77% | 11 |" in rep                # two samples; merged 800 of 810 denoised


def test_the_overhang_gate_passes_clean_reads_and_stops_every_verdict_when_bakkers_reads_carry_overhangs(tmp_path):
    import bakker_prep as BP
    import bakker_primers as bpr
    out, rc_, s = run(tmp_path / "a")
    assert "**The overhang check:**" not in s                                           # off unless asked: run 28's analysis is unchanged
    out, r, br, p26, p27 = build(tmp_path / "b")
    args = ["--out", str(out), "--references", str(r), "--bakker-references", str(br), "--previous", str(p26), "--previous-rules", str(p27), "--overhang-max", "0.01"]
    assert B.main(args) == 0 and "(at most 1%: passed)" in (out / "ITS_BAKKER_SUMMARY.md").read_text()
    tail = "".join(bpr.IUPAC[c][0] for c in BP.rc(bpr.PRIMERS["B58S3"]))[:8]
    arm(out, "bakker_pcr_its1", "std", {smp: {M1 + tail: 50, M2: 30, M3: 20} for smp in SAMPLES["bakker_pcr_its1"]})
    assert B.main(args) == 2
    summ = (out / "ITS_BAKKER_SUMMARY.md").read_text()
    assert "bakker_pcr_its1 50.00%" in summ and "FAILED (more than 1%" in summ and "- **P105** not decided" in summ


def test_a_failed_overhang_gate_names_itself_and_three_bases_pass_a_one_base_overhang(tmp_path):
    out, r, br, p26, p27 = build(tmp_path)
    arm(out, "bakker_pcr_its1", "std", {smp: {"A" + M1: 50, M2: 30, M3: 20} for smp in SAMPLES["bakker_pcr_its1"]})   # one extra base, within the 2-edit rule
    base = ["--out", str(out), "--references", str(r), "--bakker-references", str(br), "--previous", str(p26), "--previous-rules", str(p27), "--overhang-max", "0.01"]
    assert B.main(base) == 2
    s = (out / "ITS_BAKKER_SUMMARY.md").read_text()
    assert "- **P105** not decided (the overhang check failed)" in s and "consistency check failed" not in s and "exactly." in s
    assert B.main(base + ["--overhang-min-bases", "3", "--suffix", "_min3"]) == 0
    s3 = (out / "ITS_BAKKER_SUMMARY_min3.md").read_text()
    assert "of at least 3 bases: bakker_pcr_its1 0.00%" in s3 and "- **P105** MET" in s3
    assert (out / "ITS_BAKKER_SUMMARY.md").read_text() == s and (out / "its_bakker_metrics_min3.tsv").is_file()     # the first analysis is kept
