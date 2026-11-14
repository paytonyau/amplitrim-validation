"""validation/its_night_report.py: the truths, the metrics and the verdicts P83 to P93, on a synthetic night whose answers are known by construction."""
import importlib.util
import json
import math
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load():
    spec = importlib.util.spec_from_file_location("its_night_report", HERE.parent / "validation" / "its_night_report.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = load()
random.seed(3)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
A1, B1, X1, A2, B2, X2 = rnd(200), rnd(210), rnd(205), rnd(300), rnd(310), rnd(305)        # A, B: strain-exact members; X: a species-level member
K1, C1, K2, C2 = rnd(220), rnd(230), rnd(320), rnd(330)                                     # two piglet isolates


def edit(s, n):
    s = list(s)
    for i in range(n):
        s[10 + 7 * i] = "A" if s[10 + 7 * i] != "A" else "C"
    return "".join(s)


def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def test_the_banded_distance_equals_the_full_one_up_to_k():
    for _ in range(300):
        a = rnd(random.randint(5, 40))
        b = list(a)
        for _ in range(random.randint(0, 5)):
            op, i = random.choice("sid"), random.randrange(len(b) + 1)
            if op == "s" and i < len(b):
                b[i] = random.choice("ACGT")
            elif op == "i":
                b.insert(i, random.choice("ACGT"))
            elif op == "d" and i < len(b):
                del b[i]
        b = "".join(b)
        for k in (0, 1, 2, 3):
            assert R.within(a, b, k) == min(lev(a, b), k + 1), (a, b, k)


def test_spearman_handles_ties():
    assert R.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == 1.0 and R.spearman([1, 2, 3], [3, 2, 1]) == -1.0
    assert abs(R.spearman([1, 1, 2], [1, 2, 3]) - 0.866) < 0.001


# ------------------------------------------------------------------ the synthetic night
def write(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows))


def arm(out, lab, arm_, samples, seqs, flow=(1000, 800, 700, 600, 3)):
    """samples: names; seqs: {sample: {sequence: reads}}; flow: input, filtered, merged, nonchim, asvs (per sample)."""
    d = out / lab / "dada2" / arm_
    write(d / "samples.tsv", ["sample", "input", "filtered", "denoised", "merged", "nonchim", "asvs"],
          [(s, flow[0], flow[1], flow[2], flow[2], flow[3], flow[4]) for s in samples])
    write(d / "seqs.tsv", ["sample", "sequence", "abundance"], [(s, q, a) for s, v in seqs.items() for q, a in v.items()])
    (d / "DONE").write_text("")


def refs(folder):
    folder.mkdir(parents=True, exist_ok=True)
    for m, (a, b, x) in (("ITS1", (A1, B1, X1)), ("ITS2", (A2, B2, X2))):
        (folder / f"references_{m}.fasta").write_text(f">A|sp_a|GCF_1|strain_genome|copy1x1|{len(a)}nt\n{a}\n>B|sp_b|GCA_2|strain_record|copy1x1|{len(b)}nt\n{b}\n"
                                                      f">X|sp_x|GCA_3|species_genome|copy1x1|{len(x)}nt\n{x}\n")
    return folder


def build(tmp):
    out = tmp / "night"
    junk1, junk2 = rnd(200), rnd(300)
    for lab, (a, b, x, junk) in (("atcc_its1", (A1, B1, X1, junk1)), ("atcc_its2", (A2, B2, X2, junk2))):
        smp = [f"ATCC-ITS{lab[-1]}a", f"ATCC-ITS{lab[-1]}b"]
        base = {a: 50, edit(b, 1): 30, edit(x, 1): 10, junk: 10}                       # true exact, true with 1 edit, unjudged, spurious
        arm(out, lab, "std", smp, {s: dict(base) for s in smp}, flow=(1000, 800, 700, 600, 4))
        arm(out, lab, "ee23", smp, {s: dict(base) for s in smp}, flow=(1000, 900, 800, 700, 4))      # as safe, yield +10 points
        worse = dict(base)
        worse[rnd(len(junk))] = 4                                                        # spurious 14/94 against 10/90: +3.8 points
        arm(out, lab, "ee25", smp, {s: dict(worse) for s in smp}, flow=(1000, 950, 850, 750, 5))
        lost = {a: 50, junk: 10}                                                         # B lost: recall 50%
        arm(out, lab, "tq20", smp, {s: dict(lost) for s in smp}, flow=(1000, 900, 800, 700, 2))
        clean_loss = {a: 50, edit(x, 1): 10}                                             # B lost, no spurious read at all: only recall can fail
        arm(out, lab, "tq10", smp, {s: dict(clean_loss) for s in smp}, flow=(1000, 900, 800, 700, 2))
        arm(out, lab, "fwd", smp, {s: {a[:150]: 60, b[:140]: 30, junk[:150]: 10} for s in smp})
        arm(out, lab, "mo20", smp, {s: dict(base) for s in smp})
        arm(out, lab, "mo30", smp, {s: {a: 50, edit(b, 1): 30, edit(x, 1): 10} for s in smp})      # no spurious read at 30 nt: 12 nt is 11.1 points worse
        arm(out, lab, "mfirst", smp, {s: dict(base) for s in smp})
    for lab, (k, c) in (("arfken_its1", (K1, C1)), ("arfken_its2", (K2, C2))):
        arm(out, lab, "std", ["K-ITSa", "K-ITSb", "C-ITSa"], {"K-ITSa": {k: 99, rnd(50): 0.5}, "K-ITSb": {k: 100}, "C-ITSa": {c: 100}})
    for lab, (k, c) in (("arfken_mix_its1", (K1, C1)), ("arfken_mix_its2", (K2, C2))):
        arm(out, lab, "std", ["MixB-ITSa", "MixC-ITSa"], {"MixB-ITSa": {k: 60, edit(c, 2): 30, rnd(200): 10}, "MixC-ITSa": {k: 100}})
    for lab, std_flow, mf in (("its2", (1000, 800, 140, 120, 10), 400), ("its", (1000, 700, 470, 450, 10), 500)):     # ASVs: 10 at std, 11 at ee23 (+10%), 30 at ee25
        arm(out, lab, "std", ["s1", "s2"], {"s1": {rnd(300): 1}}, flow=std_flow)
        arm(out, lab, "ee23", ["s1", "s2"], {"s1": {rnd(300): 1}}, flow=(1000, 850, 160, 140, 11))
        arm(out, lab, "ee25", ["s1", "s2"], {"s1": {rnd(300): 1}}, flow=(1000, 900, 200, 180, 30))
        arm(out, lab, "mfirst", ["s1", "s2"], {"s1": {rnd(300): 1}}, flow=(1000, 800, mf, mf, 12))
    # Q1: per sample, d10 / d5 / d2 / d0 merged pairs; its2 DADA2 merged 14%: the 0-difference estimate is closest (16%)
    write(out / "its2" / "estimates.tsv", ["sample", "filtered_pairs", "merged_d10", "merged_d5", "merged_d2", "merged_d0"], [("s1", 800, 790, 600, 400, 160), ("s2", 800, 790, 600, 400, 160)])
    write(out / "its" / "estimates.tsv", ["sample", "filtered_pairs", "merged_d10", "merged_d5", "merged_d2", "merged_d0"], [("s1", 700, 690, 600, 500, 400), ("s2", 700, 690, 600, 500, 400)])
    for lab, d10, d0 in (("16s_a", 950, 700), ("16s_b", 900, 600)):
        arm(out, lab, "std", ["s1"], {"s1": {rnd(250): 1}}, flow=(1000, 960, 800, 790, 20))
        write(out / lab / "estimates.tsv", ["sample", "filtered_pairs", "merged_d10", "merged_d5", "merged_d2", "merged_d0"], [("s1", 960, d10, 900, 850, d0)])
    # AmpliTrim on two ITS groups: predicted pass and yield at the maxEE settings
    for lab, shift in (("atcc_its1", 0.002), ("its2", 0.03)):
        table = [{"maxEE": [2, 2], "pass_rate": 0.8 + shift, "estimated_yield": 0.75}, {"maxEE": [2, 3], "pass_rate": 0.9 + shift, "estimated_yield": 0.85},
                 {"maxEE": [2, 5], "pass_rate": 0.95 + shift, "estimated_yield": 0.9}]
        (out / lab / "amplitrim").mkdir(parents=True, exist_ok=True)
        (out / lab / "amplitrim" / "summary.json").write_text(json.dumps({"maxEE_table": table}))
    return out


def run(tmp):
    out = build(tmp)
    rc = R.main(["--out", str(out), "--references", str(refs(tmp / "refs"))])
    return out, rc


def verdict(out, cid):
    for line in (out / "ITS_NIGHT_SUMMARY.md").read_text().splitlines():
        if line.startswith(f"- **{cid}** "):
            return line
    raise AssertionError(f"{cid} not in the summary")


def test_the_atcc_truth_counts_exact_near_unjudged_and_spurious_reads(tmp_path):
    out = build(tmp_path)
    res = R.analyse(out, {m: R.read_references(refs(tmp_path / "r") / f"references_{m}.fasta") for m in ("ITS1", "ITS2")})
    a = res["acc"][("atcc_its1", "std")]
    assert math.isclose(a["spurious"], 10 / 90) and math.isclose(a["exact"], 50 / 90) and math.isclose(a["unjudged"], 0.10)
    assert a["recall"] == 1.0 and a["found"] == 2 and a["samples"] == 2
    assert res["acc"][("atcc_its1", "tq20")]["recall"] == 0.5
    f = res["acc"][("atcc_its1", "fwd")]
    assert math.isclose(f["spurious"], 0.10) and f["recall"] == 1.0                     # forward reads judged against the starts of the references


def test_the_isolate_reference_judges_the_mixes(tmp_path):
    out = build(tmp_path)
    res = R.analyse(out, {m: R.read_references(refs(tmp_path / "r") / f"references_{m}.fasta") for m in ("ITS1", "ITS2")})
    assert res["isolate_reference"]["arfken_mix_its1"] == {"C": 1, "K": 1}                # the rare sequence of K-ITSa (0.5%) is not in the reference
    m = res["acc"][("arfken_mix_its1", "std")]
    assert math.isclose(m["spurious"], (0.10 + 0.0) / 2) and m["recall"] == 1.0             # MixB: K and C found (C within 2 edits); MixC: K


def test_the_verdicts(tmp_path):
    out, rc = run(tmp_path)
    assert rc == 2                                                                        # some arms were not run in this small night, so some criteria are undecided
    assert verdict(out, "P83").startswith("- **P83** MET") and "(closest: 0)" in verdict(out, "P83")       # 16% against 14%: +2 points
    assert "- **P84** NOT MET" in verdict(out, "P84")                                    # on 16S the 0-difference estimate is far farther from DADA2
    assert verdict(out, "P85 ee23").startswith("- **P85 ee23** MET") and verdict(out, "P85 ee25").startswith("- **P85 ee25** NOT MET")
    assert verdict(out, "P85 tq20").startswith("- **P85 tq20** NOT MET")                 # recall falls by 50 points
    assert verdict(out, "P86 ee23").startswith("- **P86 ee23** MET") and "+10.0 points" in verdict(out, "P86 ee23")
    assert "not decided" in verdict(out, "P85 ee33")
    assert verdict(out, "P87 ee25").startswith("- **P87 ee25** NOT MET") and verdict(out, "P87 ee23").startswith("- **P87 ee23** MET")
    assert verdict(out, "P89").startswith("- **P89** MET") and verdict(out, "P90").startswith("- **P90** MET")     # its2: 40% against 12%
    assert verdict(out, "P88").startswith("- **P88** MET")                                # forward reads: spurious 10% against 11.1%, recall the same
    assert verdict(out, "P93").startswith("- **P93** NOT MET")                            # 12 nt lets 11.1 points of spurious reads through that 30 nt does not
    assert "P92" in (out / "ITS_NIGHT_SUMMARY.md").read_text() and "**P92** NOT MET" in verdict(out, "P92")             # its2's pass is 3 points off
    assert "`ee23`" in (out / "ITS_NIGHT_SUMMARY.md").read_text()                          # the advice: the only arm safe on both markers and useful
    rows = (out / "its_night_metrics.tsv").read_text().splitlines()
    assert rows[0].startswith("dataset\tarm\tinput") and any(r.startswith("atcc_its1\tstd\t") for r in rows)



def test_losing_a_member_is_unsafe_even_without_spurious_reads(tmp_path):
    out, _ = run(tmp_path)
    v = verdict(out, "P85 tq10")
    assert v.startswith("- **P85 tq10** NOT MET") and "recall -50.0 points" in v and "spurious -11.11 points" in v


def test_a_safe_arm_that_is_not_useful_is_not_the_advice(tmp_path):
    out = build(tmp_path)
    p = out / "atcc_its2" / "dada2" / "ee23" / "samples.tsv"
    p.write_text(p.read_text().replace("\t700\t4\n", "\t620\t4\n"))                      # ITS2 yield +2 points instead of +10
    assert R.main(["--out", str(out), "--references", str(refs(tmp_path / "refs"))]) == 2
    assert verdict(out, "P85 ee23").startswith("- **P85 ee23** MET") and verdict(out, "P86 ee23").startswith("- **P86 ee23** NOT MET")
    assert "keep maxEE 2,2" in (out / "ITS_NIGHT_SUMMARY.md").read_text()


def test_merging_first_is_measured_against_the_input_pairs_not_the_merged_reads(tmp_path):
    """mfirst's samples.tsv counts the reads that vsearch merged as its input; its yield is a share of the pairs, like every arm's."""
    out = build(tmp_path)
    p = out / "its2" / "dada2" / "mfirst" / "samples.tsv"
    p.write_text(p.read_text().replace("\t1000\t800\t400\t400\t400\t12", "\t500\t480\t400\t400\t400\t12"))          # 500 merged reads of 1000 pairs
    res = R.analyse(out, {m: R.read_references(refs(tmp_path / "r") / f"references_{m}.fasta") for m in ("ITS1", "ITS2")})
    f = res["flow"][("its2", "mfirst")]
    assert f["input"] == 2000 and math.isclose(f["yield"], 0.40) and math.isclose(f["pass"], 0.48)
    assert "+28.0 points" in dict((c, t) for c, _, t in R.verdicts(res))["P90"]                  # 40% against the default's 12%
