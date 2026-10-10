"""Round 4 (validation/round4_steps.py, round4.py, simulated/community_sim.py): the analyses, the simulator, the steps with their external tools replaced by fakes,
the pauses between phases and the phased runner. Nothing here needs DADA2 or vsearch."""
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
for p in (str(HERE), str(ROOT), str(ROOT / "validation"), str(ROOT / "validation" / "simulated")):
    if p not in sys.path:
        sys.path.insert(0, p)


def load(name: str, sub: str = "validation"):
    spec = importlib.util.spec_from_file_location(name, ROOT / sub / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def steps():
    import round4_steps
    return round4_steps


def helpers(ov):
    """The helpers overnight.py gives to the steps, with the tools that need R or vsearch replaced where a test says so."""
    return SimpleNamespace(verdict=ov.verdict, parse_extra=ov.parse_extra, table_lines=ov.table_lines, recommended=None, high_input_reference=lambda ctx: None, guarded=lambda n, w, fn, lines, verdicts: fn())


def write_fastq_gz(path: Path, reads: list):
    with gzip.open(path, "wt") as fh:
        for i, (seq, q) in enumerate(reads):
            fh.write(f"@r{i}\n{seq}\n+\n{q}\n")


# ------------------------------------------------------------------------------------------------ the simulator
def test_the_simulated_community_has_the_truth_it_claims(tmp_path):
    import community_sim as cs
    info = cs.write_community(tmp_path / "a", 60, 3, n_samples=2, pairs_per_sample=800, repeat_len=10, chimera_frac=0.03, indel_rate=1e-3)
    info2 = cs.write_community(tmp_path / "b", 60, 3, n_samples=2, pairs_per_sample=800, repeat_len=10, chimera_frac=0.03, indel_rate=1e-3)
    assert Path(info["truth_fasta"]).read_text() == Path(info2["truth_fasta"]).read_text() and info["counts"] == info2["counts"]       # the same seed, the same community
    seqs = [x.split("\n")[1] for x in Path(info["truth_fasta"]).read_text().split(">")[1:]]
    assert len(set(seqs)) == 60 and max(map(len, seqs)) <= 478 and min(map(len, seqs)) >= 440
    assert sum(1 for k in info["kinds"].values() if k == "repeat") == 12                                                              # a fifth of the taxa
    rep = next(x for x in Path(info["truth_fasta"]).read_text().split(">")[1:] if x.split("\n")[0].endswith("|repeat")).split("\n")[1]
    assert "ATATATATAT" in rep
    from amplitrim.subsample import find_pairs, read_fastq
    name, r1, r2 = find_pairs(tmp_path / "a")[0]
    kinds, lens = {}, set()
    for h, s, _p, q in read_fastq(r1):
        parts = h.split("|")
        kinds[parts[2]] = kinds.get(parts[2], 0) + 1
        lens.add(len(s))
        assert len(s) == len(q) and parts[1].startswith("T")
    assert kinds.keys() == {"normal", "repeat", "chimera"} and {249, 250, 251} <= lens                                                # indels make reads of 249 and 251 nt
    marked = [h for h, s, _p, q in read_fastq(r1) if h.endswith("|indel")]
    assert marked and all(len(h.split("|")) == 4 for h in marked)                                                                     # a pair with an indel says so in its header
    assert sum(info["counts"].values()) > 0 and set(info["counts"]) == set(info["lengths"])


def test_taxa_are_related_like_real_ones_and_the_community_is_uneven(tmp_path):
    import community_sim as cs
    rng = np.random.default_rng(5)
    taxa = cs.make_taxa(40, rng)
    same = [t for t in taxa if t["family"] == 0]
    a, b = same[0]["codes"], same[1]["codes"]
    assert 1 <= int((a[:min(len(a), len(b))] != b[:min(len(a), len(b))]).sum()) < 30                                                  # strains of a family differ by a few bases at most (before any indel)
    other = [t for t in taxa if t["family"] == 1][0]["codes"]
    assert len(other) > 400
    ab = cs.abundances(500, rng)
    assert abs(ab.sum() - 1) < 1e-9 and (ab < 0.001).mean() > 0.2 and cs.abundances(10, rng, even=True).std() < 1e-12
    prof = cs.sample_profiles(ab, 5, rng)
    assert len(prof) == 5 and all(abs(p.sum() - 1) < 1e-9 for p in prof) and (prof[0] == 0).any()
    assert [cs.abundance_class(n) for n in (0, 1, 2, 9, 10, 49, 50)] == ["unseen", "unseen", "rare", "rare", "medium", "medium", "abundant"]


def test_a_simulated_folder_is_never_overwritten_and_real_qualities_can_be_resampled(tmp_path):
    import community_sim as cs
    (tmp_path / "x").mkdir()
    (tmp_path / "x" / "something").write_text("x")
    with pytest.raises(SystemExit):
        cs.write_community(tmp_path / "x", 5, 1, 1, 50)
    real = tmp_path / "real"
    real.mkdir()
    write_fastq_gz(real / "S1_R1.fastq.gz", [("A" * 120, "I" * 60 + "5" * 60)] * 30)
    write_fastq_gz(real / "S1_R2.fastq.gz", [("A" * 120, "I" * 30 + "#" * 90)] * 30)
    mf, mr = cs.load_quality_pool(real, 250)
    assert mf.shape == (30, 250) and int(mf[0, 0]) == 40 and int(mf[0, 100]) == 20 and int(mf[0, 249]) == 20                           # padded with the last quality
    assert int(mr[0, 100]) == 2


# ------------------------------------------------------------------------------------------------ Q28: the split
def test_the_landscape_grid_is_read_back_exactly_as_the_landscape_step_writes_it(tmp_path):
    ov = load("overnight")
    rs = steps()
    ctx = SimpleNamespace(quick=False, cur=tmp_path, out=tmp_path, log=lambda *a: None)
    ov.recommended = lambda ctx, reads, tag, **kw: SimpleNamespace(recommendation=SimpleNamespace(trim_f=250, trim_r=240), floor=480, sub=SimpleNamespace(read_len_f=300, read_len_r=300))
    ov._thin = lambda ctx, reads, name, thin: reads
    ov.dada_on = lambda ctx, reads, tf, tr, ee, ref, tag: {"tf": tf, "tr": tr, "merged": 0.5 + 0.0005 * (tf - 250) - 0.0003 * (tr - 240) - (0.3 if tf + tr < 480 else 0.0), "variants": 1, "present": 1, "strains_found": 1, "strains_present": 1}
    res = ov.landscape_dataset(ctx, "demo", tmp_path, "2", None, 10, 10, 1000)
    grids = rs.parse_landscape_lines(res["lines"])
    assert len(grids) == 1
    g = grids[0]
    assert (g["name"], g["f0"], g["r0"], g["floor"]) == ("demo", 250, 240, 480) and len(g["cells"]) == 9
    mine = [c for c in g["cells"] if c["mine"]]
    assert len(mine) == 1 and (mine[0]["tf"], mine[0]["tr"]) == (250, 240)
    below = {(c["tf"], c["tr"]) for c in g["cells"] if not c["meets"]}
    assert below == {(240, 230)}                                                                                                        # the only cell whose total (470) is under the floor (480)
    by = {(c["tf"], c["tr"]): c["merged"] for c in g["cells"]}
    assert abs(by[(260, 230)] - (0.5 + 0.005 + 0.003)) < 0.0006                                                                         # the printed percentage, to 0.1 point
    r = rs.same_total_analysis(g)
    assert r["mine"] == "250/240" and r["best"] == "260/230" and r["df"] == 10 and r["direction"] == "forward-heavier" and r["gain"] > 0.5


def test_the_split_decision_needs_four_datasets_and_one_point():
    rs = steps()
    mk = lambda df, gain: {"df": df, "gain": gain}
    assert rs.split_decision([mk(10, 2.0), mk(10, 1.5), mk(10, 0.5), mk(10, 1.0), mk(0, 0.0), mk(-10, 0.3)])["supported"] is True
    assert rs.split_decision([mk(10, 2.0), mk(10, 1.5), mk(10, 0.5), mk(0, 0.0), mk(0, 0.0), mk(-10, 0.3)])["supported"] is False        # only three forward-heavy datasets
    assert rs.split_decision([mk(10, 0.4), mk(10, 0.5), mk(10, 0.5), mk(10, 0.3), mk(0, 0), mk(0, 0)])["supported"] is False               # four, but the mean gain is under a point
    assert rs.same_total_analysis({"name": "x", "cells": [{"tf": 5, "tr": 5, "merged": 0.5, "mine": True, "meets": True}, {"tf": 6, "tr": 6, "merged": 0.6, "mine": False, "meets": True}]}) is None


def test_the_split_step_reads_the_saved_landscape_and_gives_no_verdict_with_fewer_than_six_datasets(tmp_path):
    ov, rs = load("overnight"), steps()
    (tmp_path / "landscape").mkdir()
    lines = []
    for name, f0, r0 in (("a", 250, 240), ("b", 200, 150)):
        lines += [f"--- {name} (maxEE 2); AmpliTrim recommends {f0}/{r0} (coverage floor 100); share of pairs merged; * = AmpliTrim, . = total below the floor",
                  "  trF\\trR " + " ".join(f"{t:>7d}" for t in (r0 - 10, r0, r0 + 10))]
        for tf, vals in ((f0 - 10, (50.0, 51.0, 52.0)), (f0, (53.0, 54.0, 55.0)), (f0 + 10, (56.0, 57.0, 58.0))):
            cells = [f"{v:6.1f}{'*' if (tf, t) == (f0, r0) else ' '}" for t, v in zip((r0 - 10, r0, r0 + 10), vals)]
            lines.append(f"  {tf:>8d} " + " ".join(cells))
        lines.append("  best of the grid: ...")
    (tmp_path / "landscape" / "result.json").write_text(json.dumps({"lines": lines}))
    ctx = SimpleNamespace(out=tmp_path, cur=tmp_path)
    out = rs.step_split_tiebreak(ctx, helpers(ov))
    assert len(out["data"]) == 2 and all(d["direction"] == "forward-heavier" for d in out["data"])
    assert out["verdicts"][0]["outcome"] == "info" and "no verdict" in out["verdicts"][0]["detail"]
    (tmp_path / "landscape" / "result.json").unlink()
    with pytest.raises(RuntimeError, match="landscape"):
        rs.step_split_tiebreak(ctx, helpers(ov))


# ------------------------------------------------------------------------------------------------ Q29: the missing variant
def test_edit_positions_and_which_reads_cover_a_base():
    rs = steps()
    assert rs.edit_positions("ACGTACGT", "ACGTTCGT") == [5] and rs.edit_positions("ACGTACGT", "ACGTACGT") == []
    assert rs.edit_positions("ACGTCAGT", "ACGTAGT") == [5]                                                                              # a base of the first that the second lacks (not a repeat, so the alignment is unambiguous)
    assert rs.coverage_class(276, 500, 279, 248) == "both" and rs.coverage_class(276, 500, 274, 253) == "reverse only"
    assert rs.coverage_class(100, 500, 279, 248) == "forward only" and rs.coverage_class(260, 500, 200, 200) == "neither"
    assert rs.coverage_counts([100, 276, 400], 500, 279, 248) == {"both": 1, "forward only": 1, "reverse only": 1, "neither": 0}
    seqs = {"a": "ACGT" * 20, "b": "ACGT" * 20 + "A", "c": "TTTT" * 20}
    assert rs.nearest_sibling("a", seqs) == "b"


def test_the_missing_variant_step_finds_the_variant_one_tool_recovers_and_where_its_distinguishing_base_is_read(tmp_path, monkeypatch):
    ov, rs = load("overnight"), steps()
    import mock_check
    rng = np.random.default_rng(1)
    s1 = "".join(rng.choice(list("ACGT"), 500))
    s2 = s1[:275] + ("A" if s1[275] != "A" else "C") + s1[276:]                                                                         # differs at position 276 (1-based)
    s3 = "".join(rng.choice(list("ACGT"), 500))
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">v1|S1|x\n{s1}\n>v2|S1|x\n{s2}\n>v3|S2|x\n{s3}\n")
    folder = tmp_path / "extra_datasets" / "mockB" / "compare"
    folder.mkdir(parents=True)
    (folder / "compare_methods.tsv").write_text("name\ttruncF\ttruncR\tmaxEE\nFIGARO (amplicon 515)\t279\t248\t6,9\nAmpliTrim (FIGARO's maxEE 6,9)\t274\t253\t6,9\n")
    for key in ("FIGARO__amplicon_515_", "AmpliTrim__FIGARO_s_maxEE_6_9_"):
        (folder / f"asvs_{key}.fasta").write_text(">ASV1;size=10\nACGT\n")
    recovered = {"FIGARO__amplicon_515_": ["v1|S1|x", "v2|S1|x", "v3|S2|x"], "AmpliTrim__FIGARO_s_maxEE_6_9_": ["v1|S1|x", "v3|S2|x"]}
    monkeypatch.setattr(mock_check, "summarise", lambda asv, ref_, nid, vs: {"recovered": recovered[Path(asv).stem[len("asvs_"):]]})      # through monkeypatch: the whole suite shares one process
    ctx = SimpleNamespace(out=tmp_path, cur=tmp_path / "s", p=lambda k: None, args=SimpleNamespace(extra=[f"mockB={tmp_path / 'ready'}|{ref}"], near_id=0.99))
    out = rs.step_missing_variants(ctx, helpers(ov))
    d = [x for x in out["data"] if x["variant"] == "v2|S1|x"]
    assert len(d) == 1 and d[0]["recovered_by"] == "FIGARO" and d[0]["sibling"] == "v1|S1|x" and d[0]["positions"] == [276]
    assert d[0]["covered_recovering"]["both"] == 1 and d[0]["covered_other"]["reverse only"] == 1
    assert "1 of 1 variants" in out["verdicts"][0]["detail"] and out["verdicts"][0]["outcome"] == "info"


# ------------------------------------------------------------------------------------------------ Q30, Q23, Q16
def test_the_share_of_pairs_of_one_length_and_the_length_of_the_asvs(tmp_path):
    rs = steps()
    reads = tmp_path / "ready"
    reads.mkdir()
    pairs = [("A" * 100, "I" * 100)] * 10 + [("A" * 90, "I" * 90)] * 4 + [("A" * 110, "I" * 110)] * 2
    write_fastq_gz(reads / "S1_R1.fastq.gz", pairs)
    write_fastq_gz(reads / "S1_R2.fastq.gz", pairs)
    r = rs.pair_length_shares(reads)
    assert r["pairs"] == 16 and (r["modal_f"], r["modal_r"]) == (100, 100) and abs(r["exact"] - 10 / 16) < 1e-9 and abs(r["at_least"] - 12 / 16) < 1e-9
    ref = tmp_path / "ref.fasta"
    ref.write_text(f">a|A|x\n{'A' * 470}\n>b|B|x\n{'C' * 500}\n")
    asv = tmp_path / "asv.fasta"
    asv.write_text(f">ASV1;size=900\n{'A' * 472}\n>ASV2;size=90\n{'C' * 497}\n>ASV3;size=10\n{'G' * 520}\n")
    s = rs.asv_length_shares(asv, ref)
    assert s["reads"] == 1000 and abs(s["within"] - 0.99) < 1e-9 and abs(s["off_share"] - 0.01) < 1e-9 and s["off"][0][:2] == (10, 520)


def test_the_expected_errors_at_a_cut_and_the_rule_check(tmp_path):
    rs = steps()
    from amplitrim.trimmer import ReadPairs
    q_f, q_r = "I" * 100 + "+" * 100, "5" * 200                                                                                          # Q40 then Q10; Q20 throughout
    recs = [(f"@r{i}", "A" * 200, "+", q_f) for i in range(20)], [(f"@r{i}", "A" * 200, "+", q_r) for i in range(20)]
    pairs = ReadPairs(*recs)
    ef, er = rs.mean_ee_at_cut(pairs, 150, 150)
    assert abs(ef - (100 * 1e-4 + 50 * 0.1)) < 1e-6 and abs(er - 150 * 0.01) < 1e-6
    row = rs.rule_check_row(pairs, 150, 150, (2.0, 2.5))
    assert abs(row["pred_f"] - (ef + 1)) < 1e-9 and row["within"] is False                                                              # the forward prediction (6.0) is far from 2
    assert rs.rule_check_row(pairs, 150, 150, (round(ef + 1), round(er + 1)))["within"] is True
    short = ReadPairs([("@r", "A" * 50, "+", "I" * 50)], [("@r", "A" * 50, "+", "I" * 50)])
    assert rs.mean_ee_at_cut(short, 100, 100)[0] != rs.mean_ee_at_cut(short, 100, 100)[0]                                              # no read reaches the cut: nan, not an error


# ------------------------------------------------------------------------------------------------ the simulations (decisions and steps with fakes)
def test_the_richness_and_overmerge_decisions_follow_the_preset_rules():
    rs = steps()
    mk = lambda n, ex, ee="2": {"n_taxa": n, "ee": ee, "excess": ex}
    assert rs.richness_decision([mk(25, 0.3), mk(25, 0.5), mk(500, 14.0), mk(500, 16.0)])["outcome"] == "supported"
    assert rs.richness_decision([mk(25, 0.3), mk(500, 6.0)])["outcome"] == "partly supported"
    assert rs.richness_decision([mk(25, 0.3), mk(500, 1.0)])["outcome"] == "not supported"
    assert rs.richness_decision([mk(25, 0.3)])["outcome"] is None and rs.richness_decision([mk(25, 0.3, "1"), mk(1000, 20.0, "1")])["outcome"] is None
    assert rs.richness_decision([mk(25, 0.3), mk(250, 5.0), mk(1000, 30.0)])["high"] == 1000
    assert abs(rs.spearman([1, 2, 3], [10, 20, 30]) - 1) < 1e-9 and abs(rs.spearman([1, 2, 3], [30, 20, 10]) + 1) < 1e-9
    row = lambda ov, wrong, rec: {"arm": "realistic", "ov": ov, "wrong_clean": wrong, "wrong": 0.5, "recovered": rec}              # the share of all pairs is not what decides
    assert rs.overmerge_decision([row(12, 0.004, 20), row(20, 0.001, 20)])["ok"] is True
    assert rs.overmerge_decision([row(12, 0.02, 20), row(20, 0.001, 20)])["ok"] is False                                                # 2% wrong-length merges at 12 nt
    assert rs.overmerge_decision([row(12, 0.004, 17), row(20, 0.001, 20)])["ok"] is False                                               # three taxa fewer at 12 nt
    assert rs.overmerge_decision([row(12, 0.004, 19), row(20, 0.001, 20)])["ok"] is True                                                # one fewer is within the loss allowed
    assert rs.overmerge_decision([row(12, 0.004, 20), row(20, 0.001, 17)])["ok"] is True                                                # 12 nt recovering MORE taxa than 20 nt is not a failure (the quick run did this)
    assert rs.overmerge_decision([])["ok"] is None


def test_wrong_length_merges_are_counted_by_kind(tmp_path):
    rs = steps()
    lengths = {"T1": 100, "T2": 120}
    recs = [("simS1_1|T1|normal", "A" * 100), ("simS1_2|T1|normal", "A" * 101), ("simS1_3|T2|repeat", "A" * 99), ("simS1_4|T1xT2|chimera", "A" * 50), ("simS1_5|T2|repeat", "A" * 120)]
    r = rs.wrong_length_share(recs, lengths)
    assert r["n"] == 4 and abs(r["wrong"] - 0.5) < 1e-9 and abs(r["normal"] - 0.5) < 1e-9 and abs(r["repeat"] - 0.5) < 1e-9             # the chimera is left out
    assert abs(r["wrong_clean"] - 0.5) < 1e-9 and r["n_indel"] == 0
    r = rs.wrong_length_share(recs + [("simS1_6|T1|normal|indel", "A" * 101), ("simS1_7|T1|normal|indel", "A" * 100)], lengths)         # pairs that carry an indel are counted apart
    assert r["n"] == 6 and r["n_indel"] == 2 and abs(r["wrong_indel"] - 0.5) < 1e-9 and abs(r["wrong_clean"] - 0.5) < 1e-9 and abs(r["wrong"] - 3 / 6) < 1e-9
    f = tmp_path / "m.fasta"
    f.write_text(">simS1_1|T1|normal extra\nAAAA\n>simS1_2|T2|repeat\nCC\nCC\n")
    assert rs.merged_records(f) == [("simS1_1|T1|normal", "AAAA"), ("simS1_2|T2|repeat", "CCCC")]
    q = tmp_path / "m.fastq"
    q.write_text("@simS1_1|T1|normal x\nACGT\n+\nIIII\n")
    assert rs.merged_records(q) == [("simS1_1|T1|normal", "ACGT")]


def _fake_tools(monkeypatch, rs, tmp_path):
    import real_data
    seen = {"dada2": []}

    def fake_dada2(reads, tf, tr, ee, ov, min_len, track, asvs):
        seen["dada2"].append((tf, tr, ee, ov))
        Path(asvs).write_text(">ASV1;size=5\nACGT\n")
        return {"filtered": 0.9, "merged": 0.6 + 0.05 * float(str(ee).split(",")[0]) / 5 - 0.001 * (tf % 7), "asvs": 3, "nonchim": 0.5}
    monkeypatch.setattr(real_data, "dada2_run", fake_dada2)
    monkeypatch.setattr(rs, "score_asvs", lambda asv, truth, info, near: {"recovered": 18, "classes": {"abundant": (10, 10), "medium": (6, 8), "rare": (2, 6)}, "identical": 0.97, "off": 0.01, "asvs": 22})
    return seen


def test_the_richness_step_runs_every_level_and_maxee_writes_its_table_and_is_cautious_in_quick_mode(tmp_path, monkeypatch):
    ov, rs = load("overnight"), steps()
    seen = _fake_tools(monkeypatch, rs, tmp_path)
    H = helpers(ov)
    H.recommended = lambda ctx, reads, tag, **kw: SimpleNamespace(recommendation=SimpleNamespace(trim_f=240, trim_r=230, yield_est=0.80 + 0.01 * float(kw["max_ee"] if not isinstance(kw["max_ee"], tuple) else kw["max_ee"][0])))
    ctx = SimpleNamespace(out=tmp_path, cur=tmp_path / "step", quick=True, p=lambda k: None, log=lambda *a: None, args=SimpleNamespace(near_id=0.99))
    ctx.cur.mkdir()
    out = rs.step_richness_sim(ctx, H)
    assert len(out["data"]) == 2 * 1 * 3                                                                                                # two levels, one seed, three maxEE
    assert len(seen["dada2"]) == 6 and (ctx.cur / "richness_sim.tsv").is_file()
    tsv = (ctx.cur / "richness_sim.tsv").read_text().splitlines()
    assert tsv[0].startswith("n_taxa\tseed\tee") and len(tsv) == 7
    assert "supported" in out["verdicts"][0]["detail"] or "too few" in out["verdicts"][0]["detail"]
    assert out["verdicts"][0]["outcome"] == "info" and "quick run" in out["verdicts"][0]["detail"]                                      # no verdict from a quick run
    assert any("rank correlation" in ln for ln in out["lines"])
    again = rs.step_richness_sim(ctx, H)                                                                                                # a repeat (--fresh, or after an interruption) clears the step's scratch folders and works
    assert len(again["data"]) == len(out["data"]) and len(seen["dada2"]) == 12


def test_the_overmerge_step_sets_the_truncation_for_each_overlap_and_scores_the_merges(tmp_path, monkeypatch):
    ov, rs = load("overnight"), steps()
    seen = _fake_tools(monkeypatch, rs, tmp_path)
    calls = []

    def fake_filter(reads, out_dir, tf, tr):
        calls.append((tf, tr))
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        return Path(out_dir) / "R1.fastq", Path(out_dir) / "R2.fastq"

    def fake_merge(r1, r2, out_dir, overlap):
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        f = Path(out_dir) / "merged.fasta"
        import community_sim
        f.write_text(">simS1_1|T0001|normal\n" + "A" * 460 + "\n>simS1_2|T0001|normal\n" + "A" * 461 + "\n")
        return SimpleNamespace(merged_path=f, merge_rate=0.9)
    monkeypatch.setattr(rs, "filter_reads", fake_filter)
    monkeypatch.setattr(rs, "run_merge", fake_merge)
    import community_sim
    real_write = community_sim.write_community

    def small(out_dir, n_taxa, seed, n_samples, pairs, **kw):
        info = real_write(out_dir, n_taxa, seed, n_samples, pairs, **kw)
        info["lengths"]["T0001"] = 460
        return info
    monkeypatch.setattr(community_sim, "write_community", small)
    H = helpers(ov)
    ctx = SimpleNamespace(out=tmp_path, cur=tmp_path / "step", quick=True, p=lambda k: None, log=lambda *a: None, args=SimpleNamespace(near_id=0.99))
    ctx.cur.mkdir()
    out = rs.step_overmerge_sim(ctx, H)
    assert len(rs.step_overmerge_sim(ctx, H)["data"]) == len(out["data"])                                                               # a repeat works too
    calls.clear()
    out = rs.step_overmerge_sim(ctx, H)
    assert len(out["data"]) == 3 * 1 * 3 and len(calls) == 9                                                                            # three arms, one seed, three overlaps
    tfs = sorted({(tf + tr) for tf, tr in calls})
    assert len(tfs) == 3 and tfs[1] - tfs[0] == 4 and tfs[2] - tfs[1] == 4                                                              # the total grows with the overlap: Lmax + 12, 16, 20
    assert {d["ov"] for d in out["data"]} == {12, 16, 20} and all(abs(d["wrong"] - 0.5) < 1e-9 for d in out["data"])
    assert {w[3] for w in seen["dada2"]} == {12, 16, 20}                                                                                # DADA2 gets the minimum overlap too
    assert (ctx.cur / "overmerge_sim.tsv").is_file() and out["verdicts"][-1]["outcome"] == "info"


def test_the_new_steps_are_in_the_registry_as_tier_four():
    ov = load("overnight")
    ids = {s[0]: s for s in ov.STEPS}
    new = ("split_tiebreak", "missing_variants", "figaro_usable_pairs", "asv_length_check", "figaro_rule_check", "richness_sim", "overmerge_sim")
    assert all(i in ids and ids[i][1] == 4 and callable(ids[i][6]) for i in new)
    import round4
    assert [s for _k, _t, ss in round4.PHASES for s in ss if ids.get(s, (0, 0))[1] == 4] == list(new)                                  # every new step is in exactly one phase, in order


# ------------------------------------------------------------------------------------------------ the pauses and the phased runner
class Clock:
    """A fake clock: time passes only when the pause sleeps."""
    def __init__(self):
        self.t = 0.0

    def now(self):
        return self.t

    def sleep(self, s):
        self.t += s


def keys(*script):
    """A key reader that returns the scripted keys one per call (then nothing); each is (time of availability, key)."""
    items, clock = list(script), None

    def make(c):
        def read():
            if items and c.now() >= items[0][0]:
                return items.pop(0)[1]
            return ""
        return read
    return make


def test_the_break_goes_on_by_itself_and_obeys_the_keys():
    import round4
    for script, expect, minimum in (((), "continue", 90), (((5, "s"),), "skip", 5), (((5, "q"),), "quit", 5), (((5, "\n"),), "continue", 5), (((5, "c"),), "continue", 5),
                                    (((5, "p"), (50, "\n")), "continue", 50), (((5, "p"), (60, "s")), "skip", 60)):
        c = Clock()
        said = []
        action = round4.pause_between(90, keys(*script)(c), now=c.now, sleep=c.sleep, say=said.append)
        assert action == expect and c.t >= minimum and (c.t < minimum + 1 if script else c.t >= 90), (script, action, c.t)
    c = Clock()
    assert round4.pause_between(0, lambda: "s", now=c.now, sleep=c.sleep, say=lambda *a: None) == "continue" and c.t == 0                # no break: nothing waits
    c = Clock()
    said = []
    round4.pause_between(30, keys()(c), now=c.now, sleep=c.sleep, say=said.append)
    assert any("continuing in" in m for m in said)


def _fake_overnight(out: Path, calls: list):
    """A runner that behaves like overnight.py: it writes a result.json for each step in --only."""
    def run(cmd):
        calls.append(cmd)
        only = cmd[cmd.index("--only") + 1].split(",")
        for step in only:
            d = out / step
            d.mkdir(parents=True, exist_ok=True)
            verdicts = {"split_tiebreak": [{"criterion": "a forward-heavy split beats the balanced tie-break (best same-total cell forward-heavier in at least 4 datasets)", "outcome": "met", "detail": "5 of 6"}],
                        "figaro_rule_check": [{"criterion": "the rule 'mean expected errors at the cut + 1' reproduces FIGARO's maxEE", "outcome": "NOT MET", "detail": "2 of 6 datasets"}],
                        "richness_sim": [{"criterion": "the excess of the yield estimate grows with richness", "outcome": "met", "detail": "supported (+14.0 points)"},
                                         {"criterion": "the estimate ranks the three maxEE settings as DADA2 does", "outcome": "met", "detail": "lowest 1.00"}]}.get(step, [])
            (d / "result.json").write_text(json.dumps({"id": step, "status": "ok", "seconds": 90, "verdicts": verdicts}))
        return 0
    return run


def test_the_runner_runs_the_phases_in_order_with_breaks_and_writes_the_report(tmp_path):
    import round4
    out = tmp_path / "out"
    cfg = {"out": str(out), "extras": {}}
    calls, said = [], []
    c = Clock()
    args = SimpleNamespace(quick=False, fresh=False, break_seconds=60, phases=None)
    code = round4.run_round4(cfg, args, runner=_fake_overnight(out, calls), read_key=keys()(c), sleep=c.sleep, now=c.now, say=said.append)
    assert code == 0 and [cmd[cmd.index("--only") + 1].split(",")[0] for cmd in calls] == ["split_tiebreak", "richness_sim", "battery"]
    assert c.t >= 120                                                                                                                   # two breaks of 60 s, and nobody pressed anything
    report = (out / "ROUND4_REPORT.md").read_text()
    assert "## Phase A" in report and "## Phase B" in report and "## Phase C" in report
    assert "-> supported: add the opt-in `--split-preference forward`" in report and "-> FIGARO's rule is not reproduced" in report
    assert report.count("the excess grows with richness") == 1                                                                         # the follow-up belongs to the decision verdict only
    assert "--fresh" not in calls[0] and "--quick" not in calls[0]


def test_the_runner_can_skip_a_phase_stop_early_or_start_fresh(tmp_path):
    import round4
    out = tmp_path / "out"
    cfg = {"out": str(out), "extras": {}}
    calls = []
    c = Clock()
    code = round4.run_round4(cfg, SimpleNamespace(quick=False, fresh=True, break_seconds=60, phases="A,B,C"), runner=_fake_overnight(out, calls), read_key=keys((5, "s"))(c), sleep=c.sleep, now=c.now, say=lambda *a: None)
    assert [cmd[cmd.index("--only") + 1].split(",")[0] for cmd in calls] == ["split_tiebreak", "battery"] and "--fresh" in calls[0]       # B skipped at the first break
    assert "phase B skipped at your request" in (out / "ROUND4_REPORT.md").read_text()
    calls.clear()
    c = Clock()
    round4.run_round4(cfg, SimpleNamespace(quick=False, fresh=False, break_seconds=60, phases="A,B,C"), runner=_fake_overnight(out, calls), read_key=keys((5, "q"))(c), sleep=c.sleep, now=c.now, say=lambda *a: None)
    assert len(calls) == 1 and "stopped at your request after phase A" in (out / "ROUND4_REPORT.md").read_text()                       # the report is written all the same
    calls.clear()
    c = Clock()
    round4.run_round4({"out": str(out), "extras": {}}, SimpleNamespace(quick=True, fresh=False, break_seconds=0, phases="B"), runner=_fake_overnight(Path(str(out) + "_quick"), calls), read_key=keys()(c),
                      sleep=c.sleep, now=c.now, say=lambda *a: None)
    assert len(calls) == 1 and "--quick" in calls[0] and calls[0][calls[0].index("--out") + 1].endswith("out_quick") and c.t == 0         # --quick has its own folder; no breaks with 0 s
    with pytest.raises(SystemExit):
        round4.run_round4(cfg, SimpleNamespace(quick=False, fresh=False, break_seconds=0, phases="Z"), runner=lambda cmd: 0, say=lambda *a: None)


def test_the_steps_that_read_earlier_results_say_clearly_what_is_missing_when_nothing_is_on_hand(tmp_path):
    ov, rs = load("overnight"), steps()
    H = helpers(ov)
    ctx = SimpleNamespace(out=tmp_path, cur=tmp_path / "s", quick=False, p=lambda k: tmp_path / "nowhere" / k, log=lambda *a: None, args=SimpleNamespace(extra=None, near_id=0.99))
    ctx.cur.mkdir()
    assert rs.datasets_with_comparison(ctx, H) == []                                  # paths that point at nothing are not datasets
    for fn, text in ((rs.step_missing_variants, "no dataset"), (rs.step_asv_length_check, "no mock community"), (rs.step_figaro_rule_check, "no dataset"), (rs.step_figaro_usable_pairs, "no dataset")):
        with pytest.raises(RuntimeError, match=text):
            fn(ctx, H)
    reads = tmp_path / "ready"
    reads.mkdir()
    write_fastq_gz(reads / "S1_R1.fastq.gz", [("A" * 50, "I" * 50)] * 5)
    write_fastq_gz(reads / "S1_R2.fastq.gz", [("A" * 50, "I" * 50)] * 5)
    ctx.p = lambda k: reads if k == "v1v2" else None
    out = rs.step_figaro_usable_pairs(ctx, H)                                          # one folder on hand: the table has one row
    assert len(out["data"]) == 1 and out["data"][0]["name"] == "v1v2" and out["verdicts"][0]["outcome"] == "NOT MET"           # all 50 nt reads: 100% of the pairs have the common length, not about 47%


def test_the_ordinary_overnight_run_does_not_start_the_round_four_steps_but_naming_them_does(tmp_path, monkeypatch):
    ov = load("overnight")
    ran = []
    monkeypatch.setattr(ov, "run_steps", lambda ctx, selected: ran.append(list(selected)) or [])
    monkeypatch.setattr(ov, "write_summary", lambda ctx, results: tmp_path / "x")
    ov.main(["--out", str(tmp_path / "o1")])
    assert ran[-1] and not any(s in ran[-1] for s in ("richness_sim", "overmerge_sim", "split_tiebreak")) and "battery" in ran[-1]
    ov.main(["--out", str(tmp_path / "o2"), "--only", "split_tiebreak,richness_sim"])
    assert ran[-1] == ["split_tiebreak", "richness_sim"]


def test_the_board_lists_round_four_apart_from_what_overnight_runs(tmp_path):
    cb = load("check_build")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "out": str(tmp_path / "out"), "extras": {}}
    board = "\n".join(cb.todo(cfg, SimpleNamespace(quick=False)))
    assert "split_tiebreak" in board and "(round 4: --round4)" in board
    nxt = [ln for ln in board.splitlines() if "check_build.py --overnight" in ln]
    assert all("richness_sim" not in ln and "split_tiebreak" not in ln for ln in nxt)
    assert "--round4 --quick" in board


def test_a_quick_run_reads_the_earlier_results_of_the_full_run(tmp_path):
    ov, rs = load("overnight"), steps()
    full, quick = tmp_path / "overnight", tmp_path / "overnight_quick"
    (full / "landscape").mkdir(parents=True)
    quick.mkdir()
    (full / "landscape" / "result.json").write_text(json.dumps({"lines": []}))
    folder = full / "extra_datasets" / "mockB" / "compare"
    folder.mkdir(parents=True)
    (folder / "compare_methods.tsv").write_text("name\ttruncF\ttruncR\tmaxEE\n")
    assert rs.results_folder(SimpleNamespace(out=quick)) == full and rs.results_folder(SimpleNamespace(out=full)) == full
    assert rs.results_folder(SimpleNamespace(out=tmp_path / "x_quick")) == tmp_path / "x_quick"                                       # no full folder: its own
    ctx = SimpleNamespace(out=quick, cur=quick, quick=True, p=lambda k: None, args=SimpleNamespace(extra=[f"mockB={tmp_path / 'ready'}|{tmp_path / 'ref.fasta'}"], near_id=0.99))
    assert [d[0] for d in rs.datasets_with_comparison(ctx, helpers(ov))] == ["mockB"]
    with pytest.raises(RuntimeError, match="no landscape grid"):                                                                        # found in the full folder (it has no grids), not "has not run"
        rs.step_split_tiebreak(ctx, helpers(ov))
