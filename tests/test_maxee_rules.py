"""validation/maxee_rules.py: simple maxEE rules against FIGARO's, checked on reads whose expected errors are known (no DADA2, no vsearch)."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))


def load():
    spec = importlib.util.spec_from_file_location("maxee_rules", HERE.parent / "validation" / "maxee_rules.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pairs_of(mr, qf: int, qr: int, length: int = 250, n: int = 200, short: int = 0):
    """Read pairs whose every base has Phred quality qf (forward) or qr (reverse): expected errors = length * 10^(-Q/10). `short` reads are 100 nt."""
    from amplitrim.trimmer import ReadPairs
    def rec(i, q, ln):
        return (f"r{i}", "ACGT" * (ln // 4) + "A" * (ln % 4), "+", chr(33 + q) * ln)
    f = [rec(i, qf, 100 if i < short else length) for i in range(n)]
    r = [rec(i, qr, 100 if i < short else length) for i in range(n)]
    return ReadPairs(f, r)


def test_figaros_pad_is_truncation_plus_one_not_a_ceiling():
    mr = load()
    assert [mr.figaro_pad(x) for x in (-0.3, 0.4, 0.99, 1.0, 2.3, 2.6, 4.35)] == [1, 1, 1, 2, 3, 3, 5]
    assert all(mr.figaro_pad(x) == int(x) + 1 for x in (0.0, 0.7, 3.2, 8.9))


def test_expected_errors_are_counted_over_the_cut_and_only_for_reads_long_enough():
    mr = load()
    pairs = pairs_of(mr, qf=20, qr=30, length=250, n=100, short=40)            # 40 reads are only 100 nt long
    ef, er = mr.read_errors(pairs, 200, 150)
    assert len(ef) == len(er) == 60                                            # the short reads cannot reach the cut and are not counted
    assert np.allclose(ef, 200 * 0.01) and np.allclose(er, 150 * 0.001)       # Q20 = 0.01 expected errors per base, Q30 = 0.001


def test_the_rule_is_the_mean_expected_errors_plus_one_and_the_floor_keeps_dadas_default():
    mr = load()
    noisy = pairs_of(mr, qf=20, qr=15)                                         # 250 bases: 2.5 and 7.9 expected errors
    ef, er = mr.read_errors(noisy, 250, 250)
    assert mr.rules(ef, er)["mean"] == (3, 8) and mr.rules(ef, er)["mean_floored"] == (3, 8)
    clean = pairs_of(mr, qf=40, qr=35)                                         # 0.025 and 0.08 expected errors
    ef, er = mr.read_errors(clean, 250, 250)
    assert mr.rules(ef, er)["mean"] == (1, 1)                                  # FIGARO's 1 on clean reads ...
    assert mr.rules(ef, er)["mean_floored"] == (2, 2)                          # ... and the default where the default is enough


def test_the_median_and_the_mean_differ_on_a_skewed_set_of_reads():
    mr = load()
    ef = np.array([0.2] * 90 + [12.0] * 10)                                    # most reads clean, a few very noisy
    r = mr.rules(ef, ef)
    assert r["median"] == (1, 1) and r["mean"] == (2, 2)                       # mean 1.38, median 0.2


def test_within_one_on_each_read():
    mr = load()
    assert mr.within((3, 4), (3, 4)) and mr.within((2, 1), (1, 1)) and mr.within((6, 9), (5, 8))
    assert not mr.within((6, 9), (4, 9)) and not mr.within((2, 2), (2, 4))
    assert mr.within(2.0, (1, 3))                                              # a single maxEE means the same limit on both reads


def _fake_results(mr, mean_matches, floored_matches, clean_default, median_matches):
    """Results for all ten datasets with the given numbers of matches, to check that the verdicts count them as the register says."""
    out = {}
    for i, k in enumerate(mr.ROWS):
        out[k] = {"figaro": {"match": {"mean": i < mean_matches, "median": i < median_matches}},
                  "amplitrim": {"match": {"mean_floored": k in mr.POOR[:floored_matches]},
                                "rules": {"mean_floored": (2, 2) if k in mr.CLEAN[:clean_default] else (3, 3)}}}
    return out


def test_verdicts_follow_the_registered_thresholds():
    mr = load()
    assert len(mr.POOR) == 4 and len(mr.CLEAN) == 6                           # mock A, mock B, V4-V6, wastewater; and the six where FIGARO is at 2 or less
    met = dict((p, ok) for p, _, ok in mr.verdicts(_fake_results(mr, 8, 3, 5, 8)))
    assert met == {"P13": True, "P14": True, "P15": True, "P16": True}
    missed = dict((p, ok) for p, _, ok in mr.verdicts(_fake_results(mr, 7, 2, 4, 9)))
    assert missed == {"P13": False, "P14": False, "P15": False, "P16": False}
    partial = dict((p, ok) for p, _, ok in mr.verdicts({k: v for k, v in _fake_results(mr, 8, 3, 5, 8).items() if k != "v4"}))
    assert partial["P13"] is None                                              # not decided unless all ten datasets ran


def test_the_analysis_runs_on_simulated_reads_and_never_picks_below_the_default(reads, tmp_path):
    mr = load()
    from amplitrim.subsample import build_subsample
    from amplitrim.trimmer import ReadPairs
    info = build_subsample(reads, tmp_path, max_pairs=3000)
    pairs = ReadPairs.load(info.r1, info.r2)
    res = mr.analyse(pairs, dict(fig=(225, 200, (2, 1)), dada2_pass=0.9, amp=(231, 215)))
    for block in ("figaro", "amplitrim"):
        r = res[block]["rules"]
        assert min(r["mean_floored"]) >= 2 and min(r["median_floored"]) >= 2 and min(r["mean"]) >= 1
    assert set(res["figaro"]["q"]) == set(mr.QGRID) and set(res["amplitrim"]["q"]) == set(mr.QGRID)
    assert all(min(x["ee"]) >= 2 and 0 < x["pass"] <= 1 for x in res["amplitrim"]["q"].values())        # the floored rule never goes below 2
    pts = res["amplitrim"]["points"]
    assert set(pts) == {"default (2,2)", "mean rule, floored", "FIGARO's pick"}
    assert 0 < pts["default (2,2)"]["pass"] <= 1 and pts["FIGARO's pick"]["pass"] <= pts["default (2,2)"]["pass"] + 1e-9   # (2,1) is stricter than (2,2)
    assert "### " in "\n".join(mr.dataset_lines("mockB", res))


def test_a_dataset_that_is_not_registered_is_skipped_and_nothing_run_is_reported(tmp_path, capsys):
    mr = load()
    cfg = tmp_path / "paths.json"
    cfg.write_text("{}")
    assert mr.main(["--config", str(cfg), "--out", str(tmp_path / "o")]) == 2
    assert "nothing was run" in capsys.readouterr().out


def test_reads_too_short_for_figaros_cut_are_reported_and_skipped(reads, tmp_path, capsys):
    """A wrong reads folder (shorter than the cut FIGARO used) must be said plainly, not crash."""
    import json
    mr = load()
    cfg = tmp_path / "paths.json"
    cfg.write_text(json.dumps({"mockB": str(reads)}))                          # mock B's FIGARO cut is 279/248; the simulated reads are 231/227
    assert mr.main(["--config", str(cfg), "--only", "mockB", "--out", str(tmp_path / "o")]) == 2
    out = capsys.readouterr().out
    assert "SKIPPED" in out and "no reads are long enough" in out and "279/248" in out


def test_reads_folders_are_found_at_the_top_level_and_under_extras(tmp_path, monkeypatch):
    """Night 4's first run skipped eight datasets because only the top-level keys were read; night2.py reads both places."""
    mr = load()
    monkeypatch.setenv("DATA", str(tmp_path))
    cfg = {"mockA_high": "C:/data/a", "v1v2": "C:/data/b", "extras": {"mockB": {"ready": "%s/mockB/ready" % "C:/data"}, "v4": {"ready": "$DATA/v4/ready"}, "its": {}}}
    assert mr.reads_folder(cfg, "mockA_high") == Path("C:/data/a")
    assert mr.reads_folder(cfg, "mockB") == Path("C:/data/mockB/ready")
    assert mr.reads_folder(cfg, "v4") == tmp_path / "v4" / "ready"                # environment variables are expanded
    assert mr.reads_folder(cfg, "swabs_V3-V4") is None and mr.reads_folder({}, "mockB") is None
    assert mr.reads_folder({"mockB": "C:/top", "extras": {"mockB": {"ready": "C:/extra"}}}, "mockB") == Path("C:/top")


def test_the_implied_scale_is_the_intersection_of_what_each_pick_allows():
    mr = load()
    def res(pick, f_stat, r_stat):
        return {"figaro": {"pick": pick, "stats": {name: (f_stat, r_stat) for name in mr.STATS}}}
    both = {"mockA": res((3.0, 4.0), 1.30, 1.95), "v1v2": res((2.0, 1.0), 0.70, 0.24)}   # the first partial run's numbers
    w = mr.implied_windows(both, "mean")
    assert w["feasible"] and w["lo"] == pytest.approx(3 / 1.95) and w["hi"] == pytest.approx(4 / 1.95)    # k from 1.54 up to, not including, 2.05
    assert w["n"] == 4 and w["exact"] == 1                                       # the mean alone (k = 1) reproduces only v1v2's reverse pick
    assert not mr.implied_windows({"a": res((3.0, 4.0), 1.0, 1.0), "b": res((1.0, 1.0), 2.0, 2.0)}, "mean")["feasible"]   # one dataset wants k >= 2, the other k < 0.5
    table = "\n".join(mr.exploratory_lines(both))
    assert "| mean | 1.54 to 2.05 | 1 of 4 |" in table and "not pre-registered" in table
    assert "none (the windows do not overlap" in "\n".join(mr.exploratory_lines({"a": res((3.0, 4.0), 1.0, 1.0), "b": res((1.0, 1.0), 2.0, 2.0)}))


def test_the_percentile_rule_is_the_percentile_plus_one_and_never_below_two_when_floored():
    mr = load()
    ef = np.arange(0, 10, 0.1)                                                   # 100 reads with 0.0 ... 9.9 expected errors
    raw, floored = mr.q_rule(ef, ef * 0.5, 0.80)                                 # the 80th percentile is 7.92 (forward) and 3.96 (reverse)
    assert raw == (8, 4) and floored == (8, 4)
    raw, floored = mr.q_rule(np.full(50, 0.1), np.full(50, 0.3), 0.85)
    assert raw == (1, 1) and floored == (2, 2)                                   # clean reads: FIGARO's 1, and DADA2's default after the floor


def test_the_percentile_window_is_the_intersection_over_every_read():
    mr = load()
    res = lambda *cdf: {"figaro": {"cdf": cdf}}
    ok = mr.q_window({"a": res((0.50, 0.90), (0.70, 0.95)), "b": res((0.0, 0.88), (0.60, 0.99))})
    assert ok["feasible"] and ok["lo"] == 0.70 and ok["hi"] == 0.88 and ok["n"] == 4
    bad = mr.q_window({"a": res((0.85, 0.90), (0.0, 0.8))})                      # one read needs q above 0.85, another needs q of 0.8 or less
    assert not bad["feasible"]


def _fake_q_results(mr, matches, clean22, reach):
    """Ten datasets: `matches` of them where q = 0.85 reproduces FIGARO within 1; `clean22` of the six clean ones at (2,2) and `reach` of the four poor ones at 70% passing."""
    out = {}
    for i, k in enumerate(mr.ROWS):
        match = {q: (i < matches if q == 0.85 else i < 3) for q in mr.QGRID}                   # the other q values match on 3 datasets only
        out[k] = {"figaro": {"match": {"mean": True, "median": True}, "q": {q: {"raw": (3, 3), "match": match[q]} for q in mr.QGRID}},
                  "amplitrim": {"match": {"mean_floored": True}, "rules": {"mean_floored": (2, 2)},
                                "q": {q: {"ee": (2, 2) if (k in mr.CLEAN[:clean22]) else (3, 3), "pass": 0.75 if (k in mr.POOR[:reach] or k in mr.CLEAN) else 0.5} for q in mr.QGRID}}}
    return out


def test_p17_and_p18_follow_the_registered_thresholds_and_the_best_q_wins_ties_to_the_smaller_q():
    mr = load()
    v = {p: ok for p, _, ok in mr.verdicts(_fake_q_results(mr, matches=9, clean22=5, reach=3))}
    assert v["P17"] is True and v["P18"] is True
    v = {p: ok for p, _, ok in mr.verdicts(_fake_q_results(mr, matches=8, clean22=5, reach=3))}
    assert v["P17"] is False                                                     # 8 of 10 is below the 9 that P17 needs
    v = {p: ok for p, _, ok in mr.verdicts(_fake_q_results(mr, matches=9, clean22=4, reach=3))}
    assert v["P18"] is False                                                     # the rule relaxes on too many clean datasets
    v = {p: ok for p, _, ok in mr.verdicts(_fake_q_results(mr, matches=9, clean22=5, reach=2))}
    assert v["P18"] is False                                                     # and does not reach 70% on enough poor ones
    tie = _fake_q_results(mr, matches=3, clean22=6, reach=4)                     # every q matches on 3 datasets
    assert mr.best_q(tie)["q"] == min(mr.QGRID)


def test_the_percentile_table_is_written_with_every_q():
    mr = load()
    text = "\n".join(mr.percentile_lines({k: dict(v, figaro=dict(v["figaro"], cdf=((0.5, 0.9), (0.6, 0.95)))) for k, v in _fake_q_results(mr, 9, 6, 4).items()}))
    assert "## The percentile rule" in text and "q from 0.60 to 0.90" in text
    assert all(f"| {q:g} |" in text for q in mr.QGRID)


def test_the_cap_limits_the_floored_rule_and_never_touches_the_raw_value():
    mr = load()
    noisy = np.full(100, 12.0)                                                   # 12 expected errors on every read: raw maxEE 13
    raw, floored = mr.q_rule(noisy, noisy * 0.5, 0.75, cap=5)                    # forward 13 -> 5; reverse 7 -> 5
    assert raw == (13, 7) and floored == (5, 5)
    raw, floored = mr.q_rule(np.full(100, 3.5), np.full(100, 0.1), 0.75, cap=5)  # under the cap and above the floor: untouched
    assert raw == (4, 1) and floored == (4, 2)
    assert mr.q_rule(noisy, noisy, 0.75)[1] == (13, 13)                          # no cap: as before


def test_the_variants_table_lists_only_datasets_where_a_variant_leaves_the_default(reads, tmp_path):
    mr = load()
    from amplitrim.subsample import build_subsample
    from amplitrim.trimmer import ReadPairs
    info = build_subsample(reads, tmp_path, max_pairs=3000)
    pairs = ReadPairs.load(info.r1, info.r2)
    res = mr.analyse(pairs, dict(fig=(225, 200, (2, 1)), dada2_pass=0.9, amp=(231, 215)))
    v = res["amplitrim"]["variants"]
    assert set(v) == {"0.75|0", "0.75|5", "0.8|0", "0.8|5"} and all(min(x["ee"]) >= 2 and max(x["ee"]) <= 13 for x in v.values())
    assert all(max(x["ee"]) <= 5 for key, x in v.items() if key.endswith("|5"))
    quiet = {"clean": {"amplitrim": {"variants": {k: dict(x, ee=(2.0, 2.0)) for k, x in v.items()}, "points": res["amplitrim"]["points"]}}}
    assert "| clean |" not in "\n".join(mr.variants_lines(quiet))                  # nothing to show where every variant is (2,2)
    loud = {"mockB": {"amplitrim": {"variants": {k: dict(x, ee=(4.0, 5.0), **{"pass": 0.6}) for k, x in v.items()}, "points": res["amplitrim"]["points"]}}}
    text = "\n".join(mr.variants_lines(loud))
    assert "## Variants with a cap (exploratory" in text and "| mockB |" in text and "4,5: 60.0%" in text

