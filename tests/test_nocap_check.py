"""validation/nocap_check.py: the rule without a cap (P39), the floor-only cut (P43), reads lost to length (P45), and the evaluation of the mock B table (P40, P41). The package and the reads are stand-ins."""
import csv
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "validation"))
import test_release_check as trc                       # noqa: E402  (the stand-ins of the release check: the nights, the config, what was tested)


def load():
    spec = importlib.util.spec_from_file_location("nocap_check", HERE.parent / "validation" / "nocap_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


MOCKB_NEW = (270, 256, (4.0, 8.0))                     # the released package without the cap on mock B: a longer reverse limit, a different cut


def released_settings(mockb=MOCKB_NEW, **override):
    got = trc.as_released(trc.what_was_tested())
    got["mockB"] = mockb
    got.update(override)
    return got


def stand_ins(released, floors=None, floor_only=None, short=None, noise=None):
    """A package and a reads loader that know the datasets by the name of their folder. floor_only: the best split of the floor, by dataset (default: the released cut)."""
    floors, floor_only, short = floors or {}, floor_only or {}, short or {}

    def run(settings):
        name = Path(settings.input_dir).name
        got = released[name]
        return SimpleNamespace(recommendation=None if got is None else SimpleNamespace(trim_f=got[0], trim_r=got[1]), cfg=SimpleNamespace(max_ee=got[2] if got else 2.0),
                               floor=floors.get(name, 500), short_lost=short.get(name, 0.0), sub=SimpleNamespace(r1=name, r2=name))

    def loader(r1, r2):
        name = r1
        best = floor_only.get(name, released[name][:2] if released[name] else None)
        return SimpleNamespace(pass_table=lambda ee: SimpleNamespace(best_split=lambda total: None if best is None else (best[0], best[1], 0.9)),
                               error_load=lambda tf, tr, ee: (noise or {}).get(name))
    return run, loader


def check(tmp_path, released, **kw):
    nc = load()
    n3, n4 = trc.write_nights(tmp_path)
    cfg = trc.write_config(tmp_path)
    run, loader = stand_ins(released, **kw)
    out = tmp_path / "out"
    code = nc.main(["--config", str(cfg), "--night3", str(n3), "--night4", str(n4), "--out", str(out)], pipeline_run=run, load_pairs=loader)
    return nc, code, out


def test_the_rule_is_as_tested_everywhere_except_on_mock_b_where_the_cap_is_gone(tmp_path):
    nc = load()
    n3, n4 = trc.write_nights(tmp_path)
    by = lambda name: nc.tested_settings(n3, n4, name)                                                                      # noqa: E731
    rows = [nc.inspect(k, SimpleNamespace(recommendation=SimpleNamespace(trim_f=v[0], trim_r=v[1]), cfg=SimpleNamespace(max_ee=v[2]), floor=1, short_lost=0, sub=SimpleNamespace(r1=k, r2=k)),
                       lambda r1, r2: SimpleNamespace(pass_table=lambda ee: SimpleNamespace(best_split=lambda t: (1, 1, 1)))) for k, v in released_settings().items() if k in nc.DATASETS]
    verdict, status = nc.p39(rows, by)
    assert verdict is True and status["mockB"].startswith("expected change") and status["v4"] == "same" and status["wastewater_V3-V4"] == "same"
    unchanged = [dict(r, released=(277, 250, (4.0, 5.0))) if r["dataset"] == "mockB" else r for r in rows]                                   # the cap still there: mock B is as tested
    v, st = nc.p39(unchanged, by)
    assert v is False and "is the cap still in this package?" in st["mockB"]                                                              # the prediction was a change, so this is not met
    other = [dict(r, released=(r["released"][0], r["released"][1], (3.0, 3.0))) if r["dataset"] == "wastewater_V3-V4" else r for r in rows]  # a dataset that should not have changed did
    assert nc.p39(other, by)[0] is False and nc.p39(other, by)[1]["wastewater_V3-V4"] == "**NO**"
    wrong = [dict(r, released=(r["released"][0], r["released"][1], (4.0, 6.0))) if r["dataset"] == "mockB" else r for r in rows]            # changed, but not to the predicted (4,8)
    v, st = nc.p39(wrong, by)
    assert v is False and "predicted maxEE 4,8" in st["mockB"]
    assert nc.p39(rows[:11], by)[0] is None                                                                                 # not every dataset: not decided
    assert nc.p39(rows, lambda name: None)[0] is None                                                                       # no tested results found: not decided


def test_the_floor_only_cut_is_counted_against_the_released_cut(tmp_path):
    nc = load()
    mk = lambda name, rel, fo: {"dataset": name, "released": rel, "floor_only": fo, "short_lost": 0.0, "floor": 1, "noise": None}          # noqa: E731
    rows = [mk(k, (100, 100, (2.0, 2.0)), (100, 100, (2.0, 2.0))) for k in nc.DATASETS if k != "its"] + [mk("its", (0, 0, (2.0, 2.0)), None)]
    assert nc.p43(rows) == (True, 11, 11)
    for i in range(3):                                                                                                      # three differ: 8 of 11 is still enough
        rows[i] = mk(rows[i]["dataset"], (100, 100, (2.0, 2.0)), (101, 99, (2.0, 2.0)))
    assert nc.p43(rows) == (True, 8, 11)
    rows[3] = mk(rows[3]["dataset"], (100, 100, (2.0, 2.0)), (101, 99, (2.0, 2.0)))                                          # a fourth: 7 of 11
    assert nc.p43(rows) == (False, 7, 11)
    assert nc.p43(rows[:5])[0] is None                                                                                      # not every dataset: not decided


def test_the_table_of_the_mock_b_run_is_read_with_the_margins_fixed_for_the_battery(tmp_path):
    nc = load()
    def table(merged_h, merged_n, ident_h=0.738, ident_n=0.76, found_n=6):
        rows = [{"name": "AmpliTrim (maxEE 2)", "truncF": "279", "truncR": "248", "maxEE": "2", "merged": "0.199", "found": "6", "strains": "3/3", "identical": "0.687"},
                {"name": "AmpliTrim (maxEE 4,5)", "truncF": "277", "truncR": "250", "maxEE": "4,5", "merged": str(merged_h), "found": "6", "strains": "3/3", "identical": str(ident_h)},
                {"name": "AmpliTrim (maxEE 4,8)", "truncF": "270", "truncR": "256", "maxEE": "4,8", "merged": str(merged_n), "found": str(found_n), "strains": "3/3", "identical": str(ident_n)}]
        path = tmp_path / "compare_methods.tsv"
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) + ["error"], delimiter="\t"); w.writeheader(); [w.writerow(dict(r, error="")) for r in rows]
        return nc.evaluate_table(nc.read_table(path))
    good = table(0.550, 0.690)
    assert good["p40"] and good["p40_in_range"] and good["p41"] and good["accuracy"] == "more accurate" and abs(good["gain"] - 14.0) < 1e-6
    assert not table(0.550, 0.560)["p40"]                                                                                    # +1 point is inside the 2-point margin
    assert table(0.550, 0.800)["p40"] and not table(0.550, 0.800)["p40_in_range"]                                           # a gain, but outside the range predicted
    assert not table(0.550, 0.690, ident_n=0.700)["p41"] and table(0.550, 0.690, ident_n=0.700)["accuracy"] == "less accurate"          # accuracy lost beyond 0.1 points
    assert not table(0.550, 0.690, found_n=4)["p41"]                                                                        # fewer variants than the margin of 1 allows
    assert table(0.550, 0.690, ident_n=0.7385)["accuracy"] == "as accurate"                                                  # within the 0.1 point margin
    txt = "\n".join(nc.table_lines(good))
    assert "P40** MET" in txt and "P41** MET" in txt and "+14.0 points" in txt
    assert "run compare_methods.py with --amplitrim-max-ee 2 2,5 4,5 4,8" in nc.evaluate_table({})["error"]
    assert nc.PASSING_ACCURACY == ("as accurate", "more accurate")


def test_the_check_runs_end_to_end_and_prints_the_dada2_commands(tmp_path):
    differ = {"v1v2": (194, 169)}                                                                                          # the floor-only cut differs on one dataset
    nc, code, out = check(tmp_path, released_settings(), floor_only=differ, short={"wastewater_V3-V4": 0.35, "v4": 0.004}, noise={"mockB": (1.02, 3.31)}, floors={"mockB": 526})
    assert code == 0
    text = (out / "NOCAP_CHECK_REPORT.md").read_text()
    assert "P39** MET" in text and "P43** MET" in text and "10 of 11" in text                                                  # (its has no cut)
    assert "expected change: the cap is gone" in text
    assert "35.0%  (the package warns from 1%)" in text and "| v4 | 0.4% |" in text                                           # length loss is reported, and flagged from 1%
    assert "--amplitrim-max-ee 2 2,5 4,5 4,8 --pairs 5000" in text and "--reference" in text and "night5_mockB" in text      # the mock B run
    assert "--setting floor_only=194,169,2" in text and "night5_v1v2" in text                                                  # the run for the dataset where the floor-only cut differs
    assert "--setting released=" in text and "--baseline released" in text                                                     # the released cut is a row, and the baseline
    assert "take the commands from that report" in text                                                                         # on the package's own sample the cuts are not the nights' pairs
    assert "1.02 forward, 3.31 reverse" in text                                                                                  # P42, descriptive
    tsv = list(csv.reader(open(out / "nocap_check.tsv", encoding="utf-8"), delimiter="\t"))
    assert tsv[0][:4] == ["dataset", "released_truncF", "released_truncR", "released_maxEE"] and len(tsv) == 13
    row = next(r for r in tsv if r[0] == "mockB")
    assert row[1:4] == ["270", "256", "4,8"] and row[4] == "526"


def test_when_the_floor_only_cut_is_always_the_released_cut_nothing_is_left_to_run(tmp_path):
    nc, code, out = check(tmp_path, released_settings())
    text = (out / "NOCAP_CHECK_REPORT.md").read_text()
    assert "P44: the floor-only cut is the released cut on every dataset" in text and "--setting floor_only" not in text


def test_a_rule_that_changed_where_it_should_not_have_is_reported_as_not_met(tmp_path):
    nc, code, out = check(tmp_path, released_settings(**{"wastewater_V3-V4": (261, 183, (5.0, 5.0))}))
    text = (out / "NOCAP_CHECK_REPORT.md").read_text()
    assert "P39** NOT MET" in text and "**NO**" in text


def test_some_datasets_only_leaves_the_verdicts_undecided_and_a_package_with_the_cap_is_refused(tmp_path, capsys):
    nc = load()
    n3, n4 = trc.write_nights(tmp_path)
    cfg = trc.write_config(tmp_path)
    run, loader = stand_ins(released_settings())
    assert nc.main(["--config", str(cfg), "--night3", str(n3), "--night4", str(n4), "--only", "mockB", "--out", str(tmp_path / "o")], pipeline_run=run, load_pairs=loader) == 0
    assert "not decided" in (tmp_path / "o" / "NOCAP_CHECK_REPORT.md").read_text()
    assert nc.main(["--config", str(cfg), "--out", str(tmp_path / "o2"), "--only", "zzz"], pipeline_run=run, load_pairs=loader) == 2            # nothing was run


def test_a_missing_config_is_reported_by_name_and_a_dataset_without_reads_says_where_it_looked(tmp_path, capsys):
    nc = load()
    run, loader = stand_ins(released_settings())
    missing = tmp_path / "nowhere" / "amplitrim_paths.json"
    assert nc.main(["--config", str(missing), "--out", str(tmp_path / "o")], pipeline_run=run, load_pairs=loader) == 2
    said = capsys.readouterr().out
    assert "No amplitrim_paths.json was found at" in said and str(missing) in said and "where /r C:" in said      # it names the file and says how to find the real one
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"extras": {}}), encoding="utf-8")                                                    # a config that has no entry for mock B
    assert nc.main(["--config", str(cfg), "--only", "mockB", "--out", str(tmp_path / "o2")], pipeline_run=run, load_pairs=loader) == 2
    said = capsys.readouterr().out
    assert "mockB: no reads folder for it in" in said and "'extras' > 'mockB' > 'ready'" in said and "Using the datasets in" in said


def test_like_the_nights_gives_the_package_the_thinned_pairs_and_says_so(tmp_path):
    import night2
    nc = load()
    assert nc.PAIRS_PER_SAMPLE == night2.PAIRS and nc.DEFAULT_PAIRS_PER_SAMPLE == night2.DEFAULT_PAIRS          # the same thinning as the nights
    n3, n4 = trc.write_nights(tmp_path)
    cfg = trc.write_config(tmp_path)
    run, loader = stand_ins(released_settings())
    calls = []
    def thin(folder, out_dir, per):
        calls.append((Path(folder).name, per)); return folder                                                     # the stand-in package reads the dataset from the folder name
    out = tmp_path / "o"
    assert nc.main(["--config", str(cfg), "--night3", str(n3), "--night4", str(n4), "--like-the-nights", "--out", str(out)], pipeline_run=run, load_pairs=loader, thin_reads=thin) == 0
    thinned = dict((name, per) for name, per in calls)
    assert "mockA_high" not in thinned and set(thinned.values()) == {5000} and len(thinned) == 11                    # mock A is given whole; every other dataset 5,000 pairs per sample
    text = (out / "NOCAP_CHECK_REPORT.md").read_text()
    assert "this is P39b" in text and "thinned to 5,000 pairs per sample" in text
    calls.clear()
    out2 = tmp_path / "o2"
    assert nc.main(["--config", str(cfg), "--night3", str(n3), "--night4", str(n4), "--out", str(out2)], pipeline_run=run, load_pairs=loader, thin_reads=thin) == 0
    assert calls == [] and "this is P39." in (out2 / "NOCAP_CHECK_REPORT.md").read_text()                            # without the option nothing is thinned


def test_dada2_is_given_what_it_got_in_the_nights_all_the_reads_of_mock_a_and_5000_pairs_per_sample_elsewhere(tmp_path):
    nc, code, out = check(tmp_path, released_settings(), floor_only={"mockA_high": (274, 246), "wastewater_V3-V4": (254, 186)})
    lines = [l for l in (out / "NOCAP_CHECK_REPORT.md").read_text().splitlines() if l.startswith("python validation\\compare_methods.py")]
    mock_a = next(l for l in lines if "night5_mockA_high" in l)
    waste = next(l for l in lines if "night5_wastewater_V3-V4" in l)
    mock_b = next(l for l in lines if "night5_mockB" in l)
    assert "--pairs" not in mock_a and "--pairs 5000" in waste and "--pairs 5000" in mock_b                                   # mock A is not thinned
    assert "--setting released=273,247,2,3 --setting floor_only=274,246,2,3 --baseline released" in mock_a                  # released first, floor-only second, released as the baseline
