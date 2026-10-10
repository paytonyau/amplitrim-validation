"""validation/maxee_frontier.py: where FIGARO's maxEE sits against what a rule could pick, from the exact quality table (no DADA2, no vsearch)."""
import csv
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))


def load():
    spec = importlib.util.spec_from_file_location("maxee_frontier", HERE.parent / "validation" / "maxee_frontier.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _point(mf, ee, pass_rate, ef, er):
    return {"ee": ee, "pass": pass_rate, "ee_f": ef, "ee_r": er, "score": mf.fig_score(pass_rate, ee)}


def test_the_rules_pick_what_they_say_on_a_made_up_grid():
    mf = load()
    grid = sorted([_point(mf, (1, 1), 0.06, 0.29, 0.66), _point(mf, (2, 2), 0.20, 0.46, 1.24), _point(mf, (2, 5), 0.49, 0.63, 2.42), _point(mf, (3, 5), 0.54, 0.77, 2.49),
                   _point(mf, (5, 5), 0.57, 0.93, 2.54), _point(mf, (6, 9), 0.74, 1.30, 3.60)], key=lambda x: mf.looseness(x["ee"]))
    assert mf.strictest_reaching(grid, 0.5)["ee"] == (3, 5)                      # (2,5) keeps 49%: just short
    assert mf.strictest_reaching(grid, 0.7)["ee"] == (6, 9) and mf.strictest_reaching(grid, 0.9) is None
    assert mf.quality_first(grid)["ee"] == (3, 5)                                # the ladder rule: the strictest of the ladder that keeps half
    assert mf.best_by_score(grid)["ee"] == (3, 5)                                # retention minus (F-1)^2 minus (R-1)^2: 34 at (3,5), 32 at (2,5), 25 at (5,5)
    assert mf.fig_score(0.74, (6, 9)) == 74 - 25 - 64 < 0                        # a loose maxEE scores low; FIGARO ranks cuts, each with the maxEE it derives from the reads, so a poor-quality cut can win with a negative score


def test_kept_reads_are_called_cleaner_noisier_or_mixed():
    mf = load()
    fig = _point(mf, (3, 4), 0.76, 0.60, 1.80)
    assert mf.cleaner_or_noisier(_point(mf, (3, 3), 0.7, 0.55, 1.50), fig) == "cleaner"
    assert mf.cleaner_or_noisier(_point(mf, (5, 5), 0.8, 0.90, 2.50), fig) == "noisier"
    assert mf.cleaner_or_noisier(_point(mf, (2, 6), 0.8, 0.40, 2.50), fig) == "mixed"            # cleaner forward, noisier reverse
    assert mf.cleaner_or_noisier(_point(mf, (3, 4), 0.76, 0.61, 1.79), fig) == "about the same"


def test_the_comparison_table_is_read_without_errors_or_untruncated_rows(tmp_path):
    mf = load()
    with open(tmp_path / "compare_methods.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "truncF", "truncR", "maxEE", "error"], delimiter="\t")
        w.writeheader()
        w.writerows([{"name": "AmpliTrim (maxEE 2)", "truncF": "257", "truncR": "260", "maxEE": "2,6"}, {"name": "FIGARO (amplicon 505)", "truncF": "273", "truncR": "244", "maxEE": "3,4"},
                     {"name": "full", "truncF": "0", "truncR": "0", "maxEE": "2"}, {"name": "broken", "truncF": "", "truncR": "", "maxEE": "2", "error": "x"}])
    rows = mf.read_methods(tmp_path)
    assert [(r["name"], r["tf"], r["tr"], r["ee"]) for r in rows] == [("AmpliTrim (maxEE 2)", 257, 260, (2.0, 6.0)), ("FIGARO (amplicon 505)", 273, 244, (3.0, 4.0))]


def test_the_analysis_runs_on_simulated_reads_and_writes_its_table(tmp_path):
    import simulate
    mf = load()
    simulate.write_dataset(tmp_path / "sim", n_samples=2, pairs_per_sample=1500)
    reads = tmp_path / "sim" / "primer_trimmed"
    (tmp_path / "cmp").mkdir()
    with open(tmp_path / "cmp" / "compare_methods.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "truncF", "truncR", "maxEE"], delimiter="\t")
        w.writeheader()
        w.writerows([{"name": "AmpliTrim (maxEE 2)", "truncF": 215, "truncR": 215, "maxEE": "2"}, {"name": "FIGARO (amplicon 440)", "truncF": 220, "truncR": 212, "maxEE": "3,5"}])
    assert mf.main(["--reads", str(reads), "--compare", str(tmp_path / "cmp"), "--max-pairs", "3000", "--name", "sim"]) == 0
    rows = list(csv.DictReader(open(tmp_path / "cmp" / "maxee_frontier.tsv", encoding="utf-8"), delimiter="\t"))
    kinds = [r["what"] for r in rows]
    assert kinds.count("setting") == 2 and kinds.count("rule") == 6 and "reaches FIGARO's retention" in kinds
    fig = next(r for r in rows if r["name"].startswith("FIGARO"))
    assert 0.5 < float(fig["pass"]) <= 1.0 and fig["maxEE"] == "3,5"
    default = next(r for r in rows if r["name"] == "default (2,2)")
    assert default["maxEE"] == "2,2" and float(default["pass"]) < float(fig["pass"])             # a looser maxEE keeps more pairs
