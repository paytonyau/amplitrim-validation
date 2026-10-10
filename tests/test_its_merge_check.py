"""validation/its_merge_check.py: the exact flow of the ITS pairs through DADA2, and what an overhang explains. R is a stand-in that writes a table."""
import csv
import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("its_merge_check", HERE.parent / "validation" / "its_merge_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# two samples; reach = merged + rejected in every row; sample b recovers 80 more pairs with trimOverhang
ROWS = [("a", 1000, 850, 600, 400, 60, 40, 100, 480), ("b", 1000, 750, 500, 300, 50, 30, 120, 380)]
HEADER = ["sample", "input", "filtered", "reach", "merged", "rej_short", "rej_mismatch", "rej_other", "merged_trim_overhang"]


def write_table(path, rows=ROWS, header=HEADER):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def write_night3(root, merged_by_sample):
    from compare_methods import method_key
    folder = root / "its" / "compare"
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / f"track_{method_key('AmpliTrim (maxEE 2)')}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "input", "filtered", "denoised_f", "denoised_r", "merged", "nonchim"])
        for s, m in zip("ab", merged_by_sample):
            w.writerow([s, 1000, 800, 600, 600, m, m])
    return root


def test_the_table_is_read_and_a_missing_column_is_refused(tmp_path):
    im = load()
    write_table(tmp_path / "t.tsv")
    rows = im.read_table(tmp_path / "t.tsv")
    assert rows[0]["input"] == 1000.0 and rows[1]["sample"] == "b"
    write_table(tmp_path / "bad.tsv", rows=[("a", 1, 1, 1, 1, 0, 0, 0)], header=HEADER[:-1])
    with pytest.raises(RuntimeError, match="merged_trim_overhang"):
        im.read_table(tmp_path / "bad.tsv")


def test_the_flow_is_exact_and_adds_up(tmp_path):
    im = load()
    write_table(tmp_path / "t.tsv")
    f = im.flow(im.read_table(tmp_path / "t.tsv"))
    assert f["input"] == 2000 and f["samples"] == 2
    assert f["filtered"] == pytest.approx(0.80) and f["reach"] == pytest.approx(0.55) and f["merged"] == pytest.approx(0.35)
    assert f["lost_denoising"] == pytest.approx(0.25) and f["lost_merging"] == pytest.approx(0.20)               # (1600 - 1100) / 2000 and (1100 - 700) / 2000
    assert f["rej_short"] == pytest.approx(0.055) and f["rej_mismatch"] == pytest.approx(0.035) and f["rej_other"] == pytest.approx(0.11)
    assert f["rej_short"] + f["rej_mismatch"] + f["rej_other"] == pytest.approx(f["lost_merging"])               # every rejected pair has exactly one reason
    assert f["lost_denoising"] + f["lost_merging"] == pytest.approx(f["filtered"] - f["merged"])
    assert f["merged_trim_overhang"] == pytest.approx(0.43) and f["overhang_gain"] == pytest.approx(0.08)         # (860 - 700) / 2000
    assert im.flow([]) == {}


def test_the_verdicts_follow_the_registered_criteria():
    im = load()
    base = {"merged": 0.35, "lost_denoising": 0.25, "lost_merging": 0.20, "overhang_gain": 0.08}
    v = {p: ok for p, _, ok in im.verdicts(base, 0.351)}
    assert v == {"P26": True, "P27": True, "P28": True, "P29": True}
    v = {p: ok for p, _, ok in im.verdicts(dict(base, lost_denoising=0.05, overhang_gain=0.0199, lost_merging=0.065), 0.356)}
    assert v == {"P26": False, "P27": False, "P28": False, "P29": False}                  # 0.6 points off night 3; denoising smaller than merging; 1.99 points is under 2; 6.5 - 1.89 = 4.61 is under 5
    v = {p: ok for p, _, ok in im.verdicts(dict(base, overhang_gain=0.02, lost_merging=0.070, lost_denoising=0.30), 0.3529)}
    assert v["P28"] is True and v["P29"] is True and v["P26"] is True and v["P27"] is True   # exactly 2 points is enough; 7.0 - 1.89 = 5.11 is over 5; 0.29 points is within 0.3
    assert {p: ok for p, _, ok in im.verdicts(base, None)}["P26"] is None                                         # no night 3 track: not decided
    assert all(ok is None for _, _, ok in im.verdicts({}, 0.35))


def test_the_check_runs_end_to_end_and_writes_its_report(tmp_path):
    im = load()
    cfg = tmp_path / "paths.json"
    reads = tmp_path / "reads"
    reads.mkdir()
    cfg.write_text('{"its": "%s"}' % str(reads).replace("\\", "/"))
    night3 = write_night3(tmp_path / "night3", [350, 350])
    seen = {}

    def thin(folder, out_dir, per_sample):
        seen["thin"] = (Path(folder), per_sample)
        return out_dir

    def run(cmd, **kw):
        seen["cmd"] = [str(c) for c in cmd]
        write_table(Path(cmd[-1]))
        return SimpleNamespace(returncode=0, stdout="", stderr="")
    rc = im.main(["--config", str(cfg), "--night3", str(night3), "--out", str(tmp_path / "o")], run=run, find_rscript=lambda: "Rscript", thin=thin)
    from night2 import DEFAULT_PAIRS, PAIRS
    assert rc == 0 and seen["thin"] == (reads, PAIRS.get("its", DEFAULT_PAIRS))                                  # the same cap per sample as night 3
    assert seen["cmd"][0] == "Rscript" and seen["cmd"][1].endswith("its_merge_check.R") and seen["cmd"][3:6] == ["50", "2", "12"] and seen["cmd"][6].endswith("its_merge.tsv")
    rep = (tmp_path / "o" / "ITS_MERGE_REPORT.md").read_text()
    for p in ("P26", "P27", "P28", "P29"):
        assert f"- **{p}** MET" in rep, p
    assert "20.0 points to merging" in rep or "20.0 points to merging (rejected pairs)" in rep
    assert "25.0 points to denoising" in rep and "+8.00 points" in rep and "rejected: a clean overlap that is still refused (an overhang) | 11.00%" in rep
    assert "Night 3 merged 35.00% on the same reads." in rep


def test_missing_things_are_said_plainly(tmp_path, capsys):
    im = load()
    cfg = tmp_path / "paths.json"
    (tmp_path / "reads").mkdir()
    cfg.write_text('{"its": "%s"}' % str(tmp_path / "reads").replace("\\", "/"))
    assert im.main(["--config", str(tmp_path / "absent.json"), "--out", str(tmp_path / "o")]) == 2 and "No reads folder" in capsys.readouterr().out
    assert im.main(["--config", str(cfg), "--out", str(tmp_path / "o")], find_rscript=lambda: None) == 2 and "Rscript was not found" in capsys.readouterr().out
    boom = lambda cmd, **kw: SimpleNamespace(returncode=1, stdout="", stderr="Error in library(dada2): there is no package called 'dada2'")      # noqa: E731
    assert im.main(["--config", str(cfg), "--night3", str(tmp_path / "none"), "--out", str(tmp_path / "o")], run=boom, find_rscript=lambda: "Rscript", thin=lambda f, o, n: o) == 2
    assert "DADA2 could not run" in capsys.readouterr().out
