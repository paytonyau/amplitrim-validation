"""The helper scripts in validation/ fail plainly when R is missing."""
import importlib.util
import sys
from pathlib import Path

import pytest

VALIDATION = Path(__file__).parent.parent / "validation"


def load(name):
    sys.path.insert(0, str(VALIDATION))
    spec = importlib.util.spec_from_file_location(name, VALIDATION / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def no_rscript(*args, **kwargs):
    raise FileNotFoundError(2, "The system cannot find the file specified")


def test_example_data_says_what_to_install_when_rscript_is_missing(monkeypatch, capsys, tmp_path):
    module = load("example_data")
    monkeypatch.setattr(module.subprocess, "run", no_rscript)
    monkeypatch.setattr(module.sys, "argv", ["example_data.py", str(tmp_path / "out")])
    assert module.main() == 1
    assert "Rscript was not found" in capsys.readouterr().out and not (tmp_path / "out").exists()


def test_the_dada2_step_of_the_runner_says_what_to_install_when_rscript_is_missing(monkeypatch, tmp_path):
    module = load("real_data")
    monkeypatch.setattr(module.subprocess, "run", no_rscript)
    with pytest.raises(SystemExit, match="Rscript was not found"):
        module.dada2_run(tmp_path, 100, 100, 2.0, 12)


def test_thinning_keeps_every_sample_and_the_pairs_in_step(tmp_path):
    import gzip
    rd = load("real_data")
    raw = tmp_path / "raw"
    raw.mkdir()
    for sample, n in (("a", 100), ("b", 40)):
        for mate in ("R1", "R2"):
            with gzip.open(raw / f"{sample}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                for i in range(n):
                    fh.write(f"@{sample}{i}/{mate}\nACGT\n+\nIIII\n")
    out = rd.thin_reads(raw, tmp_path / "thin", 20)
    for sample, expected in (("a", 20), ("b", 20)):
        names = {m: [l.strip() for i, l in enumerate(gzip.open(out / f"{sample}_{m}.fastq.gz", "rt")) if i % 4 == 0] for m in ("R1", "R2")}
        assert len(names["R1"]) == expected and [x.split("/")[0] for x in names["R1"]] == [x.split("/")[0] for x in names["R2"]]


def test_primers_and_read_through_are_both_removed_from_the_same_read(tmp_path):
    """cutadapt removes one adapter per read by default; reads that run through the amplicon have two to remove."""
    import gzip
    import random
    import subprocess
    if subprocess.run([sys.executable, "-m", "cutadapt", "--version"], capture_output=True).returncode != 0:
        pytest.skip("cutadapt is not installed")
    rd = load("real_data")
    rng = random.Random(2)
    comp = str.maketrans("ACGT", "TGCA")
    rc = lambda s: s.translate(comp)[::-1]
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    fwd, rev = "ACTGGGCTTAAAGAGC"[:15], "TACGTGGGTATCTAATCC"
    raw = tmp_path / "raw"
    raw.mkdir()
    inserts = [junk(rng.randint(90, 120)) for _ in range(60)]
    with gzip.open(raw / "s_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(raw / "s_R2.fastq.gz", "wt", newline="\n") as o2:
        for i, ins in enumerate(inserts):
            r1 = (fwd + ins + rc(rev) + junk(40))[:170]                 # runs through the end into the other primer
            r2 = (rev + rc(ins) + rc(fwd) + junk(40))[:170]
            o1.write(f"@r{i}\n{r1}\n+\n{'I' * len(r1)}\n")
            o2.write(f"@r{i}\n{r2}\n+\n{'I' * len(r2)}\n")
    out = tmp_path / "out"
    rd.trim_reads(raw, out, fwd, rev)
    got1 = [l.strip() for i, l in enumerate(gzip.open(out / "s_R1.fastq.gz", "rt")) if i % 4 == 1]
    got2 = [l.strip() for i, l in enumerate(gzip.open(out / "s_R2.fastq.gz", "rt")) if i % 4 == 1]
    assert got1 == inserts and got2 == [rc(i) for i in inserts]


def test_the_error_rate_decides_whether_read_through_with_two_errors_is_removed(tmp_path):
    import gzip
    import random
    import subprocess
    if subprocess.run([sys.executable, "-m", "cutadapt", "--version"], capture_output=True).returncode != 0:
        pytest.skip("cutadapt is not installed")
    rd = load("real_data")
    rng = random.Random(5)
    comp = str.maketrans("ACGT", "TGCA")
    rc = lambda s: s.translate(comp)[::-1]
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    fwd, rev = "ACTGGGCTTAAAGAG", "TACGTGGGTATCTAATCC"
    two_errors = lambda s: "".join(c if i not in (4, 11) else ("A" if c != "A" else "C") for i, c in enumerate(s))
    raw = tmp_path / "raw"
    raw.mkdir()
    inserts = [junk(100) for _ in range(40)]
    with gzip.open(raw / "s_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(raw / "s_R2.fastq.gz", "wt", newline="\n") as o2:
        for i, ins in enumerate(inserts):
            r1 = (fwd + ins + two_errors(rc(rev)) + junk(40))[:170]     # read-through whose primer has two errors
            r2 = (rev + rc(ins) + rc(fwd) + junk(40))[:170]
            o1.write(f"@r{i}\n{r1}\n+\n{'I' * len(r1)}\n")
            o2.write(f"@r{i}\n{r2}\n+\n{'I' * len(r2)}\n")
    lengths = lambda folder: {len(l.strip()) for i, l in enumerate(gzip.open(folder / "s_R1.fastq.gz", "rt")) if i % 4 == 1}
    rd.trim_reads(raw, tmp_path / "default", fwd, rev)
    rd.trim_reads(raw, tmp_path / "tolerant", fwd, rev, error_rate=0.15)
    assert lengths(tmp_path / "tolerant") == {100}                                   # all trimmed back to the insert
    assert lengths(tmp_path / "default") != {100}                                    # one error allowed: the rest is left behind


def test_the_runner_prints_the_dada2_numbers_it_produced(tmp_path):
    rd = load("real_data")
    row = {"dataset": "study", "truncLenF": 123, "truncLenR": 115, "estimated_yield": 0.9759, "dada2_merged": 0.9612,
           "dada2_nonchim": 0.9301, "dada2_asvs": 412, "compare_trunc": "200,200", "compare_dada2_merged": 0.9105, "compare_asvs": 398}
    lines = rd.summary_lines(row, tmp_path / "results.tsv")
    assert lines[0] == "study: AmpliTrim 123/115, estimated yield 0.9759"
    assert lines[1] == "  DADA2 at 123/115: 96.1% of read pairs merged (93.0% after chimera removal), 412 ASVs"
    assert lines[2] == "  DADA2 at 200/200: 91.0% of read pairs merged, 398 ASVs" or lines[2].startswith("  DADA2 at 200/200: 91.")
    assert lines[-1].startswith("  row added to")
    none = rd.summary_lines({"dataset": "x", "truncLenF": None}, tmp_path / "r.tsv")
    assert none[0] == "x: AmpliTrim made no recommendation." and len(none) == 2


def test_mock_check_classifies_asvs_against_the_reference(tmp_path):
    import random
    import shutil
    if not shutil.which("vsearch"):
        pytest.skip("vsearch is not installed")
    mc = load("mock_check")
    rng = random.Random(8)
    refs = {f"strain{i}": "".join(rng.choice("ACGT") for _ in range(300)) for i in range(1, 5)}
    (tmp_path / "ref.fasta").write_text("".join(f">{n}\n{s}\n" for n, s in refs.items()))
    one_off = refs["strain2"][50:290]
    one_off = one_off[:100] + ("A" if one_off[100] != "A" else "C") + one_off[101:]
    asvs = [(">ASV1;size=800", refs["strain1"][40:290]), (">ASV2;size=150", refs["strain2"][50:290]), (">ASV3;size=10", one_off),
            (">ASV4;size=5", "".join(rng.choice("ACGT") for _ in range(240)))]
    (tmp_path / "asvs.fasta").write_text("".join(f"{h}\n{s}\n" for h, s in asvs))
    s = mc.summarise(tmp_path / "asvs.fasta", tmp_path / "ref.fasta")
    assert s["categories"]["exact"] == [2, 950] and s["categories"]["near"] == [1, 10] and s["categories"]["none"] == [1, 5]
    assert s["recovered"] == ["strain1", "strain2"] and s["missing"] == ["strain3", "strain4"] and s["split"] == {"strain2": 2}


def test_placeholders_and_flat_quality_are_stopped_early_with_a_plain_message():
    rd = load("real_data")
    rd.check_sequence("--primer-f", "GTGYCAGCMGCCGCGGTAA")
    rd.check_sequence("--primer-f", None)
    with pytest.raises(SystemExit) as e:
        rd.check_sequence("--primer-f", "FWD")
    assert "'FWD' is not a DNA sequence" in str(e.value) and "not a placeholder" in str(e.value)
    assert rd.dada2_blocker(2).startswith("the reads have only 2 distinct quality values") and rd.dada2_blocker(3) is None and rd.dada2_blocker(32) is None
    lines = rd.summary_lines({"dataset": "x", "truncLenF": 120, "truncLenR": 118, "estimated_yield": 0.98, "dada2_skipped": rd.dada2_blocker(2)}, Path("r.tsv"))
    assert any(l.startswith("  DADA2 not run: the reads have only 2") for l in lines)


def test_a_dada2_failure_is_explained_and_does_not_abort_the_run(tmp_path):
    import gzip
    import random
    import shutil
    import subprocess
    if not shutil.which("Rscript") or subprocess.run(["Rscript", "-e", "library(dada2)"], capture_output=True).returncode != 0:
        pytest.skip("R with DADA2 is not installed")
    rd = load("real_data")
    rng = random.Random(6)
    for sample in ("a", "b"):
        for mate in ("R1", "R2"):
            with gzip.open(tmp_path / f"{sample}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                for i in range(300):
                    s = "".join(rng.choice("ACGT") for _ in range(120))
                    fh.write(f"@{sample}{i}\n{s}\n+\n{''.join('I' if rng.random() < 0.985 else rng.choice('5+') for _ in s)}\n")     # mostly Q40, three values in all
    with pytest.raises(rd.Dada2Error) as e:
        rd.dada2_run(tmp_path, 500, 500, 2.0, 12)                      # longer than every read: the filter removes them all
    assert str(e.value).startswith("DADA2 failed:")
    flat = rd.explain_dada2_failure("Error in getErrors(err, enforce = TRUE) : Error matrix is NULL.\nCalls: learnErrors")
    assert "could not learn an error model" in flat and "AmpliTrim's own result above is unaffected" in flat
    assert rd.explain_dada2_failure("boom").startswith("DADA2 failed: boom")


def test_the_python_runner_prints_help_and_still_writes_an_outcome_file_when_it_fails(tmp_path):
    import gzip
    import subprocess
    script = Path(__file__).resolve().parents[1] / "validation" / "run_dataset.py"
    help_text = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True).stdout
    assert "run_dataset.py NAME RAW_DIR" in help_text and "--reference" in help_text
    raw = tmp_path / "raw"
    raw.mkdir()
    for mate in ("R1", "R2"):
        with gzip.open(raw / f"s_{mate}.fastq.gz", "wt", newline="\n") as fh:
            fh.write("@r\nACGTACGTAC\n+\nIIIIIIIIII\n")
    res = subprocess.run([sys.executable, str(script), "x", str(raw), "--primer-f", "FWD", "--out", str(tmp_path / "o")], capture_output=True, text=True)
    assert res.returncode == 1
    text = (tmp_path / "o" / "x_outcome.txt").read_text(encoding="utf-8")
    assert "is not a DNA sequence" in text and "No AmpliTrim result was written" in text


def test_the_readmes_refer_only_to_scripts_that_exist_and_not_to_removed_files():
    import re
    root = Path(__file__).resolve().parents[1]
    for readme in (root / "README.md", root / "validation" / "README.md"):
        text = readme.read_text(encoding="utf-8")
        for m in re.findall(r"validation[/\\]((?:simulated[/\\])?\w+\.(?:py|R))", text):
            assert (root / "validation" / m.replace("\\", "/")).exists(), f"{readme.name} mentions {m}, which does not exist"
        assert "TESTING_NEW_DATASETS" not in text and "run_dataset.sh" not in text


def test_samples_with_too_few_pairs_are_left_out_of_the_working_copy(tmp_path):
    import gzip
    rd = load("real_data")
    for sample, n in (("big", 12), ("small", 3), ("control", 1)):
        for mate in ("R1", "R2"):
            with gzip.open(tmp_path / f"{sample}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                for i in range(n):
                    fh.write(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    dropped = rd.drop_small_samples(tmp_path, 5)
    assert sorted(dropped) == [("control", 1), ("small", 3)]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["big_R1.fastq.gz", "big_R2.fastq.gz"]
    assert rd.drop_small_samples(tmp_path, 5) == []


def test_the_runner_reports_a_stable_build_id_that_changes_with_the_code(tmp_path):
    import subprocess
    script = Path(__file__).resolve().parents[1] / "validation" / "run_dataset.py"
    run = lambda: subprocess.run([sys.executable, str(script), "--version"], capture_output=True, text=True).stdout.strip()
    first = run()
    assert first.startswith("AmpliTrim 0.1.0 build ") and len(first.split()[-1]) == 10 and first == run()


def test_runner_max_ee_values_are_validated_and_normalised():
    import argparse
    rd = load("real_data")
    assert rd.ee_text("2") == "2" and rd.ee_text("2,6") == "2,6" and rd.ee_text("3,3") == "3"
    with pytest.raises(argparse.ArgumentTypeError):
        rd.ee_text("two")


def test_the_maxEE_sweep_formats_its_rows_and_explains_itself():
    import subprocess
    sw = load("max_ee_sweep")
    d = {"filtered": 0.918, "merged": 0.917, "nonchim": 0.9171, "asvs": 15}
    score = {"categories": {"exact": [15, 950], "near": [0, 0], "none": [1, 50]}, "reads": 1000, "recovered": ["a"] * 14, "references": 15}
    row = sw.format_row("2,5", d, score)
    assert "91.8%" in row and "95.00%" in row and "5.00%" in row and "14 of 15" in row and row.rstrip().endswith("1")
    assert "identical" in sw.header(True) and "identical" not in sw.header(False)
    script = Path(__file__).resolve().parents[1] / "validation" / "max_ee_sweep.py"
    assert "maxEE" in subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True).stdout


def test_the_runner_shows_the_quality_first_suggestion_when_it_was_asked_for():
    rd = load("real_data")
    row = {"dataset": "x", "truncLenF": 257, "truncLenR": 260, "estimated_yield": 0.775, "suggested_max_ee": "2,4",
           "suggested_trunc": "258/259", "suggested_yield": 0.91, "suggested_reached": True, "min_retention": 0.5}
    lines = rd.summary_lines(row, Path("r.tsv"))
    assert any("maxEE, quality first (keeping at least 50% of pairs): 2,4 with truncLen 258/259, estimated yield 91.0%; yours 77.5%" in l for l in lines)


def test_the_runners_accept_the_coverage_options():
    import subprocess
    for name in ("run_dataset.py", "real_data.py"):
        script = Path(__file__).resolve().parents[1] / "validation" / name
        text = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True).stdout
        assert "--explore-coverage" in text and "--allow-unmergeable" in text and "--only" in text


def test_only_keeps_the_named_samples_and_refuses_a_name_that_does_not_exist(tmp_path):
    import gzip
    rd = load("real_data")
    def make(folder):
        folder.mkdir()
        for name in ("A", "B", "C"):
            for mate in ("R1", "R2"):
                with gzip.open(folder / f"{name}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                    fh.write("@r\nACGT\n+\nIIII\n")
        return folder
    one = make(tmp_path / "one")
    assert sorted(rd.keep_only(one, "A,C")) == ["B"] and sorted(p.name for p in one.iterdir()) == ["A_R1.fastq.gz", "A_R2.fastq.gz", "C_R1.fastq.gz", "C_R2.fastq.gz"]
    two = make(tmp_path / "two")
    names = tmp_path / "names.txt"
    names.write_text("# the high-input samples\nB\n\nC\n", encoding="utf-8")
    assert rd.keep_only(two, str(names)) == ["A"]
    three = make(tmp_path / "three")
    with pytest.raises(SystemExit) as e:
        rd.keep_only(three, "A,Z")
    assert "not in the reads: Z" in str(e.value) and len(list(three.iterdir())) == 6        # nothing was deleted


def test_the_sweep_summarises_abundance_without_a_reference(tmp_path):
    sw = load("max_ee_sweep")
    sizes = [5000, 4000, 3000, 2000, 1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 1, 1, 1]
    (tmp_path / "a.fasta").write_text("".join(f">ASV{i};size={n}\nACGT\n" for i, n in enumerate(sizes)))
    s = sw.abundance_summary(tmp_path / "a.fasta")
    total = sum(sizes)
    assert abs(s["top15"] - sum(sizes[:15]) / total) < 1e-12 and s["n_01"] == sum(1 for n in sizes if n / total >= 0.001)
    assert abs(s["tail"] - sum(n for n in sizes if n / total < 0.001) / total) < 1e-12
    row = sw.format_row("2", {"filtered": 0.9, "merged": 0.9, "nonchim": 0.9, "asvs": 17}, None, s)
    assert f"{s['top15']:.2%}" in row and "top 15" in sw.header(False)


def test_a_setting_that_merges_nothing_is_scored_as_finding_nothing_instead_of_failing(tmp_path):
    mc, sw = load("mock_check"), load("max_ee_sweep")
    ref = tmp_path / "ref.fasta"
    ref.write_text(">a\nACGTACGTACGTACGTACGTACGT\n>b\nTTGACCAGTTGACCAGTTGACCAG\n")
    missing = tmp_path / "no_asvs_were_written.fasta"                       # DADA2 writes no file when no pair merged
    s = mc.summarise(missing, ref)
    assert s["asvs"] == 0 and s["recovered"] == [] and s["missing"] == ["a", "b"] and s["categories"]["exact"] == (0, 0)
    assert sw.abundance_summary(missing) == {"top15": 0.0, "n_01": 0, "tail": 0.0}


def test_compare_methods_parses_settings_and_refuses_a_malformed_one():
    cm = load("compare_methods")
    assert cm.parse_setting("authors=279,263,2,6") == ("authors", 279, 263, "2,6") and cm.parse_setting("full=280,283,2") == ("full", 280, 283, "2")
    assert cm.parse_setting("full=0,0,2") == ("full", 0, 0, "2")                  # 0 means no truncation, as in DADA2
    for bad in ("279,263,2", "x=279,263", "x=279,263,2,6,7", "x=a,263,2", "=279,263,2", "x=279,263,0", "x=-1,263,2", "full=LENF,LENR,2"):
        with pytest.raises(SystemExit, match="NAME=truncLenF"):
            cm.parse_setting(bad)


def test_the_paired_difference_is_exact_for_a_constant_gap_and_brackets_a_noisy_one():
    cm = load("compare_methods")
    mean, lo, hi, up, down = cm.paired_difference([0.80, 0.70, 0.90, 0.60, 0.75, 0.85], [0.78, 0.68, 0.88, 0.58, 0.73, 0.83])
    assert abs(mean - 0.02) < 1e-12 and abs(lo - 0.02) < 1e-12 and abs(hi - 0.02) < 1e-12 and (up, down) == (6, 0)
    mean, lo, hi, up, down = cm.paired_difference([0.9, 0.5, 0.8, 0.4, 0.7, 0.6], [0.6, 0.8, 0.5, 0.7, 0.4, 0.9])
    assert lo <= mean <= hi and abs(mean) < 1e-12 and (up, down) == (3, 3)
    assert cm.paired_difference([0.9, 0.5, 0.8], [0.6, 0.8, 0.5]) == cm.paired_difference([0.9, 0.5, 0.8], [0.6, 0.8, 0.5])      # seeded: repeatable


def test_figaro_is_given_only_pairs_of_the_most_common_lengths_and_the_share_is_reported(tmp_path):
    import gzip
    cm = load("compare_methods")
    reads = tmp_path / "reads"
    reads.mkdir()
    for mate, length in (("R1", 50), ("R2", 60)):
        with gzip.open(reads / f"s_{mate}.fastq.gz", "wt", newline="\n") as fh:
            for i in range(100):
                n = length - 5 if i >= 80 else length                         # 20 pairs have shorter reads
                fh.write(f"@p{i}\n{'A' * n}\n+\n{'I' * n}\n")
    folder, share, lf, lr = cm.figaro_input(reads, tmp_path / "f", max_pairs=0)
    assert (lf, lr) == (50, 60) and abs(share - 0.8) < 1e-12
    kept = [r for r in cm._read_fastq_gz(folder / "subsample_R1.fastq.gz")]
    assert len(kept) == 80 and all(len(r[1].strip()) == 50 for r in kept)
    assert sum(1 for _ in cm._read_fastq_gz(folder / "subsample_R2.fastq.gz")) == 80


def test_the_comparison_table_shows_reference_columns_only_with_a_reference_and_explains_a_failure():
    cm = load("compare_methods")
    row = {"name": "AmpliTrim (maxEE 2)", "truncF": 277, "truncR": 237, "maxEE": "2", "own": 0.566, "own_what": "yield", "est_yield": 0.568, "pass": 0.572,
           "merged": 0.559, "nonchim": 0.559, "asvs": 104, "top15": 0.5115, "n01": 95, "tail": 0.0063, "found": 104, "identical": 1.0, "extra": 0}
    failed = {"name": "x", "truncF": 240, "truncR": 160, "maxEE": "2", "error": "DADA2 could not run: nothing"}
    plain, withref = cm.format_table([row, failed], False), cm.format_table([row, failed], True)
    assert "identical" not in plain[0] and "identical" in withref[0] and "100.00%" in withref[1] and "55.9%" in plain[1]
    assert "DADA2 could not run: nothing" in plain[2]


def test_the_per_sample_merged_share_is_read_from_the_dada2_table(tmp_path):
    cm = load("compare_methods")
    track = tmp_path / "t.tsv"
    track.write_text("sample\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\nA\t1000\t900\t880\t880\t800\t780\nB\t0\t0\t0\t0\t0\t0\nC\t500\t450\t440\t440\t400\t390\n", encoding="utf-8")
    assert cm.sample_merged_share(track) == {"A": 0.8, "C": 0.8}                  # a sample with no input is left out, not divided by zero


def test_compare_methods_refuses_a_bad_argument_before_doing_any_work(tmp_path, monkeypatch):
    cm = load("compare_methods")
    reads = tmp_path / "reads"
    reads.mkdir()
    def refused(*extra):
        monkeypatch.setattr(sys, "argv", ["compare_methods.py", str(reads), "--out", str(tmp_path / "o"), *extra])
        with pytest.raises(SystemExit) as e:
            cm.main()
        assert not (tmp_path / "o").exists()                                     # nothing was started
        return str(e.value)
    assert "NAME=truncLenF" in refused("--setting", "full=LENF,LENR,2")
    assert "no such file" in refused("--figaro-script", str(tmp_path / "C:figaro" / "figaro.py"))
    assert "--reference" in refused("--reference", str(tmp_path / "missing.fasta"))
    monkeypatch.setattr(sys, "argv", ["compare_methods.py", str(tmp_path / "nothere"), "--out", str(tmp_path / "o")])
    with pytest.raises(SystemExit, match="not a folder of reads"):
        cm.main()


def test_the_table_shows_no_truncation_as_none_and_leaves_the_estimate_blank():
    cm = load("compare_methods")
    row = {"name": "full", "truncF": 0, "truncR": 0, "maxEE": "2", "own": None, "own_what": "", "est_yield": None, "pass": 0.9, "merged": 0.8,
           "nonchim": 0.75, "asvs": 10, "top15": 0.99, "n01": 9, "tail": 0.001}
    line = cm.format_table([row], False)[1]
    assert " none " in line and "80.0%" in line and "None" not in line


def _calibration_rows():
    base = [("AmpliTrim (maxEE 2)", 0.906, 0.914, 0.842, 0.915), ("FIGARO (amplicon 343)", 0.879, 0.884, 0.822, 0.885),
            ("AmpliTrim (FIGARO's maxEE 2,1)", 0.884, 0.889, 0.823, 0.891), ("more_overlap", 0.848, 0.848, 0.771, 0.848), ("much_more", 0.546, 0.548, 0.483, 0.548)]
    return [{"name": n, "est_yield": y, "est_pass": ep, "merged": m, "pass": p, "own": 0.881 if n.startswith("FIGARO") else None, "est_strict": None}
            for n, y, ep, m, p in base]


def test_the_calibration_lines_give_the_offset_and_the_rank_agreement_and_say_what_figaros_number_measures():
    cm = load("compare_methods")
    lines = cm.calibration_lines(_calibration_rows())
    assert "Across 5 settings" in lines[0] and "yield against merged share +6.4 points (range +5.7 to +7.7)" in lines[0] and "Rank correlation of the estimated yield with DADA2's merged share: +1.00" in lines[0] and "As a fraction of the estimated yield, DADA2's merged share is 0.92 (range 0.88 to 0.94)" in lines[0]
    assert "FIGARO's predicted retention minus DADA2's pass rate: -0.4 points" in lines[1] and "filter only" in lines[1]
    assert cm.calibration_lines(_calibration_rows()[:2]) == [lines[1]]                 # two settings are too few for an offset or a rank; FIGARO's line remains
    reversed_rows = _calibration_rows()                                                  # estimates fall in the order of rows 0, 2, 1, 3, 4
    for i, m in zip((0, 2, 1, 3, 4), (0.5, 0.6, 0.7, 0.8, 0.9)):
        reversed_rows[i]["merged"] = m                                                   # so DADA2's shares rise in that same order: exactly opposite
    assert "-1.00" in cm.calibration_lines(reversed_rows)[0]                            # a ranking that is wrong is reported as wrong


def test_the_baseline_defaults_to_the_first_method_and_can_be_named_by_a_unique_piece_of_its_name():
    cm = load("compare_methods")
    names = ["AmpliTrim (maxEE 2)", "FIGARO (amplicon 343)", "more_overlap"]
    assert cm.pick_baseline(None, names) == names[0] and cm.pick_baseline("figaro", names) == names[1] and cm.pick_baseline("more_overlap", names) == names[2]
    tricky = names + ["AmpliTrim (FIGARO's maxEE 2,1)"]
    assert cm.pick_baseline("FIGARO", tricky) == "FIGARO (amplicon 343)"                  # a name that starts with the text wins over one that contains it
    assert cm.pick_baseline("maxEE 2,1", tricky) == "AmpliTrim (FIGARO's maxEE 2,1)"      # a piece of one name only
    with pytest.raises(SystemExit, match="more than one"):
        cm.pick_baseline("mp", names)                                                  # inside two names ("AmpliTrim", "amplicon") and the start of none
    with pytest.raises(SystemExit, match="no method"):
        cm.pick_baseline("zzz", names)


def test_the_pairwise_comparison_is_redone_from_saved_tables_against_any_baseline(tmp_path):
    cm = load("compare_methods")
    out = tmp_path / "run"
    out.mkdir()
    names = ["AmpliTrim (maxEE 2)", "FIGARO (amplicon 343)", "more_overlap"]
    (out / "compare_methods.tsv").write_text("name\tmerged\n" + "".join(f"{n}\t0.8\n" for n in names), encoding="utf-8")
    level = {"AmpliTrim (maxEE 2)": 0.84, "FIGARO (amplicon 343)": 0.82, "more_overlap": 0.77}
    for n in names:
        rows = "".join(f"S{i}\t1000\t900\t880\t880\t{int(1000 * (level[n] + 0.001 * (i % 3)))}\t0\n" for i in range(8))
        (out / f"track_{cm.method_key(n)}.tsv").write_text("sample\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n" + rows, encoding="utf-8")
    text = "\n".join(cm.pairwise_report(out, "FIGARO"))
    assert "against FIGARO (amplicon 343)" in text
    assert "AmpliTrim (maxEE 2): +2.0 points (95% interval +2.0 to +2.0); better in 8 of 8 samples, worse in 0" in text
    assert "more_overlap: -5.0 points" in text and "worse in 8" in text
    assert "AmpliTrim (maxEE 2): +0.0" not in "\n".join(cm.pairwise_report(out, "more_overlap"))
    with pytest.raises(SystemExit, match="was not found"):
        cm.pairwise_report(tmp_path / "nothing")


def test_figaro_is_refused_for_its_before_any_work_starts(tmp_path, monkeypatch):
    cm = load("compare_methods")
    reads = tmp_path / "reads"
    reads.mkdir()
    monkeypatch.setattr(sys, "argv", ["compare_methods.py", str(reads), "--out", str(tmp_path / "o"), "--marker", "its", "--figaro-script", str(tmp_path / "figaro.py")])
    with pytest.raises(SystemExit, match="cannot be used with --marker its"):
        cm.main()
    assert not (tmp_path / "o").exists()


def test_the_strict_estimate_is_left_out_of_the_calibration_when_a_row_has_none():
    cm = load("compare_methods")
    rows = [{"name": f"AmpliTrim (maxEE {i})", "est_yield": y, "est_pass": y, "merged": m, "pass": y, "own": y, "est_strict": None}
            for i, (y, m) in enumerate(((0.60, 0.596), (0.81, 0.808), (0.90, 0.895), (0.96, 0.959)), start=1)]
    text = "\n".join(cm.calibration_lines(rows))
    assert "Rank correlation of the estimated yield with DADA2's merged share: +1.00" in text and "strict" not in text


def test_peek_reads_reports_the_commonest_starts_lengths_and_examples(tmp_path):
    import gzip
    pk = load("peek_reads")
    folder = tmp_path / "f"
    folder.mkdir()
    for mate, start in (("R1", "GAACCTGCGG"), ("R2", "GCTGCGTTCT")):
        with gzip.open(folder / f"S1_{mate}.fastq.gz", "wt", newline="\n") as fh:
            for i in range(100):
                s = (start if i < 70 else "ACGTACGTAC") + "T" * (40 if i % 2 else 30)
                fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    text = "\n".join(pk.peek(folder, "S1", scan=100, n_examples=2, k=10))
    assert "R1 starts (first 10 bases): GAACCTGCGG 70.0% | ACGTACGTAC 30.0%" in text and "R2 starts (first 10 bases): GCTGCGTTCT 70.0%" in text
    assert "min 40, median 50, max 50; most common 40 nt (50%)" in text and text.count("    R1 ") == 2
    assert "not found" in "\n".join(pk.peek(folder, "NOPE"))


def test_compare_methods_can_use_only_some_samples_and_refuses_a_name_that_is_not_there(tmp_path, monkeypatch):
    import gzip
    cm = load("compare_methods")
    reads = tmp_path / "reads"
    reads.mkdir()
    for name in ("A", "B", "C"):
        for mate in ("R1", "R2"):
            with gzip.open(reads / f"{name}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                fh.write("@r\nACGT\n+\nIIII\n")
    out = cm.subset_reads(reads, tmp_path / "sub", "A,C")
    assert sorted(p.name for p in out.iterdir()) == ["A_R1.fastq.gz", "A_R2.fastq.gz", "C_R1.fastq.gz", "C_R2.fastq.gz"]
    names = tmp_path / "names.txt"
    names.write_text("# a few samples\nB\n", encoding="utf-8")
    assert sorted(p.name for p in cm.subset_reads(reads, tmp_path / "sub2", str(names)).iterdir()) == ["B_R1.fastq.gz", "B_R2.fastq.gz"]
    monkeypatch.setattr(sys, "argv", ["compare_methods.py", str(reads), "--out", str(tmp_path / "o"), "--only", "A,Z"])
    with pytest.raises(SystemExit, match="not in .*: Z"):
        cm.main()
    assert not (tmp_path / "o").exists()


def _cutadapt_or_skip():
    import subprocess
    if subprocess.run([sys.executable, "-m", "cutadapt", "--version"], capture_output=True).returncode != 0:
        pytest.skip("cutadapt is not installed")


def _pair_through_trimming(tmp_path, r1, r2, primer_f="GAACCWGCGGARGGATCA", primer_r="GCTGCGTTCTTCATCGATGC"):
    import gzip
    rd = load("real_data")
    for name, s in (("i1", r1), ("i2", r2)):
        with gzip.open(tmp_path / f"{name}.fastq.gz", "wt", newline="\n") as fh:
            for i in range(5):
                fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    rd._cutadapt(tmp_path / "i1.fastq.gz", tmp_path / "i2.fastq.gz", tmp_path / "o1.fastq.gz", tmp_path / "o2.fastq.gz", primer_f, primer_r, None, None, 0.15)
    try:
        with gzip.open(tmp_path / "o1.fastq.gz", "rt") as a, gzip.open(tmp_path / "o2.fastq.gz", "rt") as b:
            return a.readlines()[1].strip(), b.readlines()[1].strip()
    except (IndexError, OSError):
        return None


def test_chance_matches_at_the_read_ends_cannot_pass_a_pair_that_has_no_primers_in_the_right_reads(tmp_path):
    _cutadapt_or_skip()
    body = "ACGTTGCAAGTCCATGCAAGGTCCTAAGTTCCGGAATTCCGGTTAACCGGTTCAAGTCCAGTCAGTCATAGGAT" * 2
    # the wrong way round (R1 starts with the REVERSE primer), and each read happens to end with the first bases of the 3' adapter that is applied to it
    assert _pair_through_trimming(tmp_path, "GCTGCGTTCTTCATCGATGC" + body + "GCAT", "GAACCTGCGGAGGGATCA" + body + "TGAT") is None


def test_read_through_is_still_removed_down_to_a_residue_of_three_bases(tmp_path):
    _cutadapt_or_skip()
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    rc = lambda s: "".join(comp[c] for c in s[::-1])
    body = "ACGTTGCAAGTCCATGCAAGGTCCTAAGTTCCGGAATTCCGGTTAACCGGTTCAAGTCCAGTCAGTCATAGGAT" * 2
    f18, rev = "GAACCTGCGGAGGGATCA", "GCTGCGTTCTTCATCGATGC"
    assert _pair_through_trimming(tmp_path, f18 + body + rc(rev), rev + rc(body) + rc(f18)) == (body, rc(body))       # full read-through
    assert _pair_through_trimming(tmp_path, f18 + body + "GCA", rev + rc(body) + "TGA") == (body, rc(body))              # a 3-base residue
    assert _pair_through_trimming(tmp_path, f18[1:] + body, rev + rc(body), "GAACCWGCGGARGGATCA") == (body, rc(body))      # a primer missing its first base is still found


def test_rscript_is_found_on_path_else_in_the_newest_standard_install_folder(tmp_path, monkeypatch):
    rt = load("rtools")
    for k in ("R_HOME", "ProgramFiles", "ProgramW6432", "ProgramFiles(x86)", "LOCALAPPDATA"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setattr(rt.shutil, "which", lambda name: None)
    assert rt.find_rscript() is None                                                  # nothing anywhere
    pf = tmp_path / "Program Files"
    for v in ("R-4.3.2", "R-4.5.1", "R-4.10.0"):
        exe = pf / "R" / v / "bin" / "x64" / "Rscript.exe"
        exe.parent.mkdir(parents=True)
        exe.write_text("")
    monkeypatch.setenv("ProgramFiles", str(pf))
    assert rt.find_rscript() == str(pf / "R" / "R-4.10.0" / "bin" / "x64" / "Rscript.exe")        # the newest by number (4.10 is newer than 4.5)
    home = tmp_path / "myR"
    (home / "bin").mkdir(parents=True)
    (home / "bin" / "Rscript.exe").write_text("")
    monkeypatch.setenv("R_HOME", str(home))
    assert rt.find_rscript() in (str(home / "bin" / "Rscript.exe"), str(pf / "R" / "R-4.10.0" / "bin" / "x64" / "Rscript.exe"))   # either is a real R
    monkeypatch.setattr(rt.shutil, "which", lambda name: "/usr/bin/Rscript")
    assert rt.find_rscript() == "/usr/bin/Rscript"                                    # PATH wins when R is on it


def test_the_dada2_step_uses_the_r_that_was_found_and_says_so_when_there_is_none(tmp_path, monkeypatch):
    rd = load("real_data")
    seen = {}
    class Done:
        returncode, stderr, stdout = 1, "boom", ""
    def fake_run(cmd, **kw):
        seen["cmd"] = cmd
        return Done()
    monkeypatch.setattr(rd, "find_rscript", lambda: "C:/Program Files/R/R-4.5.1/bin/x64/Rscript.exe")
    monkeypatch.setattr(rd.subprocess, "run", fake_run)
    with pytest.raises(rd.Dada2Error):
        rd.dada2_run(tmp_path, 200, 200, "2", 12)
    assert seen["cmd"][0] == "C:/Program Files/R/R-4.5.1/bin/x64/Rscript.exe"          # the found R, not a bare "Rscript"
    monkeypatch.setattr(rd, "find_rscript", lambda: None)
    with pytest.raises(SystemExit, match="neither on PATH nor in the usual install folders"):
        rd.dada2_run(tmp_path, 200, 200, "2", 12)


def _variant_toy(tmp_path):
    import random
    rng = random.Random(21)
    a = "".join(rng.choice("ACGT") for _ in range(200))
    b = "".join(rng.choice("ACGT") for _ in range(200))
    def sub(seq, positions):
        out = list(seq)
        for p in positions:
            out[p] = "A" if out[p] != "A" else "C"
        return "".join(out)
    ref = tmp_path / "ref.fasta"                                   # strain S1 has two variants one base apart, S2 has one sequence
    ref.write_text(f">v1|S1|A\n{a}\n>v2|S1|A\n{sub(a, [100])}\n>v3|S2|B\n{b}\n")
    asvs = tmp_path / "asvs.fasta"
    asvs.write_text(f">a1;size=100\n{a}\n>a2;size=10\n{b}\n>a3;size=5\n{sub(a, [50])}\n>a4;size=3\n{sub(a, [10, 40, 80, 120, 160])}\n")
    return ref, asvs, a, sub


def test_variants_are_scored_at_variant_and_strain_level_and_near_means_a_few_bases(tmp_path):
    import shutil
    if not shutil.which("vsearch"):
        pytest.skip("vsearch is not installed")
    mc = load("mock_check")
    ref, asvs, a, sub = _variant_toy(tmp_path)
    s = mc.summarise(asvs, ref)                                    # near = identity of at least 0.99 (about 2 bases in 200)
    c = s["categories"]
    assert (c["exact"][0], c["near"][0], c["none"][0]) == (2, 1, 1)                        # a1, a2 identical; a3 one base off; a4 five bases off
    assert s["recovered"] == ["v1|S1|A", "v3|S2|B"] and s["references"] == 3               # two of the three variants
    assert (s["strains_found"], s["strains_total"]) == (2, 2)                              # but both strains are detected
    only_one = tmp_path / "one.fasta"
    only_one.write_text(f">a1;size=100\n{a}\n")
    s1 = mc.summarise(only_one, ref)
    assert (s1["strains_found"], s1["strains_total"]) == (1, 2) and len(s1["recovered"]) == 1
    assert mc.summarise(asvs, ref, min_id=0.97)["categories"]["near"][0] == 2             # the wider 3% rule also calls a4 near: variants and errors blur
    flat = tmp_path / "flat.fasta"
    flat.write_text(f">x\n{a}\n>y\n{sub(a, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])}\n")
    assert mc.strain_counts(["x", "y"], ["x"]) == (None, None)                            # one sequence per strain: no strain line


def test_the_strain_column_appears_only_when_the_reference_has_variants():
    sw, cm = load("max_ee_sweep"), load("compare_methods")
    d = {"filtered": 0.9, "merged": 0.8, "nonchim": 0.8, "asvs": 60}
    score = {"categories": {"exact": (50, 900), "near": (3, 20), "none": (7, 80)}, "reads": 1000, "recovered": ["r"] * 49, "references": 53, "strains_found": 14, "strains_total": 14}
    assert " 14/14" in sw.format_row("2", d, score) and "strains" in sw.header(True, False, True) and "strains" not in sw.header(True, False, False)
    assert " 14/14" not in sw.format_row("2", d, {**score, "strains_total": None})
    row = {"name": "m", "truncF": 257, "truncR": 260, "maxEE": "2", "own": 0.8, "own_what": "yield", "est_yield": 0.8, "est_pass": 0.9, "pass": 0.9, "merged": 0.8, "nonchim": 0.8,
           "asvs": 60, "top15": 0.9, "n01": 40, "tail": 0.01, "found": 49, "identical": 0.97, "extra": 10, "strains": "14/14"}
    table = cm.format_table([row], True)
    assert "strains" in table[0] and table[1].rstrip().endswith("14/14")
    assert "strains" not in cm.format_table([{**row, "strains": None}], True)[0]


def test_the_mock_a_reference_is_built_from_regions_with_wildcards_and_close_asvs():
    pytest.importorskip("rapidfuzz")
    pytest.importorskip("numpy")
    b = load("build_mockA_reference")
    region1, region2 = "ACGTACGTTGCAAGTCCATGCA", "TTGACCGGTTAAGCTCGATTCA"
    asvs = {
        "ASV_1": "GGGG" + region1 + "CCCCAAAA",                                    # contains region 1 exactly
        "ASV_2": "TT" + region1[:10] + "A" + region1[11:] + "GGTT",               # one mismatch from region 1
        "ASV_3": "CCCCAAAA" + "GATTACAGATTACAGATTACAGATTACAGATTACA" + "TTTT",     # unrelated
        "ASV_4": "TT" + region2[:5] + "T" + region2[6:] + "CCGG",                 # one mismatch from region 2 (region 2 has no exact ASV)
        "ASV_5": "GGGG" + region1 + "CCCCAAAT",                                   # one edit from ASV_1: found only by edit distance
    }
    assert b.best_mismatches(asvs["ASV_1"], region1) == 0 and b.best_mismatches(asvs["ASV_2"], region1) == 1 and b.best_mismatches("ACG", region1) is None
    wild = region1[:4] + "Y" + region1[5:]                                         # an ambiguity code in the region allows C or T at that position
    assert b.best_mismatches("TT" + region1[:4] + "C" + region1[5:] + "GG", wild) == 0 and b.best_mismatches("TT" + region1[:4] + "G" + region1[5:] + "GG", wild) == 1
    hits = b.match_asvs(asvs, {"S1": region1, "S2": region2}, tolerance=0)
    assert set(hits["S1"]) == {"ASV_1", "ASV_5"} and hits["S2"] == {}
    fasta, rows, unmatched = b.build(asvs, {"S1": region1, "S2": region2}, {"S1": "Strain one", "S2": "Strain two"}, tolerance=1, max_edits=4)
    names = [r[0] for r in rows]
    assert names == ["ASV_1", "ASV_2", "ASV_4", "ASV_5"] and "ASV_3" not in fasta and unmatched == []
    assert ">ASV_4|S2|Strain_two" in fasta and ">ASV_1|S1|Strain_one" in fasta                   # variants carry the strain in the second field
    assert b.build(asvs, {"S1": region1, "S9": "A" * 22}, {}, tolerance=0)[2] == ["S9"]          # a strain with no ASV is reported


def test_the_reference_can_be_limited_to_the_variants_present_in_chosen_samples():
    pd = pytest.importorskip("pandas")
    b = load("build_mockA_reference")
    counts = pd.DataFrame({"103-KitA-1-1-8": [100, 0, 5, 0], "104-KitA-1-1-9": [50, 0, 0, 0], "105-KitA-1-2-2": [0, 70, 0, 9]}, index=["ASV_1", "ASV_2", "ASV_3", "ASV_4"])
    ids = ["ASV_1", "ASV_2", "ASV_3", "ASV_4"]
    assert b.present_in(counts, ids, r"^(103|104)-KitA") == ["ASV_1", "ASV_3"]                  # ASV_2 and ASV_4 exist only in another dilution
    assert b.present_in(counts, ids, r"^(103|104)-KitA", min_reads=10) == ["ASV_1"]
    assert b.present_in(counts, ids, r"-KitA-1-2-2$") == ["ASV_2", "ASV_4"]
    with pytest.raises(SystemExit, match="matches no sample"):
        b.present_in(counts, ids, "^nothing")


def test_mock_a_runs_are_grouped_by_dilution_and_the_reference_follows_the_presence_in_each():
    pd = pytest.importorskip("pandas")
    md = load("mockA_dilutions")
    design = {"ERR1": ("103-KitA-1-1-8_S103_L001", 1000), "ERR2": ("104-KitA-1-1-9_S104_L001", 1000), "ERR3": ("111-KitA-1-2-8_S111_L001", 1000),
              "ERR4": ("705-pipeline-NC", 50), "ERR5": ("113-KitA-1-3-8_S113_L001", 1000)}
    groups = md.group_samples(design, {"ERR1", "ERR2", "ERR3", "ERR4"})                    # ERR5 is not among the reads; ERR4 is a control
    assert groups == {"8": ["ERR1", "ERR3"], "9": ["ERR2"]}
    counts = pd.DataFrame({"103-KitA-1-1-8": [10, 0, 0], "104-KitA-1-1-9": [0, 0, 7], "111-KitA-1-2-8": [5, 0, 0], "105-KitA-1-2-2": [1, 9, 1]}, index=["ASV_1", "ASV_2", "ASV_3"])
    ref = [("ASV_1|S1|A", "AAAA"), ("ASV_2|S2|B", "CCCC"), ("ASV_3|S3|C", "GGGG")]
    assert [n for n, _ in md.restricted_reference(ref, counts, "8")] == ["ASV_1|S1|A"]
    assert [n for n, _ in md.restricted_reference(ref, counts, "9")] == ["ASV_3|S3|C"]
    assert [n for n, _ in md.restricted_reference(ref, counts, "2")] == ["ASV_1|S1|A", "ASV_2|S2|B", "ASV_3|S3|C"]
    row = {"dilution": "8", "samples": 3, "pass": 0.81, "merged": 0.78, "asvs": 50, "present": 39, "found": 35, "strains_present": 13, "strains_found": 11,
           "identical": 0.9976, "near": 0.0003, "none": 0.0021, "outside": 14}
    table = md.format_table([row, {"dilution": "2", "samples": 1, "error": "DADA2 could not run: nothing"}])
    assert " 1e8 " in table[1] and "35/39" in table[1] and "11/13" in table[1] and "99.76%" in table[1] and "DADA2 could not run" in table[2]


def test_the_dilution_series_runs_end_to_end_with_a_stand_in_for_dada2(tmp_path, monkeypatch, capsys):
    import gzip
    import random
    import shutil
    pytest.importorskip("pandas")
    if not shutil.which("vsearch"):
        pytest.skip("vsearch is not installed")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
    import real_data as rd
    md = load("mockA_dilutions")
    rng = random.Random(5)
    seqs = {f"ASV_{i}": "".join(rng.choice("ACGT") for _ in range(240)) for i in range(1, 5)}
    reference = tmp_path / "ref.fasta"                                         # four variants of two strains
    reference.write_text("".join(f">ASV_{i}|S{1 + (i - 1) // 2}|Strain_{1 + (i - 1) // 2}\n{seqs[f'ASV_{i}']}\n" for i in range(1, 5)))
    table = tmp_path / "asvs.csv"                                               # ASV_1, ASV_2 are in the 1e9 sample; ASV_2, ASV_3 in the 1e8 ones; ASV_4 only at 1e2
    table.write_text("OTU_ID;103-KitA-1-1-8;104-KitA-1-1-9;105-KitA-1-1-2\nASV_1;0;50;0\nASV_2;20;30;0\nASV_3;40;0;0\nASV_4;0;0;9\n")
    runinfo = tmp_path / "runs.tsv"
    runinfo.write_text("run_accession\tsample_alias\tread_count\nERR1\t103-KitA-1-1-8_S103_L001\t2000\nERR2\t104-KitA-1-1-9_S104_L001\t2000\n", encoding="utf-8")
    reads = tmp_path / "reads"
    reads.mkdir()
    for run in ("ERR1", "ERR2"):
        for mate in ("R1", "R2"):
            with gzip.open(reads / f"{run}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                fh.write("@r\nACGT\n+\nIIII\n")
    calls = []
    def fake_dada2(reads_dir, tf, tr, ee, min_overlap, min_len, track, asv):
        calls.append(sorted(p.name for p in Path(reads_dir).iterdir()))
        found = ["ASV_2", "ASV_3"] if calls[-1][0].startswith("ERR1") else ["ASV_1"]            # at 1e8 DADA2 finds both present variants; at 1e9 it misses ASV_2
        lines = [f">x{k};size=100\n{seqs[a]}\n" for k, a in enumerate(found)]
        lines.append(">junk;size=2\n" + "".join(rng.choice("ACGT") for _ in range(240)) + "\n")
        Path(asv).write_text("".join(lines))
        return {"filtered": 0.8, "merged": 0.7, "nonchim": 0.7, "asvs": len(found) + 1}
    monkeypatch.setattr(rd, "dada2_run", fake_dada2)
    monkeypatch.setattr(sys, "argv", ["mockA_dilutions.py", str(reads), "--runinfo", str(runinfo), "--asv-table", str(table), "--reference", str(reference),
                                      "--trunc", "257,260", "--out", str(tmp_path / "out")])
    assert md.main() == 0
    out = capsys.readouterr().out
    assert calls == [["ERR2_R1.fastq.gz", "ERR2_R2.fastq.gz"], ["ERR1_R1.fastq.gz", "ERR1_R2.fastq.gz"]]       # each dilution gets only its own samples, highest input first
    lines = [ln for ln in out.splitlines() if ln.lstrip().startswith("1e")]
    assert len(lines) >= 2
    nine = next(ln for ln in lines if ln.split()[0] == "1e9"); eight = next(ln for ln in lines if ln.split()[0] == "1e8")
    assert "1/2" in nine and "1/1" in nine                                    # present: ASV_1, ASV_2 (two variants of S1); found ASV_1; the strain is found
    assert eight.count("2/2") == 2 and eight.split("|")[-1].strip() == "1"   # present: ASV_2 (strain 1) and ASV_3 (strain 2), both found; one ASV outside the reference (the junk)
    assert "present but not recovered" in (tmp_path / "out" / "mockA_dilutions.txt").read_text(encoding="utf-8")


# ------------------------------------------------------------------------------------------------ overnight.py
def _overnight():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
    return load("overnight")


def test_overnight_helpers_clip_pick_and_parse_real_outputs(tmp_path):
    import gzip
    ov = _overnight()
    assert ov.pick_even(["d", "a", "c", "b"], 2) == ["a", "c"] and ov.pick_even(["a", "b"], 5) == ["a", "b"]
    with gzip.open(tmp_path / "s_R1.fastq.gz", "wt", newline="\n") as f1, gzip.open(tmp_path / "s_R2.fastq.gz", "wt", newline="\n") as f2:
        for i in range(3):
            f1.write(f"@r{i}\n{'A' * 300}\n+\n{'I' * 300}\n")
            f2.write(f"@r{i}\n{'C' * 280}\n+\n{'I' * 280}\n")
    assert ov.clip_reads(tmp_path, tmp_path / "clip", 250) == 3
    with gzip.open(tmp_path / "clip" / "s_R1.fastq.gz", "rt") as f:
        lines = f.read().splitlines()
    with gzip.open(tmp_path / "clip" / "s_R2.fastq.gz", "rt") as f:
        lines2 = f.read().splitlines()
    assert len(lines[1]) == 250 and len(lines[3]) == 250 and len(lines2[1]) == 250 and ov.longest_read(tmp_path) == 300
    # the real output of mockA_dilutions.py
    text = """   1e9       3  82.8%   78.5%    44 |   35/38    11/13   |    99.84%  0.09%  0.06% |                          9
   1e5       3  79.9%   78.9%    49 |   20/21     8/9    |    24.32%  0.00% 75.68% |                         29
1e9: present but not recovered: ASV_049|Enterococcus_faecalis, ASV_078|Staphylococcus_epidermidis, ASV_282|Lactobacillus_reuteri
1e5: present but not recovered: ASV_049|Enterococcus_faecalis
"""
    rows, misses = ov.parse_dilution_text(text)
    assert [(r["d"], r["found"], r["present"], r["sf"], r["sp"]) for r in rows] == [("9", "35", "38", "11", "13"), ("5", "20", "21", "8", "9")]
    assert misses["9"] == ["ASV_049|Enterococcus_faecalis", "ASV_078|Staphylococcus_epidermidis", "ASV_282|Lactobacillus_reuteri"] and misses["5"] == ["ASV_049|Enterococcus_faecalis"]
    pairs = ov.parse_pairwise("  AmpliTrim (maxEE 2): +11.4 points (95% interval +10.7 to +12.2); better in 8 of 8 samples, worse in 0\n  trunc200: -23.5 points (95% interval -25.2 to -21.8); better in 0 of 37 samples, worse in 37\n")
    assert pairs == [{"name": "AmpliTrim (maxEE 2)", "mean": 11.4, "lo": 10.7, "hi": 12.2}, {"name": "trunc200", "mean": -23.5, "lo": -25.2, "hi": -21.8}]
    a = {"variants": 35, "strains_found": 11, "identical": 0.9976, "extra": 14}
    assert ov.noninferior(a, {"variants": 36, "strains_found": 11, "identical": 0.9979, "extra": 14})[0]
    assert not ov.noninferior({**a, "variants": 33}, {"variants": 36, "strains_found": 11, "identical": 0.9979, "extra": 14})[0]
    assert not ov.noninferior({**a, "identical": 0.99}, {"variants": 35, "strains_found": 11, "identical": 0.9979, "extra": 14})[0]
    assert "AMPLITRIM_DADA2_POOL" in (HERE_VALIDATION / "dada2_real.R").read_text() and "AMPLITRIM_DADA2_MINFOLD" in (HERE_VALIDATION / "dada2_real.R").read_text()


HERE_VALIDATION = Path(__file__).resolve().parents[1] / "validation"


def test_overnight_isolates_a_failing_step_and_resumes_finished_ones(tmp_path, monkeypatch):
    ov = _overnight()
    calls = []
    def boom(ctx):
        calls.append("boom")
        raise SystemExit("an input is wrong")
    def fine(ctx):
        calls.append("fine")
        return {"title": "fine", "lines": ["a line"], "verdicts": [ov.verdict("c", True, "d")], "files": []}
    monkeypatch.setattr(ov, "STEPS", [("one", 1, 1, 1, [], [], boom), ("two", 2, 1, 1, [], [], fine)])
    argv = ["--out", str(tmp_path / "o"), "--quick"]
    assert ov.main(argv) == 1                                              # one step failed
    assert calls == ["boom", "fine"]                                       # and the other still ran
    summary = (tmp_path / "o" / "OVERNIGHT_SUMMARY.md").read_text(encoding="utf-8")
    assert "| one | 1 | failed |" in summary and "| two | 2 | ok |" in summary and "an input is wrong" in summary and "Not testable" in summary
    assert ov.main(argv) == 1 and calls == ["boom", "fine", "boom"]        # the finished step is not repeated, the failed one is retried
    assert ov.main(argv + ["--fresh"]) == 1 and calls[-2:] == ["boom", "fine"]
    assert ov.main(argv + ["--only", "two"]) == 0
    with pytest.raises(SystemExit, match="unknown step"):
        ov.main(argv + ["--only", "nope"])


def test_overnight_mock_grid_scores_truncations_against_the_expected_sequences(tmp_path, monkeypatch, capsys):
    import gzip
    import json
    import random
    import shutil
    pytest.importorskip("pandas")
    if not shutil.which("vsearch"):
        pytest.skip("vsearch is not installed")
    ov = _overnight()
    import real_data as rd
    rng = random.Random(3)
    seqs = {i: "".join(rng.choice("ACGT") for _ in range(n)) for i, n in ((1, 480), (2, 510), (3, 535), (4, 470))}      # lengths against the 505 and 530 nt merge limits
    reference = tmp_path / "ref.fasta"
    reference.write_text("".join(f">ASV_{i}|S{(i + 1) // 2}|Strain_{(i + 1) // 2}\n{s}\n" for i, s in seqs.items()))
    table = tmp_path / "asvs.csv"
    table.write_text("OTU_ID;103-KitA-1-1-8;104-KitA-1-1-9;137-KitA-NK-1;Sequence\n" + "".join(f"ASV_{i};10;20;0;{seqs[i]}\n" for i in seqs))
    reads = tmp_path / "ready"
    reads.mkdir()
    for run in ("ERR1", "ERR2"):
        for mate in ("R1", "R2"):
            with gzip.open(reads / f"{run}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                fh.write("@r\nACGT\n+\nIIII\n" * 40)
    calls = []
    def fake_dada2(reads_dir, tf, tr, ee, min_overlap, min_len, track, asv):
        calls.append((tf, tr, ee))
        Path(asv).write_text("".join(f">x{i};size=100\n{s}\n" for i, s in seqs.items() if len(s) <= tf + tr - 12))      # DADA2 can only merge what the overlap allows
        return {"filtered": 0.8, "merged": 0.7, "nonchim": 0.7, "asvs": 3}
    monkeypatch.setattr(rd, "dada2_run", fake_dada2)
    monkeypatch.setattr(ov, "tool_missing", lambda tool: False)
    argv = ["--out", str(tmp_path / "o"), "--quick", "--only", "mock_grid", "--mockA-high", str(reads), "--asv-table", str(table), "--reference", str(reference)]
    assert ov.main(argv) == 0
    assert calls == [(257, 260, "2,6"), (268, 249, "2,6"), (279, 263, "2,6")]
    result = json.loads((tmp_path / "o" / "mock_grid" / "result.json").read_text(encoding="utf-8"))
    assert result["status"] == "ok"
    verdicts = {v["criterion"]: v for v in result["verdicts"]}
    acc = [v for k, v in verdicts.items() if k.startswith("accuracy")]
    assert acc and all(v["outcome"] == "met" for v in acc) and "variants -1" in acc[0]["detail"]                    # one variant fewer than the authors', within the margin
    assert verdicts["coverage rule: no variant longer than (total - 12) nt is ever recovered"]["outcome"] == "met"
    ref_file = tmp_path / "o" / "_reference" / "mockA_high_input_reference.fasta"
    assert ref_file.read_text().count(">") == 4
    n = len(calls)
    assert ov.main(argv) == 0 and len(calls) == n                                                                   # resumed: nothing is run twice
    capsys.readouterr()
    pytest.importorskip("rapidfuzz")
    assert ov.main([("asv_classify" if a == "mock_grid" else a) for a in argv]) == 0                                # the classes of the ASVs of that grid
    ac = json.loads((tmp_path / "o" / "asv_classify" / "result.json").read_text(encoding="utf-8"))
    assert ac["status"] == "ok" and any("outside the reference" in v["criterion"] and "+0 in all" in v["detail"] for v in ac["verdicts"])
    assert (tmp_path / "o" / "asv_classify" / "classes_257_260_ee2_6.tsv").exists()


# ------------------------------------------------------------------------------------------------ classify_asvs.py
def test_classify_asvs_sorts_every_asv_into_the_class_with_the_evidence(tmp_path):
    import random
    import shutil
    pytest.importorskip("rapidfuzz")
    pytest.importorskip("pandas")
    if not shutil.which("vsearch"):
        pytest.skip("vsearch is not installed")
    ca = load("classify_asvs")
    rng = random.Random(11)
    rnd = lambda n=300: "".join(rng.choice("ACGT") for _ in range(n))
    A, B, C, D, E = rnd(), rnd(), rnd(), rnd(), rnd()
    def mutate(s, positions):
        out = list(s)
        for p in positions:
            out[p] = {"A": "C", "C": "G", "G": "T", "T": "A"}[out[p]]
        return "".join(out)
    (tmp_path / "ref.fasta").write_text(f">ASV_1|S1|A\n{A}\n>ASV_2|S2|B\n{B}\n>ASV_3|S3|C\n{C}\n")
    (tmp_path / "present.fasta").write_text(f">ASV_1|S1|A\n{A}\n>ASV_2|S2|B\n{B}\n")
    (tmp_path / "table.csv").write_text(f"OTU_ID;Sequence;103-KitA-1-1-8;137-KitA-NK-1;Genus\nASV_9;{D};0;5;Ralstonia\nASV_8;{A};900;0;Strain\n")
    asvs = {"x1": (A, 1000), "x2": (B, 800), "x3": (C, 50), "x4": (mutate(A, [100]), 20), "x5": (A[:150] + B[150:], 15), "x6": (D, 12),
            "x7": (mutate(A, [10, 50, 90, 130, 170, 210, 250]), 10), "x8": (E, 5)}
    (tmp_path / "asvs.fasta").write_text("".join(f">{k};size={n}\n{s}\n" for k, (s, n) in asvs.items()))
    rows = {r["id"]: r for r in ca.classify(tmp_path / "asvs.fasta", tmp_path / "ref.fasta", tmp_path / "present.fasta", tmp_path / "table.csv")}
    got = {k: rows[k]["class"] for k in asvs}
    assert got == {"x1": "expected", "x2": "expected", "x3": "listed elsewhere", "x4": "error variant", "x5": "chimera", "x6": "contaminant", "x7": "unlisted variant", "x8": "unexplained"}, got
    assert "Ralstonia" in rows["x6"]["detail"] and "1 edit" in rows["x4"]["detail"] and "ASV of 1000 reads" in rows["x4"]["detail"]
    s = ca.summarise(list(rows.values()))
    assert s["classes"]["chimera"] == (1, 15) and s["classes"]["expected"] == (2, 1800) and s["asvs"] == 8
    # without the table no contaminant call is possible: the same ASV is then unexplained
    rows2 = {r["id"]: r for r in ca.classify(tmp_path / "asvs.fasta", tmp_path / "ref.fasta", tmp_path / "present.fasta")}
    assert rows2["x6"]["class"] == "unexplained"
    with pytest.raises(SystemExit, match="no column"):
        ca.load_controls(tmp_path / "table.csv", "nothing-matches")


def test_overnight_extra_datasets_readiness_cost_benefit_and_the_landscape(tmp_path, monkeypatch):
    from types import SimpleNamespace
    ov = _overnight()
    assert ov.parse_extra(["mockB=C:\\x\\ready|C:\\x\\ref.fasta", "v4v5=/data/v4v5/ready"]) == {"mockB": (Path("C:\\x\\ready"), Path("C:\\x\\ref.fasta")), "v4v5": (Path("/data/v4v5/ready"), None)}
    with pytest.raises(SystemExit, match="use NAME=READY"):
        ov.parse_extra(["nonsense"])
    ctx = SimpleNamespace(args=SimpleNamespace(extra=None), p=lambda name: None)
    assert ov.readiness(ctx, ["extra"], []) == ["no --extra dataset given"] and "needs mock A" in ov.readiness(ctx, ["anydata"], [])[0]
    ctx2 = SimpleNamespace(args=SimpleNamespace(extra=[f"new={tmp_path}"]), p=lambda name: None)
    assert ov.readiness(ctx2, ["extra", "anydata"], []) == []
    rows = [{"level": 0.01, "tf": 257, "tr": 260, "limit": 505, "d": {"merged": 0.784, "variants": 35, "strains_found": 11, "extra": 15}, "beyond": (4, 0.0031)},
            {"level": 0.001, "tf": 265, "tr": 270, "limit": 523, "d": {"merged": 0.762, "variants": 36, "strains_found": 13, "extra": 12}}]
    text = ov.cost_benefit({"name": "mockA", "rows": rows})
    assert "257/260 to 265/270" in text and "-2.2 points" in text and "variants +1, strains +2, ASVs outside the reference -3" in text and "4 ASVs (0.31% of reads)" in text
    # the landscape: merged share rises with the forward length, so the recommendation is below the best setting that covers the amplicons
    rec = SimpleNamespace(trim_f=250, trim_r=250, yield_est=0.8)
    monkeypatch.setattr(ov, "recommended", lambda ctx, reads, tag, **kw: SimpleNamespace(recommendation=rec, floor=500, sub=SimpleNamespace(read_len_f=300, read_len_r=300), inserts=[480]))
    monkeypatch.setattr(ov, "dada_on", lambda ctx, reads, tf, tr, ee, ref, tag: {"tag": tag, "tf": tf, "tr": tr, "merged": 0.5 + 0.001 * tf if tf + tr >= 500 else 0.2, "variants": None})
    monkeypatch.setattr(ov, "_thin", lambda ctx, reads, name, pairs: reads)
    ctx3 = SimpleNamespace(quick=False, cur=tmp_path)
    res = ov.landscape_dataset(ctx3, "toy", tmp_path, "2", None, 20, 10, None)
    assert res["n"] == 25 and abs(res["regret"] - 2.0) < 1e-9                                                       # 270/x against 250/x: 2 points
    assert any("75.0*" in ln for ln in res["lines"]) and any("20.0." in ln for ln in res["lines"])                   # * marks AmpliTrim's cell, . a total below the coverage floor


# ------------------------------------------------------------------------------------------------ check_build.py
class _Resp:
    def __init__(self, text):
        self.text = text
    def __enter__(self):
        return self
    def __exit__(self, *a):
        return False
    def read(self):
        return self.text.encode()


def test_check_build_compares_manifests_hints_conda_and_registers_downloads(tmp_path, monkeypatch):
    import json
    import types
    cb = load("check_build")
    # the manifest: identical, differing, missing and extra files are told apart
    files = cb.tracked_files()
    lines = [f"{cb.file_hash(f)}  {f.relative_to(cb.ROOT).as_posix()}" for f in files]
    lines[0] = "0000000000" + lines[0][10:]
    lines.append("1111111111  validation/gone.py")
    manifest = tmp_path / "m.txt"
    manifest.write_text("build abc: files\n" + "\n".join(lines[:-2] + [lines[-1]]) + "\n", encoding="utf-8")
    c = cb.compare_manifest(manifest)
    assert c["differ"] == [files[0].relative_to(cb.ROOT).as_posix()] and c["missing"] == ["validation/gone.py"] and c["extra"] == [files[-1].relative_to(cb.ROOT).as_posix()] and c["same"] == len(files) - 2
    # a missing package gets the conda command for this environment
    real = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda name, *a: None if name == "rapidfuzz" else real(name, *a))
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "amplitrim")
    miss = next(m for m in cb.check_machine() if m[0] == "rapidfuzz")
    assert miss[1] is False and miss[3] == "conda install -n amplitrim -c conda-forge rapidfuzz"
    # the run table of a project: only paired runs are kept
    tsv = "run_accession\tsample_alias\tlibrary_layout\tinstrument_model\tread_count\nERR1\ta\tPAIRED\tMiSeq\t10\nERR2\tb\tSINGLE\tMiSeq\t10\nERR3\tc\tPAIRED\tMiSeq\t10\n"
    got = cb.fetch_run_table("PRJEB1", tmp_path / "t", opener=lambda url, timeout=0: _Resp(tsv))
    assert got == ["ERR1", "ERR3"] and (tmp_path / "t" / "runs.txt").read_text() == "ERR1\nERR3\n"
    with pytest.raises(SystemExit, match="no run table"):
        cb.fetch_run_table("PRJEB1", tmp_path / "t2", opener=lambda url, timeout=0: _Resp("nothing\n"))
    # download, prepare and register with the commands that would be run
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "extras": {}}
    cmds = []
    def fake_run(cmd, *a, **k):
        cmds.append([str(x) for x in cmd])
        if "run_dataset.py" in cmds[-1][1]:
            (tmp_path / "data" / "new" / "mockB" / "prepared" / "mockB" / "ready").mkdir(parents=True)
        return types.SimpleNamespace(returncode=0)
    args = types.SimpleNamespace(runs_file=None, max_runs=3, dry_run=False, primer_f=None, primer_r=None, select=".")                 # everything: the catalog's own selection of mock B names samples this table does not have
    assert cb.download(cfg, "mockB", args, run=fake_run, opener=lambda url, timeout=0: _Resp(tsv)) == 0
    assert "download_ena.py" in cmds[0][1] and cmds[0][3:] == ["--out", str(tmp_path / "data" / "new" / "mockB" / "raw"), "--max-runs", "3"] and "describe_reads.py" in cmds[1][1]
    with pytest.raises(SystemExit, match="runs-file"):
        cb.download(cfg, "v3v4", args, run=fake_run)
    with pytest.raises(SystemExit, match="no reads"):
        cb.prepare(cfg, "mockB", args, tmp_path / "cfg.json", run=fake_run)
    (tmp_path / "data" / "new" / "mockB" / "raw" / "ERR1_R1.fastq.gz").write_bytes(b"x")
    cfg_path = tmp_path / "cfg.json"
    assert cb.prepare(cfg, "mockB", args, cfg_path, run=fake_run) == 0
    assert "--primer-f" in cmds[-1] and "AGAGTTTGATYMTGGCTCAG" in cmds[-1] and json.loads(cfg_path.read_text())["extras"]["mockB"]["ready"].endswith("ready")
    (tmp_path / "data" / "new" / "v4v5" / "raw").mkdir(parents=True)
    (tmp_path / "data" / "new" / "v4v5" / "raw" / "R_R1.fastq.gz").write_bytes(b"x")
    with pytest.raises(SystemExit, match="primers of v4v5 are not known"):
        cb.prepare({**cfg, "extras": {}}, "v4v5", args, cfg_path, run=fake_run)
    ref = tmp_path / "zymo.fasta"
    ref.write_text(">a\nACGT\n")
    assert cb.set_reference(cfg, f"mockB={ref}", cfg_path) == 0
    # the paths reach overnight.py as the options it knows
    oa = cb.overnight_args(cfg, types.SimpleNamespace(quick=True, only="mock_grid", skip=None, fresh=False))
    assert oa[:2] == ["--out", cfg["out"] + "_quick"] and "--extra" in oa and oa[oa.index("--extra") + 1].endswith(f"|{ref}") and "--quick" in oa and oa[oa.index("--only") + 1] == "mock_grid"


def test_check_build_reports_a_missing_manifest_plainly_and_uses_the_shipped_one(tmp_path, monkeypatch, capsys):
    import types
    cb = load("check_build")
    lines = cb.check_build(types.SimpleNamespace(manifest=str(tmp_path / "nope.txt")))
    assert any("manifest file not found" in ln and "nope.txt" in ln for ln in lines)                       # a message, not a traceback
    shipped = tmp_path / "build_manifest.txt"
    files = cb.file_lines()
    shipped.write_text("build abc1234567: files\n" + "\n".join(files) + "\n", encoding="utf-8")
    monkeypatch.setattr(cb, "SHIPPED_MANIFEST", shipped)
    ok = cb.check_build(types.SimpleNamespace(manifest=None))
    assert any("(build abc1234567)" in ln and f"{len(files)} files identical" in ln for ln in ok) and any("exactly that build" in ln for ln in ok)
    shipped.write_text("build abc1234567: files\n" + "\n".join(["0000000000" + files[0][10:]] + files[1:-1]) + "\n", encoding="utf-8")   # one edited, one new file
    bad = "\n".join(cb.check_build(types.SimpleNamespace(manifest=None)))
    assert "differs from that build" in bad and "is here but not in that build" in bad
    assert cb.main(["--files"]) == 0 and capsys.readouterr().out.splitlines()[0].count("  ") == 1


def test_the_build_id_is_the_same_on_windows_and_linux_orderings():
    import hashlib
    from pathlib import PureWindowsPath
    rd = load("run_dataset")
    files = list((VALIDATION.parent / "amplitrim").glob("*.py")) + list(VALIDATION.glob("*.py")) + list(VALIDATION.glob("*.R"))
    names = [f.name for f in files]
    assert "mock_check.py" in names and "mockA_dilutions.py" in names                  # the pair that Windows orders differently from Linux
    def digest(order):
        h = hashlib.sha256()
        for p in order:
            h.update(p.name.encode())
            h.update(p.read_bytes().replace(b"\r\n", b"\n"))
        return h.hexdigest()[:10]
    windows_order = sorted(files, key=lambda p: PureWindowsPath(str(p).replace("/", "\\")))
    case_sensitive_order = sorted(files, key=lambda p: p.as_posix())                       # what Linux gives; written out so that the test means the same on every system
    assert [f.name for f in windows_order] != [f.name for f in case_sensitive_order]      # the two orderings really differ here
    assert rd.build_id() == digest(sorted(files, key=lambda p: (p.parent.name, p.name))) and digest(windows_order) != rd.build_id()
    cb = load("check_build")
    assert [f.name for f in cb.tracked_files()] == [f.name for f in sorted(files, key=lambda p: (p.parent.name, p.name))]


def test_check_build_lists_runs_by_kind_selects_them_and_passes_the_sample_choice_on(tmp_path):
    import json
    import types
    cb = load("check_build")
    tsv = ("run_accession\tsample_alias\tlibrary_layout\tinstrument_model\tread_count\n"
           "ERR1\t103-KitA-1-1-8_S103_L001\tPAIRED\tMiSeq\t30000\nERR2\t104-KitA-1-1-9_S104_L001\tPAIRED\tMiSeq\t40000\n"
           "ERR3\tZymo-D6321-5_S5_L001\tPAIRED\tMiSeq\t500\nERR4\tZymo-D6321-6_S6_L001\tPAIRED\tMiSeq\t9000\nERR5\tx\tSINGLE\tMiSeq\t9\n").replace("\\\t", "\t").replace("\\\n", "\n")
    opener = lambda url, timeout=0: _Resp(tsv)
    assert cb.sample_kind("103-KitA-1-1-8_S103_L001") == "#-KitA-#-#-#" and cb.sample_kind("Zymo-D6321-5_S5_L001") == "Zymo-D#-#"
    assert cb.fetch_run_table("P", tmp_path / "a", opener, select="zymo") == ["ERR3", "ERR4"]                              # the sample name, in any case
    assert cb.fetch_run_table("P", tmp_path / "b", opener, select="zymo", min_reads=1000) == ["ERR4"]
    assert (tmp_path / "b" / "runs.tsv").read_text().count("ERR") == 5                                                        # the table keeps every run
    with pytest.raises(SystemExit, match="matches --select 'nothing'|no run of P matches"):
        cb.fetch_run_table("P", tmp_path / "c", opener, select="nothing")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "extras": {}}
    d = cb.dataset_dirs(cfg, "mockB")
    d["raw"].mkdir(parents=True)
    (d["raw"] / "ERR3_R1.fastq.gz").write_bytes(b"x")
    lines = cb.list_runs(cfg, "mockB", opener)
    text = "\n".join(lines)
    assert cb.downloaded_aliases(d) == ["  ERR3  Zymo-D6321-5_S5_L001                       500 read pairs"]
    assert "4 paired runs in 2 kinds of sample name" in lines[0] and "#-KitA-#-#-#" in text and "Zymo-D#-#" in text and "103-KitA-1-1-8_S103_L001 (ERR1)" in text
    zymo = next(ln for ln in lines if ln.startswith("Zymo-D#-#"))
    assert zymo.split()[1] == "2" and "500 / 9,000 / 9,000" in zymo and zymo.split("9,000")[-1].split()[0] == "1"              # two runs, one of them downloaded
    # the choice of samples, the orientation and the read minimum go on to run_dataset.py
    cmds = []
    def fake_run(cmd, *a, **k):
        cmds.append([str(x) for x in cmd])
        d["ready"].mkdir(parents=True, exist_ok=True)
        return types.SimpleNamespace(returncode=0)
    args = types.SimpleNamespace(primer_f=None, primer_r=None, samples="ERR3,ERR4", mixed_orientation=True, min_reads=500)
    cfg_path = tmp_path / "cfg.json"
    ref = tmp_path / "zymo.fasta"
    ref.write_text(">a\nACGT\n")
    assert cb.set_reference(cfg, f"mockB={ref}", cfg_path) == 0                                                                # before the reads are prepared
    assert cb.overnight_args(cfg, types.SimpleNamespace(quick=False, only=None, skip=None, fresh=False)).count("--extra") == 0   # so it is not used yet
    assert any(item[0] == "extra mockB" and item[1] is False for item in cb.check_data(cfg))
    assert cb.prepare(cfg, "mockB", args, cfg_path, run=fake_run) == 0
    assert cmds[-1][cmds[-1].index("--only") + 1] == "ERR3,ERR4" and "--mixed-orientation" in cmds[-1] and cmds[-1][cmds[-1].index("--min-reads") + 1] == "500"
    saved = json.loads(cfg_path.read_text())["extras"]["mockB"]
    assert saved["reference"] == str(ref) and saved["ready"].endswith("ready")                                                 # the reference kept, the reads registered
    assert cb.overnight_args(cfg, types.SimpleNamespace(quick=False, only=None, skip=None, fresh=False)).count("--extra") == 1
    with pytest.raises(SystemExit, match="not one of"):
        cb.set_reference(cfg, f"nope={ref}", cfg_path)


def test_check_build_discovers_folders_by_sample_count_and_figaro(tmp_path, monkeypatch):
    cb = load("check_build")
    def ready(path, n):
        path.mkdir(parents=True)
        for i in range(n):
            (path / f"S{i}_R1.fastq.gz").write_bytes(b"x")
            (path / f"S{i}_R2.fastq.gz").write_bytes(b"x")
    data = tmp_path / "data"
    ready(data / "results3" / "mockA_high" / "ready", 6)
    ready(data / "results2" / "mockA_ee26" / "ready", 23)
    ready(data / "results" / "prjeb42528" / "ready", 33)
    ready(data / "results4" / "cryobank_its3" / "ready", 37)
    ready(data / "results4" / "other" / "ready", 12)                                           # a size that belongs to no known dataset
    (data / "mockA").mkdir()
    (data / "mockA" / "runs_info.tsv").write_text("run_accession\n")
    (data / "figaro" / "figaro").mkdir(parents=True)
    (data / "figaro" / "figaro" / "figaro.py").write_text("")
    monkeypatch.setenv("HOME", str(tmp_path / "home")); monkeypatch.setenv("USERPROFILE", str(tmp_path / "home"))
    cfg = {**cb.DEFAULT_CONFIG, "data": str(data), "mockA_high": "", "mockA_all": "", "its": str(data / "results4" / "cryobank_its3" / "ready"), "v1v2": "", "figaro": "",
           "runinfo": "", "asv_table": "", "reference": str(tmp_path / "nope.fasta")}
    found = cb.discover(cfg)
    assert found["mockA_high"][0].parent.name == "mockA_high" and found["mockA_all"][0].parent.name == "mockA_ee26" and found["v1v2"][0].parent.name == "prjeb42528"
    assert "its" not in found                                                                  # already set and present: left alone
    assert found["runinfo"][0].name == "runs_info.tsv" and "asv_table" not in found
    assert found["figaro"][0].name == "figaro.py"
    cfg = cb.apply_discovery(cfg, found)
    assert cfg["v1v2"].endswith("ready") and "v1v2" not in cb.discover(cfg) and "figaro" not in cb.discover(cfg)       # a second pass finds nothing the first one filled


def test_check_build_decides_the_coverage_default_by_the_rule(tmp_path):
    cb = load("check_build")
    def rows(a, b, truth, a_strains=11, b_strains=12, a_var=35, b_var=36):
        r1 = {"level": 0.01, "tf": 257, "tr": 260, "limit": 505, "left": 0.0031, "est": 0.8, "merged": a, "variants": a_var if truth else None, "present": 39, "strains_found": a_strains if truth else None, "strains_present": 13, "outside": 15}
        r2 = {**r1, "level": 0.002, "tf": 257, "tr": 263, "merged": b, "variants": b_var if truth else None, "strains_found": b_strains if truth else None, "outside": 11}
        return [r1, r2]
    mock_a = {"name": "mockA_high", "has_truth": True, "rows": rows(0.784, 0.782, True)}
    d = cb.decide_coverage([mock_a])
    assert d["verdict"] == "UNDETERMINED" and "second dataset with expected sequences" in d["why"] and "cost only" not in d["why"]
    v1v2 = {"name": "v1v2", "has_truth": False, "rows": rows(0.842, 0.838, False)}
    assert cb.decide_coverage([mock_a, v1v2])["verdict"] == "UNDETERMINED"                       # still one dataset with expected sequences
    mock_b = {"name": "mockB", "has_truth": True, "rows": rows(0.70, 0.695, True, 3, 3, 10, 10)}   # nothing lost, nothing gained there
    adopt = cb.decide_coverage([mock_a, mock_b, v1v2])
    assert adopt["verdict"] == "ADOPT 0.2% as the default" and "benefit on the datasets with expected sequences: yes" in adopt["why"]
    worse = {"name": "mockB", "has_truth": True, "rows": rows(0.70, 0.695, True, 4, 3, 10, 10)}   # a strain lost
    assert cb.decide_coverage([mock_a, worse, v1v2])["verdict"] == "KEEP 1%"
    within = {"name": "v1v2", "has_truth": False, "rows": rows(0.842, 0.80, False)}            # 4.2 points: inside both limits
    assert cb.decide_coverage([mock_a, mock_b, within])["verdict"] == "ADOPT 0.2% as the default"
    too_much = {"name": "v1v2", "has_truth": False, "rows": rows(0.842, 0.78, False)}          # 6.2 points: over the limit of 5 on one dataset
    d6 = cb.decide_coverage([mock_a, mock_b, too_much])
    assert d6["verdict"] == "KEEP 1%" and "largest cost 6.2 points (limit 5.0): NO" in d6["why"]
    slow = [{"name": f"d{i}", "has_truth": False, "rows": rows(0.80, 0.77, False)} for i in range(3)]          # 3 points on three datasets: the median is over 2
    d3 = cb.decide_coverage([mock_a, mock_b] + slow)
    assert d3["verdict"] == "KEEP 1%" and "median cost 3.0 points (limit 2.0): NO" in d3["why"]
    # the saved numbers are read back from the step results
    import json
    (tmp_path / "coverage_levels").mkdir()
    (tmp_path / "coverage_levels" / "result.json").write_text(json.dumps({"status": "ok", "data": [mock_a, v1v2], "verdicts": []}))
    assert [x["name"] for x in cb.collect_coverage_data(tmp_path)] == ["mockA_high", "v1v2"]


def test_check_build_todo_board_and_version_lines(tmp_path, monkeypatch, capsys):
    import json
    import types
    cb = load("check_build")
    out = tmp_path / "overnight"
    (out / "mock_grid").mkdir(parents=True)
    (out / "mock_grid" / "result.json").write_text(json.dumps({"status": "ok", "verdicts": [{"outcome": "NOT MET"}, {"outcome": "met"}], "title": "x", "lines": []}))
    (out / "coverage_levels").mkdir()
    (out / "coverage_levels" / "result.json").write_text(json.dumps({"status": "ok", "verdicts": [], "title": "x", "lines": []}))   # saved before the numbers were kept
    cfg = {**cb.DEFAULT_CONFIG, "out": str(out), "data": str(tmp_path / "nodata"), "mockA_high": "", "mockA_all": "", "its": "", "v1v2": "", "figaro": "", "extras": {}}
    monkeypatch.setattr(cb, "discover", lambda c: {})
    text = "\n".join(cb.todo(cfg, types.SimpleNamespace()))
    assert "DONE    mock_grid" in text and "1 criterion not met" in text and "saved before the decision rule existed" in text and "BLOCKED figaro_v1v2" in text
    assert "Default of --allow-unmergeable (1% -> 0.2%): UNDETERMINED" in text and "[ ] Title block: affiliation" in text and "[ ] Archive the release for a DOI" in text
    assert "--list-runs mockB" in text and "Next, in this order" in text
    # a version line is never a citation
    class R:
        returncode = 0
        stdout = ""
        stderr = "vsearch v2.30.0_win_x86_64, 15.9GB RAM, 12 cores\nhttps://github.com/torognes/vsearch\n\nRognes T, Flouri T (2016) VSEARCH PeerJ 4:e2584 doi: 10.7717/peerj.2584\n"
    monkeypatch.setattr(subprocess_module(), "run", lambda *a, **k: R())
    assert cb.tool_version(["vsearch", "--version"]).startswith("vsearch v2.30.0")
    R.stderr = "Rognes T (2016) PeerJ 4:e2584 doi: 10.7717/peerj.2584 https://doi.org/10.7717/peerj.2584\n"
    assert "doi" not in (cb.tool_version(["vsearch", "--version"]) or "")


def subprocess_module():
    import subprocess
    return subprocess


# ------------------------------------------------------------------------------------------------ the test-stage battery and table
def test_battery_rows_give_the_points_and_the_intervals_with_the_right_sign():
    ov = _overnight()
    rows = [{"name": "AmpliTrim (maxEE 2)", "truncF": "198", "truncR": "157", "maxEE": "2", "est_pass": "0.915", "pass": "0.913", "est_yield": "0.906", "merged": "0.842"},
            {"name": "full", "truncF": "0", "truncR": "0", "maxEE": "2", "est_pass": "", "pass": "0.60", "est_yield": "", "merged": "0.50"}]
    text = "  full: -34.0 points (95% interval -36.0 to -32.0); better in 0 of 33 samples, worse in 33\\n"      # full minus AmpliTrim, so AmpliTrim minus full is the opposite
    r = ov.battery_row("v1v2", "16s", "V1-V2", 33, rows, text, {"like_for_like": {"mean": 0.1, "lo": -0.4, "hi": 0.6}})
    assert r["setting"] == "198/157" and abs(r["filter_diff"] - 0.2) < 1e-9 and abs(r["yield_offset"] - 6.4) < 1e-9 and abs(r["ratio"] - 0.842 / 0.906) < 1e-9
    assert r["vs_other"] == {"label": "full", "mean": 34.0, "lo": 32.0, "hi": 36.0} and r["vs_figaro"] == {"mean": 0.1, "lo": -0.4, "hi": 0.6}
    its = [{"name": f"AmpliTrim (maxEE {e})", "truncF": "0", "truncR": "0", "maxEE": str(e), "est_pass": str(ep), "pass": str(dp), "est_yield": str(ey), "merged": str(m)}
           for e, ep, dp, ey, m in ((1, 0.579, 0.577, 0.564, 0.517), (2, 0.754, 0.755, 0.731, 0.633), (5, 0.958, 0.957, 0.842, 0.730))]
    r2 = ov.battery_row("its", "its", "ITS1", 8, its + [{"name": "trunc200", "maxEE": "2"}], "  trunc200: -23.5 points (95% interval -25.2 to -21.8); better in 0 of 8\\n")
    assert r2["setting"] == "no truncation" and r2["maxEE"] == "2" and abs(r2["filter_diff"] - (-0.1)) < 1e-9 and abs(r2["yield_offset"] - 9.8) < 1e-9       # the maxEE 2 row
    assert [round(x["yield_offset"], 1) for x in r2["its_rows"]] == [4.7, 9.8, 11.2] and r2["vs_other"]["mean"] == 23.5 and r2["vs_other"]["lo"] == 21.8
    with pytest.raises(RuntimeError, match="no AmpliTrim row"):
        ov.battery_row("x", "16s", "V4", 3, [{"name": "full"}], "")


def test_check_build_renders_the_test_tables_per_region(tmp_path):
    import json
    cb = load("check_build")
    row16 = {"name": "v1v2", "marker": "16s", "region": "V1-V2", "samples": 33, "setting": "198/157", "maxEE": "2", "filter_diff": 0.2, "yield_offset": 6.4, "ratio": 0.929, "merged": 0.842,
             "vs_other": {"label": "full", "mean": 34.0, "lo": 32.0, "hi": 36.0}, "vs_figaro": {"mean": 0.1, "lo": -0.4, "hi": 0.6}}
    rowits = {"name": "its", "marker": "its", "region": "ITS1", "samples": 37, "setting": "no truncation", "maxEE": "2", "filter_diff": -0.1, "yield_offset": 9.8, "ratio": 0.866, "merged": 0.633,
              "vs_other": {"label": "trunc200", "mean": 23.5, "lo": 21.8, "hi": 25.2}, "vs_figaro": None,
              "its_rows": [{"maxEE": "1", "filter_diff": 0.2, "yield_offset": 4.7, "merged": 0.517, "ratio": 0.917}, {"maxEE": "2", "filter_diff": -0.1, "yield_offset": 9.8, "merged": 0.633, "ratio": 0.866}]}
    cov = {"name": "v1v2", "has_truth": False, "rows": [{"level": 0.01, "merged": 0.842}, {"level": 0.002, "merged": 0.840}]}
    lines = cb.render_tables([row16, rowits], {"v1v2": cov}, {"v1v2": 0.91}, ["v4 (soil)"])
    text = "\\n".join(lines)
    assert "### 16S rRNA datasets" in text and "### ITS datasets" in text
    assert "| v1v2 | V1-V2 | 33 | 198/157 (maxEE 2) | +0.2 | +6.4 | 0.93 | - | - | +34.0 (+32.0 to +36.0) | +0.1 (-0.4 to +0.6) | - | -0.2 | - | - | +0.91 |" in text
    assert "| v4 (soil) | | | not yet run |" in text
    assert "| its | ITS1 | 37 | 1 | +0.2 | +4.7 | 0.92 | +23.5 (+21.8 to +25.2) |" in text and "|  |  |  | 2 | -0.1 | +9.8 | 0.87 |  |" in text
    # from saved results
    (tmp_path / "battery").mkdir()
    (tmp_path / "battery" / "result.json").write_text(json.dumps({"status": "ok", "data": [row16, rowits], "verdicts": []}))
    (tmp_path / "coverage_levels").mkdir()
    (tmp_path / "coverage_levels" / "result.json").write_text(json.dumps({"status": "ok", "data": [cov], "verdicts": []}))
    (tmp_path / "per_sample_rank").mkdir()
    (tmp_path / "per_sample_rank" / "result.json").write_text(json.dumps({"status": "ok", "data": [{"name": "v1v2", "rho": 0.91, "samples": 33, "offset_mean": 5.0}], "verdicts": []}))
    out = "\\n".join(cb.table_lines({**cb.DEFAULT_CONFIG, "out": str(tmp_path)}))
    assert "+0.91 |" in out and "| -0.2 |" in out and "| mockB (independent mock, V1-V3) | | | not yet run" in out and "| v4 (PRJNA574030, V4) | | | not yet run" in out and "| v3v4 (V3-V4, tight overlap) | | | not yet run" in out


def test_check_build_registers_the_v4_folder_by_its_36_samples_and_the_filter_grid_compares_counts(tmp_path, monkeypatch):
    cb = load("check_build")
    monkeypatch.setattr(cb, "find_file", lambda *a, **k: None)         # the files are looked for by name in the user's Downloads and in C:/figaro: not in this test (it is about the reads)
    data = tmp_path / "data" / "results" / "v4" / "ready"
    data.mkdir(parents=True)
    for i in range(36):
        (data / f"S{i}_R1.fastq.gz").write_bytes(b"x")
        (data / f"S{i}_R2.fastq.gz").write_bytes(b"x")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "mockA_high": "x", "mockA_all": "x", "v1v2": "x", "its": "x", "figaro": "x", "runinfo": "x", "asv_table": "x", "reference": "x", "extras": {}}
    found = cb.discover(cfg)
    assert list(found) == ["v4"] and found["v4"][0].name == "ready"
    cfg = cb.apply_discovery(cfg, found)
    assert cfg["extras"]["v4"]["ready"].endswith("ready") and cb.discover(cfg) == {}                # registered as a dataset that is tested like a new one
    fg = load("filter_grid_check")
    class Table:
        def rate(self, tf, tr):
            return 0.5
    class Pairs:
        def pass_table(self, ee):
            return Table()
    rows = [{"sample": "s1", "tf": "150", "tr": "130", "ee_f": "2", "ee_r": "2", "n_in": "1000", "n_out": "500"}, {"sample": "s1", "tf": "240", "tr": "220", "ee_f": "2", "ee_r": "2", "n_in": "1000", "n_out": "499"}]
    out = fg.compare(rows, lambda sample: Pairs())
    assert [r[6] - r[7] for r in out] == [0, -1]                                                  # the second setting differs by one pair


def test_the_filter_grid_check_is_identical_to_dada2_on_the_example_reads(tmp_path):
    import shutil
    import subprocess
    if not shutil.which("Rscript"):
        pytest.skip("R is not installed")
    if subprocess.run(["Rscript", "-e", "library(dada2)"], capture_output=True).returncode != 0:
        pytest.skip("the dada2 R package is not installed")
    fg = load("filter_grid_check")
    assert fg.main(["--quick", "--out", str(tmp_path / "o")]) == 0
    lines = (tmp_path / "o" / "grid_compare.csv").read_text().splitlines()
    assert len(lines) == 17 and all(ln.split(",")[-1] == ln.split(",")[-2] for ln in lines[1:])


def test_the_simulators_in_the_repository_reproduce_their_documented_communities(tmp_path):
    import subprocess
    sim = VALIDATION / "simulated"
    def run(script, *args):
        res = subprocess.run([sys.executable, str(sim / script), *args], capture_output=True, text=True)
        assert res.returncode == 0, res.stderr[-300:]
        return res.stdout.strip().splitlines()[-1]
    # the coverage rule's community (120 taxa, 439-504 nt) and the one where both reads degrade (430-510 nt), each with two samples
    assert run("coverage_sim.py", str(tmp_path / "a"), "noisy_reverse") == "120 taxa, lengths 439-504 (median 472); 2 x 14,000 pairs"
    assert run("coverage_sim.py", str(tmp_path / "b"), "both_degrade") == "120 taxa, lengths 430-510 (median 472); 2 x 14,000 pairs"
    assert sorted(p.name for p in (tmp_path / "a" / "reads").iterdir()) == ["s1_R1.fastq.gz", "s1_R2.fastq.gz", "s2_R1.fastq.gz", "s2_R2.fastq.gz"] and (tmp_path / "a" / "truth.fasta").read_text().count(">") == 120
    assert subprocess.run([sys.executable, str(sim / "coverage_sim.py"), str(tmp_path / "c"), "nonsense"], capture_output=True, text=True).returncode != 0
    # the planted chimeras (7.8% of the pairs, 12 recurring types)
    assert run("chimera_sim.py", str(tmp_path / "d")) == "36000 pairs, 2813 (7.8%) from 12 recurring chimera types; true taxa: 40"
    assert run("its_sim.py", str(tmp_path / "e")).startswith("60 taxa; fragment lengths with primers 122-402 (median 257); 29 of 60 taxa are shorter than the 250 nt reads")


def test_every_step_in_the_overnight_registry_is_a_defined_function_with_a_unique_name():
    ov = _overnight()
    ids = [s[0] for s in ov.STEPS]
    assert len(ids) == len(set(ids)) and "battery" in ids and "extra_datasets" in ids
    for sid, tier, minutes, quick, needs, tools, fn in ov.STEPS:
        assert callable(fn) and fn.__name__.startswith("step_") and tier in (1, 2, 3, 4) and minutes > 0 and quick > 0, sid
    for helper in ("battery_row", "battery_dataset", "per_sample_rank_one", "coverage_dataset", "landscape_dataset"):
        assert callable(getattr(ov, helper)), helper


def test_the_test_register_lists_every_overnight_step_and_every_development_check():
    ov = _overnight()
    register = (Path(__file__).resolve().parents[1] / "TEST_REGISTER.md").read_text(encoding="utf-8")
    for sid in (s[0] for s in ov.STEPS):
        assert f"--only {sid}" in register or f",{sid}" in register or f"{sid}," in register, f"{sid} is not in TEST_REGISTER.md"
    for dev in ("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"):
        assert f"| {dev} |" in register
    assert "filter_grid_check.py" in register and "simulated\\coverage_sim.py" in register


def test_the_standard_battery_gives_every_dataset_the_same_columns_with_figaro_end_to_end_landscape_stability_and_truth(tmp_path):
    ov = _overnight()
    cb = load("check_build")
    base = {"name": "AmpliTrim (maxEE 2,6)", "truncF": "257", "truncR": "260", "maxEE": "2,6", "est_pass": "0.814", "pass": "0.816", "est_yield": "0.807", "merged": "0.784", "found": "35", "identical": "0.9973", "extra": "15", "strains": "11/13"}
    figaro = {"name": "FIGARO (amplicon 505)", "truncF": "262", "truncR": "255", "maxEE": "2,2", "merged": "0.70", "found": "36", "identical": "0.9970", "extra": "13", "strains": "12/13"}
    full = {"name": "full", "truncF": "0", "truncR": "0", "maxEE": "2,6", "merged": "0.50", "found": "30", "identical": "0.99", "extra": "30", "strains": "9/13"}
    numbers = {"like_for_like": {"mean": 0.1, "lo": -0.4, "hi": 0.6}, "end_to_end": {"mean": 8.4, "lo": 7.0, "hi": 9.8}}
    r = ov.battery_row("mockA_high", "16s", "V1-V3", 6, [base, figaro, full], "  full: -28.0 points (95% interval -30.0 to -26.0); better in 0 of 6\n", numbers)
    assert r["vs_figaro_end"] == {"mean": 8.4, "lo": 7.0, "hi": 9.8} and r["vs_figaro"]["mean"] == 0.1
    assert [(t["method"], t["setting"], t["found"], t["strains"]) for t in r["truth"]] == [("AmpliTrim (maxEE 2,6)", "257/260", 35.0, "11/13"), ("FIGARO (amplicon 505)", "262/255", 36.0, "12/13"), ("full", "no truncation", 30.0, "9/13")]
    assert "truth" not in ov.battery_row("v1v2", "16s", "V1-V2", 33, [{"name": "AmpliTrim (maxEE 2)", "truncF": "198", "truncR": "157", "maxEE": "2", "est_pass": "0.9", "pass": "0.9", "est_yield": "0.9", "merged": "0.8"}], "")
    lines = cb.render_tables([r], {}, {"mockA_high": 0.6}, [], {"mockA_high": 0.3}, {"mockA_high": 0})
    text = "\n".join(lines)
    assert "AmpliTrim default - FIGARO as shipped" in text and "+8.4 (+7.0 to +9.8)" in text and "| 0.3 | 0 | +0.60 |" in text
    assert "### Datasets with expected sequences" in text and "| mockA_high | FIGARO (amplicon 505) | 262/255 | 2,2 | 70.0% | 36 | 12/13 | 99.70% | 13 |" in text and "| mockA_high | full | no truncation | 2,6 | 50.0% | 30 | 9/13 | 99.00% | 30 |" in text


def test_the_battery_row_summarises_the_calibration_across_the_maxEE_settings():
    ov = _overnight()
    cb = load("check_build")
    its = [{"name": f"AmpliTrim (maxEE {e})", "truncF": "0", "truncR": "0", "maxEE": str(e), "est_pass": str(ep), "pass": str(dp), "est_yield": str(ey), "merged": str(m)}
           for e, ep, dp, ey, m in ((1, 0.579, 0.577, 0.564, 0.517), (2, 0.754, 0.755, 0.731, 0.633), (5, 0.958, 0.957, 0.842, 0.730))]
    r = ov.battery_row("its", "its", "ITS1", 8, its, "")
    c = r["calibration"]
    assert c["n"] == 3 and abs(c["filter_max"] - 0.2) < 1e-9 and abs(c["offset_lo"] - 4.7) < 1e-9 and abs(c["offset_hi"] - 11.2) < 1e-9 and abs(c["offset_mean"] - 25.7 / 3) < 1e-9
    assert abs(c["ratio_lo"] - 0.866) < 1e-3 and abs(c["ratio_hi"] - 0.917) < 1e-3 and c["rank"] == 1.0
    assert "calibration" not in ov.battery_row("v", "16s", "V1-V2", 33, [its[0]], "")                               # one setting: no calibration across settings
    text = "\n".join(cb.render_tables([{**r, "marker": "16s", "setting": "198/157", "samples": 8}], {}, {}, []))
    assert "| +8.6 (+4.7 to +11.2) | +1.00 (3) |" in text


def test_the_battery_does_not_offer_figaro_to_its_reads(tmp_path):
    ov = _overnight()
    reads = tmp_path / "reads"
    reads.mkdir()
    for suffix in ("R1", "R2"):
        (reads / f"S1_{suffix}.fastq.gz").write_bytes(b"x")
    calls = []

    class FakeCtx:
        cur = tmp_path / "work"
        def p(self, key):
            return "figaro.py" if key == "figaro" else None
        def sh(self, cmd, label):
            calls.append([str(c) for c in cmd])
            out = Path(cmd[cmd.index("--out") + 1])
            out.mkdir(parents=True, exist_ok=True)
            (out / "compare_methods.tsv").write_text("name\ttruncF\ttruncR\tmaxEE\test_pass\tpass\test_yield\tmerged\nAmpliTrim (maxEE 2)\t0\t0\t2\t0.75\t0.75\t0.73\t0.63\n")
            (out / "compare_methods.txt").write_text("")
            return ""
    row, lines, verdicts = ov.battery_dataset(FakeCtx(), "its", "its", "ITS1", reads, None, 100, "2")
    assert calls and all("--figaro-script" not in c for c in calls) and row["setting"] == "no truncation" and verdicts == []


def test_the_summary_keeps_the_steps_of_earlier_runs(tmp_path):
    import json
    import types
    ov = _overnight()
    (tmp_path / "mock_grid").mkdir()
    (tmp_path / "mock_grid" / "result.json").write_text(json.dumps({"id": "mock_grid", "tier": 1, "status": "ok", "seconds": 600, "title": "grid", "lines": ["earlier run line"], "verdicts": [{"criterion": "c1", "outcome": "met", "detail": "d"}], "error": ""}))
    ctx = types.SimpleNamespace(out=tmp_path, quick=False)
    now = [{"id": "landscape", "tier": 1, "status": "ok", "seconds": 60, "title": "land", "lines": ["this run line"], "verdicts": [], "error": ""}]
    text = ov.write_summary(ctx, now).read_text(encoding="utf-8")
    assert "earlier run line" in text and "this run line" in text and text.index("| mock_grid |") < text.index("| landscape |")      # in the order of the registry
    assert "earlier run line" in ov.write_summary(ctx, []).read_text(encoding="utf-8")                                       # --summary: nothing run


# ------------------------------------------------------------------------------------------------ the one overnight command
def test_a_step_is_repeated_when_a_dataset_is_added_and_kept_otherwise(tmp_path, monkeypatch):
    import json
    import types
    ov = _overnight()
    calls = []

    def fake(ctx):
        calls.append(1)
        return {"title": "t", "lines": [], "verdicts": [], "files": [], "data": []}
    monkeypatch.setattr(ov, "STEPS", [("coverage_levels", 1, 1, 1, [], [], fake)])

    class C:
        def __init__(self, extra):
            self.out, self.quick, self.cur = tmp_path, False, None
            self.args = types.SimpleNamespace(fresh=False, extra=extra)
        def p(self, key):
            return Path("/x/v1v2") if key == "v1v2" else None
        def log(self, message):
            pass
    ov.run_steps(C([]), ["coverage_levels"])
    assert len(calls) == 1 and json.loads((tmp_path / "coverage_levels" / "result.json").read_text())["datasets"] == ["v1v2"]
    ov.run_steps(C([]), ["coverage_levels"])
    assert len(calls) == 1                                                          # the same datasets: kept
    ov.run_steps(C(["v4=/d/v4"]), ["coverage_levels"])
    assert len(calls) == 2 and json.loads((tmp_path / "coverage_levels" / "result.json").read_text())["datasets"] == ["v1v2", "v4"]       # a dataset was added: repeated
    # a result saved before the datasets were recorded: its data names stand in
    (tmp_path / "coverage_levels" / "result.json").write_text(json.dumps({"id": "coverage_levels", "tier": 1, "status": "ok", "seconds": 1, "title": "t", "lines": [], "verdicts": [], "data": [{"name": "v1v2"}]}))
    ov.run_steps(C([]), ["coverage_levels"])
    assert len(calls) == 2
    (tmp_path / "coverage_levels" / "result.json").write_text(json.dumps({"id": "coverage_levels", "tier": 1, "status": "ok", "seconds": 1, "title": "t", "lines": [], "verdicts": []}))
    ov.run_steps(C([]), ["coverage_levels"])
    assert len(calls) == 3                                                          # no record of what it covered: repeated once
    # steps that do not depend on the datasets are never repeated for that reason
    assert ov.step_datasets(C(["v4=/d/v4"]), "mock_grid") is None and ov.step_datasets(C(["v4=/d/v4"]), "its_chimeras") is None


def test_a_failing_dataset_does_not_lose_the_others(monkeypatch):
    import types
    ov = _overnight()

    def fake(ctx, name, marker, region, reads, ref, thin, ee):
        if name == "bad":
            raise RuntimeError("boom")
        return {"name": name}, [f"--- {name}"], []
    monkeypatch.setattr(ov, "battery_dataset", fake)
    ctx = types.SimpleNamespace(quick=False, args=types.SimpleNamespace(extra=["good=/g", "bad=/b", "also=/a"]))
    res = ov.step_extra_datasets(ctx)
    assert [r["name"] for r in res["data"]] == ["good", "also"]
    assert any("bad" in v["criterion"] and v["outcome"] == "NOT MET" and "boom" in v["detail"] for v in res["verdicts"]) and any("bad: battery FAILED" in ln for ln in res["lines"])


def test_the_inter_laboratory_mock_primers_trim_every_variant_in_both_orientations_and_leave_out_the_other_region(tmp_path):
    import gzip
    import random
    import shutil
    if not shutil.which("cutadapt"):
        pytest.skip("cutadapt is not installed")
    rb = load("real_data")
    cb = load("check_build")
    comp = str.maketrans("ACGT", "TGCA")
    rc = lambda s: s.translate(comp)[::-1]
    iupac = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "K": "GT", "D": "AGT"}
    variants = ["GCTCAGATTGAACGCTGGCGG", "GTTCAGATTGAACGCTGGCGG", "GCTCAGGACGAACGCTGGCGG", "GCTCAGGATTAACGCTGGCGG", "GCTCAGGATGAACGCTGGCGG", "GCTCAGAATGAACGCTGGCGG", "GCTCAGGGTGAACGCTGGCGG", "GCTCAGAGTGAACGCTGGCGG"]
    e12, e46 = cb.CATALOG["mcm_v12"], cb.CATALOG["mcm_v46"]
    assert all(len(v) == len(e12["primer_f"]) and all(b in iupac[p] for b, p in zip(v, e12["primer_f"])) for v in variants)         # the degenerate primer covers all eight of Table 3
    random.seed(1)
    raw = tmp_path / "raw"
    raw.mkdir()

    def sample(name, forward_variants, reverse_primer, n=20):
        r1, r2 = [], []
        for f in forward_variants:
            for k in range(n):
                amp = f + "".join(random.choice("ACGT") for _ in range(250)) + rc(reverse_primer)
                a, b = amp[:150], rc(amp)[:150]
                if k % 2:                                                       # a library made by adapter ligation is read from either end
                    a, b = b, a
                r1.append(a)
                r2.append(b)
        for suffix, reads in (("R1", r1), ("R2", r2)):
            with gzip.open(raw / f"{name}_{suffix}.fastq.gz", "wt") as fh:
                for i, s in enumerate(reads):
                    fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
    sample("A", variants, "ACTGCTGCCTCCCGTAGGAGT")                              # V1-V2
    sample("B", ["GTGCCAGCAGCTGCGGTAATAC"], "CACGACACGAGCTGACGACA", n=40)          # V4-V6 (the degenerate bases written out)
    for entry, kept, other, n_kept in ((e12, "A", "B", 160), (e46, "B", "A", 40)):
        out = tmp_path / f"out_{kept}"
        rb.trim_reads(raw, out, entry["primer_f"], entry["primer_r"], None, None, True, 0.15)
        counts = {s: rb.count_pairs(out, s) if hasattr(rb, "count_pairs") else None for s in (kept, other)}
        dropped = [s for s, n in rb.drop_small_samples(out, 10)]
        assert other in dropped and kept not in dropped, (entry["short"], dropped, counts)
        with gzip.open(out / f"{kept}_R1.fastq.gz", "rt") as fh:
            lengths = {len(line.strip()) for i, line in enumerate(fh) if i % 4 == 1}
        assert max(lengths) == 150 - len(entry["primer_f"]) and min(lengths) >= 100          # the primer is gone from every read; chance read-through matches trim a few more
        assert sum(1 for _ in gzip.open(out / f"{kept}_R1.fastq.gz", "rt")) // 4 == n_kept             # and every pair was kept, whichever end it was read from


def test_the_inter_laboratory_mock_regions_share_reads_and_leave_out_the_shotgun_runs(tmp_path):
    cb = load("check_build")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path)}
    a, b = cb.dataset_dirs(cfg, "mcm_v12"), cb.dataset_dirs(cfg, "mcm_v46")
    assert a["raw"] == b["raw"] and a["ready"] != b["ready"] and cb.dataset_dirs(cfg, "v4")["raw"] != b["raw"]

    class R:
        def __init__(self, text):
            self.text = text
        def __enter__(self):
            return self
        def __exit__(self, *exc):
            return False
        def read(self):
            return self.text.encode()
    tsv = ("run_accession\tsample_alias\tlibrary_layout\tinstrument_model\tread_count\nR1\ta\tPAIRED\tIllumina MiSeq\t800000\nR2\tb\tPAIRED\tIllumina NextSeq 500\t27000000\n"
           "R3\tc\tPAIRED\tIllumina MiSeq\t900000\n")
    assert cb.fetch_run_table("PRJ", tmp_path / "t", opener=lambda url, timeout=120: R(tsv), instrument="MiSeq") == ["R1", "R3"]
    assert len(cb.fetch_run_table("PRJ", tmp_path / "t2", opener=lambda url, timeout=120: R(tsv))) == 3
    # a dataset that shares reads downloads nothing when they are there
    a["raw"].mkdir(parents=True)
    (a["raw"] / "X_R1.fastq.gz").write_bytes(b"x")
    args = type("A", (), {"runs_file": None, "dry_run": False, "max_runs": 12})()
    assert cb.download(cfg, "mcm_v12", args, run=lambda cmd, **k: (_ for _ in ()).throw(AssertionError("must not run"))) == 0


def test_overnight_all_runs_every_phase_writes_the_report_and_survives_a_failed_download(tmp_path):
    import types
    cb = load("check_build")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "out": str(tmp_path / "out"), "extras": {}, "mockA_high": "", "mockA_all": "", "its": "", "v1v2": "", "figaro": "",
           "runinfo": "", "asv_table": "", "reference": ""}
    cfg_path = tmp_path / "paths.json"
    cb.save_config(cfg_path, cfg)
    calls, seen = [], {}

    def dl(cfg, name, args):
        calls.append(("download", name, args.max_runs))
        if name.startswith("mcm"):
            return 1
        d = cb.dataset_dirs(cfg, name)
        d["raw"].mkdir(parents=True, exist_ok=True)
        (d["raw"] / "S_R1.fastq.gz").write_bytes(b"x")
        return 0

    def prep(cfg, name, args, path):
        calls.append(("prepare", name, args.min_reads))
        d = cb.dataset_dirs(cfg, name)
        d["ready"].mkdir(parents=True, exist_ok=True)
        cur = cb.load_config(path)
        cur.setdefault("extras", {}).setdefault(name, {})["ready"] = str(d["ready"])
        cb.save_config(path, cur)
        return 0

    def run(cfg, args):
        seen["extras"] = dict(cfg.get("extras", {}))
        return 0
    args = types.SimpleNamespace(quick=False)
    code = cb.overnight_all(cfg, args, cfg_path, download_fn=dl, prepare_fn=prep, run_fn=run, machine_fn=lambda: [("x", True, "ok", "")], discover_fn=lambda c: {})
    assert code == 0 and ("download", "v4", 36) in calls and ("prepare", "v4", 0) in calls and not any(c[0] == "prepare" and c[1].startswith("mcm") for c in calls)
    assert list(seen["extras"]) == ["v4", "mockB"] and ("download", "mockB", 8) in calls and ("prepare", "mockB", 1000) in calls     # the datasets prepared during the run are in the run
    report = (tmp_path / "out" / "MORNING_REPORT.md").read_text(encoding="utf-8")
    assert "| dataset v4 | downloaded, trimmed and registered |" in report and "| dataset mcm_v46 | download failed" in report and "## Test tables" in report and "## What still needs you" in report
    assert (tmp_path / "out" / "overnight_all.log").exists()
    # a machine that is not ready stops before anything is fetched
    calls.clear()
    code = cb.overnight_all(cfg, args, cfg_path, download_fn=dl, prepare_fn=prep, run_fn=run, machine_fn=lambda: [("R", False, "missing", "install R")], discover_fn=lambda c: {})
    assert code == 1 and calls == [] and "NOT READY" in (tmp_path / "out" / "MORNING_REPORT.md").read_text(encoding="utf-8")


# ------------------------------------------------------------------------------------------------ references from full-length sequences, resumed downloads
def _synthetic_16s(seed, forward, reverse, amplicon_len=300):
    import random
    rnd = random.Random(seed)
    comp = str.maketrans("ACGT", "TGCA")
    bases = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    inner = bases(amplicon_len)
    return bases(200) + forward + inner + reverse.translate(comp)[::-1] + bases(400), inner


def test_a_full_length_reference_is_cut_between_the_primers_with_copies_collapsed_and_scored_by_mock_check(tmp_path):
    import zipfile
    mr = load("make_reference")
    mc = load("mock_check")
    comp = str.maketrans("ACGT", "TGCA")
    f_primer, r_primer = "GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"
    f_concrete, r_concrete = "GTGCCAGCAGCCGCGGTAA", "GGACTACAAGGGTATCTAAT"         # one way of reading the ambiguity codes
    seq_a, inner_a = _synthetic_16s(1, f_concrete, r_concrete)
    seq_b, inner_b = _synthetic_16s(2, f_concrete, r_concrete)
    copy_a = seq_a[:100] + ("A" if seq_a[100] != "A" else "C") + seq_a[101:]            # a second operon: it differs outside the amplicon only
    seq_c, inner_c = _synthetic_16s(3, f_concrete, r_concrete)
    seq_c_variant = seq_c.replace(inner_c, inner_c[:150] + ("A" if inner_c[150] != "A" else "C") + inner_c[151:])      # a third copy that differs inside it
    other_strand, inner_d = _synthetic_16s(4, f_concrete, r_concrete)
    yeast = "".join("ACGT"[(i * 7) % 4] for i in range(1700))                           # an 18S-like sequence with no primer site
    fasta = "\n".join([f">Escherichia_coli_16S_1\n{seq_a}", f">Escherichia_coli_16S_2\n{copy_a}", f">Bacillus_subtilis_16S_1\n{seq_b}",
                       f">Listeria_monocytogenes_16S_1\n{seq_c}", f">Listeria_monocytogenes_16S_2\n{seq_c_variant}",
                       f">Salmonella_enterica_16S_1\n{other_strand.translate(comp)[::-1]}", f">Saccharomyces_cerevisiae_18S_1\n{yeast}"]) + "\n"
    folder = tmp_path / "ZymoBIOMICS.STD.refseq.v2" / "ssrRNAs"
    folder.mkdir(parents=True)
    (folder / "all.fasta").write_text(fasta)
    with zipfile.ZipFile(tmp_path / "ref.zip", "w") as z:
        z.write(folder / "all.fasta", "ZymoBIOMICS.STD.refseq.v2/ssrRNAs/all.fasta")
    for source in (tmp_path / "ref.zip", tmp_path / "ZymoBIOMICS.STD.refseq.v2"):        # a zip and a folder give the same result
        text, lines = mr.build(mr.read_inputs([source]), f_primer, r_primer)
        names = [h[1:] for h in text.splitlines() if h.startswith(">")]
        assert names == ["Bacillus_subtilis_v1|Bacillus_subtilis|1_copies", "Escherichia_coli_v1|Escherichia_coli|2_copies", "Listeria_monocytogenes_v1|Listeria_monocytogenes|1_copies",
                         "Listeria_monocytogenes_v2|Listeria_monocytogenes|1_copies", "Salmonella_enterica_v1|Salmonella_enterica|1_copies"]
        assert any("Saccharomyces_cerevisiae_18S_1" in ln and "no site for both primers" in ln for ln in lines) and "6 cut between the primers" in "\n".join(lines)
    records = dict(zip(names, [s for s in text.splitlines() if not s.startswith(">")]))
    assert records["Escherichia_coli_v1|Escherichia_coli|2_copies"] == inner_a and records["Salmonella_enterica_v1|Salmonella_enterica|1_copies"] == inner_d      # a reverse-strand sequence is turned round
    # mock_check reads the names (variants and strains) and recovers an ASV that equals an expected amplicon
    ref = tmp_path / "reference.fasta"
    ref.write_text(text)
    asv = tmp_path / "asvs.fasta"
    asv.write_text(f">asv1;size=900\n{inner_a}\n>asv2;size=100\n{inner_b}\n")
    s = mc.summarise(asv, ref)
    assert s["references"] == 5 and s["strains_total"] == 4 and s["strains_found"] == 2 and s["categories"]["exact"][0] == 2


def test_check_build_makes_and_registers_a_reference_from_a_catalog_dataset(tmp_path):
    cb = load("check_build")
    seq, inner = _synthetic_16s(5, "AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG")
    src = tmp_path / "zymo.fasta"
    src.write_text(f">Pseudomonas_aeruginosa_16S_1\n{seq}\n")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "extras": {}}
    cfg_path = tmp_path / "paths.json"
    cb.save_config(cfg_path, cfg)
    args = type("A", (), {"primer_f": "AGAGTTTGATYMTGGCTCAG", "primer_r": "ATTACCGCGGCTGCTGG"})()
    assert cb.make_reference_for(cfg, f"mockB={src}", args, cfg_path) == 0
    reg = cb.load_config(cfg_path)["extras"]["mockB"]["reference"]
    assert Path(reg).read_text().splitlines()[1] == inner and reg.endswith("reference.fasta")
    with pytest.raises(SystemExit):                                                  # primers that do not fit anything: nothing is registered
        cb.make_reference_for(cfg, f"mockB={src}", type("A", (), {"primer_f": "CCCCCCCCCCCCCCCCCCC", "primer_r": "GGGGGGGGGGGGGGGGG"})(), cfg_path)


def test_a_resume_the_server_refuses_starts_the_file_again(tmp_path):
    import io
    import urllib.error
    de = load("download_ena")
    content = b"x" * 5000
    part = tmp_path / "run_R1.part"
    part.write_bytes(b"<html>error page</html>")                                     # a partial file that is not the beginning of the real one
    seen = []

    class Response:
        status = 200
        def __init__(self, data):
            self.buf = io.BytesIO(data)
        def __enter__(self):
            return self
        def __exit__(self, *exc):
            return False
        def read(self, n):
            return self.buf.read(n)

    def opener(request, timeout=60):
        seen.append(request.get_header("Range"))
        if request.get_header("Range"):
            raise urllib.error.HTTPError(request.full_url, 416, "Range Not Satisfiable", {}, None)
        return Response(content)
    de.fetch_resumable("https://x/run_R1.fastq.gz", part, len(content), opener)
    assert part.read_bytes() == content and seen[0].startswith("bytes=") and seen[1] is None


def test_an_overnight_dataset_with_some_runs_downloaded_goes_on_with_them(tmp_path):
    cb = load("check_build")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "extras": {}}
    cfg_path = tmp_path / "paths.json"
    cb.save_config(cfg_path, cfg)
    d = cb.dataset_dirs(cfg, "mcm_v46")
    d["raw"].mkdir(parents=True)
    (d["raw"] / "ERR1_R1.fastq.gz").write_bytes(b"x")
    calls = []
    out = cb.auto_dataset(cfg, "mcm_v46", cfg_path, lambda c, n, a: calls.append("download") or 1, lambda c, n, a, p: calls.append("prepare") or 0)
    assert calls == ["download", "prepare"] and out.startswith("downloaded, trimmed and registered") and "not complete" in out        # the missing runs are retried, the rest is used
    assert cb.auto_dataset(cfg, "v4", cfg_path, lambda c, n, a: 1, lambda c, n, a, p: 0).startswith("download failed")             # nothing downloaded: nothing to prepare


def test_a_reference_zip_with_genomes_and_genes_reads_only_the_genes(tmp_path):
    import zipfile
    mr = load("make_reference")
    f_primer, r_primer = "AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG"
    gene, inner = _synthetic_16s(7, f_primer, r_primer)
    genome = "ACGT" * 5000                                                           # 20,000 nt: a genome, far longer than a gene
    with zipfile.ZipFile(tmp_path / "D6321.refseq.zip", "w") as z:
        z.writestr("D6321.refseq/Genomes/Imtechella_halotolerans.fasta", f">NZ_genome\n{genome}\n")
        z.writestr("D6321.refseq/ssrRNAs/Imtechella_halotolerans_16S.fasta", f">Imtechella_halotolerans_16S_1\n{gene}\n")
    recs = mr.read_inputs([tmp_path / "D6321.refseq.zip"])
    assert [h for _, h, _ in recs] == ["Imtechella_halotolerans_16S_1"]                # the genome file is not read at all
    text, lines = mr.build(recs, f_primer, r_primer)
    assert text.splitlines()[1] == inner and text.startswith(">Imtechella_halotolerans_v1|Imtechella_halotolerans|1_copies")
    # with only a genome in the input, the record is refused by its length and not scanned
    text2, lines2 = mr.build([("g.fasta", "NZ_genome", genome)], f_primer, r_primer)
    assert text2 == "" and any("a genome or contig" in ln for ln in lines2)


def test_datasets_that_were_dropped_stay_out_of_the_table_and_the_board(tmp_path):
    import types
    cb = load("check_build")
    cfg = {**cb.DEFAULT_CONFIG, "out": str(tmp_path / "out"), "data": str(tmp_path / "data"), "mockA_high": "", "mockA_all": "", "its": "", "v1v2": "", "figaro": "", "extras": {}}
    table = "\\n".join(cb.table_lines(cfg))
    board = "\\n".join(cb.todo(cfg, types.SimpleNamespace()))
    assert "mockB" in table and "mcm_v46" in table and "sputum" not in table
    assert "mockB" in board and "sputum" not in board and "NovaSeq" not in board and "--make-reference mockB=D6321.refseq.zip" in board
    assert all(cb.CATALOG[n].get("dropped") for n in cb.CATALOG if n.startswith("sputum"))


def test_a_reference_path_that_does_not_exist_gets_a_plain_message_with_the_files_that_were_found(tmp_path, monkeypatch):
    mr = load("make_reference")
    cb = load("check_build")
    with pytest.raises(SystemExit) as err:
        mr.read_inputs([tmp_path / "nowhere" / "D6321.refseq.zip"])
    assert "does not exist" in str(err.value)                                           # not a FileNotFoundError traceback
    data = tmp_path / "data"
    (data / "downloads").mkdir(parents=True)
    (data / "downloads" / "D6321.refseq.zip").write_bytes(b"x")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("USERPROFILE", str(tmp_path / "home"))
    cfg = {**cb.DEFAULT_CONFIG, "data": str(data), "extras": {}}
    cfg_path = tmp_path / "paths.json"
    cb.save_config(cfg_path, cfg)
    with pytest.raises(SystemExit) as err:
        cb.make_reference_for(cfg, "mockB=C:\\path\\D6321.refseq.zip", type("A", (), {"primer_f": None, "primer_r": None})(), cfg_path)
    message = str(err.value)
    assert "placeholder" in message and "D6321.refseq.zip" in message and str(data / "downloads") in message


def test_compare_reference_says_which_strains_a_published_reference_shares_with_the_expected_sequences(tmp_path, capsys):
    import zipfile
    cb = load("check_build")
    f_primer, r_primer = "AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG"
    full_a, inner_a = _synthetic_16s(11, f_primer, r_primer)
    full_b, inner_b = _synthetic_16s(12, f_primer, r_primer)
    full_d, inner_d = _synthetic_16s(14, f_primer, r_primer)
    _full_c, inner_c = _synthetic_16s(13, f_primer, r_primer)
    ref = tmp_path / "mockA_reference.fasta"                                         # the expected sequences in use: strains S1 (A), S2 (B) and S3 (C, which the published file lacks)
    ref.write_text(f">v1|S1|a\n{inner_a}\n>v2|S2|b\n{inner_b}\n>v3|S3|c\n{inner_c}\n")
    (tmp_path / "mockA_reference_map.tsv").write_text("asv\tstrain_id\tspecies\nv1\tS1\tEscherichia coli\nv2\tS2\tBacillus subtilis\nv3\tS3\tCutibacterium acnes\n")
    with zipfile.ZipFile(tmp_path / "zymo.zip", "w") as z:                           # the published file has A, B and a strain D that the expected sequences lack
        z.writestr("zymo/ssrRNAs/all.fasta", f">Escherichia_coli_16S_1\n{full_a}\n>Bacillus_subtilis_16S_1\n{full_b}\n>Salmonella_enterica_16S_1\n{full_d}\n")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "reference": str(ref), "extras": {}}
    assert cb.compare_reference_for(cfg, f"mockA={tmp_path / 'zymo.zip'}", type("A", (), {"primer_f": f_primer, "primer_r": r_primer})()) == 0
    out = capsys.readouterr().out
    assert "3 sequences of 3 strains in the candidate reference; 2 of them match" in out and "3 expected variants of 3 strains" in out and "2 are matched" in out
    assert "strains of the candidate with no counterpart: Salmonella_enterica" in out
    assert "no counterpart in the candidate: S3 (Cutibacterium acnes)" in out and "different communities, or the same one only in part" in out
    assert (tmp_path / "data" / "mockA" / "reference_check" / "zymo_amplicons.fasta").is_file()
    with pytest.raises(SystemExit) as err:                                           # nothing registered to compare with
        cb.compare_reference_for({**cfg, "reference": ""}, f"mockA={tmp_path / 'zymo.zip'}", type("A", (), {"primer_f": f_primer, "primer_r": r_primer})())
    assert "no expected sequences are registered" in str(err.value)


def test_a_sequence_that_cannot_be_cut_says_which_primer_is_missing_and_odd_names_come_from_the_file(tmp_path, capsys):
    mr = load("make_reference")
    f_primer, r_primer = "AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG"
    full, inner = _synthetic_16s(21, f_primer, r_primer)
    starts_inside = full[full.index(inner):]                                              # a sequence that begins after the forward primer site
    nothing = "".join("ACGT"[(i * 5) % 4] for i in range(900))
    text, lines = mr.build([("Pseudomonas_aeruginosa_16S.fasta", "Pseudomonas_aeruginosa_16S_1", starts_inside), ("yeast.fasta", "Saccharomyces_18S_1", nothing), ("ok.fasta", "Bacillus_16S_1", full)], f_primer, r_primer)
    joined = "\n".join(lines)
    assert "Pseudomonas_aeruginosa_16S_1: no site for both primers: forward primer not found (best site" in joined and "reverse primer found" in joined and "may start or end inside the amplicon" in joined
    assert "Saccharomyces_18S_1: no site for both primers: forward primer not found" in joined and "reverse primer not found" in joined
    assert "[from Bacillus_16S_1]" in joined and text.count(">") == 1
    # a header that names no organism: the file does
    assert mr.strain_name("Consensus_Imtechella_halotolerans_16S_1", "all.fasta") == "Imtechella_halotolerans"
    assert mr.strain_name("Consensus", "Allobacillus_halotolerans_16S.fasta") == "Allobacillus_halotolerans"                    # a header that names no organism takes the cleaned name of its file, so it joins the strain of that file
    assert mr.strain_name("Consensus_16S_1", "Truepera_radiovictrix.fasta") == "Truepera_radiovictrix"
    assert mr.strain_name("Escherichia_coli_16S_3", "ssrRNA.fasta") == "Escherichia_coli"
    # --list shows what an archive holds, without primers
    src = tmp_path / "x.fasta"
    src.write_text(">Consensus_Imtechella_16S\nACGTACGT\n")
    assert mr.main([str(src), "--list"]) == 0
    assert "Consensus_Imtechella_16S\t8 nt\t-> strain Imtechella" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        mr.main([str(src)])                                                                 # without --list the primers and --out are needed


def test_mock_b_is_fetched_by_default_as_the_eight_mock_samples_of_the_paper(tmp_path):
    cb = load("check_build")
    rows = [("ERR1", "536-BiomPsy-965"), ("ERR2", "539-BiomPsy-968"), ("ERR3", "544-BiomPsy-973"), ("ERR4", "547-BiomPsy-976"), ("ERR5", "552-BiomPsy-981"), ("ERR6", "555-BiomPsy-984"),
            ("ERR7", "560-OtherMock-001"), ("ERR8", "540-BiomPsy-969")]
    tsv = "run_accession\tsample_alias\tlibrary_layout\tinstrument_model\tread_count\n" + "".join(f"{r}\t{a}\tPAIRED\tIllumina MiSeq\t50000\n" for r, a in rows)

    class R:
        def __enter__(self):
            return self
        def __exit__(self, *exc):
            return False
        def read(self):
            return tsv.encode()
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "extras": {}}
    args = type("A", (), {"runs_file": None, "dry_run": True, "max_runs": 8, "select": None, "min_reads": 1000})()
    assert cb.download(cfg, "mockB", args, run=lambda cmd, **k: type("P", (), {"returncode": 0})(), opener=lambda url, timeout=120: R()) in (0, None)
    chosen = (tmp_path / "data" / "new" / "mockB" / "runs.txt").read_text().split()
    assert chosen == ["ERR1", "ERR2", "ERR3", "ERR4"]                                   # the two controls (552, 555), another mock and the other BiomPsy samples are left out


def test_overnight_all_can_leave_a_dataset_out(tmp_path):
    import types
    cb = load("check_build")
    cfg = {**cb.DEFAULT_CONFIG, "data": str(tmp_path / "data"), "out": str(tmp_path / "out"), "extras": {}, "mockA_high": "", "mockA_all": "", "its": "", "v1v2": "", "figaro": "",
           "runinfo": "", "asv_table": "", "reference": ""}
    cfg_path = tmp_path / "paths.json"
    cb.save_config(cfg_path, cfg)
    calls = []

    def dl(cfg, name, args):
        calls.append(name)
        d = cb.dataset_dirs(cfg, name)
        d["raw"].mkdir(parents=True, exist_ok=True)
        (d["raw"] / "S_R1.fastq.gz").write_bytes(b"x")
        return 0

    def prep(cfg, name, args, path):
        d = cb.dataset_dirs(cfg, name)
        d["ready"].mkdir(parents=True, exist_ok=True)
        cur = cb.load_config(path)
        cur.setdefault("extras", {}).setdefault(name, {})["ready"] = str(d["ready"])
        cb.save_config(path, cur)
        return 0
    seen = {}
    args = types.SimpleNamespace(quick=False, without="mockB, mcm_v12")
    cb.overnight_all(cfg, args, cfg_path, download_fn=dl, prepare_fn=prep, run_fn=lambda c, a: seen.update(extras=list(c["extras"])) or 0, machine_fn=lambda: [("x", True, "ok", "")], discover_fn=lambda c: {})
    assert calls == ["v4", "mcm_v46"] and seen["extras"] == ["v4", "mcm_v46"]
    report = (tmp_path / "out" / "MORNING_REPORT.md").read_text(encoding="utf-8")
    assert "| dataset mockB | left out (--without) |" in report and "| dataset mcm_v12 | left out (--without) |" in report


def test_a_short_error_page_instead_of_the_file_is_reported_and_not_retried(tmp_path):
    import io
    de = load("download_ena")

    class Response:
        status = 200
        def __init__(self, data):
            self.buf = io.BytesIO(data)
        def __enter__(self):
            return self
        def __exit__(self, *exc):
            return False
        def read(self, n):
            return self.buf.read(n)
    calls = []

    def opener(request, timeout=60):
        calls.append(1)
        return Response(b"<html><body>404 Not Found: the file is not available</body></html>")
    part = tmp_path / "ERR9_R1.part"
    with pytest.raises(de.ServerAnswer) as err:
        de.fetch_resumable("https://x/ERR9_1.fastq.gz", part, 78477350, opener)
    assert "404 Not Found" in str(err.value) and not part.exists()
    # a partial gzip file that is short is an interruption, not a message
    class Short(Response):
        pass
    part2 = tmp_path / "ERR8_R1.part"
    with pytest.raises(OSError) as err2:
        de.fetch_resumable("https://x/ERR8_1.fastq.gz", part2, 5000, lambda request, timeout=60: Short(b"\x1f\x8b\x08" + b"x" * 100))
    assert "incomplete" in str(err2.value) and part2.exists()
    # download() stops at the first answer of that kind and says what to do
    tsv = "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\nERR9\tftp.x/ERR9_1.fastq.gz;ftp.x/ERR9_2.fastq.gz\tabc;def\t78477350;80000000\n"

    class Listing(Response):
        def read(self, n=-1):
            return tsv.encode()
    calls.clear()
    with pytest.raises(SystemExit) as stop:
        de.download("ERR9", tmp_path / "raw", retries=8, log=lambda *a: None, opener=lambda url, timeout=60: Listing(b""), file_opener=opener)
    assert len(calls) == 1 and "not retried" in str(stop.value) and "fasterq-dump ERR9" in str(stop.value) and "404 Not Found" in str(stop.value)


def test_a_step_is_repeated_when_the_reference_it_scored_against_has_changed(tmp_path, monkeypatch):
    import json
    import types
    ov = _overnight()
    calls = []

    def fake(ctx):
        calls.append(1)
        return {"title": "t", "lines": [], "verdicts": [], "files": [], "data": []}
    monkeypatch.setattr(ov, "STEPS", [("extra_datasets", 1, 1, 1, [], [], fake), ("per_sample_rank", 1, 1, 1, [], [], fake)])
    ref = tmp_path / "mockB_reference.fasta"
    ref.write_text(">a|S1|x\nACGT\n")

    class C:
        def __init__(self):
            self.out, self.quick, self.cur = tmp_path / "out", False, None
            self.args = types.SimpleNamespace(fresh=False, extra=[f"mockB=/d/mockB|{ref}"])
        def p(self, key):
            return None
        def log(self, message):
            pass
    ov.run_steps(C(), ["extra_datasets", "per_sample_rank"])
    assert len(calls) == 2 and json.loads((tmp_path / "out" / "extra_datasets" / "result.json").read_text())["datasets"][0].startswith("mockB#")
    ov.run_steps(C(), ["extra_datasets", "per_sample_rank"])
    assert len(calls) == 2                                                                  # nothing changed: both kept
    ref.write_text(">a|S1|x\nACGT\n>b|S2|y\nTTGA\n")                                      # the reference is rebuilt
    ov.run_steps(C(), ["extra_datasets", "per_sample_rank"])
    assert len(calls) == 3                                                                  # the step that scores against it is repeated; per_sample_rank, which does not, is kept


def test_a_consensus_sequence_joins_its_strain_and_is_not_an_expected_variant(tmp_path):
    mr = load("make_reference")
    f_primer, r_primer = "AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG"
    real1, inner1 = _synthetic_16s(31, f_primer, r_primer, 300)
    real2, inner2 = _synthetic_16s(32, f_primer, r_primer, 303)                          # a second operon copy
    cons, inner_c = _synthetic_16s(33, f_primer, r_primer, 303)
    lone, inner_l = _synthetic_16s(34, f_primer, r_primer, 280)
    assert mr.strain_name("Consensus", "Allobacillus.halotolerans.fasta") == "Allobacillus_halotolerans"                 # a dot is an underscore
    assert mr.strain_name("Imtechella_halotolerans_16S_2", "x.fasta") == "Imtechella_halotolerans"
    text, lines = mr.build([("Allobacillus_halotolerans.fasta", "Allobacillus_halotolerans_16S_1", real1), ("Allobacillus_halotolerans.fasta", "Allobacillus_halotolerans_16S_2", real2),
                            ("Allobacillus.halotolerans.fasta", "Consensus", cons), ("Truepera.fasta", "Consensus", lone)], f_primer, r_primer)
    names = [h[1:] for h in text.splitlines() if h.startswith(">")]
    assert names == ["Allobacillus_halotolerans_v1|Allobacillus_halotolerans|1_copies", "Allobacillus_halotolerans_v2|Allobacillus_halotolerans|1_copies", "Truepera_v1|Truepera|1_copies"]
    joined = "\n".join(lines)
    assert "Consensus: a consensus of Allobacillus_halotolerans's other sequences, not an observed variant" in joined         # left out where real copies exist
    assert inner_c not in text and inner_l in text                                                                           # kept where the consensus is all there is
