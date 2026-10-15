"""validation/night6.py: what can be tested without cutadapt, vsearch or R: the table of pairs, the single-end download, loading the released package, the DADA2
command, every verdict of the pre-registered checks, the environment check, and how a run resumes and fails."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
import download_ena  # noqa: E402
import longread_describe as lr  # noqa: E402
import night6  # noqa: E402

LISTING = "run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n{acc}\t{urls}\t{md5}\t{size}\n"


class Response:
    def __init__(self, text):
        self.text = text

    def read(self):
        return self.text.encode()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def listing_opener(acc, content=b"data", urls=None, md5=None):
    urls = urls or f"ftp.sra.ebi.ac.uk/vol1/fastq/SRR233/{acc}/{acc}.fastq.gz"
    text = LISTING.format(acc=acc, urls=urls, md5=md5 or hashlib.md5(content).hexdigest(), size=len(content))
    return lambda url, timeout=60: Response(text)


@pytest.fixture(autouse=True)
def quiet(monkeypatch):
    monkeypatch.setattr(download_ena, "PAUSE", 0)
    monkeypatch.setattr(night6, "LOG", None)
    path = list(sys.path)
    ours = {k: v for k, v in sys.modules.items() if k == "amplitrim" or k.startswith("amplitrim.")}        # the development build that the other test files have imported
    yield
    sys.path[:] = path
    for m in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:                    # load_release replaces the package by the released one under the same names:
        del sys.modules[m]                                                                                  # take whatever the test loaded out, and put the development build back
    sys.modules.update(ours)
    night6.LOG = None


def test_the_table_has_the_26_pairs_of_nine_people_and_a_pilot_spreads_over_it():
    rows = night6.read_pairs(night6.PAIRS_TABLE)
    assert len(rows) == 26 and len({r["volunteer"] for r in rows}) == 9
    assert len({r["illumina_run"] for r in rows}) == 26 and len({r["pacbio_run"] for r in rows}) == 26
    assert all(r["illumina_run"].startswith("SRR") and r["pacbio_run"].startswith("SRR") for r in rows)
    assert {r["site"] for r in rows} == {"Faeces", "Saliva", "Subgingival"}
    pilot = night6.read_pairs(night6.PAIRS_TABLE, 2)
    assert pilot == [rows[0], rows[-1]]


def test_the_single_end_file_is_the_one_named_after_the_run_and_other_listings_are_refused():
    assert night6.single_file("SRR1", listing_opener("SRR1", b"abc"))[2] == 3
    url, md5, size = night6.single_file("SRR1", listing_opener("SRR1", b"abc", urls="ftp.x/y/SRR1.fastq.gz"))
    assert url == "ftp.x/y/SRR1.fastq.gz" and md5 == hashlib.md5(b"abc").hexdigest()
    assert night6.single_file("SRR1", listing_opener("SRR1", urls="ftp.x/other.fastq.gz"))[0] == "ftp.x/other.fastq.gz"          # a single file, whatever it is called
    with pytest.raises(SystemExit, match="expected one FASTQ file"):
        night6.single_file("SRR1", listing_opener("SRR1", urls="ftp.x/SRR1_1.fastq.gz;ftp.x/SRR1_2.fastq.gz"))
    with pytest.raises(SystemExit, match="no FASTQ file"):
        night6.ena_listing("SRR1", lambda url, timeout=60: Response("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"))


def test_an_unreachable_ena_is_reported_after_the_attempts():
    calls = []

    def down(url, timeout=60):
        calls.append(1)
        raise OSError("no route")
    with pytest.raises(SystemExit, match="could not ask ENA"):
        night6.ena_listing("SRR1", down, attempts=3, wait=0)
    assert len(calls) == 3


def test_a_single_end_download_is_verified_kept_and_not_repeated(tmp_path):
    content = b"@r\nACGT\n+\nIIII\n"
    sent = []

    def retriever(url, part):
        sent.append(url)
        part.write_bytes(content)
    line = night6.download_single("SRR9", tmp_path, log=lambda *_: None, opener=listing_opener("SRR9", content), retriever=retriever)
    assert line == f"{hashlib.md5(content).hexdigest()}  SRR9.fastq.gz" and (tmp_path / "SRR9.fastq.gz").read_bytes() == content and not list(tmp_path.glob("*.part"))
    assert sent == ["https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR233/SRR9/SRR9.fastq.gz"]
    night6.download_single("SRR9", tmp_path, log=lambda *_: None, opener=listing_opener("SRR9", content), retriever=retriever)
    assert len(sent) == 1                                                                                                   # the file was there with the right sum


def test_a_wrong_checksum_starts_again_and_a_dropped_connection_resumes(tmp_path):
    good = b"good data"
    attempts = []

    def retriever(url, part):
        attempts.append(1)
        if len(attempts) == 1:
            part.write_bytes(b"corrupt")                                       # attempt 1: arrives, wrong sum
        elif len(attempts) == 2:
            raise OSError("connection reset")                                  # attempt 2: dropped
        else:
            part.write_bytes(good)
    night6.download_single("SRR8", tmp_path, retries=5, log=lambda *_: None, opener=listing_opener("SRR8", good), retriever=retriever)
    assert len(attempts) == 3 and (tmp_path / "SRR8.fastq.gz").read_bytes() == good
    with pytest.raises(SystemExit, match="no verified download after 2 attempts"):
        night6.download_single("SRR7", tmp_path, retries=2, log=lambda *_: None, opener=listing_opener("SRR7", good), retriever=lambda u, p: p.write_bytes(b"always wrong"))
    assert not (tmp_path / "SRR7.part").exists() and not (tmp_path / "SRR7.fastq.gz").exists()                              # nothing wrong is kept


def make_package(folder: Path, with_cap_warning: bool = True, nested: bool = False) -> Path:
    pkg = folder / "amplitrim" / ("amplitrim" if nested else "")
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text('__version__ = "9.9.9"\n')
    (pkg / "pipeline.py").write_text("WARN_ABOVE = 5\n" if with_cap_warning else "RULE_CAP = 5\n")
    return folder


def test_the_released_package_is_found_in_the_usual_places_and_an_old_one_is_refused(tmp_path):
    root = make_package(tmp_path / "a")
    assert night6.resolve_package(root) == root and night6.resolve_package(root / "amplitrim") == root
    nested = make_package(tmp_path / "b", nested=True)
    assert night6.resolve_package(nested / "amplitrim") == nested / "amplitrim" and night6.resolve_package(tmp_path / "nowhere") is None
    pl = night6.load_release(root)
    assert pl.WARN_ABOVE == 5 and Path(pl.__file__).parent == root / "amplitrim"
    old = make_package(tmp_path / "c", with_cap_warning=False)
    with pytest.raises(SystemExit, match="not the release without a cap"):
        night6.load_release(old)
    with pytest.raises(SystemExit, match="has no 'amplitrim' package"):
        night6.load_release(tmp_path / "nowhere")


def amp_results(released=(274, 253, "4,8"), default=(270, 250, "2"), g25=(272, 252, "2,5")):
    def entry(t):
        f, r, e = t
        ef, er = (float(x) for x in (e.split(",") if "," in e else (e, e)))
        return {"trim_f": f, "trim_r": r, "max_ee": e, "max_ee_f": ef, "max_ee_r": er}
    return {"version": "0.1.0", "package": "x", "released": entry(released), "default": entry(default), "g25": entry(g25)}


def table(**rows):
    base = {"est_pass": "0.80", "pass": "0.80", "merged": "0.60", "truncF": "270", "truncR": "250"}
    return {name: {"name": name, **{**base, **vals}} for name, vals in rows.items()}


def test_the_dada2_settings_are_run_once_each_and_the_command_names_them():
    s = night6.dada2_settings(amp_results())
    assert [x[0] for x in s] == ["default", "g25", "released"] and s[2] == ("released", 274, 253, "4,8")
    same = night6.dada2_settings(amp_results(released=(270, 250, "2")))
    assert [x[0] for x in same] == ["default", "g25"]                                                                       # the rule kept the default: DADA2 runs it once
    cmd = night6.compare_command(Path("thin"), Path("out"), s, Path("figaro.py"))
    assert cmd[1].endswith("compare_methods.py") and "--setting" in cmd and "default=270,250,2" in cmd and "released=274,253,4,8" in cmd and "g25=272,252,2,5" in cmd
    assert cmd[cmd.index("--baseline") + 1] == "default" and cmd[cmd.index("--figaro-script") + 1] == "figaro.py" and "--figaro-script" not in night6.compare_command(Path("t"), Path("o"), s)


def verdicts(rows, amp):
    return {i: (v, t) for i, v, t in night6.evaluate_illumina(rows, amp)}


def test_p46_compares_the_filter_estimate_with_dada2_at_the_released_and_the_default_setting():
    amp = amp_results()
    ok = verdicts(table(default={"merged": "0.55", "est_pass": "0.81"}, g25={}, released={"merged": "0.70", "pass": "0.78", "est_pass": "0.79"}), amp)
    assert ok["P46"][0] == "MET" and "released +1.0 points" in ok["P46"][1] and "default +1.0 points" in ok["P46"][1]
    bad = verdicts(table(default={"merged": "0.55"}, released={"merged": "0.70", "est_pass": "0.86"}), amp)
    assert bad["P46"][0] == "NOT MET" and "+6.0 points" in bad["P46"][1]
    pct = verdicts(table(default={"merged": "55", "pass": "80", "est_pass": "81"}, released={"merged": "70", "pass": "78", "est_pass": "79"}), amp)                  # percents in the table
    assert pct["P46"][0] == "MET" and pct["P47"][0] == "MET"
    assert verdicts(table(default={}), amp)["P46"][0] == "not decided"


def test_p47_asks_whether_a_relaxed_rule_merges_more_than_the_default():
    amp = amp_results()
    met = verdicts(table(default={"merged": "0.55"}, released={"merged": "0.70"}), amp)["P47"]
    assert met[0] == "MET" and "+15.0 points" in met[1]
    assert verdicts(table(default={"merged": "0.55"}, released={"merged": "0.56"}), amp)["P47"][0] == "NOT MET"
    assert verdicts(table(default={"merged": "0.55"}, released={"merged": "0.52"}), amp)["P47"][0] == "NOT MET"            # worse than the default
    kept = amp_results(released=(270, 250, "2"))
    r = verdicts(table(default={"merged": "0.55"}), kept)                                                                  # DADA2 ran the setting once, as 'default'
    assert r["P47"][0] == "not applicable" and r["P46"][0] == "MET" and "kept the default" in r["P47"][1]
    assert verdicts(table(default={"merged": "0.55"}, released={"merged": ""}), amp)["P47"][0] == "not decided"
    assert verdicts({}, {"version": "0.1.0", "package": "x"})["P47"][0] == "not decided"


def test_p48_and_p49_are_descriptive_and_say_whether_the_two_builds_chose_the_same_cut():
    amp = amp_results()
    rows = table(default={"merged": "0.55"}, g25={"merged": "0.65"}, released={"merged": "0.70"})
    rows["AmpliTrim, maxEE 2"] = {"name": "AmpliTrim, maxEE 2", "truncF": "270", "truncR": "250"}
    r = verdicts(rows, amp)
    assert r["P48"][0] == "descriptive" and "default (maxEE 2) 55.0%, 2,5 65.0%, released 70.0%" in r["P48"][1]
    assert r["P49"][0] == "descriptive" and "the same cut" in r["P49"][1]
    rows["AmpliTrim, maxEE 2"]["truncF"] = "268"
    assert "DIFFERENT cuts" in verdicts(rows, amp)["P49"][1]
    assert verdicts(table(default={}), amp)["P49"][0] == "not decided"


def pacbio_row(**kw):
    row = {"sample": "SRRx", "reads": 1000, "both_forward": 900, "both_reverse": 50, "between_primers": 950, "len_median": 1450.0, "share_within_10pct": 0.99,
           "distinct_q": 60, "max_q": 93, "p_start": 0.001, "p_mid": 0.001, "p_end": 0.001, "share_both": 0.95}
    row.update(kw)
    return row


def test_p50_to_p53_are_judged_on_the_pooled_files_and_say_what_would_be_wrong():
    good = {i: v for i, v, _ in night6.evaluate_pacbio([pacbio_row(), pacbio_row(sample="SRRy")])}
    assert set(good.values()) == {"as expected"}
    bad_rows = [pacbio_row(both_forward=300, both_reverse=0, between_primers=300, distinct_q=2, max_q=40, len_median=1200.0, p_end=0.05)]
    r = {i: (v, t) for i, v, t in night6.evaluate_pacbio(bad_rows)}
    assert all(v == "NOT as expected" for v, _ in r.values())
    assert "simplified qualities" in r["P52"][1] and "30.0% of 1,000 reads" in r["P50"][1]
    assert [v for _, v, _ in night6.evaluate_pacbio([])] == ["not decided"] * 4
    one_bad = {i: v for i, v, _ in night6.evaluate_pacbio([pacbio_row(), pacbio_row(sample="SRRz", distinct_q=5)])}
    assert one_bad["P52"] == "NOT as expected" and one_bad["P50"] == "as expected"                                           # every file must keep its quality values


def test_the_environment_check_reports_each_need_without_stopping(tmp_path):
    ctx = night6.Ctx(argparse.Namespace(out=tmp_path / "o", pairs_table=night6.PAIRS_TABLE))
    ctx.package = make_package(tmp_path / "pkg")
    ctx.pairs = night6.read_pairs(night6.PAIRS_TABLE)
    res = night6.check_environment(ctx, opener=listing_opener("SRR23380961"))
    names = [w for w, _, _ in res]
    assert any(n.startswith("Python") for n in names) and "cutadapt" in names and "vsearch on the PATH" in names and any(n.startswith("free disk space") for n in names)
    by = {w: (ok, d) for w, ok, d in res}
    assert by["ENA reachable (www.ebi.ac.uk)"][0] is True and by["table of the 26 pairs"][0] is True and by["released package (--package)"][0] is True and "9.9.9" in by["released package (--package)"][1]
    down = night6.check_environment(ctx, opener=lambda url, timeout=60: (_ for _ in ()).throw(OSError("offline")))
    assert {w: ok for w, ok, _ in down}["ENA reachable (www.ebi.ac.uk)"] is False
    short = tmp_path / "short.csv"
    short.write_text("volunteer,site,illumina_run,pacbio_run\nV1,Faeces,SRR1,SRR2\n")
    ctx.args.pairs_table = short
    wrong = {w: (ok, d) for w, ok, d in night6.check_environment(ctx, opener=listing_opener("SRR1"))}["table of the 26 pairs"]
    assert wrong[0] is False and "1 rows" in wrong[1]                                                                        # a table that is not the 26 pairs is not accepted
    ctx.args.pairs_table = tmp_path / "missing.csv"
    assert {w: ok for w, ok, _ in night6.check_environment(ctx, opener=listing_opener("SRR1"))}["table of the 26 pairs"] is False


def test_a_run_resumes_records_failures_and_skips_and_goes_on(tmp_path, monkeypatch, capsys):
    calls = []

    def ok(c):
        calls.append("download")
        return "fine"

    def boom(c):
        calls.append("describe")
        raise SystemExit("the tool failed")

    def skip(c):
        calls.append("prepare")
        raise night6.Skip("the reads are missing")
    monkeypatch.setitem(night6.STEP_FUNCTIONS, "download", ok)
    monkeypatch.setitem(night6.STEP_FUNCTIONS, "describe", boom)
    monkeypatch.setitem(night6.STEP_FUNCTIONS, "prepare", skip)
    out = tmp_path / "run"
    argv = ["--out", str(out), "--only", "download,describe,prepare", "--package", str(tmp_path)]
    assert night6.main(argv) == 1                                                                   # one step failed
    state = json.loads((out / "night6_state.json").read_text())["steps"]
    assert state["download"]["status"] == "done" and state["describe"]["status"] == "failed" and state["prepare"]["status"] == "skipped"
    assert "the tool failed" in state["describe"]["detail"] and calls == ["download", "describe", "prepare"]
    assert "download: finished earlier" not in capsys.readouterr().out
    night6.main(argv)
    assert calls == ["download", "describe", "prepare", "describe", "prepare"]                    # the finished step was not repeated; the others were tried again
    assert "download: finished earlier, skipped" in capsys.readouterr().out
    night6.main(argv + ["--fresh"])
    assert calls.count("download") == 2
    assert "night6.log" and (out / "night6.log").is_file() and "the tool failed" in (out / "night6.log").read_text()


def test_a_crash_inside_a_step_is_recorded_and_the_next_step_still_runs(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setitem(night6.STEP_FUNCTIONS, "download", lambda c: 1 / 0)
    monkeypatch.setitem(night6.STEP_FUNCTIONS, "describe", lambda c: seen.append("next") or "ok")
    out = tmp_path / "run"
    assert night6.main(["--out", str(out), "--only", "download,describe", "--package", str(tmp_path)]) == 1
    state = json.loads((out / "night6_state.json").read_text())["steps"]
    assert state["download"]["status"] == "failed" and "ZeroDivisionError" in state["download"]["detail"] and seen == ["next"]


def test_list_check_and_unknown_steps(tmp_path, capsys):
    assert night6.main(["--out", str(tmp_path / "o"), "--list", "--package", str(tmp_path)]) == 0
    listing = capsys.readouterr().out
    assert all(s in listing for s in night6.STEPS) and "[not run]" in listing
    assert night6.main(["--out", str(tmp_path / "o"), "--only", "download,bogus", "--package", str(tmp_path)]) == 2
    assert "Unknown step" in capsys.readouterr().out
    code = night6.main(["--out", str(tmp_path / "o2"), "--check", "--package", str(tmp_path)])
    text = capsys.readouterr().out
    assert code in (0, 1) and "released package (--package)" in text and ("thing(s) to fix" in text or "Everything is ready" in text)
    assert night6.main(["--out", str(tmp_path / "o3"), "--pairs-table", str(tmp_path / "missing.csv")]) == 2


def test_the_summary_applies_the_checks_to_what_the_steps_left_behind(tmp_path):
    out = tmp_path / "run"
    (out / "amplitrim").mkdir(parents=True)
    (out / "dada2").mkdir()
    (out / "amplitrim" / "amplitrim_results.json").write_text(json.dumps(amp_results()))
    rows = table(default={"merged": "0.55", "est_pass": "0.81"}, g25={"merged": "0.65"}, released={"merged": "0.70", "est_pass": "0.79", "pass": "0.78"})
    with open(out / "dada2" / "compare_methods.tsv", "w", newline="") as fh:
        keys = ["name", "truncF", "truncR", "maxEE", "est_pass", "pass", "merged"]
        fh.write("\t".join(keys) + "\n")
        for r in rows.values():
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
    # a PacBio folder described for real (the module under test above), then summarised
    import random
    random.seed(1)
    fastq = tmp_path / "pb"
    fastq.mkdir()
    FWD, REV = night6.PACBIO_PRIMERS
    import gzip
    with gzip.open(fastq / "SRR1.fastq.gz", "wt") as fh:
        for i in range(30):
            insert = "".join(random.choice("ACGT") for _ in range(1458))                                   # as long as real full-length 16S reads between 27F and 1492R
            seq = "".join(random.choice(lr.IUPAC[c]) for c in FWD) + insert + lr.revcomp("".join(random.choice(lr.IUPAC[c]) for c in REV))
            fh.write(f"@r{i}\n{seq}\n+\n{''.join(chr(33 + random.randint(10, 93)) for _ in seq)}\n")
    lr.describe_folder(fastq, out / "longread", FWD, REV, log=lambda *_: None)
    ctx = night6.Ctx(argparse.Namespace(out=out, pairs_per_sample=5000))
    ctx.state["steps"] = {"download": {"status": "done", "seconds": 600, "detail": "52 runs downloaded and checked"}}
    night6.step_summary(ctx)
    md = (out / "NIGHT6_SUMMARY.md").read_text()
    for needle in ("P46** MET", "P47** MET", "P48** descriptive", "P49** not decided", "P50** as expected", "P51** as expected", "P52** as expected", "P53** as expected",
                   "released: truncLenF 274, truncLenR 253, maxEE 4,8", "download: done in 10 min", "prepare: not run"):
        assert needle in md, needle
    empty = tmp_path / "empty"
    night6.step_summary(night6.Ctx(argparse.Namespace(out=empty, pairs_per_sample=5000)))
    assert "Not decided: the AmpliTrim step has not finished." in (empty / "NIGHT6_SUMMARY.md").read_text() and "Not decided: the longread step has not finished." in (empty / "NIGHT6_SUMMARY.md").read_text()


def test_the_released_package_is_driven_through_its_api_and_writes_its_report(tmp_path, monkeypatch):
    """The real package next to this archive (as on the researcher's machine), with a stand-in vsearch: the settings it is given, the results read back, the report written."""
    package = night6.ROOT.parent / "amplitrim"
    if not (package / "amplitrim" / "pipeline.py").is_file() or not (package / "tests" / "fake_vsearch.py").is_file():
        pytest.skip("the clone of the software repository is not next to this archive")
    pl = night6.load_release(package)
    sys.path.insert(0, str(package / "tests"))
    from fake_vsearch import install, reads_with_a_bad_reverse_tail
    install(monkeypatch, pl)
    reads = reads_with_a_bad_reverse_tail(tmp_path / "reads", share_good=0.4)
    out = night6.run_amplitrim(pl, reads, tmp_path / "work", tmp_path / "report", say=lambda *_: None)
    assert out["version"] and Path(out["package"]) == package / "amplitrim"
    for name in ("released", "default", "g25"):
        s = out[name]
        assert isinstance(s["trim_f"], int) and isinstance(s["trim_r"], int) and s["max_ee_f"] >= 2 and s["max_ee_r"] >= 2, (name, s)
    assert out["default"]["max_ee"] == "2" and out["g25"]["max_ee"] == "2,5"                    # asked for exactly these
    assert out["released"]["max_ee_f"] >= 2 and out["released"]["max_ee_r"] >= 2
    assert (tmp_path / "report" / "REPORT.md").is_file() and (tmp_path / "report" / "summary.json").is_file()
    assert json.loads((tmp_path / "report" / "summary.json").read_text())["recommended"]["truncLenF"] == out["released"]["trim_f"]       # the report is the released run's
    settings = night6.dada2_settings(out)
    assert settings and settings[0][0] == "default" and all(len(x) == 4 for x in settings)


def write_table(path: Path, rows: list) -> None:
    keys = ["name", "truncF", "truncR", "maxEE", "est_pass", "pass", "merged", "error"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join("\t".join(str(r.get(k, "")) for k in keys) for r in [dict(zip(keys, keys))] + rows) + "\n")


def dada2_context(tmp_path, monkeypatch, rows):
    """A run folder with the thinned reads and the released package's results, and a stand-in for the compare_methods.py run that writes `rows`."""
    ctx = night6.Ctx(argparse.Namespace(out=tmp_path / "run", pairs_per_sample=5000, figaro_script=None))
    ctx.thin.mkdir(parents=True)
    (ctx.thin / "SRR1_R1.fastq.gz").write_bytes(b"x")
    (ctx.amp_dir).mkdir(parents=True)
    (ctx.amp_dir / "amplitrim_results.json").write_text(json.dumps(amp_results(released=(270, 250, "2"))))      # the rule kept the default: DADA2 runs two settings
    ctx.out.mkdir(parents=True, exist_ok=True)

    def fake_run(cmd, stdout=None, stderr=None, cwd=None, **kw):
        write_table(ctx.dada2_dir / "compare_methods.tsv", rows)
        return type("R", (), {"returncode": 0})()
    monkeypatch.setattr(night6.subprocess, "run", fake_run)
    return ctx


def test_dada2_that_ran_without_giving_results_is_a_failure_and_says_why_not_done(tmp_path, monkeypatch):
    failed = [{"name": "AmpliTrim (maxEE 2)", "truncF": 270, "truncR": 250, "maxEE": "2", "est_pass": 0.8, "error": "DADA2 could not run: filterAndTrim removed every read"},
              {"name": "default", "truncF": 270, "truncR": 250, "maxEE": "2", "est_pass": 0.8, "error": "DADA2 could not run: filterAndTrim removed every read"},
              {"name": "g25", "truncF": 272, "truncR": 252, "maxEE": "2,5", "est_pass": 0.9, "error": "DADA2 could not run: filterAndTrim removed every read"}]
    ctx = dada2_context(tmp_path, monkeypatch, failed)
    with pytest.raises(SystemExit, match="DADA2 gave no result for 2 of 2 setting"):
        night6.step_dada2(ctx)
    try:
        night6.step_dada2(ctx)
    except SystemExit as err:
        assert "filterAndTrim removed every read" in str(err) and "--redo dada2" in str(err)
    good = [dict(failed[0], error="", **{"pass": 0.8, "merged": 0.7}), dict(failed[1], error="", **{"pass": 0.8, "merged": 0.7}), dict(failed[2], error="", **{"pass": 0.9, "merged": 0.75})]
    ctx2 = dada2_context(tmp_path / "second", monkeypatch, good)
    assert "DADA2 ran at 2 setting(s)" in night6.step_dada2(ctx2)
    half = [dict(good[1]), dict(good[2], error="DADA2 could not run: out of memory", **{"pass": "", "merged": ""})]
    with pytest.raises(SystemExit, match="1 of 2 setting"):
        night6.step_dada2(dada2_context(tmp_path / "third", monkeypatch, half))


def test_the_problems_name_a_missing_row_an_error_and_a_missing_share():
    rows = {"a": {"name": "a", "pass": "0.8", "merged": "0.7", "error": ""}, "b": {"name": "b", "error": "DADA2 could not run: boom"}, "c": {"name": "c", "pass": "", "merged": "0.7", "error": ""}}
    assert night6.dada2_problems(rows, ["a"]) == []
    assert night6.dada2_problems(rows, ["a", "b", "c", "d"]) == ["b: DADA2 could not run: boom", "c: DADA2 gave no share of pairs passing or merged", "d: not in the table"]


def test_the_summary_says_dada2_gave_no_result_instead_of_printing_dashes(tmp_path):
    out = tmp_path / "run"
    (out / "amplitrim").mkdir(parents=True)
    (out / "amplitrim" / "amplitrim_results.json").write_text(json.dumps(amp_results(released=(270, 250, "2"))))
    write_table(out / "dada2" / "compare_methods.tsv", [{"name": "default", "truncF": 270, "truncR": 250, "maxEE": "2", "est_pass": 0.8, "error": "DADA2 could not run: no reads left"}])
    night6.step_summary(night6.Ctx(argparse.Namespace(out=out, pairs_per_sample=5000)))
    md = (out / "NIGHT6_SUMMARY.md").read_text()
    assert "DADA2 gave no result for 2 setting(s)" in md and "no reads left" in md and "g25: not in the table" in md


def test_redo_repeats_only_the_named_steps_and_leaves_the_others_as_they_were(tmp_path, monkeypatch, capsys):
    calls = []
    for step in ("download", "describe", "prepare", "amplitrim"):
        monkeypatch.setitem(night6.STEP_FUNCTIONS, step, lambda c, s=step: calls.append(s) or "ok")
    out = tmp_path / "run"
    base = ["--out", str(out), "--package", str(tmp_path)]
    assert night6.main(base + ["--only", "download,describe,prepare,amplitrim"]) == 0
    calls.clear()
    assert night6.main(base + ["--redo", "prepare"]) == 0
    assert calls == ["prepare"]
    assert "amplitrim" in capsys.readouterr().out.split("Note:")[1]                          # the step that uses it is named, and not repeated
    recorded = json.loads((out / "night6_state.json").read_text())["steps"]
    assert "dada2" not in recorded and "longread" not in recorded                             # steps that were never run are not started by --redo
    calls.clear()
    assert night6.main(base + ["--redo", "prepare,amplitrim"]) == 0 and calls == ["prepare", "amplitrim"]
    assert night6.main(base + ["--redo", "nonsense"]) == 2
    assert (out / "NIGHT6_SUMMARY.md").is_file()


DESCRIBE = """sample                      reads  len Qvals  R1 primer (share, offset)    R2 primer (share, offset)    reversed  R1 starts with (share)
SRR23380921               113,099  300     {a}  341F (92%, 0)                806R (97%, -1)                     0%  CCTACGGGTGGCTGCAGTGG (14%)
SRR23380922                87,263  300     {b}  341F (87%, 0)                806R (95%, -1)                     0%  CCTACGGGTGGCAGCAGTGG (12%)

- 2 of 2 samples have 2 or fewer distinct quality values: the per-base quality is missing.
"""


def test_the_quality_values_are_read_from_the_describe_file_and_simplified_ones_are_recognised(tmp_path):
    f = tmp_path / "describe_illumina.txt"
    f.write_text(DESCRIBE.format(a=2, b=2))
    assert night6.illumina_quality_values(f) == [2, 2] and night6.simplified_qualities([2, 2]) and night6.simplified_qualities([1])
    f.write_text(DESCRIBE.format(a=2, b=38))
    assert night6.illumina_quality_values(f) == [2, 38] and not night6.simplified_qualities([2, 38]) and not night6.simplified_qualities([4, 5])      # binned scores are not 'simplified'
    assert night6.illumina_quality_values(tmp_path / "missing.txt") == [] and not night6.simplified_qualities([])
    f.write_text("nothing like a table\n")
    assert night6.illumina_quality_values(f) == []


def test_dada2_is_skipped_with_the_reason_when_every_sample_has_simplified_qualities(tmp_path, monkeypatch):
    ctx = dada2_context(tmp_path, monkeypatch, [])
    (ctx.out / "describe_illumina.txt").write_text(DESCRIBE.format(a=2, b=2))
    monkeypatch.setattr(night6.subprocess, "run", lambda *a, **k: (_ for _ in ()).throw(AssertionError("DADA2 must not be started")))
    with pytest.raises(night6.Skip, match="simplified-quality copy.*SRA 'lite'"):
        night6.step_dada2(ctx)


def test_the_summary_states_the_quality_values_and_that_no_verdict_can_be_drawn_from_simplified_ones(tmp_path):
    out = tmp_path / "run"
    (out / "amplitrim").mkdir(parents=True)
    (out / "amplitrim" / "amplitrim_results.json").write_text(json.dumps(amp_results(released=(270, 250, "2"))))
    (out / "describe_illumina.txt").write_text(DESCRIBE.format(a=2, b=2))
    write_table(out / "dada2" / "compare_methods.tsv", [{"name": "default", "truncF": 270, "truncR": 250, "maxEE": "2", "est_pass": 0.8, "error": "DADA2 could not run"}])      # left by an earlier try
    night6.step_summary(night6.Ctx(argparse.Namespace(out=out, pairs_per_sample=5000)))
    md = (out / "NIGHT6_SUMMARY.md").read_text()
    assert "Distinct quality values per sample in the Illumina reads: 2 to 2 (2 samples)." in md
    assert "P46 to P49 cannot be decided on this copy of the reads" in md and "no more informative than they are" in md
    assert "**P46** not decided" not in md and "DADA2 gave no result for" not in md                    # the verdict lines are not printed for data that cannot give them
    (out / "describe_illumina.txt").write_text(DESCRIBE.format(a=38, b=39))
    night6.step_summary(night6.Ctx(argparse.Namespace(out=out, pairs_per_sample=5000)))
    md = (out / "NIGHT6_SUMMARY.md").read_text()
    assert "38 to 39 (2 samples)" in md and "cannot be decided on this copy" not in md and "DADA2 gave no result for" in md       # normal scores: the usual verdicts and problems
