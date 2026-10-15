"""validation/night7.py: the download with its checks, the steps and their dependencies, the state file, --redo, --check and the summary.
No network and no external tool is used: NCBI is a fake opener, the L4 test and the tools are stand-ins."""
import ast
import gzip
import hashlib
import io
import json
import random
import sys
import urllib.error
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
import genome_amplicons as ga  # noqa: E402
import longread_describe as lr  # noqa: E402
import longread_filter_check as lf  # noqa: E402
import night7 as n7  # noqa: E402

ROW = {"atcc_number": "ATCC 700926", "name_on_atcc_page": "Escherichia coli", "share_percent": "18", "accession": "GCF_000005845.2", "assembly_name": "ASM584v2", "tier": "A"}
BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/"
REMOTE = "GCF_000005845.2_ASM584v2_genomic.fna.gz"


class Response(io.BytesIO):
    status = 200

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def fake_ncbi(files: dict, log: list = None):
    """An opener that serves the given {url: bytes} and answers 404 for anything else."""
    def opener(request, timeout=None):
        url = getattr(request, "full_url", request)
        if log is not None:
            log.append(url)
        if url in files:
            return Response(files[url])
        raise urllib.error.URLError("404 not found")
    return opener


def gz_bytes(text: str = ">c\nACGT\n") -> bytes:
    return gzip.compress(text.encode())


def md5_text(data: bytes) -> bytes:
    return f"{hashlib.md5(data).hexdigest()}  ./{REMOTE}\n".encode()


def make_ctx(tmp_path, **kw):
    args = n7.argparse.Namespace(out=tmp_path / "out", data=tmp_path / "data", table=tmp_path / "t.tsv", reads=50000, reads_monomer=0, tolerance=0.02, retries=3, fresh=False)
    for k, v in kw.items():
        setattr(args, k, v)
    return n7.Ctx(args)


def write_table(path: Path, rows=(ROW,)) -> Path:
    cols = ["atcc_number", "name_on_atcc_page", "share_percent", "accession", "assembly_name", "tier"]
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n" for r in rows))
    return path


def test_the_folder_rule_and_the_file_names_of_the_real_table():
    assert n7.genome_folder_url("GCF_000005845.2", "ASM584v2") == BASE
    rows = n7.read_genome_rows(n7.GENOME_TABLE)
    assert len(rows) == 20 and abs(sum(float(r["share_percent"]) for r in rows) - 100) < 1e-9
    assert n7.genome_file_name([r for r in rows if r["accession"] == "GCF_000005845.2"][0]) == "Escherichia_coli_ATCC700926_GCF_000005845.2.fna.gz"
    assert all(r["assembly_name"] for r in rows) and len({r["accession"] for r in rows}) == 20


def test_a_genome_is_downloaded_checked_against_the_md5_sum_and_not_fetched_again(tmp_path):
    data = gz_bytes()
    seen = []
    opener = fake_ncbi({BASE + "md5checksums.txt": md5_text(data), BASE + REMOTE: data}, seen)
    assert n7.download_genome(ROW, tmp_path, 3, lambda *_: None, opener) == "downloaded"
    assert (tmp_path / "Escherichia_coli_ATCC700926_GCF_000005845.2.fna.gz").read_bytes() == data and not list(tmp_path.glob("*.part"))
    seen.clear()
    assert n7.download_genome(ROW, tmp_path, 3, lambda *_: None, opener) == "present" and seen == []     # a valid file is left alone, with no request at all


def test_a_wrong_checksum_is_retried_and_then_reported(tmp_path):
    data = gz_bytes()
    opener = fake_ncbi({BASE + "md5checksums.txt": b"0" * 32 + f"  ./{REMOTE}\n".encode(), BASE + REMOTE: data})
    with pytest.raises(SystemExit, match="no verified download after 2 attempts"):
        n7.download_genome(ROW, tmp_path, 2, lambda *_: None, opener)
    assert not list(tmp_path.glob("*.fna.gz")) and not list(tmp_path.glob("*.part"))


def test_without_ncbis_checksum_file_the_gzip_check_decides(tmp_path):
    good = gz_bytes()
    messages = []
    assert n7.download_genome(ROW, tmp_path, 2, messages.append, fake_ncbi({BASE + REMOTE: good})) == "downloaded"
    assert any("checksum could not be read" in m for m in messages)
    other = tmp_path / "second"
    with pytest.raises(SystemExit, match="no verified download"):                                      # an error page saved as the file does not read as gzip
        n7.download_genome(ROW, other, 2, lambda *_: None, fake_ncbi({BASE + REMOTE: b"<html>moved</html>"}))
    assert not list(other.glob("*.fna.gz"))


def test_a_truncated_file_that_is_already_there_is_downloaded_again(tmp_path):
    data = gz_bytes("".join(f">c{i}\n{'ACGT' * 500}\n" for i in range(200)))
    (tmp_path / "Escherichia_coli_ATCC700926_GCF_000005845.2.fna.gz").write_bytes(data[: len(data) // 2])
    opener = fake_ncbi({BASE + "md5checksums.txt": md5_text(data), BASE + REMOTE: data})
    assert n7.download_genome(ROW, tmp_path, 3, lambda *_: None, opener) == "downloaded"
    assert n7.gzip_ok(tmp_path / "Escherichia_coli_ATCC700926_GCF_000005845.2.fna.gz")


FWD, REV = n7.PRIMERS


def synthetic_genome() -> str:
    rng = random.Random(5)
    bg = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    amp = "".join(rng.choice("ACGT") for _ in range(1450))
    site = lambda pat: "".join(rng.choice(lr.IUPAC[c]) for c in pat)
    return bg(2000) + site(FWD) + amp + site(lr.revcomp(REV)) + bg(2000)


def test_the_amplicon_step_waits_for_the_genomes_and_then_runs(tmp_path):
    c = make_ctx(tmp_path)
    write_table(c.table)
    with pytest.raises(n7.Skip, match="run the genomes step first"):
        n7.step_amplicons(c)
    c.genomes.mkdir(parents=True)
    with gzip.open(c.genomes / n7.genome_file_name(ROW), "wt") as fh:
        fh.write(">chr\n" + synthetic_genome() + "\n")
    detail = n7.step_amplicons(c)
    assert "1 of 1 genomes give an amplicon" in detail and (c.amplicons / "EXPECTED_AMPLICONS.md").is_file()


def fake_l4(monkeypatch, calls: list, met: bool = True):
    monkeypatch.setattr(lf, "check_tools", lambda: None)

    def check_folder(folder, out, primer_f, primer_r, grid, max_reads, error_rate, tolerance, keep, log):
        calls.append({"files": sorted(p.name for p in Path(folder).iterdir()), "reads": max_reads, "tolerance": tolerance, "grid": tuple(grid)})
        out.mkdir(parents=True, exist_ok=True)
        diff = 0.001 if met else 0.05
        (out / "filter_check.tsv").write_text("sample\tmaxEE\testimate\tdada2\tdifference\twithin_tolerance\treads_in\treads_out\n"
                                              f"s1\t2\t0.8\t{0.8 + diff}\t{diff}\t{met}\t100\t80\n")
        (out / "FILTER_CHECK.md").write_text(f"# t\n\n**Gate (L4): {'MET' if met else 'NOT MET'}.** The largest difference is {100 * diff:.2f} points (s1, maxEE 2).\n\n"
                                             "| sample | maxEE | estimate % | DADA2 % | difference (points) | within |\n|---|---|---|---|---|---|\n| s1 | 2 | 80.0 | 80.1 | +0.10 | yes |\n\n"
                                             "## Reads that survive the primer step (description only, no gate)\n\n| sample | reads | a |\n|---|---|---|\n| s1 | 100 | 99 |\n")
    monkeypatch.setattr(lf, "check_folder", check_folder)


def test_the_monomer_step_copies_exactly_the_three_files_and_passes_its_settings(tmp_path, monkeypatch):
    c = make_ctx(tmp_path)
    c.monomer.mkdir(parents=True)
    names = [f"m64004e.Kinnex16S_{p}.hifi_reads.fastq.gz" for p in n7.MONOMER_PAIRS] + ["m64004e.Kinnex16S_Fwd_01--Kinnex16S_Rev_13.hifi_reads.fastq.gz", "notes.txt"]
    for nm in names:
        (c.monomer / nm).write_bytes(b"x")
    calls = []
    fake_l4(monkeypatch, calls)
    detail = n7.step_filter_monomer(c)
    assert calls[0]["files"] == sorted(names[:3]) and calls[0]["reads"] == 0 and calls[0]["tolerance"] == 0.02 and calls[0]["grid"] == (1, 2, 3, 5, 10)
    assert detail.startswith("gate L4 MET")
    assert not c.monomer_in.exists()                                                                   # the input folder is removed afterwards
    assert all((c.monomer / nm).exists() for nm in names)                                              # the data themselves are untouched


def test_the_monomer_step_skips_when_a_file_is_missing(tmp_path, monkeypatch):
    c = make_ctx(tmp_path)
    c.monomer.mkdir(parents=True)
    for p in n7.MONOMER_PAIRS[:2]:
        (c.monomer / f"x.Kinnex16S_{p}.hifi_reads.fastq.gz").write_bytes(b"x")
    fake_l4(monkeypatch, [])
    with pytest.raises(n7.Skip, match="only 2 of the 3"):
        n7.step_filter_monomer(c)
    with pytest.raises(n7.Skip, match="no FASTQ files"):
        n7.step_filter_kinnex(c)


def test_the_kinnex_step_uses_only_the_first_reads_and_reports_a_failed_gate(tmp_path, monkeypatch):
    c = make_ctx(tmp_path)
    c.kinnex.mkdir(parents=True)
    (c.kinnex / "a.fastq.gz").write_bytes(b"x")
    calls = []
    fake_l4(monkeypatch, calls, met=False)
    detail = n7.step_filter_kinnex(c)
    assert calls[0]["reads"] == 50000 and calls[0]["files"] == ["a.fastq.gz"]
    assert detail.startswith("gate L4 NOT MET") and "5.00 points" in detail


def run_main(tmp_path, *extra):
    table = write_table(tmp_path / "t.tsv")
    return n7.main(["--out", str(tmp_path / "out"), "--data", str(tmp_path / "data"), "--table", str(table), *extra])


def test_a_failed_step_does_not_stop_the_others_and_the_state_is_kept(tmp_path, monkeypatch):
    monkeypatch.setitem(n7.STEP_FUNCTIONS, "genomes", lambda c: (_ for _ in ()).throw(SystemExit("NCBI is down")))
    monkeypatch.setitem(n7.STEP_FUNCTIONS, "amplicons", lambda c: (_ for _ in ()).throw(n7.Skip("no genomes")))
    assert run_main(tmp_path, "--only", "genomes,amplicons,summary") == 1
    state = json.loads((tmp_path / "out" / "night7_state.json").read_text())["steps"]
    assert state["genomes"]["status"] == "failed" and "NCBI is down" in state["genomes"]["detail"]
    assert state["amplicons"]["status"] == "skipped" and state["summary"]["status"] == "done"
    text = (tmp_path / "out" / "NIGHT7_SUMMARY.md").read_text()
    assert "genomes: failed" in text and "amplicons: skipped" in text and "filter_monomer: not run" in text


def test_finished_steps_are_skipped_until_they_are_redone(tmp_path, monkeypatch):
    runs = []
    monkeypatch.setitem(n7.STEP_FUNCTIONS, "genomes", lambda c: runs.append("g") or "ok")
    assert run_main(tmp_path, "--only", "genomes") == 0 and runs == ["g"]
    assert run_main(tmp_path, "--only", "genomes") == 0 and runs == ["g"]                            # finished earlier: not repeated
    assert run_main(tmp_path, "--redo", "genomes") == 0 and runs == ["g", "g"]                        # --redo repeats it and writes the summary
    assert (tmp_path / "out" / "NIGHT7_SUMMARY.md").is_file()
    assert run_main(tmp_path, "--only", "nonsense") == 2


def test_the_summary_collects_the_gates_the_warnings_and_the_versions(tmp_path, monkeypatch):
    monkeypatch.setattr(n7, "versions", lambda: {"Python": "3.12", "cutadapt": "5.2", "R and DADA2": "R 4.6.1 | DADA2 1.40"})
    c = make_ctx(tmp_path)
    fake_l4(monkeypatch, [], met=True)
    c.monomer.mkdir(parents=True)
    for p in n7.MONOMER_PAIRS:
        (c.monomer / f"x.Kinnex16S_{p}.fastq.gz").write_bytes(b"x")
    n7.step_filter_monomer(c)
    c.amplicons.mkdir(parents=True)
    (c.amplicons / "EXPECTED_AMPLICONS.md").write_text("# t\n\nPrimers A and B, at most 3 mismatches.\n\n| strain | copies |\n|---|---|\n| X | 7 |\n\n## To check\n\n- X: a draft genome\n\n## What make_reference.py reported\n\n```\nlog\n```\n")
    c.state["steps"]["filter_monomer"] = {"status": "done", "seconds": 120, "detail": "gate L4 MET (largest difference 0.10 points)", "finished": "now"}
    n7.step_summary(c)
    text = (tmp_path / "out" / "NIGHT7_SUMMARY.md").read_text()
    assert "cutadapt 5.2" in text and "DADA2 1.40" in text
    assert "filter_monomer: done in 2 min (gate L4 MET" in text and "**Gate (L4): MET.**" in text and "| s1 | 2 | 80.0 |" in text
    assert "Primers A and B" in text and "- X: a draft genome" in text and "make_reference.py reported" not in text
    assert "Kinnex files" in text and "No result" in text                                              # the Kinnex step did not run, and the summary says so


def test_check_reports_what_is_missing(tmp_path, monkeypatch):
    c = make_ctx(tmp_path)
    write_table(c.table)
    monkeypatch.setattr(n7, "tool_version", lambda cmd, timeout=60: (False, "No module named cutadapt") if "cutadapt" in cmd else (True, "4.6.1"))
    monkeypatch.setattr(n7, "find_rscript", lambda: None)
    results = {w: (ok, d) for w, ok, d in n7.check_environment(c, opener=fake_ncbi({n7.NCBI: b""}))}
    assert results["cutadapt"][0] is False and results["Rscript"][0] is False
    assert results["Monomer files for L4"][0] is False and results["Kinnex files for L4"][0] is False and results["strain table"][0] is True
    assert results["NCBI reachable"][0] is True
    broken = {w: ok for w, ok, _ in n7.check_environment(c, opener=fake_ncbi({}))}
    assert broken["NCBI reachable"] is False


def test_check_is_all_clear_when_everything_is_there(tmp_path, monkeypatch):
    c = make_ctx(tmp_path)
    write_table(c.table)
    for folder in (c.monomer, c.kinnex):
        folder.mkdir(parents=True)
    for p in n7.MONOMER_PAIRS:
        (c.monomer / f"x.Kinnex16S_{p}.fastq.gz").write_bytes(b"x")
    (c.kinnex / "k.fastq.gz").write_bytes(b"x")
    c.genomes.mkdir(parents=True)
    (c.genomes / "g.fna.gz").write_bytes(b"x")
    monkeypatch.setattr(n7, "tool_version", lambda cmd, timeout=60: (True, "1.0"))
    monkeypatch.setattr(n7, "find_rscript", lambda: "Rscript")
    results = n7.check_environment(c, opener=fake_ncbi({}))
    assert all(ok for _, ok, _ in results), [r for r in results if not r[1]]
    assert dict((w, d) for w, _, d in results)["NCBI reachable"].startswith("not needed")


def test_list_shows_the_steps(tmp_path, capsys):
    run_main(tmp_path, "--list")
    out = capsys.readouterr().out
    assert all(s in out for s in n7.STEPS) and "not run" in out


def test_the_script_parses_as_python_3_9():
    ast.parse(Path(n7.__file__).read_text(), feature_version=(3, 9))


def test_link_or_copy_links_where_it_can_and_copies_where_it_cannot(tmp_path, monkeypatch):
    src = tmp_path / "a.bin"
    src.write_bytes(b"data")
    n7.link_or_copy(src, tmp_path / "linked.bin")
    assert (tmp_path / "linked.bin").read_bytes() == b"data" and src.stat().st_nlink == 2           # one file with two names
    monkeypatch.setattr(n7.os, "link", lambda *a: (_ for _ in ()).throw(OSError("different drive")))
    n7.link_or_copy(src, tmp_path / "copied.bin")
    assert (tmp_path / "copied.bin").read_bytes() == b"data"
    (tmp_path / "copied.bin").write_bytes(b"other")
    n7.link_or_copy(src, tmp_path / "copied.bin")                                                      # an existing name is left alone
    assert (tmp_path / "copied.bin").read_bytes() == b"other"


def test_the_input_folder_is_removed_even_when_the_check_fails(tmp_path, monkeypatch):
    c = make_ctx(tmp_path)
    c.monomer.mkdir(parents=True)
    for p in n7.MONOMER_PAIRS:
        (c.monomer / f"x.Kinnex16S_{p}.fastq.gz").write_bytes(b"x")
    monkeypatch.setattr(lf, "check_tools", lambda: None)
    monkeypatch.setattr(lf, "check_folder", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("cutadapt failed")))
    with pytest.raises(RuntimeError, match="cutadapt failed"):
        n7.step_filter_monomer(c)
    assert not c.monomer_in.exists()


def test_the_whole_run_with_stand_ins_gives_a_summary_with_both_gates_and_the_expected_amplicons(tmp_path, monkeypatch):
    import test_genome_amplicons as t_genome
    import test_longread_filter_check as t_l4
    data = tmp_path / "data"
    rows = [dict(ROW, accession="GCF_000000001.1", assembly_name="ASMa", atcc_number="ATCC 1", name_on_atcc_page="Strain one", share_percent="18"),
            dict(ROW, accession="GCF_000000002.1", assembly_name="ASMb", atcc_number="ATCC 2", name_on_atcc_page="Strain two", share_percent="1.8")]
    table = write_table(tmp_path / "t.tsv", rows)
    genomes = data / "ref" / "genomes"
    genomes.mkdir(parents=True)
    for i, r in enumerate(rows):                                    # the genomes are already here, so the download step makes no request
        with gzip.open(genomes / n7.genome_file_name(r), "wt") as fh:
            fh.write(">chr\n" + t_genome.genome([(t_genome.amplicon(i + 1), "+")] * (3 - i), seed=20 + i, background=1500) + "\n")
    for folder, names in ((data / "sq2_monomer_fastq", [f"m.Kinnex16S_{p}.hifi_reads.fastq.gz" for p in n7.MONOMER_PAIRS]), (data / "sq2_kinnex_fastq", ["k.fastq.gz"])):
        folder.mkdir(parents=True)
        for j, nm in enumerate(names):
            t_l4.write_fastq(folder / nm, t_l4.make_reads(60, 50, 20, 10, seed=30 + j))
    monkeypatch.setattr(lf, "run_cutadapt", t_l4.stand_in_cutadapt)
    monkeypatch.setattr(lf, "dada2_filter", t_l4.stand_in_filter)
    monkeypatch.setattr(lf, "dada2_remove_primers", t_l4.stand_in_remove_primers)
    monkeypatch.setattr(lf, "check_tools", lambda: None)
    monkeypatch.setattr(n7, "versions", lambda: {"Python": "3.12"})
    assert n7.main(["--out", str(tmp_path / "out"), "--data", str(data), "--table", str(table), "--reads", "0"]) == 0
    state = json.loads((tmp_path / "out" / "night7_state.json").read_text())["steps"]
    assert {k: v["status"] for k, v in state.items()} == {s: "done" for s in n7.STEPS}
    text = (tmp_path / "out" / "NIGHT7_SUMMARY.md").read_text()
    assert text.count("**Gate (L4): MET.**") == 2 and "Monomer files" in text and "Kinnex files" in text
    assert "Strain_one_ATCC1" in text and "Strain_two_ATCC2" in text and "expected read share" in text
    assert "| sample | reads | estimate's search %" in text                                          # the primer table of both runs is in the summary
    assert not (tmp_path / "out" / "l4_monomer_input").exists() and not (tmp_path / "out" / "filter_check_monomer" / "work").exists()
