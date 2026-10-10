"""The ENA downloader, with the network replaced by stand-ins."""
import hashlib
import importlib.util
import sys
import urllib.error
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("download_ena", Path(__file__).parent.parent / "validation" / "download_ena.py")
dl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dl)

dl.PAUSE = 0                                   # no waiting between mocked attempts
R1, R2, ORPHAN = b"read-one-data", b"read-two-data", b"unpaired"
md5 = lambda b: hashlib.md5(b).hexdigest()


class FakeResponse:
    def __init__(self, text):
        self.text = text

    def read(self):
        return self.text.encode()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def ena_answer(acc, with_orphan=True):
    files = [("ftp.sra.ebi.ac.uk/vol1/x/%s.fastq.gz" % acc, ORPHAN)] if with_orphan else []
    files += [(f"ftp.sra.ebi.ac.uk/vol1/x/{acc}_1.fastq.gz", R1), (f"ftp.sra.ebi.ac.uk/vol1/x/{acc}_2.fastq.gz", R2)]
    return FakeResponse("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
                        f"{acc}\t{';'.join(u for u, _ in files)}\t{';'.join(md5(b) for _, b in files)}\t"
                        f"{';'.join(str(len(b)) for _, b in files)}\n")


def make_retriever(payloads, calls):
    def retrieve(url, path):
        calls.append(url)
        Path(path).write_bytes(payloads[url.rsplit("_", 1)[-1][0]] if "_" in url else ORPHAN)
    return retrieve


def test_accessions_come_from_a_run_selector_csv_with_a_bom_and_quotes(tmp_path):
    f = tmp_path / "SraRunInfo.csv"
    f.write_bytes('\ufeff"Run","ReleaseDate","spots"\nSRR13969386,2021,1\n"SRR13969385",2021,2\nSRR13969386,dup,3\n'.encode("utf-8"))
    assert dl.read_accessions(str(f)) == ["SRR13969386", "SRR13969385"]


def test_accessions_from_a_text_file_and_a_single_argument(tmp_path):
    f = tmp_path / "runs.txt"
    f.write_text("# my runs\nSRR1 label\n\nERR22\n", encoding="utf-8")
    assert dl.read_accessions(str(f)) == ["SRR1", "ERR22"]
    assert dl.read_accessions("SRR13969386") == ["SRR13969386"] and dl.read_accessions("not-a-run") == []


def test_a_paired_run_is_downloaded_verified_and_named_r1_r2(tmp_path):
    calls = []
    lines = dl.download("SRR9", tmp_path, opener=lambda *a, **k: ena_answer("SRR9"), retriever=make_retriever({"1": R1, "2": R2}, calls),
                        log=lambda m: None)
    assert (tmp_path / "SRR9_R1.fastq.gz").read_bytes() == R1 and (tmp_path / "SRR9_R2.fastq.gz").read_bytes() == R2
    assert lines == [f"{md5(R1)}  SRR9_R1.fastq.gz", f"{md5(R2)}  SRR9_R2.fastq.gz"]
    assert len(calls) == 2 and not list(tmp_path.glob("*.part"))          # the unpaired file is not fetched


def test_files_already_downloaded_with_the_right_checksum_are_not_fetched_again(tmp_path):
    calls = []
    kw = dict(opener=lambda *a, **k: ena_answer("SRR9"), retriever=make_retriever({"1": R1, "2": R2}, calls), log=lambda m: None)
    dl.download("SRR9", tmp_path, **kw)
    dl.download("SRR9", tmp_path, **kw)
    assert len(calls) == 2


def test_a_corrupted_download_is_retried_and_never_kept(tmp_path):
    calls = []
    with pytest.raises(SystemExit, match="no verified download after 2 attempts"):
        dl.download("SRR9", tmp_path, retries=2, opener=lambda *a, **k: ena_answer("SRR9"),
                    retriever=make_retriever({"1": b"CORRUPT", "2": R2}, calls), log=lambda m: None)
    assert not (tmp_path / "SRR9_R1.fastq.gz").exists() and not list(tmp_path.glob("*.part"))


def test_network_failures_are_plain_errors_and_downloads_are_retried(tmp_path):
    import urllib.error

    def blocked(*a, **k):
        raise urllib.error.HTTPError("u", 403, "Forbidden", {}, None)

    with pytest.raises(SystemExit, match="ENA refused the request .*403"):
        dl.download("SRR9", tmp_path, opener=blocked, retriever=lambda *a: None, log=lambda m: None)
    attempts = []

    def flaky(url, path):
        attempts.append(url)
        if len(attempts) == 1:
            raise urllib.error.URLError("connection reset")
        Path(path).write_bytes(R1 if "_1." in url else R2)

    dl.download("SRR9", tmp_path, opener=lambda *a, **k: ena_answer("SRR9"), retriever=flaky, log=lambda m: None)
    assert (tmp_path / "SRR9_R1.fastq.gz").read_bytes() == R1 and len(attempts) == 3      # one retry, then two clean downloads


def test_a_run_without_a_paired_set_is_a_clear_error(tmp_path):
    single = FakeResponse(f"run_accession\tfastq_ftp\tfastq_md5\nSRR7\tftp.sra.ebi.ac.uk/vol1/SRR7.fastq.gz\t{md5(ORPHAN)}\n")
    with pytest.raises(SystemExit, match="does not list a paired FASTQ set"):
        dl.download("SRR7", tmp_path, opener=lambda *a, **k: single, retriever=lambda *a: None, log=lambda m: None)
    with pytest.raises(SystemExit, match="no FASTQ files"):
        dl.download("SRR8", tmp_path, opener=lambda *a, **k: FakeResponse("run_accession\tfastq_ftp\tfastq_md5\n"), retriever=lambda *a: None, log=lambda m: None)


def test_the_dry_run_reports_the_paired_size_without_downloading(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(dl.urllib.request, "urlopen", lambda *a, **k: ena_answer("SRR9"))
    monkeypatch.setattr(dl.urllib.request, "urlretrieve", lambda *a: pytest.fail("a dry run must not download"))
    monkeypatch.setattr(dl.sys, "argv", ["download_ena.py", "SRR9", "--out", str(tmp_path), "--dry-run"])
    assert dl.main() == 0
    out = capsys.readouterr().out
    assert "SRR9: 0 MB" in out and "1 run(s)" in out and "Nothing was downloaded" in out
    assert not list(tmp_path.iterdir())
    assert dl.paired_files("SRR9", lambda *a, **k: ena_answer("SRR9"))["1"][2] == len(R1)


# ---- resuming interrupted transfers and carrying on after a failed run ----------------------------------------
import io


class FileResponse:
    """A file download that can drop after `drop_after` bytes, and honours Range requests."""
    def __init__(self, payload, request, drop_after=None, honour_range=True):
        range_header = request.get_header("Range")
        start = int(range_header.split("=")[1].rstrip("-")) if range_header else 0
        self.status = 206 if (start and honour_range) else 200
        self.buffer = io.BytesIO(payload[start:] if honour_range else payload)
        self.drop_after, self.sent = drop_after, 0

    def read(self, n):
        if self.drop_after is not None and self.sent >= self.drop_after:
            raise ConnectionResetError("connection reset by peer")
        data = self.buffer.read(n if self.drop_after is None else min(n, self.drop_after - self.sent))
        self.sent += len(data)
        return data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def test_an_interrupted_transfer_resumes_from_the_partial_file(tmp_path):
    payload, requests = bytes(range(256)) * 40, []

    def server(request, timeout=None):
        requests.append(request.get_header("Range"))
        return FileResponse(payload, request, drop_after=3000 if len(requests) == 1 else None)

    part = tmp_path / "x.part"
    with pytest.raises(ConnectionResetError):
        dl.fetch_resumable("https://h/x", part, len(payload), server)
    assert part.stat().st_size == 3000                                   # what arrived is kept
    dl.fetch_resumable("https://h/x", part, len(payload), server)
    assert part.read_bytes() == payload and requests == [None, "bytes=3000-"]


def test_a_server_that_ignores_ranges_makes_the_download_start_again(tmp_path):
    payload = b"abcdefgh" * 500
    part = tmp_path / "y.part"
    part.write_bytes(b"STALE-PARTIAL")
    dl.fetch_resumable("https://h/y", part, len(payload), lambda req, timeout=None: FileResponse(payload, req, honour_range=False))
    assert part.read_bytes() == payload


def test_download_survives_repeated_drops_and_verifies_the_result(tmp_path):
    payload_by_mate, attempts = {"1": R1 * 600, "2": R2 * 600}, []

    def api(*a, **k):
        r1, r2 = payload_by_mate["1"], payload_by_mate["2"]
        return FakeResponse("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"
                            f"SRR5\th/SRR5_1.fastq.gz;h/SRR5_2.fastq.gz\t{md5(r1)};{md5(r2)}\t{len(r1)};{len(r2)}\n")

    def server(request, timeout=None):
        attempts.append(request.full_url)
        mate = request.full_url.split("_")[-1][0]
        return FileResponse(payload_by_mate[mate], request, drop_after=4000 if len(attempts) % 2 else None)

    messages = []
    lines = dl.download("SRR5", tmp_path, retries=6, opener=api, file_opener=server, log=messages.append)
    assert (tmp_path / "SRR5_R1.fastq.gz").read_bytes() == payload_by_mate["1"] and len(lines) == 2
    assert any("interrupted" in m and "resuming" in m for m in messages) and not list(tmp_path.glob("*.part"))


def test_the_batch_carries_on_after_a_failed_run_and_says_how_to_resume(tmp_path, monkeypatch, capsys):
    def urlopen(req, timeout=None):
        url = req if isinstance(req, str) else req.full_url
        if "filereport" in url:
            return ena_answer("SRR2") if "SRR2" in url else FakeResponse("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n")
        return FileResponse(R1 if "_1." in url else R2, req)

    monkeypatch.setattr(dl.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(dl.sys, "argv", ["download_ena.py", "SRR1", "SRR2", "--out", str(tmp_path)])
    assert dl.main() == 1
    out = capsys.readouterr().out
    assert "SRR1: ENA lists no FASTQ files" in out and "1 of 2 run(s) downloaded and verified" in out
    assert "Not completed: SRR1. Run the same command again" in out
    assert (tmp_path / "SRR2_R1.fastq.gz").exists() and "SRR2_R1" in (tmp_path / "fastq_checksums.md5").read_text()


def test_an_empty_answer_from_ena_is_retried_before_giving_up(tmp_path):
    answers = [FakeResponse("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n"), ena_answer("SRR3")]     # empty, then fine
    pair = dl.paired_files("SRR3", lambda *a, **k: answers.pop(0))
    assert set(pair) == {"1", "2"} and not answers
    always_empty = lambda *a, **k: FakeResponse("run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\n")
    with pytest.raises(SystemExit, match="ENA lists no FASTQ files"):
        dl.paired_files("SRR3", always_empty)


def test_a_header_line_in_an_accession_list_is_not_taken_for_a_run(tmp_path):
    f = tmp_path / "SraAccList.csv"
    f.write_text("acc\nSRR24572428\nSRR24572429\n\n# comment\nERR5,extra\n", encoding="utf-8")
    skipped = []
    assert dl.read_accessions(str(f), skipped) == ["SRR24572428", "SRR24572429", "ERR5"] and skipped == ["acc"]


def _http(code):
    return urllib.error.HTTPError("u", code, "x", {}, None)


def test_a_server_error_is_retried_and_called_ena_s_own(capsys):
    calls = []
    def flaky(*a, **k):
        calls.append(1)
        if len(calls) < 3:
            raise _http(500)
        return ena_answer("SRR5")
    assert dl.paired_files("SRR5", flaky)["1"][2] == len(R1) and len(calls) == 3
    assert "server error (HTTP 500)" in capsys.readouterr().out
    def down(*a, **k):
        raise _http(503)
    with pytest.raises(SystemExit) as e:
        dl.paired_files("SRR5", down)
    assert "server error (HTTP 503)" in str(e.value) and "on ENA's side" in str(e.value) and "run the same command again later" in str(e.value)


def test_a_refusal_is_not_retried_and_a_lost_connection_is_blamed_on_the_connection():
    calls = []
    def refused(*a, **k):
        calls.append(1)
        raise _http(404)
    with pytest.raises(SystemExit, match="refused the request .*404"):
        dl.paired_files("SRR5", refused)
    assert len(calls) == 1
    def offline(*a, **k):
        raise urllib.error.URLError("no route")
    with pytest.raises(SystemExit, match="could not reach ENA .*firewall"):
        dl.paired_files("SRR5", offline)


def test_the_dry_run_lists_what_it_can_when_one_run_cannot_be_listed(monkeypatch, capsys, tmp_path):
    def fake(acc, *a, **k):
        if acc == "SRR2":
            raise SystemExit("SRR2: ENA answered with a server error (HTTP 500).")
        return {"1": ("u1", "m1", 5_000_000), "2": ("u2", "m2", 7_000_000)}
    monkeypatch.setattr(dl, "paired_files", fake)
    monkeypatch.setattr(sys, "argv", ["download_ena.py", "SRR1", "SRR2", "SRR3", "--out", str(tmp_path), "--dry-run"])
    assert dl.main() == 1
    out = capsys.readouterr().out
    assert "SRR1: 12 MB" in out and "SRR3: 12 MB" in out and "2 of 3 run(s) listed" in out and "Could not be listed: SRR2" in out


def test_a_pilot_takes_runs_spread_evenly_over_the_list_and_ena_s_run_table_works_as_the_list(tmp_path, monkeypatch, capsys):
    runs = [f"ERR{100 + i}" for i in range(10)]
    assert dl.evenly_spaced(runs, 3) == ["ERR100", "ERR104", "ERR109"] and dl.evenly_spaced(runs, 1) == ["ERR100"]
    assert dl.evenly_spaced(runs, 10) == runs and dl.evenly_spaced(runs, 0) == runs and dl.evenly_spaced(runs, 25) == runs
    for n, k in ((5, 4), (7, 6), (30, 29), (100, 7)):                   # always k distinct runs, in order
        picked = dl.evenly_spaced(list(range(n)), k)
        assert len(picked) == k and picked == sorted(set(picked)) and picked[0] == 0 and picked[-1] == n - 1
    table = tmp_path / "runs.tsv"                                       # ENA's run table: the accession is the first column
    table.write_text("run_accession\tsample_alias\tread_count\n" + "".join(f"{r}\tS{i}\t1000\n" for i, r in enumerate(runs)), encoding="utf-8")
    monkeypatch.setattr(dl, "paired_files", lambda acc, *a, **k: {"1": ("u", "m", 1_000_000), "2": ("u", "m", 1_000_000)})
    monkeypatch.setattr(sys, "argv", ["download_ena.py", str(table), "--out", str(tmp_path / "o"), "--dry-run", "--max-runs", "2"])
    assert dl.main() == 0
    out = capsys.readouterr().out
    assert "ignored run_accession (not run accessions" in out                      # only the header word is skipped
    assert "using 2 runs" in out and "ERR100: 2 MB" in out and "ERR109: 2 MB" in out and "ERR104" not in out and "2 of 2 run(s) listed" in out


def test_a_saved_error_page_is_recognised_as_a_web_page_not_a_list(tmp_path):
    page = tmp_path / "PRJEB1.tsv"
    page.write_text("<!doctype html>\n<html><head><title>Error: 500 | EMBL\u2019s European Bionformatics Institute</title></head><body>unavailable</body></html>\n", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        dl.read_accessions(str(page))
    assert "is a web page, not a list of runs" in str(e.value) and "Error: 500" in str(e.value) and "curl -f" in str(e.value)
    ok = tmp_path / "ok.tsv"                                           # a real run table is still read
    ok.write_text("run_accession\tsample_alias\nERR1\ta\nERR2\tb\n", encoding="utf-8")
    assert dl.read_accessions(str(ok)) == ["ERR1", "ERR2"]
