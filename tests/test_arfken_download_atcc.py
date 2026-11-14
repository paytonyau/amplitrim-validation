"""validation/arfken_download.py, the ATCC standard: its runs come from ENA's run table (no accession is written in the script), are saved to a table, and are downloaded
and checked from it. ENA is replaced by fake answers."""
import gzip
import importlib.util
import io
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


def load():
    spec = importlib.util.spec_from_file_location("arfken_download_atcc", HERE.parent / "validation" / "arfken_download.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


A = load()
HEAD = "run_accession\tlibrary_name\tsample_alias\tsample_title\tread_count\tbase_count"


def table(*rows):
    return "\n".join([HEAD] + ["\t".join(map(str, r)) for r in rows]) + "\n"


GOOD = table(("SRRA1", "ATCC-ITS1a", "x", "", 4, 2000), ("SRRA2", "ATCC-ITS1b", "x", "", 5, 2500), ("SRRA3", "", "", "ATCC_ITS1C", 6, 3000),
             ("SRRB1", "ATCC-ITS2a", "x", "", 3, 1800), ("SRRB2", "ATCC-ITS2b", "x", "", 2, 1200), ("SRRB3", "ATCC-ITS2c", "x", "", 7, 4200),
             ("SRRX1", "ATCC2-ITS1a", "x", "", 9, 9), ("SRRX2", "ATCC2-ITS2a", "x", "", 9, 9), ("SRRC1", "C-ITS1a", "x", "", 1, 1), ("SRRN1", "Neg-ITS1a", "x", "", 1, 1))


class Response(io.BytesIO):
    def __init__(self, data: bytes):
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data))}


class Net:
    """ENA's run table at the report URL; a FASTQ (gzip) with the right number of reads at a run's URL."""

    def __init__(self, report_text, counts=None):
        self.report, self.counts, self.calls = report_text, counts or {}, []

    def __call__(self, url):
        self.calls.append(url)
        if url == A.ENA_REPORT:
            return Response(self.report.encode())
        run = url.rsplit("/", 1)[1].split("_")[0]
        n = self.counts[run]
        return Response(gzip.compress("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(n)).encode()))


def test_the_six_atcc_runs_are_found_by_name_and_atcc2_and_the_others_are_not():
    runs = A.find_atcc(Net(GOOD))
    assert runs["atcc_its1"] == [("SRRA1", "ATCC-ITS1a", 4, 2000), ("SRRA2", "ATCC-ITS1b", 5, 2500), ("SRRA3", "ATCC-ITS1c", 6, 3000)]      # a name in the sample title counts too
    assert [r[0] for r in runs["atcc_its2"]] == ["SRRB1", "SRRB2", "SRRB3"]
    assert not any(r[0].startswith("SRRX") for k in runs for r in runs[k])


def test_a_missing_replicate_or_duplicated_markers_are_refused():
    two = table(("SRRA1", "ATCC-ITS1a", "", "", 4, 2000), ("SRRA2", "ATCC-ITS1b", "", "", 5, 2500),
                ("SRRB1", "ATCC-ITS2a", "", "", 3, 1800), ("SRRB2", "ATCC-ITS2b", "", "", 2, 1200), ("SRRB3", "ATCC-ITS2c", "", "", 7, 4200))
    with pytest.raises(RuntimeError, match="2 ATCC libraries for atcc_its1"):
        A.find_atcc(Net(two))
    dup = table(*[(f"SRRA{i}", f"ATCC-ITS1{c}", "", "", 4 + i, 2000 + i) for i, c in enumerate("abc")],
                *[(f"SRRB{i}", f"ATCC-ITS2{c}", "", "", 4 + i, 2000 + i) for i, c in enumerate("abc")])
    with pytest.raises(RuntimeError, match="identical read and base counts"):
        A.find_atcc(Net(dup))


def test_find_saves_the_table_and_the_download_and_the_check_use_it(tmp_path, capsys):
    counts = {"SRRA1": 4, "SRRA2": 5, "SRRA3": 6, "SRRB1": 3, "SRRB2": 2, "SRRB3": 7}
    net = Net(GOOD, counts)
    assert A.main(["--base", str(tmp_path), "--find-atcc"], opener=net, sleep=lambda s: None) == 0
    saved = (tmp_path / A.ATCC_TABLE).read_text()
    assert "atcc_its1\tSRRA3\tATCC-ITS1c\t6\t3000" in saved and "SRRX1" not in saved
    assert "saved:" in capsys.readouterr().out and len(net.calls) == 1                          # nothing downloaded yet
    assert A.main(["--base", str(tmp_path), "--only", "atcc"], opener=net, sleep=lambda s: None) == 0
    assert (tmp_path / "20_mock_ITS1_arfken_atcc_PRJNA693350" / "raw" / "ATCC-ITS1c_R2.fastq.gz").is_file()
    assert (tmp_path / "21_mock_ITS2_arfken_atcc_PRJNA693350" / "raw" / "ATCC-ITS2a_R1.fastq.gz").is_file()
    assert A.main(["--base", str(tmp_path), "--only", "atcc", "--check"], opener=net) == 0
    assert "every file is there" in capsys.readouterr().out
    bad = tmp_path / "21_mock_ITS2_arfken_atcc_PRJNA693350" / "raw" / "ATCC-ITS2b_R1.fastq.gz"
    bad.write_bytes(gzip.compress(b"@r\nA\n+\nI\n"))                                            # 1 read instead of 2
    assert A.main(["--base", str(tmp_path), "--only", "atcc_its2", "--check"], opener=net) == 2
    assert "ATCC-ITS2b mate 1: 1 reads, the run table says 2" in capsys.readouterr().out


def test_without_the_saved_table_the_atcc_folders_say_what_to_run(tmp_path):
    with pytest.raises(SystemExit, match="--find-atcc"):
        A.main(["--base", str(tmp_path), "--only", "atcc", "--dry-run"], opener=Net(GOOD))


def test_the_default_is_still_the_isolates_and_the_mixes(tmp_path, capsys):
    assert A.main(["--base", str(tmp_path), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "ITS1:" in out and "MIX_ITS2:" in out and "ATCC" not in out
    assert A.ATCC_NAME.match("ATCC2-ITS1a") is None and A.ATCC_NAME.match("atcc-its2b")
