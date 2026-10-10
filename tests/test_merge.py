"""vsearch output parsing and merge behaviour."""
import math
from pathlib import Path

import pandas as pd

from amplitrim.merge_wrapper import merged_lengths, parse_vsearch_log, run_vsearch
from conftest import need_vsearch

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_real_vsearch_log():
    o = parse_vsearch_log((FIXTURES / "vsearch_mergepairs_log.txt").read_text(), "strict")
    assert (o.total_pairs, o.merged_pairs) == (10000, 4944)
    assert abs(o.mean_fragment_len - 426.35) < 1e-9
    assert o.fail_reasons["too many differences"] == 4975 and len(o.fail_reasons) == 3


def test_zero_merged_does_not_crash():
    o = parse_vsearch_log("     100  Pairs\n       0  Merged (0.0%)\n")
    assert o.merge_rate == 0.0 and math.isnan(o.mean_fragment_len)


@need_vsearch
def test_tolerant_merges_at_least_as_much_as_strict_and_never_misassembles(sim, tmp_path):
    r1, r2 = sim / "primer_trimmed/sampleA_R1.fastq.gz", sim / "primer_trimmed/sampleA_R2.fastq.gz"
    strict = run_vsearch(r1, r2, tmp_path, "t", 12, 0, "strict", keep_merged=False)
    tol = run_vsearch(r1, r2, tmp_path, "t", 12, None, "tolerant")
    assert tol.merge_rate >= strict.merge_rate and tol.merge_rate > 0.9
    truth = pd.read_csv(sim / "truth.tsv", sep="\t").set_index("read_id")["true_insert_len"]
    names = [ln[1:].split()[0] for ln in open(tol.merged_path) if ln.startswith(">")]
    assert all(truth[n] == length for n, length in zip(names, merged_lengths(tol.merged_path)))


def test_discarding_merged_reads_uses_the_platform_null_device(tmp_path, monkeypatch):
    """'/dev/null' does not exist on Windows; the portable name is os.devnull."""
    import os
    import subprocess
    from amplitrim import merge_wrapper
    seen = {}

    def fake_run(cmd, *a, **k):
        seen["cmd"] = cmd
        (tmp_path / "x.log").write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    monkeypatch.setattr(merge_wrapper.subprocess, "run", fake_run)
    try:
        merge_wrapper.run_vsearch(tmp_path / "a.fq", tmp_path / "b.fq", tmp_path, "x", 12, None, "tolerant",
                                  keep_merged=False)
    except Exception:
        pass                                                         # parsing the empty fake log is not the point
    assert os.devnull in seen["cmd"]


@need_vsearch
def test_a_vsearch_log_containing_a_non_utf8_path_is_still_read(sim, tmp_path, monkeypatch):
    """On Windows vsearch writes input paths into its log in the ANSI code page (e-acute is the single byte 0xe9)."""
    import subprocess
    from pathlib import Path
    from amplitrim import merge_wrapper
    real_run = subprocess.run

    def run_then_add_an_ansi_path(cmd, *args, **kwargs):
        res = real_run(cmd, *args, **kwargs)
        log = Path(cmd[cmd.index("--log") + 1])
        log.write_bytes(b"Reading file C:\\donn\xe9es\\work\\subsample_R1.fastq.gz\n" + log.read_bytes())
        return res

    monkeypatch.setattr(merge_wrapper.subprocess, "run", run_then_add_an_ansi_path)
    r1, r2 = sim / "primer_trimmed" / "sampleA_R1.fastq.gz", sim / "primer_trimmed" / "sampleA_R2.fastq.gz"
    outcome = merge_wrapper.run_vsearch(r1, r2, tmp_path, "t", 12, None, "tolerant", keep_merged=False)
    assert outcome.merged_pairs > 0
