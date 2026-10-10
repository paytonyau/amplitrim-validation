"""Pairing of input files and the shared subsample."""
import gzip

import pytest

from amplitrim.subsample import build_subsample, count_pairs, find_pairs


def write_fq(path, recs):
    with gzip.open(path, "wt", newline="\n") as fh:
        for name, seq, qual in recs:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def read_names(path):
    with gzip.open(path, "rt") as fh:
        return [line[1:].strip() for i, line in enumerate(fh) if i % 4 == 0]


def make_input(tmp_path, samples, n=10, length=30):
    inp = tmp_path / "in"
    inp.mkdir()
    for s in samples:
        recs = [(f"{s}{i}", "ACGT" * 20, "I" * 80) for i in range(n)]
        recs = [(name, seq[:length], qual[:length]) for name, seq, qual in recs]
        write_fq(inp / f"{s}_R1.fastq.gz", recs)
        write_fq(inp / f"{s}_R2.fastq.gz", recs)
    return inp


def test_find_pairs_accepts_plain_and_illumina_names(tmp_path):
    for n in ("s1_S1_L001_R1_001.fastq.gz", "s1_S1_L001_R2_001.fastq.gz", "plain_R1.fastq.gz",
              "plain_R2.fastq.gz", "notes.txt", "lone_R1.fastq.gz"):
        (tmp_path / n).touch()
    assert [s for s, _, _ in find_pairs(tmp_path)] == ["plain", "s1_S1_L001"]


def test_find_pairs_accepts_sra_style_names_with_an_underscore_before_fastq(tmp_path):
    for n in ("SRR39528347_1_fastq.gz", "SRR39528347_2_fastq.gz", "other_R1_fq.gz", "other_R2_fq.gz"):
        (tmp_path / n).touch()
    assert [s for s, _, _ in find_pairs(tmp_path)] == ["SRR39528347", "other"]


def test_no_pairs_error_lists_the_files_seen(tmp_path):
    (tmp_path / "weird_name.fastq.gz").touch()
    with pytest.raises(ValueError, match="weird_name.fastq.gz"):
        find_pairs(tmp_path)


def test_count_pairs(tmp_path):
    write_fq(tmp_path / "a.fq.gz", [(f"r{i}", "ACGT", "IIII") for i in range(7)])
    assert count_pairs(tmp_path / "a.fq.gz") == 7


def test_every_nth_is_deterministic(tmp_path):
    info = build_subsample(make_input(tmp_path, ["s1"]), tmp_path / "sub", every_nth=3)
    assert info.n_pairs == 4                                  # reads 0, 3, 6, 9
    assert (info.read_len_f, info.read_len_r, info.uniform_length) == (30, 30, True)
    assert read_names(info.r1) == read_names(info.r2) == ["s10", "s13", "s16", "s19"]


def test_cap_raises_the_stride_and_samples_every_sample(tmp_path):
    info = build_subsample(make_input(tmp_path, ["a", "b"], n=100), tmp_path / "sub", max_pairs=40)
    assert (info.n_input_pairs, info.stride, info.n_pairs) == (200, 5, 40)
    assert {n[0] for n in read_names(info.r1)} == {"a", "b"}   # not just the first sample


def test_no_cap_keeps_everything(tmp_path):
    assert build_subsample(make_input(tmp_path, ["s"], n=30), tmp_path / "sub", max_pairs=0).n_pairs == 30


def test_reads_of_different_lengths_are_kept_and_reported(tmp_path):
    inp = tmp_path / "in"
    inp.mkdir()
    write_fq(inp / "s_R1.fastq.gz", [("a", "A" * 100, "I" * 100), ("b", "A" * 90, "I" * 90)])
    write_fq(inp / "s_R2.fastq.gz", [("a", "C" * 80, "I" * 80), ("b", "C" * 80, "I" * 80)])
    info = build_subsample(inp, tmp_path / "sub")
    assert info.n_pairs == 2 and not info.uniform_length
    assert (info.min_len_f, info.read_len_f, info.min_len_r, info.read_len_r) == (90, 100, 80, 80)


def test_the_subsample_never_contains_carriage_returns(tmp_path):
    info = build_subsample(make_input(tmp_path, ["s"]), tmp_path / "sub")
    for path in (info.r1, info.r2):
        with gzip.open(path, "rb") as fh:
            assert b"\r" not in fh.read()


def test_the_number_of_distinct_quality_values_is_reported(tmp_path):
    def make(q_pattern):
        inp = tmp_path / f"q{len(set(q_pattern))}"
        inp.mkdir()
        recs = [(f"r{i}", "ACGT" * 5, "".join(q_pattern[(i + j) % len(q_pattern)] for j in range(20))) for i in range(30)]
        write_fq(inp / "s_R1.fastq.gz", recs)
        write_fq(inp / "s_R2.fastq.gz", recs)
        return build_subsample(inp, tmp_path / f"sub{len(set(q_pattern))}")
    assert make(["?", "$"]).quality_values == 2                    # Q30 and Q3 only, as in a simplified-quality file
    assert make(["I", "?", "8", "'"]).quality_values == 4
