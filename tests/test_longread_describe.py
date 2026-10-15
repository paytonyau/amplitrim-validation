"""validation/longread_describe.py on simulated long reads whose properties are known exactly."""
import gzip
import random
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
import longread_describe as lr  # noqa: E402

FWD, REV = "AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT"          # 27F and 1492R
INSERT = 1411                                                      # a read between the primers is this long: 1,450 nt with both primers
CHOICES = lambda p: "".join(random.choice(lr.IUPAC[c]) for c in p)  # noqa: E731


def amplicon(rng_seed: int) -> str:
    rng = random.Random(rng_seed)
    return "".join(rng.choice("ACGT") for _ in range(INSERT))


def read(insert: str, mutate: int = 0, reverse: bool = False, five_only: bool = False, none: bool = False) -> str:
    if none:
        return "".join(random.choice("ACGT") for _ in range(len(insert) + 39))
    seq = (CHOICES(FWD) + insert + lr.revcomp(CHOICES(REV))) if not five_only else (CHOICES(FWD) + insert)
    if mutate:
        s = list(seq)
        for i in random.sample(range(len(FWD)), mutate):                    # errors inside the 5' primer
            s[i] = "ACGT"[("ACGT".index(s[i]) + 1) % 4] if s[i] in "ACGT" else s[i]
        seq = "".join(s)
    return lr.revcomp(seq) if reverse else seq


def write_fastq(path: Path, records) -> None:
    with gzip.open(path, "wt") as fh:
        for i, (seq, qual) in enumerate(records):
            fh.write(f"@r{i}\n{seq}\n+\n{qual}\n")


def phred(values) -> str:
    return "".join(chr(33 + int(v)) for v in values)


def hifi_quality(n: int, seed: int, ends_low: bool = False) -> str:
    rng = random.Random(seed)
    q = [rng.choice([60, 70, 80, 93, 40, 50, 30, 20]) for _ in range(n)]
    if ends_low:
        q[:60] = [10] * 60
        q[-60:] = [10] * 60
    return phred(q)


@pytest.fixture
def mixed_file(tmp_path):
    random.seed(7)
    recs = []
    for i in range(60):
        ins = amplicon(i)
        recs.append((read(ins, mutate=i % 3), hifi_quality(INSERT + 39, i)))                  # 60 forward reads, up to 2 primer errors
    for i in range(30):
        recs.append((read(amplicon(100 + i), reverse=True), hifi_quality(INSERT + 39, 100 + i)))   # 30 reverse
    for i in range(6):
        recs.append((read(amplicon(200 + i), five_only=True), hifi_quality(INSERT + 20, 200 + i)))  # 6 with one primer only
    for i in range(4):
        recs.append(("".join(random.choice("ACGT") for _ in range(1450)), hifi_quality(1450, 300 + i)))                              # 4 with none
    path = tmp_path / "SRR1.fastq.gz"
    write_fastq(path, recs)
    return path


def locator():
    return lr.Locator(FWD, REV)


def test_primers_orientation_and_the_length_between_them_are_found_exactly(mixed_file):
    row = lr.describe_file(mixed_file, locator())
    assert row["reads"] == 100 and row["both_forward"] == 60 and row["both_reverse"] == 30 and row["five_only"] == 6 and row["none"] == 4
    assert row["share_both"] == pytest.approx(0.90) and row["share_both_reverse"] == pytest.approx(30 / 90)
    assert row["between_primers"] == 90 and row["len_median"] == INSERT and row["len_min"] == INSERT and row["len_max"] == INSERT   # the primers are cut off exactly
    assert row["share_within_10pct"] == 1.0 and row["share_over_1p5x"] == 0.0 and row["len_mode_bin"] == 1410
    assert row["raw_len_median"] == INSERT + 39


def test_a_read_with_too_many_primer_errors_is_not_counted_as_having_it():
    random.seed(1)
    ins = amplicon(5)
    bad = list(read(ins))
    for i in range(8):                                             # 8 of the 20 bases wrong: over the 15% allowed
        bad[i] = "ACGT"[("ACGT".index(bad[i]) + 2) % 4] if bad[i] in "ACGT" else "A"
    orientation, category, _, _ = locator().locate("".join(bad))
    assert category in ("three_only", "none") and category != "both"


def test_expected_errors_are_the_sum_of_the_error_probabilities_and_retention_follows(tmp_path):
    random.seed(3)
    recs = [(read(amplicon(i)), phred([30] * (INSERT + 39))) for i in range(40)]                 # Q30 everywhere: 1,411 x 0.001
    recs += [(read(amplicon(50 + i)), phred([20] * (INSERT + 39))) for i in range(10)]           # Q20 everywhere: 1,411 x 0.01
    path = tmp_path / "q.fastq.gz"
    write_fastq(path, recs)
    row = lr.describe_file(path, locator())
    assert row["ee_p25"] == pytest.approx(INSERT * 0.001) and row["ee_p95"] == pytest.approx(INSERT * 0.01)
    assert row["ee_per_1000_median"] == pytest.approx(1.0)
    assert row["keep_maxEE_1"] == 0.0 and row["keep_maxEE_2"] == pytest.approx(0.8) and row["keep_maxEE_15"] == 1.0          # 1.411 passes maxEE 2; 14.11 passes 15 only
    assert row["distinct_q"] == 2 and row["max_q"] == 30


def test_qualities_at_the_ends_are_compared_with_the_middle_in_the_amplicon_not_the_read(tmp_path):
    random.seed(4)
    n = INSERT + 39
    forward = (read(amplicon(1)), phred([40] * 20 + [10] * (INSERT) + [40] * 19))               # the primers are Q40, the amplicon Q10 ...
    path = tmp_path / "e.fastq.gz"
    write_fastq(path, [forward] * 5)
    row = lr.describe_file(path, locator())
    assert row["q_start"] == pytest.approx(10) and row["q_mid"] == pytest.approx(10) and row["q_end"] == pytest.approx(10)       # ... so the primers do not count
    ends = [(read(amplicon(i), reverse=(i % 2 == 1)), None) for i in range(20)]
    recs = []
    for i, (seq, _) in enumerate(ends):
        q = [40] * n
        q[20:80] = [10] * 60                                       # the first 60 bases of the amplicon are poor ...
        qual = phred(q) if i % 2 == 0 else phred(q[::-1])          # ... and so are the same bases of a reverse read, which sits at the other end of the read
        recs.append((seq, qual))
    path2 = tmp_path / "e2.fastq.gz"
    write_fastq(path2, recs)
    row2 = lr.describe_file(path2, locator())
    assert row2["q_start"] < row2["q_mid"] - 3 and row2["q_end"] == pytest.approx(row2["q_mid"], abs=0.5)
    assert any("an end" in f for f in lr.flags(row2))


def test_concatemers_and_short_reads_are_flagged_and_a_simplified_quality_is_noticed(tmp_path):
    random.seed(5)
    recs = [(read(amplicon(i)), phred([40] * (INSERT + 39))) for i in range(90)]
    double = read(amplicon(500) + amplicon(501))                                                 # two amplicons between the primers: about 2 x insert
    recs += [(double, phred([40] * len(double))) for _ in range(5)]
    short = read(amplicon(600)[:300])
    recs += [(short, phred([40] * len(short))) for _ in range(5)]
    path = tmp_path / "c.fastq.gz"
    write_fastq(path, recs)
    row = lr.describe_file(path, locator())
    assert row["share_over_1p5x"] == pytest.approx(0.05) and row["share_under_half"] == pytest.approx(0.05)
    f = " ".join(lr.flags(row))
    assert "concatemers" in f and "under half" in f and "only 1 distinct quality" in f


def test_without_primers_every_read_is_used_as_it_is(tmp_path):
    random.seed(6)
    recs = [(amplicon(i), phred([35] * INSERT)) for i in range(10)]
    path = tmp_path / "n.fastq"
    with open(path, "w") as fh:
        for i, (s, q) in enumerate(recs):
            fh.write(f"@r{i}\n{s}\n+\n{q}\n")
    row = lr.describe_file(path, None)
    assert row["reads"] == 10 and row["between_primers"] == 10 and row["len_median"] == INSERT and row["share_both"] == 1.0


def test_the_folder_run_writes_the_table_and_the_report(tmp_path, mixed_file):
    rows = lr.describe_folder(mixed_file.parent, tmp_path / "out", FWD, REV, log=lambda *_: None)
    assert len(rows) == 1
    tsv = (tmp_path / "out" / "longread_describe.tsv").read_text().splitlines()
    assert tsv[0].split("\t")[:3] == ["sample", "reads", "unreadable"] and len(tsv) == 2
    md = (tmp_path / "out" / "LONGREAD_DESCRIBE.md").read_text()
    assert "1 files, 100 reads; 90.0% have both primers" in md and "## What is unusual" in md and "SRR1" in md
    with pytest.raises(SystemExit):
        lr.describe_folder(mixed_file.parent, tmp_path / "x", FWD, None, log=lambda *_: None)
    with pytest.raises(SystemExit):
        lr.describe_folder(tmp_path / "out", tmp_path / "y", log=lambda *_: None)                  # the output folder holds no FASTQ files


def test_the_pooled_summary_weights_files_by_their_reads():
    a = {"reads": 100, "both_forward": 90, "both_reverse": 0, "between_primers": 90, "len_median": 1400.0, "share_within_10pct": 1.0, "distinct_q": 50, "max_q": 90,
         "p_start": 0.001, "p_mid": 0.001, "p_end": 0.001}
    b = {"reads": 300, "both_forward": 100, "both_reverse": 100, "between_primers": 200, "len_median": 1500.0, "share_within_10pct": 0.5, "distinct_q": 30, "max_q": 60,
         "p_start": 0.01, "p_mid": 0.01, "p_end": 0.01}
    p = lr.pooled([a, b])
    assert p["reads"] == 400 and p["share_both"] == pytest.approx(290 / 400) and p["len_median_of_medians"] == 1450.0
    assert p["share_within_10pct"] == pytest.approx((90 * 1.0 + 200 * 0.5) / 290) and p["distinct_q_min"] == 30 and p["max_q_min"] == 60
    assert p["q_mid"] == pytest.approx(-10 * np.log10((0.001 * 90 + 0.01 * 200) / 290))


# ----------------------------------------------------------------------------- primers with insertions and deletions (the case of night 6: errors in the primer hid primers that were there)

def indel_read(insert: str, five: str = "", three: str = "", reverse: bool = False) -> str:
    """A read whose primers carry the given edits: 'sub', 'del', 'ins' or '' for the 5' primer and for the 3' primer (inside the primer, away from its ends)."""
    def edit(primer: str, kind: str) -> str:
        if kind == "del":
            return primer[:8] + primer[9:]
        if kind == "ins":
            return primer[:8] + "G" + primer[8:]
        if kind == "sub":
            return primer[:8] + ("A" if primer[8] != "A" else "C") + primer[9:]
        return primer
    f = edit("".join(random.choice(lr.IUPAC[c]) for c in FWD), five)
    r = edit(lr.revcomp("".join(random.choice(lr.IUPAC[c]) for c in REV)), three)
    seq = f + insert + r
    return lr.revcomp(seq) if reverse else seq


def test_the_edit_distance_search_gives_the_hand_checked_answers():
    primer = lr.Primer("ACGTTGCAAGCTAGCTAGCA")
    flank = "TTTTCCCCGGGG"

    def run(read):
        codes = lr.LUT[np.frombuffer(read.encode(), dtype=np.uint8)][None, :]
        d, e = lr.edit_matches(codes, primer)
        return int(d[0]), int(e[0])
    tail = "GGGGAAAA"
    assert run(flank + primer.seq + tail) == (0, 32)
    sub = list(primer.seq)
    sub[5] = "A" if sub[5] != "A" else "C"
    assert run(flank + "".join(sub) + tail) == (1, 32)
    assert run(flank + primer.seq[:7] + primer.seq[8:] + tail) == (1, 31)                       # a deletion: the stretch is one base shorter
    assert run(flank + primer.seq[:7] + "T" + primer.seq[7:] + tail) == (1, 33)                 # an insertion: one base longer
    assert run(flank + primer.seq[:4] + primer.seq[5:12] + "G" + primer.seq[12:] + tail) == (2, 32)
    last = list(primer.seq)
    last[-1] = "T" if last[-1] != "T" else "A"
    assert run(flank + "".join(last) + tail) == (1, 32)                                          # a substitution at the end is not a missing base: the end does not move
    assert run("A" * 40)[0] > 3                                                                  # nothing like the primer: far over the budget


def test_primers_with_an_insertion_or_a_deletion_are_found_and_the_cut_is_still_exact(tmp_path):
    random.seed(21)
    recs = []
    kinds = [("del", ""), ("ins", ""), ("", "del"), ("", "ins"), ("del", "ins"), ("ins", "del"), ("sub", "del"), ("", "")]
    for i in range(160):
        five, three = kinds[i % len(kinds)]
        seq = indel_read(amplicon(i), five, three, reverse=(i % 2 == 1))
        recs.append((seq, hifi_quality(len(seq), i)))                                       # a read with an indel has its own length
    path = tmp_path / "indels.fastq.gz"
    write_fastq(path, recs)
    row = lr.describe_file(path, locator())
    assert row["both_forward"] + row["both_reverse"] == 160 and row["share_both"] == 1.0 and row["five_only"] == row["three_only"] == row["none"] == 0
    assert row["len_min"] == INSERT and row["len_max"] == INSERT and row["len_median"] == INSERT       # the amplicon is cut between the primers wherever the edit was
    assert row["share_five"] == 1.0 and row["share_three"] == 1.0 and row["both_reverse"] == 80


def test_a_primer_with_more_edits_than_the_budget_is_not_found_and_the_two_primers_are_counted_separately(tmp_path):
    random.seed(22)
    recs = []
    for i in range(40):
        ins = amplicon(i)
        seq = read(ins)
        recs.append((seq, hifi_quality(len(seq), i)))                                          # both primers, intact
    for i in range(20):                                                                        # the 5' primer wrecked: 6 of its 20 bases changed
        ins = amplicon(100 + i)
        f = list("".join(random.choice(lr.IUPAC[c]) for c in FWD))
        for k in (1, 4, 7, 10, 13, 16):
            f[k] = "ACGT"[("ACGT".index(f[k]) + 1) % 4] if f[k] in "ACGT" else "A"
        seq = "".join(f) + ins + lr.revcomp("".join(random.choice(lr.IUPAC[c]) for c in REV))
        recs.append((seq, hifi_quality(len(seq), 100 + i)))
    for i in range(10):                                                                        # the 3' primer missing altogether
        ins = amplicon(200 + i)
        seq = "".join(random.choice(lr.IUPAC[c]) for c in FWD) + ins
        recs.append((seq, hifi_quality(len(seq), 200 + i)))
    path = tmp_path / "budget.fastq.gz"
    write_fastq(path, recs)
    row = lr.describe_file(path, locator())
    assert row["reads"] == 70 and row["both_forward"] == 40 and row["three_only"] == 20 and row["five_only"] == 10 and row["none"] == 0
    assert row["share_five"] == pytest.approx(50 / 70) and row["share_three"] == pytest.approx(60 / 70)          # the 5' primer was found in 40 + 10 reads, the 3' primer in 40 + 20


def test_the_result_does_not_depend_on_how_the_reads_are_split_into_chunks(tmp_path, monkeypatch):
    random.seed(23)
    recs = []
    for i in range(53):
        five, three = [("del", ""), ("", "ins"), ("", ""), ("sub", "sub")][i % 4]
        seq = indel_read(amplicon(i), five, three, reverse=(i % 3 == 0))
        recs.append((seq, hifi_quality(len(seq), i)))
    recs += [("".join(random.choice("ACGT") for _ in range(900)), hifi_quality(900, 999)) for _ in range(4)]       # none
    path = tmp_path / "chunks.fastq.gz"
    write_fastq(path, recs)
    reference = lr.describe_file(path, locator())
    for size in (1, 7, 53, 57, 5000):
        monkeypatch.setattr(lr, "CHUNK", size)
        assert lr.describe_file(path, locator()) == reference, size                                                  # also the last, shorter chunk is processed
    singly = [locator().locate(seq) for seq, _ in recs]
    assert singly == locator().locate_many([seq for seq, _ in recs])                                                 # one read at a time or all at once
    assert reference["reads"] == 57 and reference["none"] == 4


def test_a_short_read_gives_no_primer_and_bases_that_cannot_match_never_do():
    loc = lr.Locator("A" * 20, "G" * 20)                                                    # primers that would match any padding made of A or G
    assert loc.locate("C" * 30) == (None, "none", 0, 30)
    assert loc.locate_many(["C" * 5, "T" * 12]) == [(None, "none", 0, 5), (None, "none", 0, 12)]
    codes = np.full((1, 30), 4, dtype=np.uint8)
    assert int(lr.edit_matches(codes, lr.Primer("A" * 20))[0][0]) == 20                    # nothing to match at all: every primer base is an edit
