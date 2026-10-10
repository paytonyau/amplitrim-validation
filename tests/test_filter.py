"""DADA2-style pair filtering. The vectorised filter is checked against a plain per-read reference."""
import gzip
import random

import pytest

from amplitrim.subsample import read_fastq
from amplitrim.trimmer import ReadPairs

ERR = [10 ** (-q / 10) for q in range(94)]


def write_fq(path, recs):
    with gzip.open(path, "wt", newline="\n") as fh:
        for i, (seq, qual) in enumerate(recs):
            fh.write(f"@p{i}\n{seq}\n+\n{qual}\n")


def reference(recs1, recs2, tf, tr, max_ee):
    """Straightforward per-read implementation of the filter."""
    counts = dict(short=0, low_q=0, has_n=0, ee=0)
    kept = []
    for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(recs1, recs2):
        if len(s1) < tf or len(s2) < tr:
            counts["short"] += 1
            continue
        s1, q1, s2, q2 = s1[:tf], q1[:tf], s2[:tr], q2[:tr]
        if min(map(ord, q1)) - 33 <= 2 or min(map(ord, q2)) - 33 <= 2:
            counts["low_q"] += 1
        elif "N" in s1 or "N" in s2:
            counts["has_n"] += 1
        elif sum(ERR[ord(c) - 33] for c in q1) > max_ee or sum(ERR[ord(c) - 33] for c in q2) > max_ee:
            counts["ee"] += 1
        else:
            kept.append(f"{h1}\n{s1}\n{p1}\n{q1}\n")
    return counts, kept


def random_reads(rng, n, lengths):
    recs = []
    for _ in range(n):
        length = rng.choice(lengths)
        seq = "".join("N" if rng.random() < 0.003 else rng.choice("ACGT") for _ in range(length))
        qual = "".join(chr(33 + (2 if rng.random() < 0.004 else rng.randint(8, 40))) for _ in range(length))
        recs.append((seq, qual))
    return recs


def test_each_failure_type_is_counted_and_pairs_stay_in_sync(tmp_path):
    good, q10 = "I" * 20, chr(33 + 10) * 20            # Q10 = 0.1 expected errors per base
    r1 = [("A" * 20, good), ("A" * 20, good), ("A" * 20, good), ("A" * 10, "I" * 10), ("A" * 20, good)]
    r2 = [("C" * 20, good), ("C" * 20, "I" * 10 + "#" * 10), ("N" + "C" * 19, good), ("C" * 20, good),
          ("C" * 20, q10)]
    write_fq(tmp_path / "a_R1.fq.gz", r1)
    write_fq(tmp_path / "a_R2.fq.gz", r2)
    rp = ReadPairs.load(tmp_path / "a_R1.fq.gz", tmp_path / "a_R2.fq.gz")
    st = rp.filter(20, 20, tmp_path / "o1.fq", tmp_path / "o2.fq", max_ee=1.5)
    assert (st.n_in, st.n_low_q, st.n_has_n, st.n_too_short, st.n_ee_fail, st.n_out) == (5, 1, 1, 1, 1, 1)
    assert [r[0] for r in read_fastq(tmp_path / "o1.fq")] == [r[0] for r in read_fastq(tmp_path / "o2.fq")] == ["@p0"]


@pytest.mark.parametrize("seed", range(10))
def test_matches_reference_on_random_reads(tmp_path, seed):
    rng = random.Random(seed)
    n = rng.randint(20, 80)
    write_fq(tmp_path / "a_R1.fq.gz", random_reads(rng, n, [40, 38, 25] if seed % 2 else [40]))
    write_fq(tmp_path / "a_R2.fq.gz", random_reads(rng, n, [40, 35] if seed % 3 else [40]))
    recs1, recs2 = list(read_fastq(tmp_path / "a_R1.fq.gz")), list(read_fastq(tmp_path / "a_R2.fq.gz"))
    rp = ReadPairs(recs1, recs2)
    for tf, tr, max_ee in [(30, 28, 2.0), (40, 40, 0.5), (12, 55, 3.5)]:
        counts, kept = reference(recs1, recs2, tf, tr, max_ee)
        st = rp.filter(tf, tr, tmp_path / "o1.fq", tmp_path / "o2.fq", max_ee)
        assert (st.n_too_short, st.n_low_q, st.n_has_n, st.n_ee_fail) == \
            (counts["short"], counts["low_q"], counts["has_n"], counts["ee"])
        assert st.n_out == len(kept)
        assert (tmp_path / "o1.fq").read_text() == "".join(kept)


def _dataset(tmp_path, seed, n=60):
    rng = random.Random(seed)
    write_fq(tmp_path / "a_R1.fq.gz", random_reads(rng, n, [40, 38, 25]))
    write_fq(tmp_path / "a_R2.fq.gz", random_reads(rng, n, [40, 35]))
    return ReadPairs.load(tmp_path / "a_R1.fq.gz", tmp_path / "a_R2.fq.gz")


@pytest.mark.parametrize("max_ee", [0.5, 2.0])
def test_pass_table_equals_the_filter_at_every_truncation(tmp_path, max_ee):
    rp = _dataset(tmp_path, seed=3)
    table = rp.pass_table(max_ee)
    for tf in range(1, 44, 3):
        for tr in range(1, 44, 3):
            st = rp.filter(tf, tr, tmp_path / "o1.fq", tmp_path / "o2.fq", max_ee)
            assert table.rate(tf, tr) == pytest.approx(st.n_out / st.n_in, abs=1e-12)


def test_pass_rate_only_falls_as_truncation_grows(tmp_path):
    table = _dataset(tmp_path, seed=4).pass_table(2.0)
    for tf in range(1, 40):
        for tr in range(1, 40):
            assert table.rate(tf + 1, tr) <= table.rate(tf, tr) + 1e-12
            assert table.rate(tf, tr + 1) <= table.rate(tf, tr) + 1e-12


def test_best_split_is_the_highest_pass_rate_and_ties_go_to_the_preferred_split(tmp_path):
    table = _dataset(tmp_path, seed=5).pass_table(2.0)
    for total in (30, 45, 60, 75):
        tf, tr, rate = table.best_split(total)
        brute = max(table.rate(f, total - f) for f in range(1, 41) if 0 < total - f <= 40)
        assert (tf + tr, rate) == (total, pytest.approx(brute))
        tied = [f for f in range(1, 41) if 0 < total - f <= 40 and table.rate(f, total - f) >= brute - 1e-12]
        assert table.best_split(total, prefer_f=tied[0])[0] == tied[0]
    assert table.best_split(500) is None


def test_near_best_splits_span_the_range_when_pass_rate_cannot_tell_them_apart(tmp_path):
    clean = tmp_path / "c"
    clean.mkdir()
    write_fq(clean / "a_R1.fq.gz", [("A" * 40, "I" * 40)] * 20)       # every split passes equally
    write_fq(clean / "a_R2.fq.gz", [("C" * 40, "I" * 40)] * 20)
    table = ReadPairs.load(clean / "a_R1.fq.gz", clean / "a_R2.fq.gz").pass_table(2.0)
    splits = table.near_best_splits(60, k=5)
    assert len(splits) == 5 and all(f + r == 60 for f, r in splits)
    assert splits[0][0] == 20 and splits[-1][0] == 40                  # from the shortest to the longest forward
    assert table.near_best_splits(500) == []


def reference_untruncated(recs1, recs2, max_ee, min_len):
    """Plain per-read version of DADA2's filterAndTrim with no truncLen (truncQ = 2, maxN = 0)."""
    kept = []
    for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(recs1, recs2):
        cut = []
        for s, q in ((s1, q1), (s2, q2)):
            n = next((i for i, c in enumerate(q) if ord(c) - 33 <= 2), len(q))
            cut.append((s[:n], q[:n]))
        if any(len(s) < min_len for s, _ in cut):
            continue
        if any("N" in s for s, _ in cut):
            continue
        if any(sum(ERR[ord(c) - 33] for c in q) > max_ee for _, q in cut):
            continue
        kept.append(f"{h1}\n{cut[0][0]}\n{p1}\n{cut[0][1]}\n")
    return kept


@pytest.mark.parametrize("seed", range(8))
def test_no_truncation_filter_matches_the_reference(tmp_path, seed):
    rp = _dataset(tmp_path, seed=20 + seed, n=70)
    recs1, recs2 = rp.recs1, rp.recs2
    for max_ee, min_len in [(2.0, 10), (0.5, 20), (3.5, 1)]:
        kept = reference_untruncated(recs1, recs2, max_ee, min_len)
        st = rp.filter_untruncated(tmp_path / "o1.fq", tmp_path / "o2.fq", max_ee, min_len)
        assert st.n_out == len(kept) and (tmp_path / "o1.fq").read_text() == "".join(kept)
        assert st.n_in == st.n_out + st.n_too_short + st.n_has_n + st.n_ee_fail


def test_written_fastq_never_contains_carriage_returns(tmp_path):
    """On Windows, text-mode writing turns newlines into CRLF unless newline="\\n" is given."""
    rp = _dataset(tmp_path, seed=9)
    rp.filter(10, 10, tmp_path / "o1.fq", tmp_path / "o2.fq", 2.0)
    rp.filter_untruncated(tmp_path / "u1.fq", tmp_path / "u2.fq", 2.0, 5)
    for name in ("o1.fq", "o2.fq", "u1.fq", "u2.fq"):
        assert b"\r" not in (tmp_path / name).read_bytes()


def test_forward_and_reverse_can_have_different_maxEE(tmp_path):
    """maxEE = c(F, R): each read is held to its own limit. Checked against a plain per-read implementation."""
    rng = random.Random(11)
    def rec(n, bad):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        qual = "".join(chr(33 + (rng.choice([5, 8, 12]) if (i > 40 and bad) else rng.choice([30, 35, 38]))) for i in range(60))
        return seq, qual
    recs1 = [rec(i, bad=(i % 3 == 0)) for i in range(300)]
    recs2 = [rec(i, bad=(i % 2 == 0)) for i in range(300)]
    write_fq(tmp_path / "a_R1.fastq.gz", recs1)
    write_fq(tmp_path / "a_R2.fastq.gz", recs2)
    pairs = ReadPairs.load(tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz")
    for tf, tr, ef, er in ((60, 60, 2.0, 6.0), (55, 50, 1.0, 3.0), (60, 45, 4.0, 1.0), (50, 50, 2.0, 2.0)):
        want = sum(1 for (s1, q1), (s2, q2) in zip(recs1, recs2)
                   if min(ord(c) - 33 for c in q1[:tf]) > 2 and min(ord(c) - 33 for c in q2[:tr]) > 2
                   and sum(ERR[ord(c) - 33] for c in q1[:tf]) <= ef and sum(ERR[ord(c) - 33] for c in q2[:tr]) <= er)
        got = pairs.filter(tf, tr, tmp_path / "o1.fq", tmp_path / "o2.fq", (ef, er)).n_out
        assert got == want, (tf, tr, ef, er)
        same = pairs.filter(tf, tr, tmp_path / "p1.fq", tmp_path / "p2.fq", ef if ef == er else (ef, er)).n_out
        assert same == got
        assert abs(pairs.pass_table((ef, er)).rate(tf, tr) - want / len(recs1)) < 1e-12      # the table the recommendation uses


def test_maxEE_values_are_parsed_and_formatted():
    import argparse
    from amplitrim.config import ee_pair, parse_ee
    assert parse_ee("2") == 2.0 and parse_ee("2,6") == (2.0, 6.0) and parse_ee("3,3") == 3.0 and parse_ee("2;5") == (2.0, 5.0)
    assert ee_pair(2) == (2.0, 2.0) and ee_pair((2, 6)) == (2.0, 6.0)
    for bad in ("0", "-1", "a", "1,2,3", ""):
        with pytest.raises(argparse.ArgumentTypeError):
            parse_ee(bad)


def test_the_error_load_of_kept_reads_matches_a_plain_calculation_and_grows_as_maxEE_loosens(tmp_path):
    rng = random.Random(23)
    def rec(bad):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        qual = "".join(chr(33 + (rng.choice([6, 9, 14]) if (i > 35 and bad) else rng.choice([28, 33, 37]))) for i in range(60))
        return seq, qual
    recs1 = [rec(i % 3 == 0) for i in range(400)]
    recs2 = [rec(i % 2 == 0) for i in range(400)]
    write_fq(tmp_path / "a_R1.fastq.gz", recs1)
    write_fq(tmp_path / "a_R2.fastq.gz", recs2)
    pairs = ReadPairs.load(tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz")
    tf, tr = 55, 50
    ee = lambda q, n: sum(ERR[ord(c) - 33] for c in q[:n])
    def want_load(limits):
        kept = [(ee(q1, tf), ee(q2, tr)) for (s1, q1), (s2, q2) in zip(recs1, recs2)
                if min(ord(c) - 33 for c in q1[:tf]) > 2 and min(ord(c) - 33 for c in q2[:tr]) > 2 and ee(q1, tf) <= limits[0] and ee(q2, tr) <= limits[1]]
        return sum(a for a, _ in kept) / len(kept), sum(b for _, b in kept) / len(kept)
    for limits in ((1.0, 1.0), (2.0, 2.0), (2.0, 4.0), (2.0, 8.0), (8.0, 8.0), (1.0, 8.0)):
        got, want = pairs.error_load(tf, tr, limits), want_load(limits)
        assert abs(got[0] - want[0]) < 1e-9 and abs(got[1] - want[1]) < 1e-9, limits
    # loosening ONE read's limit never lowers that read's mean error load (the other read's mean can move either way)
    reverse = [pairs.error_load(tf, tr, (2.0, r))[1] for r in (1.0, 2.0, 4.0, 8.0)]
    forward = [pairs.error_load(tf, tr, (f, 8.0))[0] for f in (1.0, 2.0, 8.0)]
    assert reverse == sorted(reverse) and forward == sorted(forward) and reverse[-1] > reverse[0] and forward[-1] > forward[0]
    assert all(x != x for x in pairs.error_load(tf, tr, 1e-9))                          # nothing kept: nan, not a crash


def test_the_error_profile_matches_a_plain_cumulative_calculation(tmp_path):
    import numpy as np
    rng = random.Random(31)
    recs1 = [("".join(rng.choice("ACGT") for _ in range(50)), "".join(chr(33 + rng.choice([12, 25, 36])) for _ in range(50))) for _ in range(300)]
    recs2 = [("".join(rng.choice("ACGT") for _ in range(40)), "".join(chr(33 + rng.choice([8, 20, 33])) for _ in range(40))) for _ in range(300)]
    write_fq(tmp_path / "a_R1.fastq.gz", recs1)
    write_fq(tmp_path / "a_R2.fastq.gz", recs2)
    prof = ReadPairs.load(tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz").error_profile()
    for key, recs, width in (("f", recs1, 50), ("r", recs2, 40)):
        cum = np.array([np.cumsum([ERR[ord(c) - 33] for c in q]) for _, q in recs])
        assert prof[key].shape == (3, width)
        for pos in (0, 9, width - 1):
            assert np.allclose(prof[key][:, pos], np.percentile(cum[:, pos], (10, 50, 90)))
        assert (np.diff(prof[key][1]) >= 0).all()                  # errors only accumulate along a read


def test_the_reads_a_looser_maxEE_adds_are_noisier_and_their_load_adds_up(tmp_path):
    rng = random.Random(77)
    def rec(bad):
        seq = "".join(rng.choice("ACGT") for _ in range(60))
        qual = "".join(chr(33 + (rng.choice([6, 9, 14]) if (i > 30 and bad) else rng.choice([28, 33, 37]))) for i in range(60))
        return seq, qual
    recs1 = [rec(i % 4 == 0) for i in range(500)]
    recs2 = [rec(i % 2 == 0) for i in range(500)]
    write_fq(tmp_path / "a_R1.fastq.gz", recs1)
    write_fq(tmp_path / "a_R2.fastq.gz", recs2)
    pairs = ReadPairs.load(tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz")
    tf, tr = 55, 55
    base, loose = (2.0, 2.0), (2.0, 5.0)
    n_base = pairs.filter(tf, tr, tmp_path / "b1.fq", tmp_path / "b2.fq", base).n_out
    n_loose = pairs.filter(tf, tr, tmp_path / "l1.fq", tmp_path / "l2.fq", loose).n_out
    add_f, add_r = pairs.added_error_load(tf, tr, loose, base)
    assert add_r > 2.0                                                  # extra reads exceed the baseline limit on the reverse read
    kb, kl = pairs.error_load(tf, tr, base), pairs.error_load(tf, tr, loose)
    assert abs(kl[1] * n_loose - (kb[1] * n_base + add_r * (n_loose - n_base))) < 1e-6      # the loads add up exactly
    assert abs(kl[0] * n_loose - (kb[0] * n_base + add_f * (n_loose - n_base))) < 1e-6
    assert all(x != x for x in pairs.added_error_load(tf, tr, base, base)) and all(x != x for x in pairs.added_error_load(tf, tr, (1.0, 1.0), base))


def test_the_load_of_the_pairs_dadas_untruncated_filter_keeps_matches_a_plain_python_count(tmp_path):
    rng = random.Random(404)
    def rec(n=90):
        seq = "".join("N" if rng.random() < 0.003 else rng.choice("ACGT") for _ in range(n))
        q = [rng.choice([14, 22, 30, 36, 38]) for _ in range(n)]
        if rng.random() < 0.4:                                            # a quality collapse: the read is cut at its first Q <= 2 base
            start = rng.randrange(40, n)
            q[start:] = [2] * (n - start)
        for i in range(n):
            if rng.random() < 0.02:
                q[i] = rng.choice([4, 8, 10])
        return seq, "".join(chr(33 + x) for x in q)
    recs1, recs2 = [rec() for _ in range(300)], [rec() for _ in range(300)]
    write_fq(tmp_path / "a_R1.fastq.gz", recs1)
    write_fq(tmp_path / "a_R2.fastq.gz", recs2)
    pairs = ReadPairs.load(tmp_path / "a_R1.fastq.gz", tmp_path / "a_R2.fastq.gz")

    def plain(ee, min_len=50):
        out = {}
        for i, ((s1, q1), (s2, q2)) in enumerate(zip(recs1, recs2)):
            cut = lambda q: next((j for j, c in enumerate(q) if ord(c) - 33 <= 2), len(q))
            l1, l2 = cut(q1), cut(q2)
            if l1 < min_len or l2 < min_len or "N" in s1[:l1] or "N" in s2[:l2]:
                continue
            e1, e2 = (sum(10 ** (-(ord(c) - 33) / 10) for c in q[:l]) for q, l in ((q1, l1), (q2, l2)))
            if e1 <= ee[0] and e2 <= ee[1]:
                out[i] = (l1, l2, e1, e2)
        return out
    for ee in ((1.0, 1.0), (2.0, 2.0), (2.0, 5.0), (5.0, 5.0)):
        kept, base = plain(ee), plain((2.0, 2.0))
        got = pairs.untruncated_load(ee, 50, baseline=(2.0, 2.0))
        assert got["n_kept"] == len(kept) and len(kept) > 0
        assert abs(got["ee_f"] - sum(v[2] for v in kept.values()) / len(kept)) < 1e-9 and abs(got["ee_r"] - sum(v[3] for v in kept.values()) / len(kept)) < 1e-9
        assert abs(got["len_f"] - sum(v[0] for v in kept.values()) / len(kept)) < 1e-9 and abs(got["len_r"] - sum(v[1] for v in kept.values()) / len(kept)) < 1e-9
        new = [kept[i] for i in kept if i not in base]
        if new:
            assert abs(got["add_f"] - sum(v[2] for v in new) / len(new)) < 1e-9 and abs(got["add_r"] - sum(v[3] for v in new) / len(new)) < 1e-9
        else:
            assert got["add_f"] != got["add_f"] and got["add_r"] != got["add_r"]              # nothing added: nan, not zero
