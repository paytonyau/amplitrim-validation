"""validation/overmerge_real.py: false joins counted against a full-length reference, and the expected errors in the overlap. vsearch is a stand-in whose mistakes the test controls."""
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("overmerge_real", HERE.parent / "validation" / "overmerge_real.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fastq_records(path):
    lines = Path(path).read_text().splitlines()
    return [(lines[i][1:], len(lines[i + 1])) for i in range(0, len(lines), 4)]


def make_merge(true_l, wrong=(), spurious=(), strict_drops_multiples_of=10):
    """vsearch as the test wants it: a pair merges when its reads span the amplicon plus the minimum overlap, at the true length; `wrong` pairs are merged at the wrong length when they
    merge at all (a false join), `spurious` pairs are merged although the reads cannot span the amplicon, and with no difference allowed every pair whose number is a multiple of
    `strict_drops_multiples_of` is refused. The reference merge (tag 'ref') makes no mistake."""
    calls = []

    def merge(r1, r2, out_dir, tag, min_overlap, max_diffs, setting, keep_merged=True, threads=1):
        calls.append((tag, max_diffs, setting))
        out = Path(out_dir) / f"{tag}.{setting}.merged.fasta"
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as fh:
            for (name, l1), (_, l2) in zip(fastq_records(r1), fastq_records(r2)):
                i = int(name[1:])
                L = true_l[i]
                spans = l1 + l2 - L >= min_overlap
                if max_diffs == 0 and i % strict_drops_multiples_of == 0:
                    continue
                if spans:
                    length = L + 7 if (tag != "ref" and i in wrong) else L
                elif tag != "ref" and i in spurious:
                    length = L + 3
                else:
                    continue
                seq = "A" * length
                fh.write(f">{name}\n" + "\n".join(seq[k:k + 60] for k in range(0, length, 60)) + "\n")          # wrapped: the parser must join the lines
        return SimpleNamespace(merged_path=out)
    merge.calls = calls
    return merge


@pytest.fixture(scope="module")
def sim_pairs(tmp_path_factory):
    from simulate import write_dataset
    from amplitrim.subsample import read_fastq
    from amplitrim.trimmer import ReadPairs
    d = tmp_path_factory.mktemp("om")
    write_dataset(d, n_samples=2, pairs_per_sample=400, seed=3)
    r1, r2 = [], []
    for s in ("sampleA", "sampleB"):
        r1 += list(read_fastq(d / "primer_trimmed" / f"{s}_R1.fastq.gz"))
        r2 += list(read_fastq(d / "primer_trimmed" / f"{s}_R2.fastq.gz"))
    truth = [int(x.split("\t")[2]) for x in (d / "truth.tsv").read_text().splitlines()[1:]]
    for name, recs in (("subsample_R1.fastq.gz", r1), ("subsample_R2.fastq.gz", r2)):                                  # the pooled files that build_subsample would write
        with gzip.open(d / name, "wt", newline="\n") as fh:
            fh.write("".join(f"{h}\n{s}\n{p}\n{q}\n" for h, s, p, q in recs))
    return SimpleNamespace(pairs=ReadPairs(r1, r2), truth=truth, dir=d)


def test_the_fasta_is_read_across_wrapped_lines_and_cells_are_parsed(tmp_path):
    om = load()
    (tmp_path / "m.fasta").write_text(">p3\nACGT\nAC\n>p10 extra words\nAAAAAAAAAA\n>p11\nA\n")
    assert om.read_lengths(tmp_path / "m.fasta") == {3: 6, 10: 10, 11: 1}
    assert om.parse_cell("120/118") == (120, 118) and om.cell_text((130, 108)) == "130/108"
    g = om.grid((120, 118))
    assert len(g) == 9 and (110, 108) in g and (130, 128) in g and (120, 118) in g and len(set(g)) == 9


def test_the_reference_is_the_full_length_merge_of_the_pairs_with_a_long_overlap(sim_pairs, tmp_path):
    om = load()
    t = sim_pairs.truth
    merge = make_merge(t)
    ref = om.build_reference(sim_pairs.pairs, merge, tmp_path, ref_overlap=20)
    assert set(ref) == set(range(len(t))) and all(ref[i] == t[i] for i in ref)                                       # reads of 231 + 227 overlap by 24 to 35 nt
    lengths = sorted(set(t))
    mid = next(L for L in lengths if L + 1 in lengths)                                                                # a length whose next length also occurs: its pairs overlap by exactly `edge` nt, the next ones by one less
    edge = 231 + 227 - mid
    refe = om.build_reference(sim_pairs.pairs, merge, tmp_path / "b", ref_overlap=edge)
    want = {i for i, L in enumerate(t) if 231 + 227 - L >= edge}
    assert set(refe) == want and 0 < len(want) < len(t)
    assert any(t[i] == mid for i in refe)                                                                             # a pair that overlaps by exactly the least counts (at least, not more than)
    assert not any(t[i] == mid + 1 for i in refe) and all(231 + 227 - t[i] >= edge for i in refe)                      # one nucleotide less overlap does not
    assert merge.calls[0][0] == "ref" and merge.calls[0][1] is None                                                  # the reference is a tolerant merge of the full reads


def test_correct_wrong_and_unmerged_pairs_are_counted_against_the_reference(sim_pairs, tmp_path):
    om = load()
    t, pairs, cell = sim_pairs.truth, sim_pairs.pairs, (225, 220)                                                    # 445 nt: the amplicons of 434 nt cannot merge (445 - 434 < 12)
    keep = pairs._masks(*cell, 2.0)[4]
    ref = {i: L for i, L in enumerate(t)}
    geometric = [i for i in range(len(t)) if keep[i] and sum(cell) - t[i] >= 12]
    short = [i for i in range(len(t)) if keep[i] and sum(cell) - t[i] < 12]
    assert geometric and short                                                                                       # the cut separates the pairs that can merge from those that cannot
    wrong = set(geometric[::7])
    spurious = set(short[:5])
    row = om.analyse_cell(pairs, ref, cell, make_merge(t, wrong=wrong, spurious=spurious), tmp_path)
    n = len(ref)
    assert row["pass_ref"] == pytest.approx(int(keep.sum()) / n) and row["total"] == 445
    assert row["correct_tol"] == pytest.approx((len(geometric) - len(wrong)) / n)                                    # merged at the reference length
    assert row["wrong_tol"] == pytest.approx((len(wrong) + len(spurious)) / n)                                       # a wrong offset, and a merge that the reads cannot span: both false joins
    assert row["wrong_of_merged_tol"] == pytest.approx((len(wrong) + len(spurious)) / (len(geometric) - len(wrong) + len(wrong) + len(spurious)))
    ok_strict = [i for i in geometric if i % 10 != 0 and i not in wrong]
    wrong_strict = [i for i in list(wrong) + list(spurious) if i % 10 != 0]
    assert row["correct_strict"] == pytest.approx(len(ok_strict) / n) and row["wrong_strict"] == pytest.approx(len(wrong_strict) / n)    # no difference allowed: a different set merges


def test_pairs_outside_the_reference_and_pairs_that_fail_the_filter_do_not_count(sim_pairs, tmp_path):
    om = load()
    t, pairs, cell = sim_pairs.truth, sim_pairs.pairs, (225, 220)
    ref = {i: t[i] for i in range(0, len(t), 2)}                                                                      # half of the pairs have no reliable reference
    keep = pairs._masks(*cell, 2.0)[4]
    row = om.analyse_cell(pairs, ref, cell, make_merge(t), tmp_path)
    want = sum(1 for i in ref if keep[i] and sum(cell) - t[i] >= 12)
    assert row["reference_pairs"] == len(ref) and row["correct_tol"] == pytest.approx(want / len(ref)) and row["wrong_tol"] == 0
    assert 0 < row["pass_ref"] < 1                                                                                    # some reference pairs fail the maxEE 2 filter at this cut


def test_the_errors_in_the_overlap_are_those_of_the_last_bases_of_both_reads(tmp_path):
    om = load()
    from amplitrim.trimmer import ReadPairs
    def rec(i, quals): return (f"p{i}", "A" * len(quals), "+", "".join(chr(33 + q) for q in quals))
    # cut 60/60 on amplicons of 100 nt: the overlap is 20 nt, the last 20 bases of each cut read. Qualities: 40 bases of Q40, then 20 of the quality given.
    def pair(q_f, q_r): return rec(0, [40] * 40 + [q_f] * 20 + [40] * 20), rec(0, [40] * 40 + [q_r] * 20 + [40] * 20)
    specs = [(20, 13), (40, 40), (30, 30)]                                                                            # forward and reverse quality of the overlap, per pair
    recs = [pair(*s) for s in specs]
    recs1, recs2 = [r[0] for r in recs], [r[1] for r in recs]
    recs1.append(rec(3, [40] * 60 + [40] * 40)); recs2.append(rec(3, [40] * 100))                                    # a pair whose amplicon (50 nt) is shorter than the cut: read-through
    pairs = ReadPairs(recs1, recs2)
    ref = {0: 100, 1: 100, 2: 100, 3: 50}
    row = om.analyse_cell(pairs, ref, (60, 60), make_merge([100, 100, 100, 50]), tmp_path)
    e = lambda q, n: n * 10 ** (-q / 10)                                                                              # noqa: E731
    want = [e(20, 20) + e(13, 20), e(40, 20) * 2, e(30, 20) * 2]
    assert row["overlap_pairs"] == 3                                                                                  # the read-through pair (overlap 70 > 60) is left out
    assert row["overlap_ee_mean"] == pytest.approx(np.mean(want)) and row["overlap_ee_ge1"] == pytest.approx(1 / 3)  # only the first pair carries one expected error or more
    assert want[0] > 1 > want[2] > want[1]


def test_the_verdicts_follow_the_registered_criteria():
    om = load()
    def row(wrong_of_merged, wrong, ee, n=1000): return {"reference_pairs": n, "wrong_of_merged_tol": wrong_of_merged, "wrong_tol": wrong, "overlap_pairs": 50, "overlap_ee_mean": ee}
    a, b, now = om.EARLIER_RECOMMENDED, om.EARLIER_BEST, (124, 115)
    rows = {a: row(0.004, 0.003, 0.50), b: row(0.006, 0.0035, 0.30), now: row(0.010, 0.004, 0.40)}
    v = {p: ok for p, _, ok in om.verdicts(rows, now)}
    assert v == {"P31": True, "P32": True, "P33": True}                                                                # 1.00% is at the limit and passes; 0.05 points apart; a cleaner overlap
    rows[now] = row(0.0101, 0.004, 0.40)
    assert {p: ok for p, _, ok in om.verdicts(rows, now)}["P31"] is False                                              # 1.01%
    rows[now] = row(0.002, 0.004, 0.40); rows[b] = row(0.006, 0.0129, 0.30)
    assert {p: ok for p, _, ok in om.verdicts(rows, now)}["P32"] is True                                               # 0.99 points apart
    rows[b] = row(0.006, 0.0131, 0.30)
    assert {p: ok for p, _, ok in om.verdicts(rows, now)}["P32"] is False                                              # 1.01 points apart
    rows[b] = row(0.006, 0.003, 0.60)
    assert {p: ok for p, _, ok in om.verdicts(rows, now)}["P33"] is False                                              # the overlap of the best cell is noisier
    assert {p: ok for p, _, ok in om.verdicts({a: rows[a]}, None)} == {"P31": True, "P32": None, "P33": None}         # a cell is missing: not decided
    assert all(ok is None for _, _, ok in om.verdicts({}, None)) and "wrong_of_merged_tol" not in om.verdicts({}, None)[0][1]


def test_the_check_runs_end_to_end_and_writes_its_report(sim_pairs, tmp_path):
    om = load()
    reads = tmp_path / "reads"
    reads.mkdir()
    cfg = tmp_path / "paths.json"
    cfg.write_text(json.dumps({"extras": {"v4": {"ready": str(reads).replace("\\", "/")}}}))
    merge = make_merge(sim_pairs.truth, wrong={5, 6})
    seen = {}

    def build(folder, out_dir):
        seen["build"] = Path(folder)
        return SimpleNamespace(r1=sim_pairs.dir / "subsample_R1.fastq.gz", r2=sim_pairs.dir / "subsample_R2.fastq.gz")
    rc = om.main(["--config", str(cfg), "--ref-overlap", "20", "--out", str(tmp_path / "o")], merge=merge, build=build, recommend=lambda folder: (225, 220))
    assert rc == 0 and seen["build"] == reads
    rep = (tmp_path / "o" / "OVERMERGE_REPORT.md").read_text()
    for need in ("- **P31**", "- **P32**", "- **P33**", "## The cuts (800 pairs analysed; 800 reference pairs", "120/118 (round B6 recommendation)", "130/108 (best cell of the B6 grid)", "225/220 (recommended now)",
                 "What vsearch sees of the gap of round B6"):
        assert need in rep, need
    rows = [ln.split("\t") for ln in (tmp_path / "o" / "overmerge.tsv").read_text().splitlines()]
    cells = {r[0] for r in rows[1:]}
    assert rows[0][0] == "cell" and {"120/118", "130/108", "225/220", "110/108", "130/128"} <= cells                      # the grid of the study and the grid around the recommendation
    assert "235/230" not in cells and "225/230" not in cells                                                             # a cut longer than the reads is not run (231 and 227)


def test_missing_things_are_said_plainly(sim_pairs, tmp_path, capsys, monkeypatch):
    om = load()
    assert om.main(["--config", str(tmp_path / "absent.json"), "--out", str(tmp_path / "o")]) == 2 and "No reads folder for 'v4'" in capsys.readouterr().out
    reads = tmp_path / "reads"
    reads.mkdir()
    cfg = tmp_path / "paths.json"
    cfg.write_text(json.dumps({"extras": {"v4": {"ready": str(reads).replace("\\", "/")}}}))
    monkeypatch.setattr("shutil.which", lambda name: None)
    assert om.main(["--config", str(cfg), "--out", str(tmp_path / "o")]) == 2 and "vsearch was not found" in capsys.readouterr().out
    def build(folder, out_dir): return SimpleNamespace(r1=sim_pairs.dir / "subsample_R1.fastq.gz", r2=sim_pairs.dir / "subsample_R2.fastq.gz")
    assert om.main(["--config", str(cfg), "--ref-overlap", "500", "--out", str(tmp_path / "o")], merge=make_merge(sim_pairs.truth), build=build, recommend=lambda f: (225, 220)) == 2
    assert "No pair overlaps by 500 nt at full length" in capsys.readouterr().out
