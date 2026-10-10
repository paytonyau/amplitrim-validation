"""Coverage of the amplicon range and the selection rule."""
import math

import pytest

from amplitrim.ranking import (Candidate, coverage_floor, frac_unmergeable, select_recommendation,
                               to_dataframe)


def cand(trim_f=222, trim_r=226, yield_est=0.99, overlap=18.0, diffs=0.2, covers=True):
    return Candidate(trim_f=trim_f, trim_r=trim_r, ee_pass_rate=1.0, merge_rate_strict=0.8,
                     merge_rate_tolerant=yield_est, mean_overlap=overlap, mean_overlap_diffs=diffs,
                     max_mergeable_insert=trim_f + trim_r - 12, covers_longest=covers,
                     yield_est=yield_est)


def test_frac_unmergeable():
    inserts = [420] * 90 + [434] * 10
    assert frac_unmergeable(inserts, 433) == pytest.approx(0.10)
    assert frac_unmergeable(inserts, 434) == 0.0
    assert math.isnan(frac_unmergeable([], 434))


@pytest.mark.parametrize("inserts, max_frac, expected_longest", [
    ([420] * 90 + [434] * 10, 0.0, 434),        # must cover everything
    ([420] * 90 + [434] * 10, 0.10, 420),       # 10% may be longer than the longest merged
    ([420] * 995 + [440] * 5, 0.01, 420),       # a 0.5% tail is allowed to be lost
    ([430], 0.01, 430),
])
def test_coverage_floor_is_the_shortest_total_that_covers(inserts, max_frac, expected_longest):
    floor, longest = coverage_floor(inserts, 12, max_frac)
    assert (floor, longest) == (expected_longest + 12, expected_longest)
    assert frac_unmergeable(inserts, floor - 12) <= max_frac
    assert frac_unmergeable(inserts, floor - 13) > max_frac or floor - 13 < inserts[0]


def test_candidate_that_does_not_cover_is_never_chosen():
    assert select_recommendation([cand(covers=False, yield_est=1.0), cand(yield_est=0.5)]).yield_est == 0.5


def test_shorter_overlap_wins_within_the_yield_tolerance():
    long_overlap, short_overlap = cand(yield_est=1.0, overlap=24.0), cand(yield_est=0.995, overlap=14.0)
    assert select_recommendation([long_overlap, short_overlap], 0.01) is short_overlap
    assert select_recommendation([long_overlap, short_overlap], 0.001) is long_overlap


def test_higher_yield_wins_an_exact_tie_on_overlap():
    low, high = cand(yield_est=0.96, overlap=14.0), cand(yield_est=1.0, overlap=14.0)
    assert select_recommendation([low, high], 0.05) is high


def test_overlaps_within_two_nt_count_as_equal_so_the_higher_yield_wins():
    degenerate, balanced = cand(trim_f=206, trim_r=60, yield_est=0.994, overlap=13.9), cand(trim_f=133, trim_r=133, yield_est=0.999, overlap=13.0)
    assert select_recommendation([degenerate, balanced], 0.01) is balanced
    clearly_longer = cand(yield_est=1.0, overlap=16.0)                  # 2+ nt more overlap is not a tie
    assert select_recommendation([clearly_longer, balanced], 0.01) is balanced


def test_equal_yield_and_overlap_prefer_the_more_balanced_split():
    skewed, balanced = cand(yield_est=0.99, overlap=13.0), cand(yield_est=0.99, overlap=13.0)
    skewed.balance, balanced.balance = 0.6, 0.0
    assert select_recommendation([skewed, balanced], 0.01) is balanced


def test_nothing_eligible_gives_none_and_unmeasurable_candidates_are_skipped():
    assert select_recommendation([]) is None
    assert select_recommendation([cand(covers=False)]) is None
    assert select_recommendation([cand(overlap=float("nan")), cand(yield_est=0.0)]) is None


def test_table_has_plain_column_names():
    df = to_dataframe([cand()])
    assert list(df.columns) == ["truncLenF", "truncLenR", "total_length", "estimated_yield", "pass_filter",
                                "merge_rate", "merge_rate_strict", "mean_overlap",
                                "longest_mergeable_insert", "covers_longest_amplicon", "suggested_by", "recommended"]
    assert df.loc[0, "total_length"] == 448


def _row(ee, yield_est, ee_f, ee_r):
    return {"ee": ee, "yield_est": yield_est, "pass_rate": yield_est + 0.003, "ee_f": ee_f, "ee_r": ee_r}


def test_quality_first_keeps_the_default_unless_too_few_pairs_remain():
    from amplitrim.ranking import quality_first_suggestion
    rows = [_row((1, 1), 0.46, 0.22, 0.40), _row((2, 2), 0.65, 0.40, 0.66), _row((2, 3), 0.72, 0.40, 0.88), _row((3, 3), 0.76, 0.54, 0.88),
            _row((2, 6), 0.81, 0.38, 1.54), _row((3, 5), 0.83, 0.53, 1.27), _row((5, 5), 0.87, 0.75, 1.26)]
    pick, reached = quality_first_suggestion(rows, 0.5)
    assert pick["ee"] == (2, 2) and reached                              # 65% is above the 50% floor: the default stays
    pick, reached = quality_first_suggestion(rows, 0.70)
    assert pick["ee"] == (2, 3) and reached                              # the least noisy row that reaches 70%
    pick, _ = quality_first_suggestion(rows, 0.80)
    assert pick["ee"] == (3, 5)                                          # (2,6) also reaches 80% but is beaten by (3,5): never suggested
    pick, reached = quality_first_suggestion(rows, 0.95)
    assert not reached and pick["ee"] == (5, 5)                          # nothing reaches 95%: the row that keeps the most
    assert quality_first_suggestion([_row((1, 1), 0.9, 0.1, 0.2)], 0.5) == (None, False)      # never stricter than 2
    assert quality_first_suggestion([], 0.5) == (None, False)


def test_a_row_is_beaten_only_when_another_is_better_on_both_counts():
    from amplitrim.ranking import dominators
    rows = [{"pass_rate": 0.66, "ee_f": 0.4, "ee_r": 0.66}, {"pass_rate": 0.81, "ee_f": 0.38, "ee_r": 1.54},
            {"pass_rate": 0.84, "ee_f": 0.53, "ee_r": 1.27}, {"pass_rate": 0.88, "ee_f": 0.75, "ee_r": 1.26},
            {"pass_rate": 0.70, "ee_f": 0.4, "ee_r": 0.66}, {"pass_rate": None, "ee_f": None, "ee_r": None}]
    beaten = dominators(rows)
    assert [len(b) for b in beaten] == [1, 1, 0, 0, 0, 0]          # row 0 is beaten by row 4 (same noise, more pairs); row 1 by row 2
    assert beaten[0][0] is rows[4] and beaten[1][0] is rows[2]
    assert dominators([{"pass_rate": 0.8, "ee_f": 1.0, "ee_r": 1.0}] * 2) == [[], []]       # identical rows do not beat each other
