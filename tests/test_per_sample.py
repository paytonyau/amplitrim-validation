"""per_sample.py: one row per sample, and plain findings."""
import importlib.util
from pathlib import Path

from conftest import need_vsearch

VALIDATION = Path(__file__).parent.parent / "validation"
spec = importlib.util.spec_from_file_location("per_sample", VALIDATION / "per_sample.py")
ps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ps)


@need_vsearch
def test_every_sample_gets_a_row_with_sensible_numbers(reads):
    rows = ps.per_sample(reads, 220, 220, pairs_per_sample=400)
    assert len(rows) == 2 and {r["sample"] for r in rows} == {"sampleA", "sampleB"}
    for r in rows:
        assert 0 < r["estimated_yield"] <= 1 and r["pairs_used"] <= 400 + 1
        assert abs(r["estimated_yield"] - r["pass_filter"] * r["merge_rate"]) < 2e-4


def test_weak_samples_are_named_and_the_pooled_estimate_is_compared():
    rows = [{"sample": f"s{i}", "pairs_used": 100, "estimated_yield": y} for i, y in enumerate([0.9, 0.91, 0.89, 0.6])]
    text = " ".join(ps.summarise(rows, pooled=0.88))
    assert "4 samples: estimated yield median 89.5%" in text and "pooled estimate: 88.0%" in text
    assert "1 sample more than 10 points below the median: s3 (60%)" in text
    assert "No sample is more than 10 points below" in " ".join(ps.summarise(rows[:3]))


def test_a_folder_that_real_data_has_not_filled_gets_a_plain_message(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(ps.sys, "argv", ["per_sample.py", str(tmp_path / "study")])
    assert ps.main() == 1
    assert "Run real_data.py first" in capsys.readouterr().out


def test_the_design_file_names_the_samples_and_groups_them(tmp_path):
    csv_file = tmp_path / "SraRunInfo.csv"
    csv_file.write_text("Run,SampleName,spots\nSRR1,G1L1,100\nSRR2,G1L2,200\nSRR3,G1\'\'L1,300\nSRR4,other,5\n", encoding="utf-8")
    design = ps.load_design(csv_file)
    assert design["SRR2"] == ("G1L2", 200)
    assert [ps.group_of(n) for n in ("G1L1", "G1\'\'L3", "other", "L2")] == ["G1", "G1\'\'", "other", "L2"]
    rows = [{"sample": "SRR1", "estimated_yield": 0.9}, {"sample": "SRR2", "estimated_yield": 0.8},
            {"sample": "SRR3", "estimated_yield": 0.7}, {"sample": "SRR9", "estimated_yield": 0.5}]
    ps.add_design(rows, design)
    assert [r["group"] for r in rows] == ["G1", "G1", "G1\'\'", ""] and rows[1]["archive_reads"] == 200
    assert ps.summarise_groups(rows) == ["G1: 2 samples, median 85.0%, range 80.0%-90.0%", "G1\'\': 1 samples, median 70.0%, range 70.0%-70.0%"]


def test_rank_correlation_handles_order_reversal_and_ties():
    assert ps.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == 1.0
    assert ps.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0
    assert abs(ps.spearman([1, 2, 2, 3], [1, 2, 2, 3]) - 1.0) < 1e-12 and ps.spearman([1, 1, 1], [1, 2, 3]) is None


def test_dada2_per_sample_results_are_joined_and_compared_with_the_estimate_and_depth(tmp_path):
    track = tmp_path / "t.tsv"
    # ten samples: more pairs given -> higher merged share; estimates fall as DADA2's share rises
    lines = ["sample\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim"]
    for i in range(10):
        n = 1000 * (i + 1)
        merged = round(n * (0.40 + 0.03 * i))
        lines.append(f"S{i}\t{n}\t{n}\t{n}\t{n}\t{merged}\t{merged}")
    track.write_text("\n".join(lines) + "\n", encoding="utf-8")
    rows = [{"sample": f"S{i}", "estimated_yield": 0.99 - 0.01 * i, "strict_yield": 0.5 + 0.05 * i} for i in range(10)] + [{"sample": "unlisted", "estimated_yield": 0.5, "strict_yield": 0.1}]
    ps.add_dada2(rows, track)
    assert [round(r["dada2_merged"], 2) for r in rows[:10]] == [round(0.40 + 0.03 * i, 2) for i in range(10)] and rows[10]["dada2_merged"] is None
    text = " ".join(ps.summarise_dada2(rows))
    assert "DADA2 merged 40.0% to 67.0% of the pairs it was given across 10 samples" in text
    assert "sample by sample (n = 10): -1.00" in text and "(no mismatch allowed in the overlap, on raw reads) and DADA2's merged share: +1.00" in text and "pairs DADA2 was given and its merged share: +1.00" in text
    assert "at least 10" in " ".join(ps.summarise_dada2(rows[:6]))                 # six samples are too few for a rank correlation



def test_the_ena_report_works_as_a_design_file_and_its_names_group_by_dilution(tmp_path):
    f = tmp_path / "runs_info.tsv"
    f.write_text("run_accession\tsample_alias\tlibrary_name\tread_count\n"
                 "ERR1\t113-KitA-1-3-2_S113_L001\tmock\t140822\nERR2\t099-KitA-1-1-4_S99_L001\tmock\t35944\nERR3\t137-KitA-NK-1_S137_L001\tneg\t8328\n", encoding="utf-8")
    design = ps.load_design(f)
    assert design["ERR1"] == ("113-KitA-1-3-2_S113_L001", 70411) and design["ERR2"][1] == 17972      # ENA counts both mates: pairs are half
    pattern = r"^\d+-KitA-1-\d-(\d)_"
    assert [ps.group_of(design[r][0], pattern) for r in ("ERR1", "ERR2")] == ["2", "4"]            # the dilution digit (log10 of the cells)
    assert ps.group_of(design["ERR3"][0], pattern) == "137-KitA-NK-1_S137_L001"                       # no match: the whole name

