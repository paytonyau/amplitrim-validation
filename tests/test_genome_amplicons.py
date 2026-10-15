"""validation/genome_amplicons.py on simulated genomes whose 16S copies, strands, primer mismatches and distances are known exactly."""
import ast
import csv
import gzip
import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
import genome_amplicons as ga  # noqa: E402
import longread_describe as lr  # noqa: E402
import make_reference as mr  # noqa: E402

FWD, REV = "AGRGTTYGATYMTGGCTCAG", "RGYTACCTTGTTACGACTT"            # 27F and 1492R
REV_SITE = lr.revcomp(REV)                                          # as it appears on the strand of the gene


def amplicon(seed: int, length: int = 1450) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def site(pattern: str, mismatches: int, rng: random.Random) -> str:
    """A concrete sequence matching an IUPAC pattern, with exactly `mismatches` positions changed to a base the pattern does not allow."""
    s = [rng.choice(lr.IUPAC[c]) for c in pattern]
    for i in rng.sample(range(len(pattern)), mismatches):
        s[i] = rng.choice([b for b in "ACGT" if b not in lr.IUPAC[pattern[i]]])
    return "".join(s)


def segment(amp: str, rng, mm_f: int = 0, mm_r: int = 0) -> str:
    return site(FWD, mm_f, rng) + amp + site(REV_SITE, mm_r, rng)


def genome(copies: list, seed: int = 1, background: int = 2500, spacer=None) -> str:
    """Random background with the given copies in between: each is (amplicon, strand) or (amplicon, strand, forward mismatches, reverse mismatches)."""
    rng = random.Random(seed)
    parts = ["".join(rng.choice("ACGT") for _ in range(background))]
    for c in copies:
        amp, strand = c[0], c[1]
        seg = segment(amp, rng, c[2] if len(c) > 2 else 0, c[3] if len(c) > 3 else 0)
        parts.append(seg if strand == "+" else lr.revcomp(seg))
        parts.append("".join(rng.choice("ACGT") for _ in range(spacer or background)))
    return "".join(parts)


def write_genome(path: Path, sequences, gz: bool = False) -> Path:
    text = "".join(f">contig{i} some description\n" + "\n".join(s[k:k + 70] for k in range(0, len(s), 70)) + "\n" for i, s in enumerate(sequences, 1))
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def write_table(path: Path, rows: list) -> Path:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["atcc_number", "name_on_atcc_page", "share_percent", "accession", "tier"])
        w.writerows(rows)
    return path


def run(tmp_path, files: dict, table_rows=None, **kw):
    folder = tmp_path / "genomes"
    folder.mkdir(exist_ok=True)
    for name, seqs in files.items():
        write_genome(folder / name, seqs, gz=name.endswith(".gz"))
    table = write_table(tmp_path / "table.tsv", table_rows) if table_rows is not None else None
    return ga.run(folder, table, FWD, REV, tmp_path / "out", log=lambda *_: None, **kw)


def fasta_records(text: str) -> dict:
    out, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:]
            out[name] = ""
        else:
            out[name] += line
    return out


A, A_VARIANT, B = amplicon(1), None, amplicon(2)
A_VARIANT = A[:700] + ("G" if A[700] != "G" else "C") + A[701:]


def test_copies_on_both_strands_collapse_into_variants(tmp_path):
    g = genome([(A, "+"), (A, "-"), (A_VARIANT, "+")], seed=3)
    r = run(tmp_path, {"Escherichia_coli_GCF_000000001.1.fna": [g]}, [("ATCC 1", "Escherichia coli", "18", "GCF_000000001.1", "A")])
    recs = fasta_records(r["fasta"])
    label = "Escherichia_coli_ATCC1"
    assert recs == {f"{label}_v1|{label}|2_copies": A, f"{label}_v2|{label}|1_copies": A_VARIANT}
    s = r["strains"][label]
    assert s["copies"] == 3 and len(s["variants"]) == 2 and s["mm_f"] == 0 and s["mm_r"] == 0


def test_a_site_with_three_mismatches_is_found_and_one_with_four_is_not(tmp_path):
    g = genome([(A, "+", 3, 3), (B, "+", 4, 0)], seed=4)
    r = run(tmp_path, {"x_GCF_000000002.1.fna": [g]}, [("ATCC 2", "Strain two", "10", "GCF_000000002.1", "A")])
    s = r["strains"]["Strain_two_ATCC2"]
    assert s["copies"] == 1 and s["mm_f"] == 3 and s["mm_r"] == 3
    assert list(fasta_records(r["fasta"]).values()) == [A]
    r4 = run(tmp_path, {"x_GCF_000000002.1.fna": [genome([(B, "+", 4, 0)], seed=5)]}, [("ATCC 2", "Strain two", "10", "GCF_000000002.1", "A")])
    assert r4["strains"]["Strain_two_ATCC2"]["copies"] == 0 and any("NO amplicon found" in w for w in r4["warnings"])
    assert r4["fasta"] == ""


def test_pairs_too_close_or_too_far_apart_are_not_amplicons(tmp_path):
    g = genome([(amplicon(7, 600), "+"), (amplicon(8, 3000), "+"), (A, "+")], seed=6)
    r = run(tmp_path, {"x_GCF_000000003.1.fna": [g]}, [("ATCC 3", "Three", "1", "GCF_000000003.1", "A")])
    assert list(fasta_records(r["fasta"]).values()) == [A]                                   # 600 and 3,000 nt lie outside 1,000 to 2,000
    r2 = run(tmp_path, {"x_GCF_000000003.1.fna": [g]}, [("ATCC 3", "Three", "1", "GCF_000000003.1", "A")], min_length=500, max_length=3500)
    assert r2["strains"]["Three_ATCC3"]["copies"] == 3


def test_expected_read_shares_come_from_the_copy_numbers(tmp_path):
    ga_ = genome([(A, "+")] * 7, seed=8, background=1200)
    gb = genome([(B, "-")] * 2, seed=9, background=1200)
    r = run(tmp_path, {"a_GCF_000000011.1.fna": [ga_], "b_GCF_000000012.1.fna": [gb]},
            [("ATCC 11", "Strain a", "18", "GCF_000000011.1", "A"), ("ATCC 12", "Strain b", "1.8", "GCF_000000012.1", "A")], depth=10000)
    a, b = r["strains"]["Strain_a_ATCC11"], r["strains"]["Strain_b_ATCC12"]
    assert (a["copies"], b["copies"]) == (7, 2)
    wa, wb = 18 * 7, 1.8 * 2
    assert a["read_share"] == pytest.approx(wa / (wa + wb)) and b["read_share"] == pytest.approx(wb / (wa + wb))
    rows = list(csv.DictReader(open(tmp_path / "out" / "genome_amplicons.tsv"), delimiter="\t"))
    by = {x["strain"]: x for x in rows}
    assert float(by["Strain_b_ATCC12"]["expected_read_share_percent"]) == pytest.approx(100 * wb / (wa + wb), abs=1e-3)
    assert float(by["Strain_b_ATCC12"]["expected_reads_per_10000"]) == pytest.approx(10000 * wb / (wa + wb), abs=0.06)
    variants = list(csv.DictReader(open(tmp_path / "out" / "expected_variants.tsv"), delimiter="\t"))
    assert sum(float(v["expected_read_share_percent"]) for v in variants) == pytest.approx(100, abs=0.01)


def test_a_variant_gets_its_share_of_the_strains_reads(tmp_path):
    g = genome([(A, "+"), (A, "+"), (A, "+"), (A_VARIANT, "-")], seed=10, background=1200)
    r = run(tmp_path, {"a_GCF_000000021.1.fna": [g]}, [("ATCC 21", "Solo", "100", "GCF_000000021.1", "A")])
    v = {x["variant"]: float(x["expected_read_share_percent"]) for x in csv.DictReader(open(tmp_path / "out" / "expected_variants.tsv"), delimiter="\t")}
    assert v["Solo_ATCC21_v1"] == pytest.approx(75, abs=1e-3) and v["Solo_ATCC21_v2"] == pytest.approx(25, abs=1e-3)


def test_files_are_matched_to_the_table_by_accession_and_the_gaps_are_reported(tmp_path):
    g = genome([(A, "+")], seed=11, background=1200)
    r = run(tmp_path, {"anything_GCF_000000031.1.fna": [g], "stranger.fna": [g]},
            [("ATCC 31", "Known strain", "5", "GCF_000000031.1", "A"), ("ATCC 32", "Absent strain", "5", "GCF_000000032.1", "A")])
    text = " ".join(r["warnings"])
    assert "stranger.fna: no row of the table" in text and "Absent_strain_ATCC32" in text and "no genome file" in text
    assert "Known_strain_ATCC31" in r["strains"] and r["strains"]["stranger"]["share"] is None and r["strains"]["stranger"]["read_share"] is None


def test_warnings_for_draft_genomes_and_weak_primer_sites(tmp_path):
    g = genome([(A, "+", 2, 0)], seed=12, background=1200)
    r = run(tmp_path, {"d_GCF_000000041.1.fna": [g]}, [("ATCC 41", "Draftish", "1", "GCF_000000041.1", "A (draft)")])
    text = " ".join(r["warnings"])
    assert "a draft genome" in text and "up to 2 mismatches" in text


def test_gzipped_and_multi_contig_genomes(tmp_path):
    c1, c2 = genome([(A, "+")], seed=13, background=900), genome([(A, "-"), (B, "+")], seed=14, background=900)
    r = run(tmp_path, {"m_GCF_000000051.1.fna.gz": [c1, c2]}, [("ATCC 51", "Multi", "2", "GCF_000000051.1", "A")])
    s = r["strains"]["Multi_ATCC51"]
    assert s["contigs"] == 2 and s["copies"] == 3 and sorted(c for _, c, _ in s["variants"]) == [1, 2]


def test_the_output_is_the_form_mock_check_reads(tmp_path):
    g = genome([(A, "+")], seed=15, background=1200)
    r = run(tmp_path, {"a_GCF_000000061.1.fna": [g]}, [("ATCC BAA-61", "Genus species", "1", "GCF_000000061.1", "A")])
    (header,) = fasta_records(r["fasta"]).keys()
    assert header == "Genus_species_ATCCBAA-61_v1|Genus_species_ATCCBAA-61|1_copies"
    assert mr.strain_of_name(header) == "Genus_species_ATCCBAA-61"


def test_strains_with_the_same_amplicon_are_flagged_by_the_reference_builder(tmp_path):
    ga1, gb = genome([(A, "+")], seed=16, background=1200), genome([(A, "-")], seed=17, background=1200)
    r = run(tmp_path, {"a_GCF_000000071.1.fna": [ga1], "b_GCF_000000072.1.fna": [gb]},
            [("ATCC 71", "One", "1", "GCF_000000071.1", "A"), ("ATCC 72", "Two", "1", "GCF_000000072.1", "A")])
    assert any("identical to another strain's" in line for line in r["report"])


def test_the_command_line_writes_the_four_files(tmp_path):
    folder = tmp_path / "genomes"
    folder.mkdir()
    write_genome(folder / "a_GCF_000000081.1.fna", [genome([(A, "+")], seed=18, background=1200)])
    table = write_table(tmp_path / "t.tsv", [("ATCC 81", "Cli strain", "3", "GCF_000000081.1", "A")])
    assert ga.main([str(folder), "--table", str(table), "--primer-f", FWD, "--primer-r", REV, "--out", str(tmp_path / "o")]) == 0
    assert {p.name for p in (tmp_path / "o").iterdir()} == {"expected_amplicons.fasta", "genome_amplicons.tsv", "expected_variants.tsv", "EXPECTED_AMPLICONS.md"}
    assert "Cli_strain_ATCC81" in (tmp_path / "o" / "EXPECTED_AMPLICONS.md").read_text()


def test_an_empty_folder_is_refused(tmp_path):
    (tmp_path / "none").mkdir()
    with pytest.raises(SystemExit, match="No genome files"):
        ga.run(tmp_path / "none", None, FWD, REV, tmp_path / "o", log=lambda *_: None)


def test_the_script_parses_as_python_3_9():
    ast.parse(Path(ga.__file__).read_text(), feature_version=(3, 9))
