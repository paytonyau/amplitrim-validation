"""describe_reads.py: what it reports about a folder of samples."""
import gzip
import importlib.util
import random
from pathlib import Path

VALIDATION = Path(__file__).parent.parent / "validation"
spec = importlib.util.spec_from_file_location("describe_reads", VALIDATION / "describe_reads.py")
dr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dr)

F338, R806 = "ACTCCTACGGGAGGCAGCAG", "GGACTACCAGGGTATCTAAT"       # 806R with its degenerate bases resolved


def write_sample(folder, name, n=300, barcode=7, q_alphabet="I?5+&", reverse_share=0.0, seed=1):
    rng = random.Random(seed)
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    with gzip.open(folder / f"{name}_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(folder / f"{name}_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(n):
            fwd, rev = junk(barcode) + F338 + junk(60), R806 + junk(70)
            if rng.random() < reverse_share:
                fwd, rev = R806 + junk(70), junk(barcode) + F338 + junk(60)
            for fh, s in ((o1, fwd), (o2, rev)):
                fh.write(f"@{name}_{i}\n{s}\n+\n{''.join(rng.choice(q_alphabet) for _ in s)}\n")


def test_barcodes_primers_quality_and_orientation_are_reported_per_sample(tmp_path):
    write_sample(tmp_path, "good", barcode=7)
    write_sample(tmp_path, "flat", q_alphabet="?$", reverse_share=0.4, seed=2)
    rows = {r["sample"]: r for r in dr.describe(tmp_path, 300)}
    good, flat = rows["good"], rows["flat"]
    assert (good["r1_primer"], good["r1_offset"], good["r2_primer"], good["r2_offset"]) == ("338F", 7, "806R", 0)
    assert good["r1_share"] == 1.0 and good["quality_values"] == 5 and good["reversed_share"] == 0.0 and good["reads"] == 300
    assert flat["quality_values"] == 2 and 0.25 < flat["reversed_share"] < 0.55


def test_the_summary_names_each_problem_in_plain_words(tmp_path):
    write_sample(tmp_path, "good")
    write_sample(tmp_path, "flat", q_alphabet="?$", reverse_share=0.4, seed=2)
    text = " ".join(dr.summarise(dr.describe(tmp_path, 300)))
    assert "1 of 2 samples have 2 or fewer distinct quality values (flat)" in text and "DADA2 cannot learn" in text
    assert "reverse orientation" in text and "--mixed-orientation" in text
    assert "Reads per sample: min 300, median 300, max 300" in text


def test_reads_that_all_start_at_the_same_place_are_reported_as_conserved(tmp_path):
    rng = random.Random(6)
    with gzip.open(tmp_path / "s_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(tmp_path / "s_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(200):
            body = "".join(rng.choice("ACGT") for _ in range(100))
            for fh in (o1, o2):
                fh.write(f"@r{i}\nTACGGAGGGTGCAAGCGTTA{body}\n+\n{'I' * 120}\n")
    row = dr.describe(tmp_path, 200)[0]
    assert row["start_motif"] == "TACGGAGGGTGCAAGCGTTA" and row["start_motif_share"] == 1.0 and row["start_conservation"] == 1.0
    text = " ".join(dr.summarise([row]))
    assert "begins with TACGGAGGGTGCAAGCGTTA in 100% of reads" in text and "run AmpliTrim directly" in text


def test_a_named_primer_is_found_with_mismatches_at_any_offset_and_read_through_is_detected():
    primer = "ACTCCTACGGGAGGCAGCA"
    rng = random.Random(9)
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    mutated = primer[:5] + ("A" if primer[5] != "A" else "C") + primer[6:]                     # one mismatch
    reads = [(junk(7) + primer + junk(80), "I" * 106)] * 40 + [(junk(25) + mutated + junk(60), "I" * 104)] * 40 \
        + [(junk(100), "I" * 100)] * 20 + [(junk(50) + dr.reverse_complement(primer) + junk(10), "I" * 79)] * 10
    found = dr.search_primer(reads, primer, max_mismatches=2)
    assert found["share"] == 80 / 110 and found["offset"] in (7, 25) and found["rc_share"] == 10 / 110
    assert dr.search_primer(reads, primer, max_mismatches=0)["share"] == 40 / 110      # exact matches only
    assert dr.mismatches("ACGT", "ANGT") == 0 and dr.mismatches("ACGT", "AGGT") == 1       # degenerate letters match


def test_a_degenerate_primer_is_read_off_the_start_of_the_reads_and_then_found_by_search():
    rng = random.Random(12)
    alphabet = {"R": "AG", "Y": "CT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
    f563, r802 = "AYTGGGYDTAAAGNG", "TACNVGGGTATCTAATCC"
    realise = lambda p: "".join(rng.choice(alphabet.get(c, c)) for c in p)
    tail = lambda: "".join(rng.choice("ACGT") for _ in range(120))
    r1 = [(realise(f563) + tail(), "") for _ in range(800)]
    r2 = [(realise(r802) + tail(), "") for _ in range(800)]
    assert dr.start_consensus(r1, 15) == f563
    assert dr.start_consensus(r2, 18) == r802
    assert dr.search_primer(r1, f563, 0)["share"] == 1.0 and dr.search_primer(r1, f563, 0)["offset"] == 0
    assert dr.find_primer(r1[0][0])[0] == "563F" and dr.find_primer(r2[0][0])[0] == "802R"       # now in the built-in list


def test_leftover_read_through_clusters_where_a_match_by_chance_does_not():
    rng = random.Random(21)
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    primer = "GGATTAGATACCCTGTAG"
    rc = dr.reverse_complement(primer)
    leftover = [(junk(rng.randint(185, 205)) + rc + junk(10), "") for _ in range(80)]      # amplicon ends at about the same place
    by_chance = [(junk(i % 150) + rc + junk(200 - i % 150), "") for i in range(80)]         # anywhere along the read
    a, b = dr.search_primer(leftover, primer, 0), dr.search_primer(by_chance, primer, 0)
    assert a["rc_share"] == 1.0 and a["rc_clustered"] > 0.9 and 185 <= a["rc_median_offset"] <= 205
    assert b["rc_share"] == 1.0 and b["rc_clustered"] < 0.4


def test_the_position_table_separates_conserved_starts_from_unrelated_sequence():
    rng = random.Random(3)
    conserved = [("TGGGGAATATTG" + "".join(rng.choice("ACGT") for _ in range(30)), "") for _ in range(200)]
    random_reads = [("".join(rng.choice("ACGT") for _ in range(42)), "") for _ in range(200)]
    t = dr.position_table(conserved, 12)
    assert [r["top"] for r in t] == list("TGGGGAATATTG") and all(r[r["top"]] == 1.0 for r in t)
    assert all(max(r[b] for b in "ACGT") < 0.5 for r in dr.position_table(random_reads, 12))


def test_reads_without_a_known_primer_and_an_empty_folder(tmp_path):
    rng = random.Random(4)
    with gzip.open(tmp_path / "s_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(tmp_path / "s_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(100):
            for fh in (o1, o2):
                s = "".join(rng.choice("ACGT") for _ in range(120))
                fh.write(f"@r{i}\n{s}\n+\n{'I' * 120}\n")
    assert any("no known primer" in n for n in dr.summarise(dr.describe(tmp_path, 100)))
    assert dr.summarise([]) == ["No paired FASTQ files were found."]


def test_length_bins_show_a_short_tail():
    bins = dict(dr.length_bins([207] * 90 + [195] * 7 + [120] * 3))
    assert bins["206-209"] == 0.9 and bins["190-199"] == 0.07 and bins["100-149"] == 0.03 and bins[">=230"] == 0.0


def test_a_primer_site_inside_the_gene_is_not_reported_as_a_primer():
    rng = random.Random(3)
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    site_563f = "ATTGGGTTTAAAGGG"
    inside_amplicon = junk(29) + site_563f + junk(200)                  # a read that starts just after 515F: 563F is 29 bases in
    assert dr.find_primer(inside_amplicon) is None
    with_barcode = junk(6) + "GTGCCAGCAGCCGCGGTAA" + junk(100)
    assert dr.find_primer(with_barcode) == ("515F", 6)


def test_a_degenerate_primer_variant_is_still_found():
    """27F-YM has Y (C or T) at position 11; a list entry with only C found 57% of reads and missed the T variants."""
    for variant in ("AGAGTTTGATCATGGCTCAG", "AGAGTTTGATTATGGCTCAG", "AGAGTTTGATCCTGGCTCAG", "AGAGTTTGATTCTGGCTCAG"):
        assert dr.find_primer(variant + "TTGAACGCTGGCGGCAGGCC") == ("27F", 0)


def test_a_few_chance_matches_are_not_reported_as_a_primer_or_as_reversed_pairs(tmp_path):
    import gzip
    rng = random.Random(9)
    junk = lambda k: "".join(rng.choice("ACGT") for _ in range(k))
    f27, r534 = "AGAGTTTGATCATGGCTCAG", "ATTACCGCGGCTGCTGG"
    for mate in ("R1", "R2"):
        with gzip.open(tmp_path / f"s_{mate}.fastq.gz", "wt", newline="\n") as fh:
            for i in range(100):
                if i < 2:                                                    # two pairs that look reversed: 534R on R1, 27F on R2
                    seq = (r534 if mate == "R1" else f27) + junk(80)
                else:
                    seq = junk(100)
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    row = dr.describe(tmp_path)[0]
    assert row["r1_primer"] == "" and row["r2_primer"] == "" and row["reversed_share"] == ""


def test_the_length_messages_are_quantitative_and_call_a_constant_offset_normal():
    assert dr.pair_length_summary([280] * 100, [283] * 100) == "R2 is 3 nt longer than R1 in 100% of pairs; this is normal when the two primers differ in length."
    assert dr.pair_length_summary([250] * 50, [250] * 50) == "R1 and R2 have the same length in 100% of pairs."
    mixed = dr.pair_length_summary([280] * 60 + [270] * 40, [283] * 100)
    assert "R2 is 3 nt longer than R1 in 60% of pairs" in mixed and "the rest differ by other amounts" in mixed


def test_read_lengths_count_as_varying_only_when_more_than_five_percent_are_off_the_common_length(tmp_path):
    import gzip
    rng = random.Random(5)
    for name, odd in (("few", 2), ("many", 20)):
        for mate in ("R1", "R2"):
            with gzip.open(tmp_path / f"{name}_{mate}.fastq.gz", "wt", newline="\n") as fh:
                for i in range(100):
                    n = 120 if i >= odd else 100
                    fh.write(f"@r{i}\n{''.join(rng.choice('ACGT') for _ in range(n))}\n+\n{'I' * n}\n")
    rows = {r["sample"]: r for r in dr.describe(tmp_path)}
    assert rows["few"]["lengths_differ"] is False and rows["many"]["lengths_differ"] is True and rows["many"]["off_modal_share"] == 0.2


def test_the_kyo_its_primers_are_recognised_in_all_their_degenerate_variants():
    for variant in ("CTAGGTCATTTAGAGGAACTAA", "CTCGGTCATTTAGAGGAAGTAA", "CTTGGTCATTTAGAGGAACTAA"):          # H = A/C/T, S = C/G
        assert dr.find_primer(variant + "GTAGAAATGTCTAACAGG") == ("ITS1F-KYO1", 0)
    for variant in ("TTCGCTGCGTTCTTCATC", "TTTACTACGTTCTTCATC", "TTCACTGCGTTCTTCATC"):                        # Y = C/T, R = A/G
        assert dr.find_primer(variant + "GATGCCAGAACCAAGAGATC") == ("ITS2-KYO2", 0)


def test_a_primer_whose_first_bases_are_missing_is_still_recognised_and_gets_a_negative_offset():
    gene = "GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGG"
    p27 = "AGAGTTTGATCATGGCTCAG"
    assert dr.find_primer(p27 + gene) == ("27F", 0)
    assert dr.find_primer(p27[1:] + gene) == ("27F", -1)                  # the first base is missing
    assert dr.find_primer(p27[2:] + gene) == ("27F", -2)
    assert dr.find_primer("TT" + p27 + gene) == ("27F", 2)                # a spacer in front still works
    assert dr.find_primer(gene) is None and dr.find_primer("ACGT" * 15) is None
    assert dr.find_primer("T" * 8 + p27[1:] + gene) is None or dr.find_primer("T" * 8 + p27[1:] + gene)[1] >= 0     # a late partial match is not taken for a shifted one
    reads = [(p27[1:] + gene, "I" * 60)] * 20 + [("ACGT" * 15, "I" * 60)]
    name, share, offset = dr.top_primer(reads)
    assert (name, offset) == ("27F", -1) and share > 0.9


def test_a_plain_primer_is_not_called_by_its_degenerate_superset_but_a_read_that_only_the_degenerate_one_fits_is():
    tail = "GTAGAAATGTCTAACAGGAAC"
    assert dr.find_primer("CTTGGTCATTTAGAGGAAGTAA" + tail) == ("ITS1F", 0)               # fits both; the more specific one is named
    assert dr.find_primer("CTAGGTCATTTAGAGGAACTAA" + tail) == ("ITS1F-KYO1", 0)          # A at the H position: only KYO1 fits
    assert dr.find_primer("CTTGGTCATTTAGAGGAACTAA" + tail) == ("ITS1F-KYO1", 0)          # differs from ITS1F only beyond the 14-base seed
    assert dr.find_primer("GCTGCGTTCTTCATCGATGC" + tail) == ("ITS2", 0) and dr.find_primer("TTCGCTGCGTTCTTCATC" + tail) == ("ITS2-KYO2", 0)
    assert dr.find_primer("AGAGTTTGATCATGGCTCAG" + tail) == ("27F", 0) and dr.find_primer("GAGTTTGATCATGGCTCAG" + tail) == ("27F", -1)


def test_its1_fl2_is_recognised_at_the_start_of_the_reads_it_begins():
    tail = "TTACAGAGTTGCAAACTCCC"
    for variant in ("GAACCTGCGGAGGGATCA", "GAACCAGCGGAAGGATCA", "GAACCTGCGGAAGGATCA"):         # W = A/T, R = A/G
        assert dr.find_primer(variant + tail) == ("ITS1-Fl2", 0)
    assert dr.find_primer("TCCGTAGGTGAACCTGCGG" + tail)[0] == "ITS1"                         # the classic ITS1 primer is not mistaken for it
    assert dr.find_primer("CTTGGTCATTTAGAGGAAGTAA" + tail) == ("ITS1F", 0) and dr.find_primer("GCTGCGTTCTTCATCGATGC" + tail) == ("ITS2", 0)
