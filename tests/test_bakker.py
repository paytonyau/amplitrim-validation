"""validation/bakker_download.py and validation/bakker_primers.py: the run table of PRJNA377530 (the real one, tests/data), ENA's addresses, the download and
the count check with ENA replaced, and the primer census on synthetic reads."""
import gzip
import importlib.util
import io
import random
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
RUNINFO = HERE / "data" / "SraRunInfo_PRJNA377530.csv"


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE.parent / "validation" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


D, P = load("bakker_download"), load("bakker_primers")


def test_enas_address_follows_the_number_of_digits():
    assert D.ena_url("SRR123456", 1) == "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR123456/SRR123456_1.fastq.gz"
    assert D.ena_url("SRR5314315", 2) == "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR531/005/SRR5314315/SRR5314315_2.fastq.gz"          # Bakker: 7 digits
    assert D.ena_url("SRR13489922", 1) == "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR134/022/SRR13489922/SRR13489922_1.fastq.gz"      # Arfken: 8 digits
    assert D.ena_url("SRR123456789", 1).endswith("/SRR123/789/SRR123456789/SRR123456789_1.fastq.gz")


def test_the_real_run_table_gives_27_mock_runs_in_each_sequencing_run():
    runs = D.read_runinfo(RUNINFO)
    assert {b: len(v) for b, v in runs.items()} == {"pcr": 27, "polymerase": 27}
    assert all(lib.startswith("FMock") for _, lib, _, _ in runs["pcr"]) and all(lib.startswith("Mock") for _, lib, _, _ in runs["polymerase"])
    with_neg = D.read_runinfo(RUNINFO, with_negatives=True)
    assert {b: len(v) for b, v in with_neg.items()} == {"pcr": 31, "polymerase": 30}
    assert ("SRR5838505", "MockB_PreAmp_1", 44854, 22442318) in runs["polymerase"]


def test_a_foreign_or_single_end_table_is_refused(tmp_path):
    text = RUNINFO.read_text(encoding="utf-8-sig")
    (tmp_path / "other.csv").write_text(text.replace("PRJNA377530", "PRJNA1", 1))
    with pytest.raises(SystemExit, match="not PRJNA377530"):
        D.read_runinfo(tmp_path / "other.csv")
    (tmp_path / "single.csv").write_text(text.replace("PAIRED", "SINGLE", 1))
    with pytest.raises(SystemExit, match="not paired"):
        D.read_runinfo(tmp_path / "single.csv")


class Response(io.BytesIO):
    def __init__(self, data):
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data))}


def test_download_and_check_with_ena_replaced(tmp_path, capsys):
    table = tmp_path / "runinfo.csv"
    lines = RUNINFO.read_text(encoding="utf-8-sig").splitlines()
    keep = [lines[0]] + [l for l in lines[1:] if ",MockB_PreAmp_1," in l or ",FMockA1," in l]
    table.write_text("\n".join(keep) + "\n")
    counts = {"SRR5838505": 3}
    import csv
    for r in csv.DictReader(open(table)):
        counts[r["Run"]] = 3
    small = table.read_text().replace(",44854,", ",3,")
    for r in csv.DictReader(open(table)):
        small = small.replace(f",{r['spots']},", ",3,")
    table.write_text(small)
    calls = []

    def opener(url):
        calls.append(url)
        return Response(gzip.compress("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(3)).encode()))
    assert D.main(["--base", str(tmp_path), "--runinfo", str(table)], opener=opener, sleep=lambda s: None) == 0
    assert any("/SRR583/005/SRR5838505/SRR5838505_1.fastq.gz" in c for c in calls)
    assert any("/SRR531/007/SRR5314357/SRR5314357_1.fastq.gz" in c for c in calls)                  # FMockA1 (SRR5314357): the old rule would give 057
    assert (tmp_path / "23_mock_ITSmix_bakker_polymerase_PRJNA377530" / "raw" / "MockB_PreAmp_1_R2.fastq.gz").is_file()
    assert (tmp_path / "22_mock_ITS1_bakker_pcr_PRJNA377530" / "raw" / "FMockA1_R1.fastq.gz").is_file()
    assert "MockB_PreAmp_1" in (tmp_path / "bakker_runs.tsv").read_text()
    assert D.main(["--base", str(tmp_path), "--runinfo", str(table), "--check"]) == 0
    p = tmp_path / "22_mock_ITS1_bakker_pcr_PRJNA377530" / "raw" / "FMockA1_R2.fastq.gz"
    p.write_bytes(gzip.compress(b"@r\nA\n+\nI\n"))
    assert D.main(["--base", str(tmp_path), "--runinfo", str(table), "--check"]) == 2
    assert "FMockA1 mate 2: 1 reads, the run table says 3" in capsys.readouterr().out


# ------------------------------------------------------------------ the primer census
random.seed(4)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))


def change(s, *pos):
    s = list(s)
    for p in pos:
        s[p] = "A" if s[p] != "A" else "C"
    return "".join(s)


def test_a_primer_is_found_after_a_spacer_with_mismatches_and_degenerate_bases():
    assert P.best_primer("ACG" + change(P.PRIMERS["ITS1F"], 5) + rnd(50)) == ("ITS1F", 3, 1)
    assert P.best_primer("GATGAAGAACGCAGCGAAATGC" + rnd(40))[0] == "ITS3_KYO2"            # Bakker's ITS2 references start like this: Y=C, Y=C, R=G
    assert P.best_primer(P.PRIMERS["ITS4"] + rnd(40)) == ("ITS4", 0, 0)
    assert P.best_primer(rnd(60), max_mm=0) is None


def fq(path, seqs):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(seqs)))


def test_the_census_counts_primer_pairs_and_shows_unknown_starts(tmp_path):
    raw = tmp_path / "base" / "23_mock_ITSmix_bakker_polymerase_PRJNA377530" / "raw"
    its1 = [("TT" + P.PRIMERS["ITS1F"] + rnd(80), P.PRIMERS["ITS2"] + rnd(80)) for _ in range(6)]
    its2 = [("GATGAAGAACGCAGCGAAATGC" + rnd(80), P.PRIMERS["ITS4"] + rnd(80)) for _ in range(3)]
    unknown = [("GGGGGGGGGGGGGGGGGGGGGG" + rnd(60), P.PRIMERS["ITS4"] + rnd(60))]
    pairs = its1 + its2 + unknown
    fq(raw / "MockE_GoTaq_1_R1.fastq.gz", [a for a, _ in pairs])
    fq(raw / "MockE_GoTaq_1_R2.fastq.gz", [b for _, b in pairs])
    assert P.main(["--base", str(tmp_path / "base"), "--folders", "23_mock_ITSmix_bakker_polymerase_PRJNA377530"]) == 0
    rep = (tmp_path / "base" / "PRIMER_CENSUS.md").read_text()
    assert "| ITS1F | ITS2 | 6 | 60.0% |" in rep and "| ITS3_KYO2 | ITS4 | 3 | 30.0% |" in rep and "| - | ITS4 | 1 | 10.0% |" in rep
    assert "ITS1F 2: 6" in rep and "R1 GGGGGGGGGGGGGGGGGGGG (1)" in rep


# ------------------------------------------------------------------ the preparation (bakker_prep.py)
PREP = load("bakker_prep")
APPX = HERE / "data" / "AppendixS2_Bakker2018.txt"
F1, R1_ = P.PRIMERS["ITS1F"], P.PRIMERS["ITS2"]


def test_the_real_appendix_gives_amplicons_for_all_19_isolates_with_each_datasets_primers(tmp_path):
    rep = PREP.references(PREP.appendix(APPX), PREP.AMPLICONS, tmp_path)
    assert PREP.AMPLICONS == {"bakker_pcr_its1": ("ITS1", "BITS", "B58S3"), "bakker_pol_its1": ("ITS1", "ITS1F", "ITS2"), "bakker_pol_its2": ("ITS2", "ITS3_KYO2", "ITS4_KYO3")}
    assert {k: len(r["members"]) for k, r in rep.items()} == {"bakker_pcr_its1": 18, "bakker_pol_its1": 18, "bakker_pol_its2": 19}    # F. oxysporum and F. verticillioides as one on ITS1
    assert "Fusarium.oxysporum+verticillioides" in rep["bakker_pol_its1"]["members"] and "Fusarium.verticillioides" in rep["bakker_pol_its2"]["members"]
    assert [rep[k]["written"] for k in ("bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2")] == [47, 47, 23]
    assert rep["bakker_pol_its2"]["missing"] == ["Rhizophagus.irregularis.A1.e"] and rep["bakker_pcr_its1"]["missing"] == []
    lengths = {k: (min(r["lengths"]), max(r["lengths"])) for k, r in rep.items()}
    assert lengths == {"bakker_pcr_its1": (105, 366), "bakker_pol_its1": (177, 438), "bakker_pol_its2": (296, 514)}
    for label in PREP.AMPLICONS:
        body = "".join(l for l in (tmp_path / f"references_{label}.fasta").read_text().splitlines() if not l.startswith(">"))
        assert set(body) <= set("ACGT")
    first = (tmp_path / "references_bakker_pcr_its1.fasta").read_text().splitlines()[1]
    assert not first.startswith("ACCTGCGGA") and (tmp_path / "references_bakker_pol_its1.fasta").read_text().startswith(">Alternaria.alternata|")


def test_a_pair_is_assigned_trimmed_and_swapped_when_its_forward_primer_is_on_r2():
    body1, body2 = rnd(60), rnd(60)
    r1 = ("@a", "AC" + F1 + body1, "+", "I" * (2 + len(F1) + 60))
    r2 = ("@a", R1_ + body2, "+", "J" * (len(R1_) + 60))
    x, y, swapped = PREP.assign(r1, r2, F1, R1_, 2, 10)
    assert (x[1], y[1], swapped) == (body1, body2, False) and len(x[3]) == 60
    x, y, swapped = PREP.assign(r2, r1, F1, R1_, 2, 10)                                  # the mates the other way round
    assert (x[1], y[1], swapped) == (body1, body2, True)
    assert PREP.assign(("@b", rnd(80), "+", "I" * 80), r2, F1, R1_, 2, 10) is None


def test_prepare_splits_a_mixed_run_by_primer_pair(tmp_path):
    raw = tmp_path / "raw"
    its1 = [("T" + F1 + rnd(70), R1_ + rnd(70))] * 4
    its2 = [("GATGAAGAACGCAGCGAAATGC" + rnd(70), P.PRIMERS["ITS4"] + rnd(70))] * 3
    junk = [(rnd(90), rnd(90))]
    fq(raw / "MockE_GoTaq_1_R1.fastq.gz", [a for a, _ in its1 + its2 + junk])
    fq(raw / "MockE_GoTaq_1_R2.fastq.gz", [b for _, b in its1 + its2 + junk])
    rows1 = PREP.prepare(raw, tmp_path / "its1", F1, R1_)
    rows2 = PREP.prepare(raw, tmp_path / "its2", P.PRIMERS["ITS3_KYO2"], P.PRIMERS["ITS4"])
    assert rows1 == [("MockE_GoTaq_1", 8, 4, 0)] and rows2 == [("MockE_GoTaq_1", 8, 3, 0)]
    assert len(list(PREP.fastq(tmp_path / "its2" / "MockE_GoTaq_1_R2.fastq.gz"))) == 3


def test_ambiguity_codes_are_expanded_into_the_variants_they_allow():
    assert sorted(PREP.expand("ARY")) == ["AAC", "AAT", "AGC", "AGT"] and PREP.expand("ACGT") == ["ACGT"]
    with pytest.raises(ValueError, match="more than 64"):
        PREP.expand("N" * 4)
