"""validation/ncbi_references.py: in-silico PCR, the choice of assembly, and the run, with NCBI replaced by fake answers (no network)."""
import importlib.util
import io
import json
import random
import zipfile
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


def load():
    spec = importlib.util.spec_from_file_location("ncbi_references", HERE.parent / "validation" / "ncbi_references.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


N = load()
F1, R1 = N.PRIMERS["ITS1"]
F2, R2 = N.PRIMERS["ITS2"]
random.seed(7)
rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
SPACER1, SPACER2 = rnd(180), rnd(230)
# an rDNA unit as it reads on the genome: ITS1F, the ITS1 spacer, the 5.8S site (where ITS2 binds as the reverse primer and ITS3 as the forward one), the ITS2 spacer, the ITS4 site
assert N.rc(R1) == F2
UNIT = F1 + SPACER1 + F2 + SPACER2 + N.rc(R2)
GENOME = rnd(3000) + UNIT + rnd(2000)


def test_pcr_finds_both_amplicons_without_their_primers():
    p1 = N.pcr(GENOME, F1, R1)
    p2 = N.pcr(GENOME, F2, R2)
    assert [x[3] for x in p1] == [SPACER1] and [x[3] for x in p2] == [SPACER2]
    assert p1[0][0] == "+"


def test_pcr_reads_the_other_strand_too():
    p = N.pcr(N.rc(GENOME), F1, R1)
    assert len(p) == 1 and p[0][0] == "-" and p[0][3] == SPACER1


def change(s, *positions):
    s = list(s)
    for p in positions:
        s[p] = "A" if s[p] != "A" else "C"
    return "".join(s)


def test_a_primer_site_with_three_mismatches_is_accepted_and_four_are_not():
    three, four = change(F1, 2, 4, 6), change(F1, 2, 4, 6, 8)
    assert N.mismatches(three, F1) == 3 and N.mismatches(four, F1) == 4
    assert N.pcr(rnd(500) + three + SPACER1 + F2 + rnd(300), F1, R1)
    assert not N.pcr(rnd(500) + four + SPACER1 + F2 + rnd(300), F1, R1)


def test_a_mismatch_at_the_3_end_of_the_primer_does_not_hide_the_site():
    """ITS1F has a 3'-terminal mismatch against Saccharomycetes, which Arfken et al. 2023 still amplified."""
    end = change(F1, len(F1) - 1)
    both_ends = change(F1, 0, len(F1) - 1)
    assert N.pcr(rnd(500) + end + SPACER1 + F2 + rnd(300), F1, R1)[0][3] == SPACER1
    assert N.pcr(rnd(500) + both_ends + SPACER1 + F2 + rnd(300), F1, R1)[0][3] == SPACER1
    assert N.pcr(rnd(500) + F1 + SPACER1 + change(F2, 0) + rnd(300), F1, R1)[0][3] == SPACER1      # the reverse primer's 3' end


def test_the_length_limits_hold():
    assert not N.pcr(GENOME, F1, R1, max_len=100)                       # the ITS1 spacer is 180 nt
    assert not N.pcr(GENOME, F1, R1, min_len=200)


def test_the_choice_prefers_the_strain_then_the_reference_genome():
    reports = [{"accession": "GCA_1", "strain": "X12", "assembly_name": "", "level": "Contig", "category": "", "n50": 10},
               {"accession": "GCF_2", "strain": "Af293", "assembly_name": "ASM265v1", "level": "Chromosome", "category": "reference genome", "n50": 100},
               {"accession": "GCA_3", "strain": "ATCC MYA-4609", "assembly_name": "", "level": "Scaffold", "category": "", "n50": 50}]
    pick, how = N.choose(reports, ["MYA-4609", "Af293"])
    assert how == "strain" and pick["accession"] == "GCF_2"            # both match the strain; the reference genome wins
    pick, how = N.choose(reports[:1] + [dict(reports[1], strain="other")], ["MYA-4609"])
    assert how == "species" and pick["accession"] == "GCF_2"
    assert N.choose([], ["x"]) == (None, "none")
    assert N.choose([{"accession": "GCA_9", "strain": "a-b c", "assembly_name": "", "level": "", "category": "", "n50": 0}], ["AB C"])[1] == "strain"   # spaces and hyphens do not count


def zipped(records):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("ncbi_dataset/data/GCF_2/GCF_2_genomic.fna", "".join(f">{n} some description\n" + "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n" for n, s in records))
        z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}")
    return buf.getvalue()


class FakeNCBI:
    """Answers the URLs the script asks for. assemblies: taxon name -> list of reports; genomes: accession -> list of (name, sequence); its: taxon name -> list of (name, sequence)."""

    def __init__(self, assemblies, genomes, its=None, its_strain=None):
        self.assemblies, self.genomes, self.its, self.its_strain, self.calls, self.pending = assemblies, genomes, its or {}, its_strain or {}, [], []

    def __call__(self, url):
        self.calls.append(url)
        if "/dataset_report" in url:
            name = N.urllib.parse.unquote(url.split("/genome/taxon/")[1].split("/dataset_report")[0])
            reports = self.assemblies.get(name, [])
            if "reference_only=true" in url:
                reports = [r for r in reports if r.get("assembly_info", {}).get("refseq_category")]
            return json.dumps({"reports": reports}).encode()
        if "/download?" in url:
            acc = url.split("/genome/accession/")[1].split("/download")[0]
            return zipped(self.genomes[acc])
        if "esearch.fcgi" in url:
            quoted = lambda x: N.urllib.parse.quote(f'"{x}"')
            strain = [k for k in self.its_strain if quoted(k[0]) in url and quoted(k[1]) in url]
            if strain:
                self.pending = self.its_strain[strain[0]]
            elif "%5BAll%20Fields%5D" in url:                                   # a strain was asked for, and there is no record of it
                self.pending = []
            else:
                self.pending = next((v for n, v in self.its.items() if quoted(n) in url), [])
            return json.dumps({"esearchresult": {"idlist": [f"id{i}" for i in range(len(self.pending))]}}).encode()
        if "efetch.fcgi" in url:
            return "".join(f">{n}\n{s}\n" for n, s in self.pending).encode()
        raise AssertionError(f"unexpected URL {url}")


def report(acc, strain, level="Chromosome", category="reference genome", n50=100):
    return {"accession": acc, "organism": {"organism_name": "x", "infraspecific_names": {"strain": strain}},
            "assembly_info": {"assembly_level": level, "refseq_category": category, "assembly_name": ""}, "assembly_stats": {"contig_n50": n50}}


def test_assemblies_are_read_from_both_queries_once_each():
    fake = FakeNCBI({"Aspergillus fumigatus": [report("GCF_2", "Af293"), report("GCA_7", "other", "Contig", "")]}, {})
    found = N.assemblies("Aspergillus fumigatus", get=fake)
    assert [f["accession"] for f in found] == ["GCF_2", "GCA_7"] and found[0]["strain"] == "Af293"
    assert sum("dataset_report" in c for c in fake.calls) == 2


def test_the_genome_is_downloaded_once_and_read_from_the_zip(tmp_path):
    fake = FakeNCBI({}, {"GCF_2": [("chr1", GENOME), ("chr2", N.rc(UNIT))]})
    z = N.download("GCF_2", tmp_path, get=fake)
    assert N.download("GCF_2", tmp_path, get=fake) == z and sum("/download?" in c for c in fake.calls) == 1
    names = [n for n, _ in N.zip_records(z)]
    assert names == ["chr1", "chr2"]
    found = N.products(N.zip_records(z), 2, 50, 1500)
    assert found["ITS1"] == {SPACER1: 2} and found["ITS2"] == {SPACER2: 2}          # one copy on each strand


def seqs_tsv(path, isolates):
    lines = ["sample\tk\tsequence\tabundance"]
    for iso, seq in isolates.items():
        lines += [f"{iso}-ITS1a\t0\t{seq}\t90", f"{iso}-ITS1a\t0\t{seq[:-3]}\t5", f"{iso}-ITS1b\t2\tAAAA\t99"]
    path.write_text("\n".join(lines) + "\n")


def full_world():
    return FakeNCBI({"Aspergillus fumigatus": [report("GCF_2", "Af293")], "Kazachstania slooffiae": [report("GCA_5", "piglet isolate", "Scaffold", "")]},
                    {"GCF_2": [("chr1", GENOME)], "GCA_5": [("s1", rnd(4000))]},                      # the K genome holds no rDNA
                    its={"Kazachstania slooffiae": [("KX1", rnd(100) + UNIT + rnd(100))]})


def test_list_shows_the_choice_and_downloads_nothing(tmp_path, capsys):
    fake = full_world()
    assert N.main(["--out", str(tmp_path), "--list", "--only", "Af,K"], get=fake) == 0
    out = capsys.readouterr().out
    assert "chosen: GCF_2 (Af293, Chromosome, reference genome) [strain]" in out and "GCA_5" in out and "Nothing was downloaded" in out
    assert not any("/download?" in c for c in fake.calls)


def test_the_run_writes_the_references_falls_back_to_its_records_and_compares_with_the_isolates(tmp_path, capsys):
    fake = full_world()
    run25 = tmp_path / "run25"
    run25.mkdir()
    other = lambda b: "A" if b != "A" else "C"
    seqs_tsv(run25 / "arfken_its1_seqs.tsv", {"K": SPACER1[:-2] + other(SPACER1[-2]) + other(SPACER1[-1])})     # two bases surely different
    assert N.main(["--out", str(tmp_path / "out"), "--only", "Af,K", "--compare-isolates", str(run25)], get=fake) == 0
    f1 = (tmp_path / "out" / "references_ITS1.fasta").read_text()
    assert f1.count(">") == 2 and f">Af|Aspergillus_fumigatus|GCF_2|strain_genome|copy1x1|{len(SPACER1)}nt\n{SPACER1}\n" in f1
    rep = (tmp_path / "out" / "REFERENCES_REPORT.md").read_text()
    assert "| K | piglet | Kazachstania slooffiae | GCA_5 | piglet isolate | species (record) | species (record) |" in rep
    assert "| Af | ATCC | Aspergillus fumigatus | GCF_2 | Af293 | strain (genome) | strain (genome) |" in rep
    assert "ITS1 from NCBI ITS records of the species" in rep and "Members without any product: none." in rep
    assert "ATCC members with a reference of the very strain (from its genome or its own ITS records): ITS1 1 of 1 (Af); ITS2 1 of 1 (Af)." in rep
    assert "2 (isolate 180 nt)" in rep                                     # the isolate's sequence differs from the reference at its last two bases


def test_a_member_without_any_product_is_named_and_the_exit_code_says_so(tmp_path):
    fake = FakeNCBI({"Fusarium keratoplasticum": [report("GCA_8", "ATCC 36031", "Contig", "")]}, {"GCA_8": [("c", rnd(3000))]})
    assert N.main(["--out", str(tmp_path), "--only", "Fk"], get=fake) == 2
    assert "Members without any product: Fk." in (tmp_path / "REFERENCES_REPORT.md").read_text()


def test_use_overrides_the_choice_and_unknown_members_are_refused(tmp_path, capsys):
    fake = FakeNCBI({"Aspergillus fumigatus": [report("GCF_2", "Af293"), report("GCA_3", "ATCC MYA-4609", "Scaffold", "")]}, {"GCA_3": [("c", GENOME)]})
    assert N.main(["--out", str(tmp_path), "--only", "Af", "--use", "Af=GCA_3"], get=fake) == 0
    assert "chosen: GCA_3 (ATCC MYA-4609, Scaffold, -) [strain]" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="unknown member"):
        N.main(["--out", str(tmp_path), "--only", "Zz"], get=fake)
    with pytest.raises(SystemExit, match="ID=ACCESSION"):
        N.main(["--out", str(tmp_path), "--use", "Af"], get=fake)


def test_the_members_are_the_ten_atcc_strains_and_the_five_piglet_isolates():
    assert [t[0] for t in N.TARGETS if t[1] == "ATCC"] == ["Af", "Cn", "Ti", "Pc", "Fk", "Ca", "Cg", "Mg", "Sc", "Cd"]
    assert [t[0] for t in N.TARGETS if t[1] == "piglet"] == ["C", "K", "M", "P", "T"]
    assert all(t[3] for t in N.TARGETS if t[1] == "ATCC")                   # every ATCC strain has its strain numbers to look for


def test_an_atcc_strain_without_rdna_in_its_genome_gets_its_own_records_first_and_only_those_count_as_the_strain(tmp_path):
    no_rdna = [("s1", rnd(4000))]
    fake = FakeNCBI({"Candida albicans": [report("GCA_C", "ATCC 10231", "Scaffold", "")]}, {"GCA_C": no_rdna},
                    its={"Candida albicans": [("SPECIES1", rnd(80) + UNIT + rnd(80))]},
                    its_strain={("Candida albicans", "10231"): [("STRAIN1", rnd(80) + UNIT + rnd(80))]})
    assert N.main(["--out", str(tmp_path / "a"), "--only", "Ca"], get=fake) == 0
    rep = (tmp_path / "a" / "REFERENCES_REPORT.md").read_text()
    assert "| Ca | ATCC | Candida albicans | GCA_C | ATCC 10231 | strain (record) | strain (record) |" in rep and "records of 10231" in rep
    assert "ITS1 1 of 1 (Ca); ITS2 1 of 1 (Ca)" in rep
    only_species = FakeNCBI({"Candida albicans": [report("GCA_C", "ATCC 10231", "Scaffold", "")]}, {"GCA_C": no_rdna},
                            its={"Candida albicans": [("SPECIES1", rnd(80) + UNIT + rnd(80))]})
    assert N.main(["--out", str(tmp_path / "b"), "--only", "Ca"], get=only_species) == 0
    rep = (tmp_path / "b" / "REFERENCES_REPORT.md").read_text()
    assert "| strain (record) |" not in rep and "species (record) | species (record) |" in rep
    assert "ITS1 0 of 1 (none); ITS2 0 of 1 (none)" in rep                              # a record of the species does not count as the strain


def test_every_site_with_up_to_three_mismatches_is_found_wherever_they_are():
    random.seed(11)
    for _ in range(300):
        pos = random.sample(range(len(F1)), random.randint(0, 3))
        assert N.sites(rnd(200) + change(F1, *pos) + rnd(200), F1, 3) == [200], pos
