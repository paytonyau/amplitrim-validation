"""The test results of the validation archive, one record per dataset, extracted from docs/FINDINGS.md and docs/validation.md.
Nothing here is typed in by hand except the names of the datasets and the short notes; every number is read from the tables."""
import pathlib
import re

# id, display name, accessions, aliases (lower case, as the tables write them)
DATASETS = [
    ("mock-a", "Mock community A", ["PRJEB67797"], {"mock community a", "mock a"}),
    ("mock-b", "Mock community B", ["PRJEB67827"], {"mock community b", "mock b"}),
    ("interlab-v1v2", "Inter-laboratory mock, V1-V2", ["PRJEB34919"], {"inter-laboratory v1-v2", "inter-laboratory mock, v1-v2"}),
    ("interlab-v4v6", "Inter-laboratory mock, V4-V6", ["PRJEB34919"], {"inter-laboratory v4-v6", "inter-laboratory mock, v4-v6"}),
    ("d6300-v3v4", "ZymoBIOMICS D6300, V3-V4", ["PRJNA1261676"], {"d6300 v3-v4", "zymobiomics d6300, v3-v4"}),
    ("d6300-v4v5", "ZymoBIOMICS D6300, V4-V5", ["SRR16961581", "SRR16961582", "SRR16961583"], {"d6300 v4-v5", "zymobiomics d6300, v4-v5"}),
    ("v1v2", "Real samples, V1-V2", ["PRJEB42528"], {"v1-v2", "prjeb42528"}),
    ("soil-v4", "Soil, V4", ["PRJNA574030"], {"soil v4", "prjna574030"}),
    ("its", "Spring wheat, ITS", ["PRJEB58189"], {"its", "prjeb58189"}),
    ("wastewater", "Wastewater, V3-V4", ["PRJEB75155"], {"wastewater v3-v4", "wastewater", "prjeb75155"}),
    ("spike-in", "Mouse gut with spike-ins, V3-V4", ["PRJNA1243306"], {"spike-in v3-v4", "spike-in", "prjna1243306"}),
    ("swabs", "Meat-processing swabs, V3-V4", ["PRJNA1261676"], {"swabs v3-v4", "swabs", "swabs of meat-processing surfaces"}),
]
IDS = [d[0] for d in DATASETS]
ALIAS = {a: d[0] for d in DATASETS for a in d[3]}

# short notes: each restates something the findings say, with the numbers as the findings give them (a test checks that every number is in the findings)
NOTES = {
    "mock-a": "We used this dataset while designing AmpliTrim, so it is not independent evidence. It has six high-input samples. At every maxEE setting we tried, the same variants and strains were recovered, so a looser maxEE cost no accuracy here. Skipping truncation altogether merged 23.3 points fewer pairs.",
    "mock-b": "The poorest reads in the set: only 20% of the pairs pass DADA2's default filter. A looser maxEE gave a higher share of reads that match an expected sequence exactly (68.7% at maxEE 2, 73.6% to 75.2% for AmpliTrim's looser settings, 75.2% at maxEE 4,8 which is the released rule's setting, and 78.3% for the other tool's maxEE 6,9, which also found one more variant: 7 against 6). This is the one dataset where AmpliTrim's cut differed from the other tool's by more than the 2-point margin, by +2.6 points. In the first round of tests, which used a stricter margin of 1.0 point, it was also beyond the margin in AmpliTrim's favour.",
    "interlab-v1v2": "Six samples with no known truth (15 bacteria; only digital PCR abundances of families). The reads are clean: 96% of the pairs pass DADA2's default filter. When the other tool was given the typical amplicon length instead of the longest, it merged 54.4% of the pairs against 93.6%.",
    "interlab-v4v6": "A real case of tight overlap: amplicons of about 564 bp. Only 43% of the pairs pass DADA2's default filter, and relaxing to maxEE 2,5 merged 18.8 points more.",
    "d6300-v3v4": "NovaSeq reads with binned quality scores; the submitter had already removed the primers. The reads are clean: 95% of the pairs pass DADA2's default filter. Every setting recovered the same 10 variants and 8 of 8 strains, with 97.9% of the reads matching an expected sequence exactly.",
    "d6300-v4v5": "About 1% of the amplicons are 550 nt long against a median of 374 (a tail we suspect is eukaryotic 18S amplified by the same primers; not checked). With the default setting, which leaves out only the longest 0.2%, AmpliTrim had to keep nearly full-length reads (282/280) and DADA2 merged 57.2% of the pairs; leaving out the longest 1% gave 268/118 and 89.8%. The other tool gave no result at the longest amplicon length. The source study is not identified.",
    "v1v2": "Thirty-three samples with no known truth. The reads are clean: 91% of the pairs pass DADA2's default filter. At the same maxEE, AmpliTrim's cut merged 1.1 points more than the other tool's, and the split between the two reads differed by 14 nt (forward 225 against 211). Given the typical amplicon length, the other tool merged 43.3% against 81.2% at the longest.",
    "soil-v4": "We used this dataset while designing AmpliTrim. The merge estimate is furthest off here (+32.1 points too high) and we do not know why: a strict estimate still overshot by +25.0. Repeated runs move the cut by 13 and 11 nt along a flat stretch where the yield barely changes, not to a different outcome. Skipping truncation merged 12.7 points fewer pairs.",
    "its": "The ITS case study (the --its option). The other tool cannot be run, because it needs reads of one length. The filter estimate was within 0.1 points. The merge estimate was +21.1 points too high (with 5,000 pairs per sample). Raising maxEE to 3 gained 3.7 points on 37 of 37 samples.",
    "wastewater": "The only real dataset in these tests with poor read quality: its reads are noisy (an estimated 54% of its pairs pass DADA2's default filter), and the loss is to read quality, not to short reads (0.6% of its pairs are shorter than the cut). Its reads were also trimmed to mixed lengths and every pair was in reverse orientation. It is the one dataset where the other tool was ahead at the same maxEE: AmpliTrim chose a cut 6 nt longer (446 against 440 nt) because its estimate rewarded the extra overlap, and the extra length cost 0.9 points at the filter. AmpliTrim's cut was the most stable of any dataset across sample sets, which we had predicted would not be so. In the first round of tests, which used a stricter margin of 1.0 point, the difference was beyond the margin, against AmpliTrim.",
    "spike-in": "Mouse gut with synthetic spike-ins; the spike-in sequences were not available, so accuracy was not assessed. NovaSeq reads with binned quality scores. At the same maxEE, AmpliTrim's cut merged 1.3 points more than the other tool's (+0.4 to +2.4). In the first round of tests, which used a stricter margin of 1.0 point, this was beyond the margin in AmpliTrim's favour.",
    "swabs": "Thirty-six of 371 eligible runs, chosen by a fixed rule. The reads are clean: 85% of the pairs pass DADA2's default filter. The two tools chose the same cut. The other tool's estimate of how many pairs pass was off by +7.7 points because it saw only 9% of the pairs; AmpliTrim's was within 0.1.",
}
NAMED_REFS = {"mock-a": "7", "mock-b": "7", "its": "8", "interlab-v1v2": "10", "interlab-v4v6": "10"}     # the numbered reference in the findings, where it has one (reference 9 and 13 are placeholders)


def anon(s):
    """The pages do not name the other tool (as the README does not); the findings do."""
    return s.replace("FIGARO's", "the other tool's").replace("FIGARO", "the other tool")


def cell(s):
    s = re.sub(r"\\\[(\d+)\\\]", r"[\1]", s).replace("\\_", "_").replace("&#91;", "[")
    return anon(s.strip())


def key(name):
    n = re.sub(r"\\\[.*?\\\]", "", name)
    n = re.sub(r",\s*PRJ\w+.*$", "", n)
    n = re.sub(r"\s*\(.*?\)", "", n)
    return n.strip().lower()


def ds(name):
    k = key(name)
    return ALIAS.get(k)


def table(text, marker, stop=None):
    """(header, rows) of the first markdown table after `marker`; continuation rows (empty first cell) repeat the dataset above."""
    i = text.index(marker); seg = text[i:] if stop is None else text[i:text.index(stop, i)]
    lines = []
    for l in seg.splitlines()[1:]:
        if l.startswith("|"): lines.append(l)
        elif lines: break
    cells = [[c.strip() for c in l.strip().strip("|").split("|")] for l in lines]
    header, rows = cells[0], [r for r in cells[2:]]
    return header, rows


def num(s):
    m = re.match(r"\s*([+\-\u2212]?\d+(?:\.\d+)?)", s.replace("\u2212", "-"))
    return float(m.group(1)) if m else None


def load(archive):
    archive = pathlib.Path(archive)
    f = (archive / "docs" / "FINDINGS.md").read_text(encoding="utf-8"); v = (archive / "docs" / "validation.md").read_text(encoding="utf-8")
    D = {d[0]: {"id": d[0], "name": d[1], "accessions": d[2], "notes": NOTES[d[0]], "ref": NAMED_REFS.get(d[0])} for d in DATASETS}
    def fill(marker, stop, slot, multi=False):
        header, rows = table(f, marker, stop); header = [anon(h) for h in header]; last = None; seen = []
        if slot == "t9": header = ["Dataset", "AmpliTrim", "AmpliTrim holds", "the other tool", "the other tool holds"]      # the table has two columns called Holds
        for r in rows:
            name = r[0]
            if multi and not name: i = last
            else: i = ds(name)
            if i is None: raise ValueError(f"{slot}: dataset not recognised: {name!r}")
            last = i; rec = {anon(h): cell(c) for h, c in zip(header, r)}
            if "(" in name and re.search(r"\((\d+)%\)", name): rec["pass_default"] = int(re.search(r"\((\d+)%\)", name).group(1))
            if multi: D[i].setdefault(slot, []).append(rec)
            else: D[i][slot] = rec
            seen.append(i)
        return header, seen
    meta = {}
    h1, _ = fill("**Table 1.", "**DADA2 runs.**", "t1"); meta["t1_header"] = h1
    h3, _ = fill("**Table 3.", "**Table 4.", "t3"); meta["t3_header"] = h3
    h4, _ = fill("**Table 4.", "### 4.3", "t4"); meta["t4_header"] = h4
    h5, _ = fill("**Table 5.", "On D6300 V4-V5", "t5"); meta["t5_header"] = h5
    h6, _ = fill("**Table 6.", "### 4.5", "t6", multi=True); meta["t6_header"] = h6
    h7, _ = fill("**Table 7.", "AmpliTrim measures the amplicon", "t7"); meta["t7_header"] = h7
    h8, _ = fill("**Table 8.", "### 4.6", "t8"); meta["t8_header"] = h8
    h9, _ = fill("**Table 9.", "### 4.8", "t9"); meta["t9_header"] = h9
    h11, _ = fill("**Table 11.", "### 5", "t11") if "### 5" in f[f.index("**Table 11."):] else fill("**Table 11.", "## 5. Discussion", "t11"); meta["t11_header"] = h11
    # the final build's rule (validation page): five datasets where it differs from the default
    hv, rows = table(v, "DADA2's merged share of the input pairs on the five datasets", "It merges more pairs than the default")
    for r in rows:
        i = ds(r[0])
        if i: D[i]["rule"] = {h: c for h, c in zip(["dataset", "maxEE 2", "2,5", "the rule", "the rule's maxEE"], [x.replace("`", "") for x in r])}
        elif r[0].startswith("mean"): meta["rule_mean"] = [x.replace("`", "") for x in r]
    # the predictions of night 2 and the arms
    _, pr = table(f, "**Table 10.", "### 4.9"); meta["predictions"] = [[cell(c) for c in r] for r in pr]
    _, arms = table(f, "**Table 2.", "**Measures.**"); meta["arms"] = [[cell(c) for c in r] for r in arms]
    return D, meta


def pct(s):
    m = re.search(r"(\d+(?:\.\d+)?)%", s or "")
    return float(m.group(1)) if m else None


def groups(D):
    """Real samples, mock communities and wastewater, as the README groups the comparison of Table 4."""
    kind = {i: D[i]["t3"]["Kind"] for i in D}
    out = {"real": [], "mock": [], "wastewater": []}
    for i, d in D.items():
        if "t4" not in d: continue
        diff = num(d["t4"]["Difference"])
        out["wastewater" if i == "wastewater" else kind[i]].append((i, diff))
    stats = {}
    for g, rows in out.items():
        vals = [x for _, x in rows]
        stats[g] = {"n": len(vals), "ahead": sum(x > 0 for x in vals), "tied": sum(x == 0 for x in vals), "behind": sum(x < 0 for x in vals), "mean": sum(vals) / len(vals)}
    allv = [x for rows in out.values() for _, x in rows]
    stats["all"] = {"n": len(allv), "ahead": sum(x > 0 for x in allv), "tied": sum(x == 0 for x in allv), "behind": sum(x < 0 for x in allv), "mean": sum(allv) / len(allv)}
    return stats


def derived(D):
    ORDER = IDS
    st = groups(D)
    fe = [num(D[i]["t3"]["Filter error"]) for i in ORDER]
    ranks = [(i, num(D[i]["t3"]["Rank correlation"])) for i in ORDER if num(D[i]["t3"]["Rank correlation"]) is not None]
    ye = [num(D[i]["t3"]["Yield excess"]) for i in ORDER]
    t9 = [D[i]["t9"] for i in ORDER if "t9" in D[i]]
    return {"stats": st, "max_abs_fe": max(abs(x) for x in fe), "mean_abs_fe": sum(abs(x) for x in fe) / len(fe), "n_rank": len(ranks), "n_rank_one": sum(1 for _, r in ranks if r == 1.0),
            "low_rank": [(i, r) for i, r in ranks if r < 1.0], "ye_min": min(ye), "ye_max": max(ye), "hold_a": sum(t["AmpliTrim holds"] == "yes" for t in t9), "hold_o": sum(t["the other tool holds"] == "yes" for t in t9), "n_t9": len(t9)}


# ---- plain words for the labels of the findings (display only: the data keep the findings' own cells)
ROLE_PLAIN = {"development": ("dev", "Used in the design", "used while designing AmpliTrim, so not independent evidence"),
              "held out": ("heldout", "Independent test", "kept aside and not used in the design"),
              "case study": ("case", "ITS case study", "a case study of the --its option")}
ARM_PLAIN = {"A (maxEE 2)": "AmpliTrim, maxEE 2", "B (3,5)": "AmpliTrim, maxEE 3,5", "B2 (5)": "AmpliTrim, maxEE 5", "B4 (4,6)": "AmpliTrim, maxEE 4,6",
             "the other tool (6,9)": "The other tool, maxEE 6,9", "E (no truncation)": "No truncation, maxEE 2",
             "A, B, B2, B4, the other tool": "AmpliTrim at each maxEE tried, and the other tool", "every arm": "Every setting tried",
             "A, B, B2, B4, E": "AmpliTrim at each maxEE tried, and no truncation", "the other tool (typical length)": "The other tool, given the typical amplicon length",
             "A with the longest 1% left out": "AmpliTrim with the longest 1% of amplicons left out"}
NO_COMPARISON = {"cannot be run": "The other tool cannot be run on this dataset, because it needs reads of one length.",
                 "no result": "The other tool gave no result on this dataset, so there is no comparison at the same maxEE."}
VALUE_PLAIN = {"no truth": "no expected sequences", "unknown": "not known", "mixed": "mixed", "noisier": "noisier", "more accurate": "more accurate", "as accurate": "as accurate"}


def role_plain(role):
    for k, v in ROLE_PLAIN.items():
        if role.startswith(k): return v
    raise ValueError(role)


# the one dataset on which the cap changed the rule (TEST_REGISTER.md, section 24): with the cap of 5 the rule set 4,5 and DADA2 merged 55.0%; the released rule has no cap
CAP_NOTE = {"mock-b": {"capped_max_ee": "4,5", "capped_merged": "55.0%"}}
