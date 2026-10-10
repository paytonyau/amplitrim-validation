"""
Simulates a paired-end 16S V3-V4 mock community (2x250 -> cropped to 248, as
in the original SOP) with Illumina-like quality decay and sequencing errors
drawn from the quality scores. Every read pair's TRUE insert length and taxon
is recorded, so tool output can be judged against ground truth.

This is a smoke-test/benchmark helper, not a substitute for real
mock-community data (see README).
"""

import gzip
import random
from pathlib import Path

import numpy as np

COMP = str.maketrans("ACGT", "TGCA")
FWD_PRIMER = "CCTACGGGAGGCAGCAG"        # 341F, N/W resolved to concrete bases (17 nt)
REV_PRIMER = "GACTACCAGGGTATCTAATCC"     # 805R, H/V resolved to concrete bases (21 nt)
READ_LEN = 248                           # after Trimmomatic CROP:248 in the SOP


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def make_taxa(n_taxa=8, base_insert=428, seed=1):
    """Related-but-distinct inserts: shared ancestor, ~10% substitutions,
    a few indels -> natural length variation around base_insert."""
    rng = random.Random(seed)
    ancestor = "".join(rng.choice("ACGT") for _ in range(base_insert + 12))
    taxa = []
    for _ in range(n_taxa):
        seq = list(ancestor)
        for i in range(len(seq)):
            if rng.random() < 0.10:
                seq[i] = rng.choice("ACGT")
        for _ in range(rng.randint(0, 3)):          # indels
            pos = rng.randint(40, len(seq) - 40)
            if rng.random() < 0.5:
                del seq[pos]
            else:
                seq.insert(pos, rng.choice("ACGT"))
        insert = "".join(seq)[:base_insert + rng.randint(-6, 6)]
        taxa.append(insert)
    return taxa


def quality_profile(read_len, start, mid, end):
    x = np.arange(read_len)
    return np.interp(x, [0, read_len * 0.6, read_len - 1], [start, mid, end])


def simulate_reads(templates, abundances, n_pairs, rng):
    """Returns r1, r2 (lists of str), q1, q2 (uint8 arrays), truth rows."""
    n_taxa = len(templates)
    taxon_idx = rng.choice(n_taxa, size=n_pairs, p=abundances)
    r1 = np.empty((n_pairs, READ_LEN), dtype="<U1")
    r2 = np.empty((n_pairs, READ_LEN), dtype="<U1")
    full_templates = [FWD_PRIMER + t + revcomp(REV_PRIMER) for t in templates]
    for i, t in enumerate(full_templates):
        idx = np.where(taxon_idx == i)[0]
        a = np.array(list(t[:READ_LEN]))
        b = np.array(list(revcomp(t)[:READ_LEN]))
        r1[idx] = a
        r2[idx] = b

    prof1 = quality_profile(READ_LEN, 37, 33, 27)
    prof2 = quality_profile(READ_LEN, 36, 28, 17)       # reverse reads degrade faster
    q = []
    for prof in (prof1, prof2):
        offset = rng.normal(0, 2.0, size=(n_pairs, 1))
        noise = rng.normal(0, 2.5, size=(n_pairs, READ_LEN))
        qq = np.clip(np.rint(prof + offset + noise), 2, 41).astype(np.uint8)
        q.append(qq)
    # a subset of reverse reads "crash" partway through (common on Illumina)
    crash = rng.random(n_pairs) < 0.08
    crash_pos = rng.integers(150, 230, size=n_pairs)
    for i in np.where(crash)[0]:
        q[1][i, crash_pos[i]:] = rng.integers(5, 13, size=READ_LEN - crash_pos[i])

    bases = np.array(list("ACGT"))
    for reads, qq in ((r1, q[0]), (r2, q[1])):
        # NOTE: cast to float BEFORE negating - qq is uint8, and -uint8 wraps
        # around (e.g. -37 -> 219), which made every base an "error".
        p_err = 10.0 ** (-(qq.astype(np.float64)) / 10.0)
        err = rng.random(qq.shape) < p_err
        shift = rng.integers(1, 4, size=qq.shape)
        cur = np.searchsorted(bases, reads)
        new = bases[(cur + shift) % 4]
        reads[err] = new[err]

    truth = [(taxon_idx[i], len(templates[taxon_idx[i]])) for i in range(n_pairs)]
    return r1, r2, q[0], q[1], truth


def _fastq_lines(prefix, reads, quals, mate):
    out = []
    for i in range(len(reads)):
        seq = "".join(reads[i])
        qs = "".join(chr(33 + int(v)) for v in quals[i])
        out.append(f"@{prefix}{i} {mate}:N:0:1\n{seq}\n+\n{qs}\n")
    return out


def write_dataset(out_dir: Path, n_samples=3, pairs_per_sample=10000, seed=7):
    """Writes:
         figaro_input/<sample>_R{1,2}.fastq.gz   248 nt, primers present
         primer_trimmed/<sample>_R{1,2}.fastq.gz primers removed (231 / 227 nt)
         truth.tsv                                read id, taxon, true insert length
    """
    out_dir = Path(out_dir)
    (out_dir / "figaro_input").mkdir(parents=True, exist_ok=True)
    (out_dir / "primer_trimmed").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    taxa = make_taxa()
    abundances = np.array([0.30, 0.20, 0.15, 0.10, 0.10, 0.07, 0.05, 0.03])
    truth_rows = ["read_id\ttaxon\ttrue_insert_len"]
    for s in range(n_samples):
        name = f"sample{chr(65 + s)}"
        r1, r2, q1, q2, truth = simulate_reads(taxa, abundances, pairs_per_sample, rng)
        prefix = f"{name}_"
        for mate, reads, quals, trim, in ((1, r1, q1, len(FWD_PRIMER)), (2, r2, q2, len(REV_PRIMER))):
            with gzip.open(out_dir / "figaro_input" / f"{name}_R{mate}.fastq.gz", "wt", newline="\n") as fh:
                fh.writelines(_fastq_lines(prefix, reads, quals, mate))
            with gzip.open(out_dir / "primer_trimmed" / f"{name}_R{mate}.fastq.gz", "wt", newline="\n") as fh:
                fh.writelines(_fastq_lines(prefix, reads[:, trim:], quals[:, trim:], mate))
        for i, (tx, ins) in enumerate(truth):
            truth_rows.append(f"{prefix}{i}\t{tx}\t{ins}")
    (out_dir / "truth.tsv").write_text("\n".join(truth_rows) + "\n")
    return taxa


if __name__ == "__main__":
    import sys
    target = Path(sys.argv[1] if len(sys.argv) > 1 else "sim_data")
    taxa = write_dataset(target)
    print("insert lengths per taxon:", [len(t) for t in taxa])
    print("written to", target)
