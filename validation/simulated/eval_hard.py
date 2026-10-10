"""Independent evaluation on the stress dataset. Imports no AmpliTrim code.
For each candidate trim pair: truncate, maxEE<=2 filter, vsearch merge (minovlen 12,
stagger allowed, default diffs), then compare every merged read with its TRUE template."""
import gzip, json, subprocess, sys
from collections import Counter
from pathlib import Path
import pandas as pd

H = Path(sys.argv[1])   # base dir holding reads/, work/, out/
truth = json.load(open(H / "reads/ground_truth.json"))
cand = pd.read_csv(H / "out/candidates.csv").rename(columns={
    "truncLenF": "trim_f", "truncLenR": "trim_r", "estimated_yield": "yield_est", "recommended": "is_recommended"})
P = H / "work/subsample"


def fq(p):
    with gzip.open(p, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().strip(); f.readline(); q = f.readline().strip()
            yield h.strip(), s, q


ee = lambda q: sum(10 ** (-(ord(c) - 33) / 10) for c in q)
pairs = list(zip(fq(P / "subsample_R1.fastq.gz"), fq(P / "subsample_R2.fastq.gz")))
N = len(pairs)
kind_n = Counter(h1[1:].split()[0].split("|")[2] for (h1, _, _), _ in pairs)
tmp = H / "tmp"; tmp.mkdir(exist_ok=True)
rows = []
for _, c in cand.drop_duplicates(["trim_f", "trim_r"]).iterrows():
    tf, tr = int(c.trim_f), int(c.trim_r)
    with open(tmp / "r1.fq", "w", newline="\n") as o1, open(tmp / "r2.fq", "w", newline="\n") as o2:
        for (h1, s1, q1), (h2, s2, q2) in pairs:
            s1, q1, s2, q2 = s1[:tf], q1[:tf], s2[:tr], q2[:tr]
            if len(s1) < tf or len(s2) < tr or ee(q1) > 2 or ee(q2) > 2:
                continue
            nm = h1[1:].split()[0]
            o1.write(f"@{nm}\n{s1}\n+\n{q1}\n"); o2.write(f"@{nm}\n{s2}\n+\n{q2}\n")
    subprocess.run(["vsearch", "--fastq_mergepairs", str(tmp / "r1.fq"), "--reverse", str(tmp / "r2.fq"),
                    "--fastaout", str(tmp / "m.fa"), "--fastq_minovlen", "12",
                    "--fastq_allowmergestagger", "--quiet"], capture_output=True)
    merged = Counter(); ok_len = Counter(); wrong = Counter()
    name, seq = None, []

    def flush():
        if name is None:
            return
        kind = name.split("|")[2]; t = truth[name.split("|")[1]]["seq"]
        merged[kind] += 1
        if len("".join(seq)) == len(t):
            ok_len[kind] += 1
        else:
            wrong[kind] += 1

    for line in open(tmp / "m.fa"):
        if line.startswith(">"):
            flush(); name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    flush()
    tot_m, tot_ok = sum(merged.values()), sum(ok_len.values())
    rows.append(dict(
        trim_f=tf, trim_r=tr, tool_est=tot_m / N,
        correct_len_yield=tot_ok / N,
        valid_yield=(ok_len["normal"] + ok_len["repeat"]) / N,        # excludes chimeras
        wrong_len_normal=wrong["normal"] / kind_n["normal"],
        wrong_len_repeat=wrong["repeat"] / kind_n["repeat"],
        recall_normal=ok_len["normal"] / kind_n["normal"],
        recall_repeat=ok_len["repeat"] / kind_n["repeat"],
        recall_chimera=ok_len["chimera"] / kind_n["chimera"]))
res = pd.DataFrame(rows)
out = cand.drop_duplicates(["trim_f", "trim_r"])[["trim_f", "trim_r", "yield_est", "is_recommended"]].merge(
    res, on=["trim_f", "trim_r"])
out.to_csv(H / "hard_truth.csv", index=False)
print(f"pairs {N}; by kind {dict(kind_n)}", file=sys.stderr)
pd.set_option("display.width", 220)
print(out.sort_values("valid_yield", ascending=False).round(4).to_string(index=False))
