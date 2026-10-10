"""Convert a simulator ground_truth.json to the label/kind/seq TSV that dada2_check.R reads.
Handles both formats: {"templates": {label: seq}} (clean simulator) and {label: {seq, kind}} (stress)."""
import json, sys
g = json.load(open(sys.argv[1]))
g = g.get("templates", g)
with open(sys.argv[2], "w") as f:
    f.write("label\tkind\tseq\n")
    for k, v in g.items():
        seq, kind = (v["seq"], v["kind"]) if isinstance(v, dict) else (v, "normal")
        f.write(f"{k}\t{kind}\t{seq}\n")
