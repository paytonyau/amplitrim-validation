"""Build docs/index.html of the software repository from its README (so the two cannot differ), with the page's layout around it.
Usage: python build_site.py [path to the software repository]   (default: the folder next to this archive)"""
import html, pathlib, re, sys
HERE = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent.parent / "amplitrim"     # the software repository, next to this archive
README = (ROOT / "README.md").read_text(encoding="utf-8")
REPO = "https://github.com/paytonyau/amplitrim"
HOW = "https://paytonyau.github.io/amplitrim-validation/how-it-works/"      # the how-it-works page lives on the validation archive's site
VAL = "https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md"

def link(url):
    if url.startswith(("http://", "https://")): return url
    return f"{REPO}/blob/main/{url}"

def inline(md):
    s = html.escape(md.replace("\\*", "\u0001"), quote=False)             # a backslash-escaped asterisk is a literal one
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", lambda m: f'<a href="{link(html.unescape(m.group(2)))}">{m.group(1)}</a>', s)
    return s.replace("&lt;br&gt;", "<br>").replace("\u0001", "*")

def section(title):
    m = re.search(rf"^## {re.escape(title)}\n(.*?)(?=^## |\Z)", README, flags=re.S | re.M)
    assert m, title
    return m.group(1).strip("\n")

def bullets(text): return [re.sub(r"^- ", "", b).strip() for b in re.findall(r"^- .*(?:\n(?!- |\n).*)*", text, flags=re.M)]
def paragraphs(text):
    text = re.sub(r"```.*?```", "", text, flags=re.S); text = re.sub(r"<details>.*?</details>", "", text, flags=re.S)
    return [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip() and not p.lstrip().startswith(("|", "- ", "#"))]
def code_blocks(text): return [m.strip("\n") for m in re.findall(r"```[a-z]*\n(.*?)```", text, flags=re.S)]
def table_rows(text):
    rows = []
    for l in text.splitlines():
        if l.startswith("|") and not re.match(r"^\|[-| ]+\|$", l): rows.append([c.strip() for c in l.strip("|").split("|")])
    return rows
def table_html(rows, wide=False):
    head, body = rows[0], rows[1:]
    th = "".join(f"<th scope=\"col\">{inline(c)}</th>" if c else '<th scope="col"><span class="sr">Measure</span></th>' for c in head)
    tr = "".join("<tr>" + "".join((f'<th scope="row">{inline(c)}</th>' if i == 0 and not wide else f"<td>{inline(c)}</td>") for i, c in enumerate(r)) + "</tr>" for r in body)
    return f'<div class="scroll"><table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>'

# ---- pieces taken from the README
top = README.split("## Why use it")[0]
tagline = re.search(r"^\*\*(.+)\*\*$", top, flags=re.M).group(1)
intro = [p for p in paragraphs(top) if p.startswith("AmpliTrim is for")][0]
flow_block = [b for b in code_blocks(top) if "→" in b][0]
flow_steps = [s.strip() for s in flow_block.split("→")]
after_flow = [p for p in paragraphs(top) if p.startswith("Run AmpliTrim")][0]
sentences = re.split(r"(?<=[.!?])\s+(?=[A-Z])", intro)
hero_sents = [s for s in sentences if s.startswith(("AmpliTrim is for", "AmpliTrim tests", "It does not change"))]
problem_sents = [s for s in sentences if s not in hero_sents]
assert len(hero_sents) == 3 and len(problem_sents) == 4, (len(hero_sents), len(problem_sents))
why = bullets(section("Why use it"))
cmp_sec = section("How it compares"); cmp_paras = paragraphs(cmp_sec)
cmp_tables = [table_rows(t) for t in re.split(r"<details>", cmp_sec)[0:1]][0]
det = re.search(r"<details>\s*<summary>(.*?)</summary>(.*?)</details>", cmp_sec, flags=re.S)
qs = section("Quick start"); qs_blocks = code_blocks(qs); qs_paras = paragraphs(qs)
need = bullets(section("What you need"))
opts = paragraphs(section("Options"))[0]
how = section("How the settings are chosen"); how_rows = table_rows(how); how_para = paragraphs(how)[0]
limits = bullets(section("Limits"))
cite = paragraphs(section("Citation, licence and contributing"))[0]
assert len(why) == 4 and len(qs_blocks) == 3 and len(need) == 3 and len(limits) == 4

def card(b):
    m = re.match(r"\*\*(.+?)\*\*\s*(.*)", b, flags=re.S); return f'<li class="card"><h3>{inline(m.group(1))}</h3><p>{inline(m.group(2))}</p></li>'
def pre(block, copy=False): return f'<div class="codeblock{" copyable" if copy else ""}"><pre><code>{html.escape(block, quote=False)}</code></pre></div>'
flow_html = "".join(f'<li><span class="step{" me" if s == "AmpliTrim" else ""}">{html.escape(s)}</span>' + ('<span class="arrow" aria-hidden="true">→</span>' if i < len(flow_steps) - 1 else "") + "</li>" for i, s in enumerate(flow_steps))
example_intro = [p for p in qs_paras if p.startswith("Example output")][0]
after_example = [p for p in qs_paras if p.startswith("`report/REPORT.md`")][0]
sim_note = [p for p in qs_paras if p.startswith("These reads are simulated")][0]
footnote = [p for p in cmp_paras if p.startswith("\\*")][0]
first_cmp_para, merging_para, held_para = cmp_paras[0], [p for p in cmp_paras if p.startswith("Merging is a tie")][0], [p for p in cmp_paras if p.startswith("What held")][0]
cmp_table = [t for t in [table_rows(x) for x in re.findall(r"(?:^\|.*\n)+", re.split(r"<details>", cmp_sec)[0], flags=re.M)]][0]
det_rows = table_rows(det.group(2))

CSS = (HERE / "style.css").read_text(encoding="utf-8")
JS = (HERE / "copy.js").read_text(encoding="utf-8")
TEMPLATE = (HERE / "template.html").read_text(encoding="utf-8")
subs = {
 "CSS": CSS, "JS": JS, "THEMEHEAD": (HERE / "theme_head.js").read_text(encoding="utf-8"), "THEME": (HERE / "theme.js").read_text(encoding="utf-8"), "REPO": REPO, "HOW": HOW, "VAL": VAL,
 "TAGLINE": inline(tagline), "HERO": " ".join(inline(s) for s in hero_sents), "PROBLEM": " ".join(inline(s) for s in problem_sents),
 "FLOW": flow_html, "AFTERFLOW": inline(after_flow), "CARDS": "\n".join(card(b) for b in why),
 "CMPINTRO": inline(first_cmp_para), "CMPTABLE": table_html(cmp_table, wide=False), "FOOTNOTE": inline(footnote), "MERGING": inline(merging_para),
 "DETSUMMARY": inline(det.group(1).strip()), "DETTABLE": table_html(det_rows, wide=True), "HELD": inline(held_para),
 "BASH1": pre(qs_blocks[0], True), "EXAMPLEINTRO": inline(example_intro), "EXAMPLE": pre(qs_blocks[1]), "AFTEREXAMPLE": inline(after_example),
 "BASH2": pre(qs_blocks[2], True), "SIMNOTE": inline(sim_note),
 "NEED": "\n".join(f"<li>{inline(b)}</li>" for b in need), "OPTS": inline(opts),
 "HOWTABLE": table_html(how_rows), "HOWPARA": inline(how_para), "LIMITS": "\n".join(f"<li>{inline(b)}</li>" for b in limits),
 "CITE": inline(cite).replace("<code>LICENSE</code>", f'<a href="{REPO}/blob/main/LICENSE"><code>LICENSE</code></a>').replace("<code>CITATION.cff</code>", f'<a href="{REPO}/blob/main/CITATION.cff"><code>CITATION.cff</code></a>'),
}
out = TEMPLATE
for k, v in subs.items(): out = out.replace(f"@@{k}@@", v)
assert "@@" not in out, re.findall(r"@@\w+@@", out)
dest = ROOT / "docs" / "index.html"; dest.write_text(out, encoding="utf-8"); (ROOT / "docs" / ".nojekyll").write_text("", encoding="utf-8")
print("written", dest, len(out.encode()) // 1024, "KB")
