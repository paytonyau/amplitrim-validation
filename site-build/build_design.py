"""Build docs/how-it-works/index.html of THIS archive (its site: paytonyau.github.io/amplitrim-validation/how-it-works/). Numbers come from the code of the software repository, and the
module table from its source files, so the page cannot drift.
Usage: python build_design.py [path to the software repository [folder to write into]]   (defaults: the repository next to this archive; this archive's docs/)"""
import ast, pathlib, sys
HERE = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent.parent / "amplitrim"     # the software repository, next to this archive
OUT = pathlib.Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else HERE.parent / "docs"                       # where the page is written: this archive's site
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(HERE))
import design_svgs as S
import amplitrim
from amplitrim.subsample import DEFAULT_MAX_PAIRS
from amplitrim.pipeline import Settings, RULE_Q, WARN_ABOVE, MIN_PAIRS_FOR_ESTIMATE
from amplitrim.config import Config

REPO, BLOB = "https://github.com/paytonyau/amplitrim", "https://github.com/paytonyau/amplitrim/blob/main/"
VAL = "https://github.com/paytonyau/amplitrim-validation/blob/main/docs/validation.md"
HOME = "https://paytonyau.github.io/amplitrim/"                                    # the software's home page, on its own site
cfg = Config(); step = {f.name: f.default for f in Settings.__dataclass_fields__.values()}["step"]
F = dict(pairs=f"{DEFAULT_MAX_PAIRS:,}", overlap=cfg.dada2_min_overlap, unmerge=f"{cfg.max_unmergeable_frac * 100:g}%", step=step,
         tol=f"{cfg.yield_tolerance * 100:g} point", q=int(RULE_Q * 100), warn=WARN_ABOVE, minpairs=MIN_PAIRS_FOR_ESTIMATE, version=amplitrim.__version__)

modules = ["cli", "config", "subsample", "trimmer", "merge_wrapper", "ranking", "pipeline", "report"]
roles = {"cli": "The command line: reads the arguments, runs the pipeline, prints the summary and writes the report.",
         "config": "Defaults that mirror the DADA2 call you will run (minimum overlap, mismatches, <code>maxEE</code>), and how a <code>maxEE</code> value is read.",
         "subsample": "Finds the read pairs and takes the shared, repeatable sample.",
         "trimmer": "Holds the sample in memory and replicates DADA2's filter: the quality table, filtering at a cut, and the error counts.",
         "merge_wrapper": "Runs vsearch to merge reads, and reads its log.",
         "ranking": "Defines a candidate cut, the floor, and the rule that picks the recommendation.",
         "pipeline": "The whole run: the steps in order, the <code>maxEE</code> rule, the ITS mode and the warnings.",
         "report": "Writes <code>REPORT.md</code>, <code>summary.json</code> and <code>summary.png</code>, and the summary for the terminal."}
def uses(mod):
    t = ast.parse((ROOT / "amplitrim" / f"{mod}.py").read_text(encoding="utf-8")); s = set()
    for n in ast.walk(t):
        if isinstance(n, ast.ImportFrom) and n.level == 1:
            if n.module: s.add(n.module.split(".")[0])
            else: s.update(a.name for a in n.names)
    return sorted(s & set(modules))
code_rows = "".join(f"<tr><th scope=\"row\"><code>{m}.py</code></th><td>{roles[m]}</td><td>{', '.join(f'<code>{u}</code>' for u in uses(m)) or 'nothing else in the package'}</td></tr>" for m in modules)

def fig(svg, cap="Schematic, not your data."): return f"<figure>{svg}<figcaption>{cap}</figcaption></figure>"
def step(n, title, text, where, svg):
    return (f'<li class="step"><div><h3><span class="num" aria-hidden="true">{n}</span><span>{title}</span></h3>{text}<p class="where">In the code: {where}</p></div>{fig(svg)}</li>')

steps = "\n".join([
 step(1, "Take a sample of your reads",
      f"<p>AmpliTrim takes every Nth read pair across all your samples, up to {F['pairs']} pairs, and uses the same pairs for everything below. The sample is the same every time you run it.</p>",
      "<code>build_subsample</code>", S.sample()),
 step(2, "Measure your amplicons",
      f"<p>It merges the full-length reads with vsearch. The lengths of the merged reads are your amplicon lengths, so you do not need to look them up.</p><p>Only pairs that overlap at full length can be seen. If fewer than {F['minpairs']} pairs do, the run stops with a message that says so.</p>",
      "<code>estimate_insert_sizes</code>", S.hist()),
 step(3, "Set the floor",
      f"<p>AmpliTrim leaves out the longest {F['unmerge']} of the amplicons and adds {F['overlap']}, DADA2's minimum overlap. The result is the shortest total cut length that can merge the rest, and AmpliTrim never tests anything shorter.</p>",
      "<code>coverage_floor</code>", S.hist(True, F["unmerge"])),
 step(4, "Try cuts above the floor",
      f"<p>From the floor up to the full read length, in steps of {F['step']} nt, AmpliTrim asks the <strong>quality table</strong> for the best split between the forward and reverse reads. That table holds the share of pairs that pass DADA2's filter at every cut.</p>"
      f"<p>It then filters the reads at that cut and merges them with vsearch. The share that passes times the share that merges is the <strong>yield</strong>. The search stops when a longer cut gains less than {F['tol']}.</p>",
      "<code>_search</code>, <code>pass_table</code>, <code>run_vsearch</code>", S.search()),
 step(5, "Pick the cut",
      f"<p>Cuts whose yield is within {F['tol']} of the best are treated as tied. Among them, AmpliTrim picks the one with the least overlap, which means the shortest reads and the fewest noisy bases. A cut that cannot merge your longest amplicons is never picked.</p>",
      "<code>select_recommendation</code>", S.pick()),
 step(6, "Set maxEE",
      f"<p>At the chosen cut, for each read, AmpliTrim finds the limit that covers three quarters of the reads: the {F['q']}th percentile of their expected errors, plus 1. It is never below 2, DADA2's default. There is no upper limit: the report notes it when a limit goes above {F['warn']}.</p>"
      "<p>If the limit is not the default for both reads, the search and the pick are repeated with it.</p>",
      "<code>_rule_max_ee</code>", S.ee()),
 step(7, "Write the report",
      "<p>You get <code>REPORT.md</code> (the settings, copy-paste snippets for DADA2 and QIIME 2, what to expect, and warnings), <code>summary.json</code> (the same result for scripts) and <code>summary.png</code> (a figure).</p>",
      "<code>write_outputs</code>", S.report()),
])
its = (f'<li class="step"><div><h3><span class="num" aria-hidden="true">+</span><span>With <code>--its</code></span></h3>'
       "<p>ITS amplicons vary a lot in length, so DADA2's ITS workflow does not truncate. AmpliTrim evaluates exactly that: each read is cut at its first base of quality 2 or below, then kept if it is long enough, free of N and within <code>maxEE</code>.</p>"
       "<p>The floor, the search and the <code>maxEE</code> rule are skipped. The report shows where the pairs go instead.</p>"
       f'<p class="where">In the code: <code>_run_its</code></p></div>{fig(S.its())}</li>')

def row(choice, why, chip, klass, how):
    return f'<tr><th scope="row">{choice}</th><td>{why}</td><td><span class="chip {klass}">{chip}</span><p>{how}</p></td></tr>'
decisions = "\n".join([
 row("Replicate DADA2's quality filter", "The share of pairs that pass is then computed, not guessed.", "Held", "held", "The predicted pass rate was within 0.8 points of DADA2's in every test dataset."),
 row(f"Measure merging with real merges, at DADA2's minimum overlap of {F['overlap']}", "Merging depends on the real overlaps and mismatches in your reads, which arithmetic cannot know.", "Held", "held", "With these cuts, DADA2 merged as many pairs as with the other tool's cuts: a tie."),
 row("Find the amplicon length from your reads", "You do not have to look it up, and it is right for your data.", "Design choice", "design", "Only pairs that overlap at full length are seen, so the longest amplicons are under-represented."),
 row(f"Use a sample of at most {F['pairs']} pairs", "The run time and the memory use stay bounded, however large your dataset is.", "Not yet timed", "timing", "The bound is by design. We have not yet timed the released version."),
 row("Count merges with a tolerant setting for the yield", "A strict count tracked DADA2 less well in our tests.", "Optimistic", "optimistic", "The yield is an estimate made before DADA2 denoises the reads. Use it to compare settings, and run DADA2 once to confirm."),
 row("Break ties by the least overlap", "Shorter reads mean fewer noisy bases.", "Design choice", "design", "Cuts within the tolerance of the best are equivalent on yield, so this is a matter of preference."),
 row(f"Set <code>maxEE</code> to cover three quarters of the reads, never below 2", "Reads with noisy ends fail DADA2's default limit, and a looser limit can recover them.", "Sensible default, not independently tested", "unvalidated",
     "It follows DADA2's author: keep maxEE 2 unless your reads need more, and relax it somewhat when they do. The 75th percentile is our own choice, made after seeing the validation data (the author gives no ceiling, so the rule has none; the report notes a limit above 5), and the rule did not reach the result we had decided in advance would make us adopt it. We keep it as a trade-off that we state openly."),

 row("Offer few options", "Three required arguments, <code>--its</code> and <code>--verbose</code> leave little room to misuse it.", "Design choice", "design", "Everything else is chosen for you, including <code>maxEE</code>."),
])

CSS = (HERE / "style.css").read_text(encoding="utf-8") + (HERE / "design_css.css").read_text(encoding="utf-8")
THEMEHEAD, THEME = (HERE / "theme_head.js").read_text(encoding="utf-8"), (HERE / "theme.js").read_text(encoding="utf-8")
page = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>How AmpliTrim works: the design in pictures</title>
<meta name="description" content="How AmpliTrim chooses DADA2 truncation lengths and maxEE: what it computes exactly, what it measures, what it does not model, and where each step is in the code.">
<meta name="theme-color" content="#1b201d">
<meta property="og:type" content="website">
<meta property="og:title" content="How AmpliTrim works">
<meta property="og:description" content="The design of AmpliTrim in pictures: what it computes, what it measures and what it does not model.">
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Crect width='32' height='32' fill='%231b201d'/%3E%3Ctext x='16' y='24' font-size='24' text-anchor='middle' fill='%23e8ebe6' font-family='Georgia,serif' font-weight='700'%3EA%3C/text%3E%3C/svg%3E">
<script>{THEMEHEAD}</script>
<style>
{CSS}
</style>
</head>
<body>
<a class="skip" href="#main">Skip to content</a>
<header class="site"><div class="wrap bar">
<a class="brand" href="{HOME}">Ampli<span class="trim">Trim</span></a>
<nav aria-label="Main"><ul>
<li><a href="{HOME}">Home</a></li>
<li><a href="../">Results</a></li>
<li><a href="#idea">The idea</a></li>
<li><a href="#cut">The cut</a></li>
<li><a href="#steps">Steps</a></li>
<li><a href="#choices">Choices</a></li>
<li><a href="#code">Code</a></li>
<li><a href="{REPO}">GitHub</a></li>
</ul></nav>
</div></header>
<main id="main">
<section class="hero" id="top"><div class="wrap">
<h1>How Ampli<span class="trim">Trim</span> works</h1>
<p class="tagline">It computes what it can compute exactly, and measures what it cannot.</p>
<p class="lead">This page shows the design one step at a time. The pictures are schematic: they explain the idea and are not drawn from your data.</p>
</div></section>

<section class="wrap" id="idea">
<h2>The idea</h2>
<p class="intro">DADA2 does three things to your reads, in this order. AmpliTrim treats each one differently.</p>
<ul class="trio">
<li><span class="chip exact">Replicated exactly</span><h3>Filter and trim</h3><p>AmpliTrim applies the same rules as <code>filterAndTrim</code>: cut each read at the chosen length, drop reads that have a base of quality 2 or below, or any N, in the part you keep, and drop reads over the expected-error limit (<code>maxEE</code>). The share that passes is computed, not guessed.</p></li>
<li><span class="chip none">Not modelled</span><h3>Denoise</h3><p>Denoising (<code>dada</code>) is DADA2's most expensive step, and a quick test cannot reproduce it. This is why AmpliTrim's yield is an estimate, and an optimistic one.</p></li>
<li><span class="chip measured">Measured on raw reads</span><h3>Merge pairs</h3><p>DADA2's <code>mergePairs</code> works on denoised reads. AmpliTrim merges your trimmed reads with vsearch instead, at DADA2's minimum overlap of {F['overlap']}, so you see how many pairs really merge.</p></li>
</ul>
</section>

<section class="wrap" id="cut">
<h2>What the cut does</h2>
<p>A forward and a reverse read each start at one end of the amplicon, the stretch of DNA that was amplified. <code>truncLenF</code> and <code>truncLenR</code> say where to cut each read. Together they have to be at least as long as the amplicon plus {F['overlap']}, so that the reads overlap, and no longer than the quality allows.</p>
<p class="legend"><span><i class="sw f"></i>forward read</span><span><i class="sw r"></i>reverse read</span><span><i class="sw o"></i>overlap</span><span><i class="sw w"></i>noisy end</span></p>
<ul class="cases">
<li><h3>Too short</h3>{fig(S.cut(100, 100, 0, "a gap: no overlap"), "Schematic, not your data.")}<p>The reads stop before they meet. There is a gap, so they cannot merge.</p></li>
<li><h3>Just enough</h3>{fig(S.cut(150, 150, 0, "they overlap"), "Schematic, not your data.")}<p>The reads overlap by at least {F['overlap']} nt, DADA2's minimum. They merge.</p></li>
<li><h3>Too long</h3>{fig(S.cut(250, 250, 40, "a long overlap, noisy ends"), "Schematic, not your data.")}<p>They overlap, but the last bases of each read are noisy. More reads fail the quality filter, and you lose them.</p></li>
</ul>
<p>AmpliTrim looks for the shortest cut that still lets your longest amplicons merge, and splits it between the two reads where the quality allows.</p>
</section>

<section class="wrap" id="steps">
<h2>Step by step</h2>
<p class="intro">AmpliTrim runs these steps on a sample of your reads.</p>
<ol class="steps">
{steps}
{its}
</ol>
</section>

<section class="wrap" id="choices">
<h2>Why these choices</h2>
<p class="intro">Every choice is listed with how well it held up in our tests. Some are design choices with no test behind them, and we say which.</p>
<div class="scroll"><table class="decisions"><thead><tr><th scope="col">Choice</th><th scope="col">Why</th><th scope="col">How it held up</th></tr></thead><tbody>
{decisions}
</tbody></table></div>
<p>The full results, including what did not hold, are on the <a href="{VAL}">validation page</a>. The user guide explains <a href="{BLOB}docs/usage.md#how-maxee-is-chosen">how <code>maxEE</code> is chosen</a>.</p>
</section>

<section class="wrap" id="limits">
<h2>What AmpliTrim does not do</h2>
<ul>
<li><strong>It does not denoise.</strong> The yield is measured before DADA2 denoises the reads, so DADA2 can keep fewer pairs.</li>
<li><strong>It checks quality and mergeability, not accuracy.</strong></li>
<li><strong>It does not recommend settings for other denoisers,</strong> such as Deblur.</li>
<li><strong>It cannot help when your reads do not overlap.</strong> If your amplicons are longer than your two reads can span, no cut can merge them.</li>
</ul>
</section>

<section class="wrap" id="code">
<h2>Where things are in the code</h2>
<p>The package is eight modules. The last column shows which others each one uses.</p>
<div class="scroll"><table class="mini mods"><thead><tr><th scope="col">Module</th><th scope="col">What it does</th><th scope="col">Uses</th></tr></thead><tbody>
{code_rows}
</tbody></table></div>
<h3>If you want to change something</h3>
<div class="scroll"><table class="mini"><thead><tr><th scope="col">To change</th><th scope="col">Look in</th></tr></thead><tbody>
<tr><th scope="row">A default, such as the minimum overlap or the share of long amplicons set aside</th><td><code>config.py</code>, and the constants at the top of <code>subsample.py</code> and <code>pipeline.py</code></td></tr>
<tr><th scope="row">How reads are filtered</th><td><code>trimmer.py</code></td></tr>
<tr><th scope="row">How the best cut is chosen</th><td><code>ranking.py</code> (<code>select_recommendation</code>)</td></tr>
<tr><th scope="row">How <code>maxEE</code> is set</th><td><code>pipeline.py</code> (<code>_rule_max_ee</code>)</td></tr>
<tr><th scope="row">What the report says</th><td><code>report.py</code></td></tr>
</tbody></table></div>
<p>Contributions are welcome: see <a href="{BLOB}CONTRIBUTING.md">CONTRIBUTING.md</a>.</p>
</section>
</main>
<footer><div class="wrap">
<p>AmpliTrim {F['version']} · <a href="{HOME}">Home</a> · <a href="../">Results</a> · <a href="{REPO}">Source on GitHub</a> · <a href="{VAL}">Validation page</a> · <a href="{REPO}/issues">Report an issue</a></p>
<p>By <a href="https://orcid.org/0000-0002-3283-0370">Payton (Tung On) Yau</a>.</p>
</div></footer>
<script>
{THEME}
</script>
</body>
</html>
'''
dest = OUT / "how-it-works" / "index.html"; dest.parent.mkdir(parents=True, exist_ok=True); dest.write_text(page, encoding="utf-8")
print("written", dest, len(page.encode()) // 1024, "KB")
