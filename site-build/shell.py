"""The shared page shell of the validation archive's website: the styles, the theme switch, the navigation and the footer. Used by build_results.py (the 12 pages of
the first validation) and by build_its.py (the ITS pages, the recovery, criteria, definitions, runs and dataset catalogue pages)."""
import html
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
ARCH, SOFT = "https://github.com/paytonyau/amplitrim-validation", "https://github.com/paytonyau/amplitrim"
BLOB, ENA = ARCH + "/blob/main/", "https://www.ebi.ac.uk/ena/browser/view/"
esc = lambda s: html.escape(str(s), quote=True)
CSS = "".join((HERE / n).read_text(encoding="utf-8") for n in ("style.css", "design_css.css", "results_css.css", "its_css.css"))
THEMEHEAD, THEME = (HERE / "theme_head.js").read_text(encoding="utf-8"), (HERE / "theme.js").read_text(encoding="utf-8")
FAVICON = "data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Crect width='32' height='32' fill='%231b201d'/%3E%3Ctext x='16' y='24' font-size='24' text-anchor='middle' fill='%23e8ebe6' font-family='Georgia,serif' font-weight='700'%3EA%3C/text%3E%3C/svg%3E"


def page_shell(title, desc, body, base, source_note=None):
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<meta name="description" content="{esc(desc)}">
<meta name="theme-color" content="#1b201d">
<meta property="og:type" content="website">
<meta property="og:title" content="{esc(title)}">
<meta property="og:description" content="{esc(desc)}">
<link rel="icon" href="{FAVICON}">
<script>{THEMEHEAD}</script>
<style>
{CSS}
</style>
</head>
<body>
<a class="skip" href="#main">Skip to content</a>
<header class="site"><div class="wrap bar">
<a class="brand" href="{base}">Ampli<span class="trim">Trim</span> results</a>
<nav aria-label="Main"><ul>
<li><a href="{base}">Overview</a></li>
<li><a href="{base}datasets/">Datasets</a></li>
<li><a href="{base}its/">ITS</a></li>
<li><a href="{base}recovery/">Recovery</a></li>
<li><a href="{base}criteria/">Criteria</a></li>
<li><a href="{base}definitions/">Definitions</a></li>
<li><a href="{base}how-it-works/">How it works</a></li>
<li><a href="{SOFT}">Software</a></li>
<li><a href="{ARCH}">Archive</a></li>
</ul></nav>
</div></header>
<main id="main">
{body}
</main>
<footer><div class="wrap">
<p><a href="{base}">Overview</a> · <a href="{SOFT}">AmpliTrim on GitHub</a> · <a href="{ARCH}">Validation archive</a> · <a href="{BLOB}docs/FINDINGS.md">The findings</a> · <a href="{BLOB}TEST_REGISTER.md">Test register</a></p>
<p>By <a href="https://orcid.org/0000-0002-3283-0370">Payton (Tung On) Yau</a>. {source_note or 'The numbers come from the findings of the validation archive. The tool we compare with is not named on these pages; the findings name it.'}</p>
</div></footer>
<script>
{THEME}
</script>
</body>
</html>
'''


