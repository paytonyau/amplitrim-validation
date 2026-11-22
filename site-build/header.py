"""The header of every page of both websites: the software's home page (build_site.py), the how-it-works page (build_design.py) and the validation results (shell.py).
One definition, so the three cannot drift apart again.

Row 1 is the same on every page: the brand, the three tabs (Software, How it works, Validation, the current one marked) and the theme button, which theme.js adds to the
`tools` box. On a narrow screen the button moves up beside the brand and the tabs take their own row.
Row 2 holds the links of the page itself (the anchors of a long page, or the sections of the results) on one row that scrolls sideways when the screen is narrow, so a long
list can never push the button or the tabs onto extra rows."""
import html

SOFTWARE_SITE = "https://paytonyau.github.io/amplitrim/"
RESULTS_SITE = "https://paytonyau.github.io/amplitrim-validation/"
TABS = (("software", "Software"), ("how", "How it works"), ("results", "Validation"))
esc = lambda s: html.escape(str(s), quote=True)


def _cur(flag):
    return ' aria-current="page"' if flag else ""


def header(active, hrefs, links, brand, here=None):
    """active: the tab of this page ('software', 'how' or 'results'); hrefs: where each tab leads from this page; links: [(label, href)] for row 2; brand: where the brand leads;
    here: the href in `links` that is the page (or its section) itself."""
    tabs = "".join(f'<li><a href="{esc(hrefs[k])}"{_cur(k == active)}>{esc(label)}</a></li>' for k, label in TABS)
    row = "".join(f'<li><a href="{esc(h)}"{_cur(h == here)}>{esc(t)}</a></li>' for t, h in links)
    return (f'<header class="site"><div class="wrap bar">\n<a class="brand" href="{esc(brand)}">Ampli<span class="trim">Trim</span></a>\n'
            f'<nav aria-label="Sites"><ul>{tabs}</ul></nav>\n<div class="tools"></div>\n</div>\n<div class="wrap subbar"><nav aria-label="Main"><ul>{row}</ul></nav></div></header>')
