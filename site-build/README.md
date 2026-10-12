# Building the website of the software

The website of `amplitrim` is two files, `docs/index.html` (the home page) and `docs/how-it-works/index.html`. They are generated, not written by hand, so that they cannot drift from the software:

- the **home page** is built from the software's `README.md` (every sentence, command and table row comes from it);
- the **how-it-works page** takes its numbers from the software's code (`Settings`, `Config`, the constants in `subsample.py` and `pipeline.py`) and its module table from the source files.

After you change the README, or a default in the code, rebuild and test:

```bat
python site-build\\build_site.py
python site-build\\build_design.py
cd ..\\amplitrim
python -m pytest tests\\test_site.py -q
```

The builders expect the software repository next to this archive (`..\\amplitrim`); give another path as the one argument. They use the standard library and import the software itself, so install its dependencies first (numpy).

The **results site of this archive** (`docs/index.html`, `docs/datasets/<id>/index.html`, `docs/reports/<id>.md`) is built the same way, from `docs/FINDINGS.md` and `docs/validation.md`:

```bat
python site-build\build_results.py
python -m pytest tests\test_results_site.py -q
```

`results_data.py` reads the tables of the findings into one record per dataset and checks that every dataset is recognised; `build_results.py` draws the pages and charts from those records. After the findings change, rebuild; the test fails until you do, and says which page is out of date.

| file | what it holds |
|---|---|
| `build_site.py`, `template.html`, `copy.js` | the home page |
| `build_design.py`, `design_svgs.py` | the how-it-works page and its schematic pictures |
| `style.css`, `design_css.css` | the styles (the second only for the how-it-works page and the results site). The colours are defined once, at the top of `style.css`, for the light and the dark theme |
| `theme_head.js`, `theme.js` | the light/dark switch: the first sets the saved theme before the page is drawn, the second adds the button |
| `build_results.py`, `results_data.py`, `results_css.css` | the results site of the archive: its data, its pages and charts, and its styles |

The pictures are schematic, drawn from nothing but the idea, and each is captioned so. To change one, edit its function in `design_svgs.py`.
