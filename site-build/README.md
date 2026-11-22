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

The **results of the validation** are built the same way: from `docs/FINDINGS.md` (the 16S findings) and from the files the runs committed to `results/`, nothing typed in. Every page of both websites has the same header (`header.py`: the brand, the tabs Software | How it works | Validation, and a second row of links).

```bat
python site-build\build_results.py
python site-build\build_16s.py
python site-build\build_its.py
python -m pytest tests\test_16s_site.py tests\test_its_site.py tests\test_site_header.py -q
```

- `build_16s.py` writes the overview of both markers (`docs/index.html`), `docs/16s/` (an overview and one page for each 16S dataset) and a small "moved" page at each address of the first version (`docs/datasets/<old name>/`). Its numbers are those of `docs/FINDINGS.md`, read by `results_data.py`; FIGARO is named, and the letters the nights gave their arms are replaced by words.
- `build_its.py` writes `docs/its/` (an overview, a page for each ITS dataset, the twelve trims, the other settings) and the pages for both markers: `docs/recovery/`, `docs/criteria/`, `docs/definitions/`, `docs/runs/` and `docs/datasets/index.html`. It reads `results/*` (the metrics tables, `V04_SUMMARY.md`, `RECOVERY.md`, the criteria of each run's summary, `predict.json`), `docs/DATASET_NAMES.tsv` and `docs/DEFINITIONS.md`.
- `build_results.py` writes only the per-dataset text reports (`docs/reports/`). `results_data.py` reads the tables of the findings into one record per dataset and checks that every dataset is recognised.
- Datasets and settings have readable names (`names.py`, from `docs/DATASET_NAMES.tsv`; the archive's codes appear only inside `<code>` paths). To add a dataset, add a row to that file; to name a new setting, add it to `SETTINGS` in `names.py`. The tests rebuild every page and fail, naming the page, until you rebuild.

| file | what it holds |
|---|---|
| `build_site.py`, `template.html`, `copy.js` | the software's home page |
| `build_design.py`, `design_svgs.py` | the how-it-works page and its schematic pictures |
| `header.py` | the header of every page of both websites; no builder writes its own |
| `style.css`, `design_css.css`, `results_css.css`, `its_css.css` | the styles. The colours are defined once, at the top of `style.css`, for the light and the dark theme |
| `theme_head.js`, `theme.js` | the light/dark switch: the first sets the saved theme before the page is drawn, the second adds the button |
| `shell.py` | the page shell of the validation results (head, footer, theme); its header comes from `header.py` |
| `build_16s.py`, `build_its.py`, `names.py` | the 16S pages and the overview, the ITS and cross-marker pages, and the names of datasets and settings |
| `build_results.py`, `results_data.py` | the text reports and the reading of the findings |
