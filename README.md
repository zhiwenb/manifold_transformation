# Manifold transformation figures

Run the figures with the original neuron collections. No SNR screening.
Code comments are in English. All required plotting inputs are included.

## Run

MATLAB R2024b with Statistics and Machine Learning Toolbox and Optimization
Toolbox was used. In the repository folder:

```matlab
run_figures()                 % All analytical panel groups
run_figures('fig5C')          % One group
run_figures('fig4BCD_checked') % Mixed geometry with per-bar p values
run_figures('fig4BCD_global')  % Pure-global alternative
```

## Simple layout

```text
data/
  s1/ ... s5/     Each session has fr/ and metrics/
  extra_var/      Original s2, variation only
  manifold/       Shared manifold result tables
panels/           Figure drawing code
results/          Generated previews (PNG and PDF)
references/       Your original Fig1–5 PDFs
provenance/       Checksums, date inventory and session mapping
```

| New session | Original session | Animal |
|---|---|---|
| s1 | s1 | M1 |
| s2 | s4 | M1 |
| s3 | s5 | M2 |
| s4 | s6 | M2 |
| s5 | s7 | M2 |
| extra_var | s2 (variation only) | M1 |

Inside each session, `fr/var` and `fr/global` contain firing-rate inputs.
`metrics/` contains saved SDI and decoding results grouped by analysis.
Duplicate FR is stored once. Input values, neurons and date groups are unchanged;
MAT-file internal historical metadata retains the original recording IDs.
Scripts and displayed session labels use the new IDs. `extra_var` is retained
without adding it to global or pooled analyses.

`run_figures` automatically prepares disposable input views inside
`outputs_reference/inputs` for the original plotting code. You do not need to edit
paths or manage these files. Generated figures go to `outputs_reference/fig*/`.
The whole checkout can be moved to another location and used independently.

## Coverage and validation

Included analytical groups: Fig2 D/E/F; Fig3 A/B/C/D/E/F; Fig4 B/C/D/E/F;
Fig5 B/C. `fig4BCD` now reproduces the historical reference combination: variation M1
and global M2. Its historical annotation bug is explicitly documented below.
`fig4BCD_global` is the pure-global alternative. Fig4E/F now read the original
prototype caches, rather than the later recomputed caches.

Fig1, Fig2 A–C/G, Fig4 A and Fig5 A have reference PDFs, but their complete drawing
source has not been identified. They are not claimed as regenerated panels.
Manifold and decoding panels plot saved estimates, rather than rerunning raw
spike analysis. Historical statistical choices are retained. Bootstrap intervals
may vary slightly because the original random sampling is retained.

```sh
python3 validate_bundle.py
```

All 13 drawing entry points passed after renumbering. Verification records are
in `provenance/simple_layout_validation.json`. The numerical comparison with
the previous layout is recorded in `provenance/numeric_comparison.json`. Git LFS stores binary inputs and
figures; install Git LFS and run `git lfs pull` after cloning.

Repository: https://github.com/zhiwenb/manifold_transformation

## Reference audit (2026-10-06)

See [side-by-side comparison](reference_check/comparison.pdf),
[interactive image index](reference_check/index.html), and
[complete coverage table](reference_check/coverage.csv).

The analytic plot content was visually checked against the supplied PDFs;
full-page layout and exact pixels have **not** been reproduced. Standalone font
size, cropping, thumbnail assets, point jitter and bootstrap intervals can differ.
Earlier run-success and layout-renaming checks did not establish reference equality.

Corrections: Fig2E/F labels; Fig3F violin rather than delta scatter; original
prototype decoding and SDI caches for Fig4E/F; historical mixed geometry for
Fig4B-D; and isolation of figures during export to avoid saving a previous panel.

**Fig4B-D significance bug:** the original nested plotting helper comments out
`p = p_vec_dir(i)`. Its shared `p` retains the last significance calculation,
so stars are not assigned to their own bars. The historical plot can look like
the reference while conveying incorrect significance. For example, mixed M1
Dimension has p around 0.094, although the reference displays a star. Use
`fig4BCD_checked` to see per-bar stars. Random sign-flip estimates near 0.05 can
vary between reruns; no seed was selected to force a match to the reference.

Original cached Fig3F uses 46 pairs and skips the dimension-mismatched original
s4 (new s2) recording pair. Original cached Fig4E/F use 70 pairs. These historical
scopes differ from later rebuilt caches (52 or 50 pairs); they are preserved for
reference reproduction, rather than described as a common corrected pipeline.

Earlier alternatives remain available as `fig3F_delta`, `fig4E_recomputed`,
`fig4F_recomputed`, and `fig4BCD_global`.

The optional comparison report can be regenerated with
`python3 tools/create_reference_audit.py` after `run_figures()`. It requires
Pillow, reportlab and the Poppler `pdftoppm` command. MATLAB drawing does not
require these Python packages.
