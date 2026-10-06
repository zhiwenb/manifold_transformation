# Manifold transformation figures

Run the figures with the original neuron collections. No SNR screening.
Code comments are in English. All required plotting inputs are included.

## Run

MATLAB R2024b with Statistics and Machine Learning Toolbox and Optimization
Toolbox was used. In the repository folder:

```matlab
run_figures()                 % All analytical panel groups
run_figures('fig5C')          % One group
run_figures('fig4BCD_mixed')  % Historical mixed geometry
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
`outputs_simple/inputs` for the original plotting code. You do not need to edit
paths or manage these files. Generated figures go to `outputs_simple/fig*/`.
The whole checkout can be moved to another location and used independently.

## Coverage and validation

Included analytical groups: Fig2 D/E/F; Fig3 A/B/C/D/E/F; Fig4 B/C/D/E/F;
Fig5 B/C. `fig4BCD` uses global data for both monkeys; `fig4BCD_mixed` preserves
the historical variation-M1/global-M2 input combination.

Fig1, Fig2 A–C, Fig4 A and Fig5 A have reference PDFs, but their complete drawing
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
