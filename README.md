# Training effect figures — original neuron collections

A standalone MATLAB figure repository containing plotting code, historical FR,
saved SDI/SVM/decoding and manifold results, and the supplied Fig1–5 reference PDFs.
All runtime inputs are inside `data/`; no original local folders are required.
English code comments. Original neuron collections and panel-specific dates are
preserved. No SNR or V2 neuron screening is performed.

## Run

MATLAB R2024b, Statistics and Machine Learning Toolbox, and Optimization Toolbox
were used. Open MATLAB in this checkout and run:

```matlab
run_figures() % All available analytical panel groups
run_figures({'fig4BCD','fig4E','fig4F','fig5B','fig5C'})
run_figures('fig4BCD_mixed') % Historical mixed variation-M1/global-M2 geometry
```

Results are written to `outputs/`. `fig4BCD` uses global data for both monkeys,
as requested previously; `fig4BCD_mixed` preserves the historical mixed input.
Neither changes the original neuron collection. No additional S2 global analysis
or new pooled statistics are introduced.

```sh
python3 validate_bundle.py
```

## Figure coverage

| Figure panels | Entry points | Packaged inputs |
|---|---|---|
| Fig2 D/E/F | fig2D, fig2E, fig2F | Original FR and historical SDI/SVM results |
| Fig3 A/B/C/D/E/F | fig3A, fig3BCD, fig3E, fig3F | Variation FR, manifold tables, decoding/SDI |
| Fig4 B/C/D/E/F | fig4BCD, fig4E, fig4F | Global manifold tables and old-neuron global decoding/SDI |
| Historical Fig4 B/C/D | fig4BCD_mixed | Variation M1 and global M2 manifold tables |
| Fig5 B/C | fig5B, fig5C | Global FR with frozen historical dates |
| Fig1; Fig2 A–C; Fig4 A; Fig5 A | Reference PDFs only | Full drawing source has not been established |

**Coverage limitation:** the complete generation chain for the schematic/example
panels in the last row is still unidentified. Their supplied PDFs are included,
but they are not claimed as independently regenerated panels. The analytical
scripts preserve historical statistical choices and inconsistencies. Manifold
plots read saved estimates; decoding plots read saved accuracy matrices. This
repository reproduces plotting from supplied inputs, not a fresh raw-spike
analysis or a validated common training-stage model.

## Layout

- `panels/`: English-commented drawing scripts.
- `data/historical_fr/`: frozen original FR snapshots.
- `data/code_manifold_combine/`: historical SDI/SVM/decoding and manifold caches.
- `data/audit_20261005/old_global_followup/`: old-neuron global inputs and FR.
- `references/`: supplied figure PDFs.
- `provenance/`: input checksums and historical panel/data inventory.
- `results/`: checked-in previews from the successful 13-entry validation run.
- `outputs/`: regenerated figures and summaries (ignored by Git).

All 13 supplied entry points were run successfully on 2026-10-06.
Results are recorded in `provenance/run_validation.json`. A separate copied
checkout also successfully generated Fig4E using only its bundled inputs; see
`provenance/relocation_validation.json`.
Bootstrap/permutation intervals may vary slightly because legacy random sampling
is retained; the data and neuron collections are unchanged.

Git LFS tracks MATLAB data and binary figures. A clone requires Git LFS and
`git lfs pull`. GitHub repository: https://github.com/zhiwenb/manifold_transformation
No redistribution license is assigned on behalf of the data owners.
