# Analysis code

`run_analysis` recomputes metrics from the supplied firing rates, retaining every stored neuron. Results are saved under `outputs_analysis/<session>/<task>/`. The plotting inputs in `data/` are preserved.

```matlab
run_analysis('sdi', 's1')          % All available recordings in this session
run_analysis('cross_svm', 's2', 1) % First recording only
```

| Task | Input within `data/<session>/fr/` | Computation | Related panels |
|---|---|---|---|
| `sdi` | `var` | Condition-pair SDI, signal distance, and trial variance | Fig. 2E–F |
| `within_svm` | `var` | Condition-pair linear SVM, leave-one-out validation | Fig. 2E |
| `category_sdi` | `variation` | Between-category SDI | Fig. 3F |
| `cross_svm` | `variation` | Between-category linear SVM, five-fold validation | Fig. 3E |
| `global_sdi` | `prototype` | Between-category SDI | Fig. 4F |
| `global_svm` | `prototype` | Between-category linear SVM, five-fold validation | Fig. 4E |
| `manifold_var` | `var` | Manifold capacity, radius, and dimension | Fig. 3B–D; M1 in default Fig. 4B–D |
| `manifold_global` | `global` | Manifold capacity, radius, and dimension | Global Fig. 4B–D |
| `cvpca` | `global` | cvPCA neural and stimulus-parameter distance matrices | Supporting representation analysis |
| `geometry` | `var` | Direct radius, pairwise distances, and category geometry | Supporting geometry analysis |
| `internal_svm` | `var` | Mean condition-pair decoding within each category | Supporting decoding analysis |

`var` and `global` retain the firing rates used by the original figure calculations. `variation` and `prototype` provide the separate processed inputs associated with the saved between-category results. Their category layouts and recording subsets can differ. The task selects its corresponding input automatically. `extra_var` has no global or prototype inputs.

Fig. 2D and Fig. 3A perform their calculations within their plotting scripts. Fig. 5B–C also compute the orthogonality and representation statistics directly from FR in their plotting scripts; these calculations are already included in `panels/`.

## Manifold summaries

Manifold tasks save a per-recording `output` and a `session_results.mat` file. To compute differences in differences, use `calculate_manifold_did(session_files, controls, baseline_indices, after_indices, output_file)`. The three index arguments are cell arrays with one entry per session. Specify the category controls and recording indices for the intended comparison; the function does not infer new early/late boundaries. Baseline and after indices refer to the order saved in `session_results.files`. The original estimator removes empty categories before manifold analysis, so category indices refer to the retained categories.

The full manifold estimator uses 1,000 random test vectors and numerical optimization. Condition-pair SVM evaluates every pair and can also take considerable time. Recomputed estimates may differ because of random sampling. Newly computed outputs are not automatically substituted for the supplied figure caches.

## Raw CDT extraction

`extract_fr_var` and `extract_fr_global` provide the historical CDT-to-FR extraction methods with explicit input and output paths. They require raw CDT files and a caller-supplied `neuron_selector(session, filename)` returning the original 2-by-N channel/unit list. This avoids applying a new neuron-selection rule. No SNR filter is applied. These low-level extractors retain the recording identifiers used in CDT and the original stimulus definitions.

```matlab
addpath('analysis')
% Supply the recording ID, local CDT root, output directory, and exact
% neuron-list lookup from the original recordings:
extract_fr_var(recording_id, cdt_root, output_dir, neuron_selector)
```

Raw CDT files and their original neuron-list lookup are not distributed here. Raw extraction has not been validated against every released FR file, and the `variation`/`prototype` preprocessing history is not established by these extractors. Use the supplied processed FR to rerun the included analyses.

Source-file provenance is recorded in `provenance/analysis_sources.json`; the additional FR inputs are listed in `provenance/analysis_inputs.json`. Validation scope and cache comparisons are recorded in `provenance/analysis_validation.json`.
