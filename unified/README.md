# Unified analysis comparison

This branch contains a retrospective sensitivity analysis of the manuscript inputs. `main` remains a figure-only release. No SNR or V2 reselection is applied.

## Inputs and candidate definitions

The 52 main recordings comprise 22 variation and 30 global recordings in s1-s5. The historical cross-category family has 50 files. The original variation-only extra session remains in the parent release but is not added to analyses requiring both stimulus families.

- `canonical`: released `data/s*/fr/var` and `global` inputs.
- `stored_cross`: released `variation` and `prototype` inputs, including their extra U labels.
- `stored_cross_common_labels`: identical daily results restricted to the labels shared with canonical inputs. This separates category-universe differences from other source-family differences.
- `union_native`: raw CDT reconstruction using the within-recording union of units already used by the historical input versions, the canonical condition definitions, and 630 ms to the recorded stimulus offset.
- `union_common240`: the same original-unit union and condition definitions, with an identical 630-870 ms response window for every recording.
- `fixed_count_common240`: sensitivity only. Three random draws of 27 original units (seeds 11, 29, 47) per recording; the primary data retain all units.

All rebuilt inputs retain the historical pooled 350-500 ms baseline subtraction and balanced first-trial selection. `release_rebuild` verifies the canonical inputs against raw CDT. The maximum FR difference across all 52 recordings was approximately 1.6e-14 Hz. The primary collection does not assert that unit-0 labels are isolated single neurons or that every unit is V2.

## Analysis and stage rules

All candidates use identical SDI, direct geometry, and representation estimators. SVM is recomputed with linear standardized five-fold CV, identical seeds for corresponding pairs, and balanced category sample sizes. The original within-category leave-one-out SVM is additionally recomputed for canonical and union_common240. Shared pair labels must be finite in all selected early and late recordings; missing pairs are not silently rematched by position.

Matched-cohort comparisons use source-family date intersections separately within each session and modality. Complete-cohort results additionally include the two global recordings absent from the prototype inputs. Stage sensitivities include first/last thirds, first/last halves, first/last recording, and calendar thirds shared across modalities. Recording order is chronological. These are retrospective descriptive stages, not verified pretraining/post-training periods. The supplied exposure helper contains manual counts without a validated complete recording-to-exposure mapping. `training_timeline_template.csv` records the missing annotation fields.

Results are summarized once per session and separately by animal. `animal_balanced_summary.json` gives equal weight to the two animals. Session t-tests and exact two-sided sign-flip tests are exploratory; they do not establish population-level inference from two animals. Multiple sensitivity choices are reported, not selected by their p-values.

Variation DiD contrasts the retained T categories against U categories using the same trained condition indices in both groups. Global direct geometry uses R/S/T versus G in s1/s3/s5 and S versus the mean of G/H/R/T in s2/s4, following the task-specific target/control configuration. Fig5 representation metrics instead measure Polar R/S/T versus Grating G in all five sessions, with Hyperbolic summaries restricted to s1/s5. These are explicitly different estimands, not interchangeable panels.

RMS radius is divided by sqrt(number of units). Signal radius removes the expected within-stimulus trial-noise contribution to the condition-mean cloud. Signal dimension is a participation ratio of the covariance after diagonal noise subtraction and PSD clipping. These estimators are distinct from MFT radius, dimension, and capacity. Historical daily MFT caches are re-staged in a separately labeled diagnostic, without claiming raw-input provenance validation. A full 1,000-vector MFT benchmark is kept separately; it is not a whole-cohort capacity validation.

## Run

In MATLAB R2024b with Statistics and Machine Learning Toolbox and Optimization Toolbox:

```matlab
addpath('unified')
run_unified() % Reuse supplied daily outputs; recompute missing outputs
% Rebuild from original CDT files if locally available:
run_unified(true, '/path/to/CDTTable')
```

To force a numerical rerun, move the applicable `unified/daily/<candidate>` folder aside before running. Rebuilding candidate FR alone does not invalidate existing daily results automatically. The raw CDT files are external; all processed candidate FR needed for the supplied comparison are included. `raw_sources.json` contains SHA-256 hashes of the original recordings.

`summary.csv`, `summary.json`, `session_effects.json`, and `stage_manifest.json` contain the full numerical comparison. `fig2_baseline_*.json` isolates baseline-date changes while holding post dates fixed. Companion figures are PDF only. The Chinese PDF report explains the conclusions and validation limits.

## Rebuild the report

The optional Chinese report builder uses Python with reportlab and Pillow, plus Poppler `pdftoppm`. It renders companion PDFs automatically into temporary files. Set `UNIFIED_REPORT_FONT` to a CJK-capable TrueType font on non-macOS systems; the generated PDF embeds the font. Run `python3 unified/build_report.py` after the MATLAB comparisons.
