# Manifold transformation

## Data and code availability

The data and MATLAB code used to generate the analytical figure panels are available in this repository: [https://github.com/zhiwenb/manifold_transformation](https://github.com/zhiwenb/manifold_transformation).

The repository contains processed firing-rate data, saved SDI and decoding results, manifold estimates, and plotting scripts. The original neuron collections are retained without additional SNR filtering. Raw spike recordings are not included. The supplied data support the following panels: Fig. 2D–F, Fig. 3A–F, Fig. 4B–F, and Fig. 5B–C. Plotting code for the remaining panels is not included.

## Repository contents

| Directory | Contents |
|---|---|
| `data/s1`–`data/s5` | Session-specific firing-rate data and saved analysis results |
| `data/extra_var` | Additional session containing variation data only |
| `data/manifold` | Shared manifold result tables |
| `panels` | MATLAB plotting scripts |
| `results` | Exported figure panels in PDF format |

Sessions are numbered sequentially for this release:

| Session | Animal |
|---|---|
| s1 | M1 |
| s2 | M1 |
| s3 | M2 |
| s4 | M2 |
| s5 | M2 |
| extra_var | M1 |

Within each session, `fr/var` and `fr/global` contain firing-rate inputs, and `metrics/` contains saved analysis results. The variation-only session is not included in global analyses.

## Reproducing the figures

The plotting code was run using MATLAB R2024b with the Statistics and Machine Learning Toolbox and the Optimization Toolbox. Binary data and figures are stored using Git LFS. After installing Git LFS, download the repository and its data:

```sh
git clone https://github.com/zhiwenb/manifold_transformation.git
cd manifold_transformation
git lfs pull
```

From the repository directory in MATLAB, run:

```matlab
run_figures()                  % Generate all included analytical panel groups
run_figures('fig5C')           % Generate a selected panel group
run_figures('fig4BCD_checked') % Assign significance to each geometry bar
run_figures('fig4BCD_global')  % Generate geometry panels using global data
```

Figures are exported as PDF and written to `outputs/fig*/`. Required input paths are prepared automatically; no manual path edits are needed.

## Analysis notes

The plotting scripts use the supplied processed data and saved estimates. They do not reconstruct firing rates, manifold estimates, or decoding results from raw spike recordings. Session membership and historical analysis subsets are retained. Fig. 3F uses 46 paired observations; Figs. 4E–F use 70 paired observations. Bootstrap intervals and point jitter may vary between runs.

The default Fig. 4B–D script retains the original combination of variation data for M1 and global data for M2. It also retains the original significance-label behavior for figure reproduction. For significance labels calculated separately for each bar, use `fig4BCD_checked`; for global data in both animals, use `fig4BCD_global`.
