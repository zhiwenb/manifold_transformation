# Reference audit

`comparison.pdf` and `index.html` compare original PDF crops (left) with
independently generated analytical plots (right). Each row in `coverage.csv`
records whether the panel was checked or is still missing its full source.

- Fig2 D/E/F, Fig3 A/E/F and Fig5 B/C: analytical content agrees visually after
  correcting the documented plot types, labels and source choices.
- Fig3 B-D: original saved variation estimates reproduce bar means and the
  displayed significance pattern; random intervals and layout are not identical.
- Fig4 E/F: original prototype caches restore the reference distributions and
  ** / *** significance patterns. The caches contain 70 paired observations.
- Fig4 B-D: original mixed bar means reproduce the reference. The original
  star assignment uses a stale shared p value, so reference significance is
  not validated. `fig4_geometry_statistics.json` records actual calculations;
  `fig4_geometry_checked_statistics.json` is from the per-bar corrected plot.
- Fig1, Fig2 A-C/G, Fig4 A and Fig5 A: full generation source not established.

Source PDF rendering permits visual comparison; it cannot recover every original
underlying data value. No exact full-page or pixel identity is claimed, and no
reference crop is substituted for generated analytical data.
