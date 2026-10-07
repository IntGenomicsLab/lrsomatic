# ReConPlot wrapper

`run_reconplot.R` and `R/` turn lrsomatic caller output (ASCAT, Wakhan, Severus, SAVANA) into
[ReConPlot](https://github.com/cortes-ciriano-lab/ReConPlot) figures and harmonised CN/SV tables.
The `RECONPLOT` module stages this directory as its wrapper input.

Vendored from [Tim-Yu/ReConPlot](https://github.com/Tim-Yu/ReConPlot) at the commit recorded in
`VERSION`. To update, copy `run_reconplot.R` and `R/` from that repository and bump `VERSION`.

Local changes on top of the vendored wrapper commit (recorded in `VERSION` as `<commit>+lrsomatic.<date>`): tables are
written without scientific notation; an SV with an undetermined orientation is dropped instead of aborting; mate records
are collapsed on the caller's mate ID (positions + orientation as the fallback, insertions never collapsed); `--min-svlen`
applies to SAVANA too; `--genes` symbols are resolved against ReConPlot's gene table (unknown or duplicated symbols are
skipped with a warning instead of losing the panel); a failing BAF-track reader skips the track instead of failing the run;
`--genome` is validated up front; thousands separators are accepted in `--regions`; malformed and duplicate regions are
dropped; the Wakhan BED header is found with or without the comment block; `--sv-format` defaults to `auto`, which reads
SAVANA's classified somatic VCF when no BEDPE is staged (SAVANA's PacBio classifier writes no BEDPE).
