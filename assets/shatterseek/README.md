# ShatterSeek wrapper

`run_shatterseek.R` turns lrsomatic caller output into ShatterSeek input, runs ShatterSeek, labels each chromosome's
cluster with the PCAWG confidence tiers and draws a figure for every called chromosome. The `SHATTERSEEK` module
stages this directory, plus `assets/reconplot/` whose `R/parsers/` read ASCAT, Wakhan, Severus and SAVANA, so
ShatterSeek and ReConPlot see the same SVs and segments.

ShatterSeek is the [robert-a-forsyth/ShatterSeek](https://github.com/robert-a-forsyth/ShatterSeek/tree/chm13-support)
fork inside the module's image (recipe in `containers/shatterseek/`): T2T-CHM13 support, the statistical criteria fixed
to the PCAWG definitions, one cluster per chromosome, and `classify_chromothripsis()`.

Conversion to ShatterSeek's format:

- SVs: insertions, single breakends and anything without two `+`/`-` strands are dropped, as are SVs with a breakend
  outside chromosomes 1-22 and X. Intra-chromosomal SVs are typed by strand (`+-` DEL, `-+` DUP, `++` h2hINV,
  `--` t2tINV), inter-chromosomal ones are TRA. Breakends are ordered so pos1 <= pos2, strands moving with them.
- Copy number: total copy number rounded to an integer (SAVANA's is fractional), chromosomes 1-22 and X, overlaps
  trimmed, and neighbouring segments with the same copy number merged (ShatterSeek counts copy-number states row by row).

Both tables are written next to the summary (`<sample>.shatterseek_{sv,cn}_input.tsv`). Tier cut-offs default to PCAWG's
and can be changed with `--use`, `--alpha`, `--min-sv-hc1`, `--min-osc-hc`, `--min-sv-hc2`, `--min-tra-hc2`,
`--min-sv-hc3`, `--min-osc-lc`, `--max-osc-lc` and `--no-hc3`; see `Rscript run_shatterseek.R --help`.
