# lrsomatic_report

Standalone reporting tool for the [LRSomatic](https://github.com/nf-core/lrsomatic) Nextflow pipeline. Generates a self-contained HTML report per sample with:

- **Summary header**: purity, ploidy, coverage, N50, variant counts
- **Circos plot**: somatic SNVs (6-class SBS colours), non-BND SVs, ASCAT copy number, translocation links
- **Breakend circos**: a second circos over just the chromosomes a breakend touches, with panel gene bodies and names; selecting a row in the SV table highlights that rearrangement's arc, and the gene-panel filter dims the arcs it hides
- **Interactive variant table**: VEP-annotated somatic small variants, optionally filtered to a gene panel, with VAF, depth and phasing
- **Interactive SV table**: Severus structural variants, one row per rearrangement with both breakend loci, sharing the same gene-panel filter (matched on breakend position)
- **QC details**: mosdepth coverage, samtools flagstat, cramino read stats, and per-chromosome WhatsHap phasing statistics (germline)

## Quick start

```bash
S=/path/to/sample-dir

Rscript bin/render_report.R \
  --sample-dir  $S \
  --sample-id   SAMPLE_ID \
  --sex         male \
  --reference   auto            # auto-detects t2t vs hg38 from VCF headers
```

The output file `SAMPLE_ID_report.html` will be written to the current directory.

Matched and tumour-only runs take the same command: the run mode and every input file
are discovered from the sample directory.

Tables render unfiltered. Add `--gene-panel lymphoid` (or a path to your own TSV) to have a
panel selected when the report opens — see [Gene panels](#gene-panels).

## All options

```
--sample-dir   Path to the sample output directory (required)
--sample-id    Sample identifier (default: directory name)
--reference    t2t | hg38 | auto  (default: auto)
--sex          male | female | XY | XX  (required)
--gene-panel   none | builtin panel name (e.g. lymphoid) | path to a custom TSV
               (default: none — tables render unfiltered)
--output       Output HTML path  (default: <sample-id>_report.html in current dir)
--title        Report title
```

> **Changed in v1.1.0:**
> - `--mode` and `--somatic-vcf` were removed. Run mode is derived from whether normal-side
>   QC is present, and the VCF supplying VAF is now discovered (see below), so neither needs
>   to be declared. Scripts passing them will fail on an unknown option.
> - `--gene-panel` now defaults to `none` instead of `lymphoid`: reports open unfiltered
>   unless a panel is asked for. Pass `--gene-panel lymphoid` to restore the old default.

> **Changed in v1.2.1:**
> - The SV table is **one row per rearrangement**, not one per breakend record: Severus's
>   `_1`/`_2` mate records are collapsed, and each row carries both loci
>   (`chrom_a`/`pos_a`, `chrom_b`/`pos_b`) plus a `svclass` separating interchromosomal
>   translocations from intra-chromosomal breakends. The SV count and the circos links
>   halve accordingly — they were double-counting. The single `gene_hits` column is
>   replaced by per-breakend `gene_a`/`gene_b`, and a `panel_hit` column names which panel
>   gene matched, on which side, and how: each entry reads `GENE (side, how)` — for
>   example `RB1 (span, direct)` when the span overlaps the gene, or `RB1 (B, 188.5 kb)`
>   when breakend B only fell inside the search window. The windows are wide (1 Mb around a
>   BND), so that second token is what separates a disrupted gene from a nearby one.
> - SVs are matched against a gene panel by **coordinate**, not gene symbol, whenever the
>   panel carries `chrom`/`start`/`end` — see [Gene panels](#gene-panels). The bundled
>   `lymphoid.tsv` is replaced by `lymphoid.hg38.tsv` and `lymphoid.t2t.tsv`; a custom
>   symbol-only TSV still works and still matches on symbols.

## Gene panels

Reports are **unfiltered by default**. `--gene-panel` only chooses which panel is selected when
the report opens; the rendered HTML always contains every variant and every builtin panel, so a
reader can switch panels (or paste a custom gene list) in the browser without re-rendering.

Built-in panels live in `assets/gene_lists/`. Each is a TSV with a `gene` column (HGNC symbols)
and, optionally, `chrom`/`start`/`end` — which changes how structural variants are matched:

| Panel columns | Small variants | Structural variants |
|---|---|---|
| `gene` only | symbol match | symbol match on the annotated breakend genes — no positional window |
| `gene, chrom, start, end` | symbol match | within **1 Mb of either breakend** of a BND, or **100 kb of the span** of any other type |

Coordinate matching is the reliable mode: whether VEP annotates a breakend with a gene symbol
at all depends on the sample's VEP invocation (1.6%–90% of breakend rows across the samples
measured), so a symbol-only panel can hide exactly the translocations it exists to find. The
note under the SV table says which mode is in force.

Because coordinates are only valid for one genome, a coordinate panel must declare its
reference (a leading `# reference: hg38` line, or a `reference` column) and a mismatch with the
rendered reference is a hard error. Builtins ship one file per reference and are offered as a
single entry, resolved against the detected one.

| Panel | Description |
|---|---|
| `lymphoid` | 72 recurrently mutated genes in B-cell lymphomas (DLBCL, FL, CLL, MCL, BL, MALT), as `lymphoid.hg38.tsv` and `lymphoid.t2t.tsv` |

```bash
--gene-panel lymphoid                # open with the builtin lymphoid panel applied
--gene-panel /path/to/my_genes.tsv   # must have a 'gene' column or be a single-column file
```

A `--gene-panel` value that is neither `none`, a builtin name, nor an existing file is an error —
a typo will not silently produce an unfiltered report. See
[`assets/gene_lists/README.md`](assets/gene_lists/README.md) for the full file format.

## Expected input layout

The `--sample-dir` must be the root of a single-sample LRSomatic output. Files are discovered
**recursively** by their distinctive filename suffix, so they can be nested in any directory
structure underneath it — for example:

```
<sample-dir>/
├── *_SOMATIC_VEP.vcf.gz                                 VEP-annotated somatic small variants
├── variants/phased/somatic_smallvariants.vcf.gz         VAF / depth / phasing source
├── severus_somatic.vcf.gz                               Severus SV calls
├── *_SV_VEP.vcf.gz                                      VEP-annotated SVs
├── *.segments_raw.txt, *.purityploidy.txt               ASCAT
├── *.mosdepth.summary.txt, *.mosdepth.global.dist.txt   mosdepth (tumour)
├── *_cramino.txt, *.flagstat, *.stats                   cramino / samtools (tumour)
├── qc/whatshap_stats/*_whatshap_stats.tsv               phasing statistics (germline)
├── wakhan/                                              Wakhan copy-number solutions
└── **/normal/**                                         same QC file set, normal side
                                                         (matched mode; e.g. qc/normal/)
```

Normal-side QC is picked up from any `normal/` directory in the tree, wherever the pipeline
nests it, and is also what determines the run mode.

**Small variants come from the VEP annotation only.** `*_SOMATIC_VEP.vcf.gz` defines the
variant set; VAF, depth, genotype and phase set are joined from the VCF that VEP annotated —
`variants/phased/somatic_smallvariants.vcf.gz`. Runs predating `variants/phased/` fall back to
`variants/clairs{,to}/somatic.vcf.gz` (then any non-germline VCF in those directories), which
yields VAF and depth but no phase set. If none is found the table still renders, without
those columns.

A footnote under the variant table names the file those columns actually came from and how
many variants they cover. **One VCF supplies them for the whole table**, so if the run
combined several somatic callers by consensus, the VAF, depth, genotype and phase set come
from whichever caller won each merge and need not match the `callers` column beside them.
Variants reported by more than one caller are highlighted in that column, and the footnote
says so. (Only the joined columns are ambiguous — the variant *set* is taken per record from
the VEP file.)

VEP writes indels at a different position and sometimes in a different allele notation than
the VCF it was given, so the join is made on a normalised key — see `variant_key()` in
`R/parse_smallvariants.R`.

`*_SOMATIC_VEP.vcf.gz` ships in two formats. Usually it is VEP *default text output* despite
the `.vcf.gz` name. If VEP was run with `--vcf` it is a genuine VCF with a `CSQ` field, and
may be a *merged* multi-caller VCF carrying germline calls (DeepVariant, Clair3) alongside
somatic ones, tagged in `INFO/CALLER`; the report then keeps only `PASS` records from a
somatic caller (ClairS, ClairS-TO, DeepSomatic). Both formats are detected automatically.

Missing files are handled gracefully: the corresponding report section shows a "not available" notice.

### Not covered: methylation

There is no methylation section. The only methylation output the pipeline publishes is
`methylation/<type>/modkit_pileup/<id>.bed.gz` — measured at 38–42 GB gzipped per sample,
unfiltered and with no tabix index, which cannot be read at render time.

`modkit pileup` already runs with `--bgzf`, so emitting a `tabix -p bed` index next to the
`.bed.gz` would be enough to unblock this: region queries on an indexed pileup measured
~0.16 s per Mb, making a binned genome-wide profile or per-locus lookup practical.

## Supported references

| `--reference` | Cytobands source | chr1 length |
|---|---|---|
| `t2t` | CHM13v2.0 | 248,387,328 bp |
| `hg38` | GRCh38 (UCSC) | 248,956,422 bp |

Auto-detection reads `##contig` lines from the VEP somatic VCF.

## R package requirements

`recipe/meta.yaml` is the source of truth for runtime dependencies — it is what the
Bioconda package and the pipeline's container are built from. The list below mirrors it;
if the two ever disagree, the recipe is right.

Install in your R environment if missing:

```r
install.packages(c("data.table", "dplyr", "DT", "htmltools", "optparse",
                   "quarto", "yaml", "ggplot2", "svglite", "knitr",
                   "R.utils", "base64enc"))
BiocManager::install("circlize")
```

Plus the `quarto` CLI itself. `R.utils` is not called directly — `data.table::fread()`
requires it to read the gzipped VCFs.

Tested with R 4.4.1 and Quarto 1.5.57.

## Repository structure

```
lrsomatic_report/
├── bin/render_report.R          CLI entrypoint
├── R/
│   ├── utils.R                  Shared helpers (gene panel, Extra-field parser)
│   ├── references.R             Cytoband + chrom-length loading, reference auto-detection
│   ├── locate_outputs.R         Discover per-tool output files in a sample directory
│   ├── parse_smallvariants.R    VEP text + raw caller VCF parsers; build variant table
│   ├── parse_severus.R          Severus VCF parsing (mate-collapsed), SV table, panel matching
│   ├── parse_ascat.R            ASCAT segments + purity/ploidy parsers
│   ├── parse_qc.R               Mosdepth, cramino, flagstat parsers
│   ├── circos.R                 draw_circos() — the genome-wide circos SVG
│   └── circos_bnd.R             draw_bnd_circos() — the breakend circos, inline and row-linked
├── templates/per_sample.qmd    Quarto template (HTML report)
├── assets/
│   ├── references/{t2t,hg38}/  Cytobands + chrom lengths (bundled, no network needed)
│   └── gene_lists/             lymphoid.{hg38,t2t}.tsv + README
└── tests/                      Unit tests (testthat)
```

## Roadmap

- **v2**: Cohort report (oncoprint, recurrence tables across multiple samples)
- **v2**: Wakhan haplotype-resolved copy-number integration
