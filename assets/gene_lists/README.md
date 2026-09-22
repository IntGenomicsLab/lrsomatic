# Gene Panel Lists

`LRSOMATICREPORT` stages this directory as `--gene-lists-dir`, replacing the set bundled in
`ghcr.io/ljwharbers/lrsomatic-report`. Drop a TSV in and it becomes a builtin that
`--report_gene_panel` accepts, with no release of
[lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) needed. The files originated in
that tool (MIT licensed).

Each file is a TSV with a required `gene` column (HGNC symbol). Coordinate columns
`chrom` (or `chr`), `start` and `end` are optional but **all-or-nothing** — a file with
some but not all three is rejected rather than quietly falling back to symbol matching.

| Panel columns             | Small-variant filter | SV filter                                                                                |
| ------------------------- | -------------------- | ---------------------------------------------------------------------------------------- |
| `gene` only               | symbol match         | direct-hit symbol match on either breakend's VEP gene — no windows                       |
| `gene, chrom, start, end` | symbol match         | coordinate match: within 1 Mb of a breakend (BND) or 100 kb of the SV span (other types) |

Coordinate matching is what makes breakend filtering reliable: whether a BND carries a
VEP gene symbol at all depends on the sample's VEP invocation (1.6%–90% of breakends
across the samples measured), so a symbol-only panel can hide the very translocations it
exists to find. Matching on coordinates needs no annotation on the row.

Optional metadata columns (`panel`, `notes`) are ignored by the loader and kept for the
reader.

## Scoping a gene to one table: `applies_to`

By default every gene filters **both** tables. An optional `applies_to` column scopes a
row to one of them:

| `applies_to`           | Small-variant table | SV table |
| ---------------------- | ------------------- | -------- |
| _blank_, `both`, `all` | ✓                   | ✓        |
| `snv`, `small`         | ✓                   | —        |
| `sv`, `structural`     | —                   | ✓        |

Values are case-insensitive; **anything else is a hard error**, so a typo cannot quietly
change what is filtered. A file with no `applies_to` column behaves exactly as before.

The two tables want different genes: rearrangement partners (`IGH`, `IGK`, `IGL`, `TRA/D`, `TRB`,
`TRG`) belong in an SV panel but add noise to a small-variant one, and coding-mutation targets
such as `MYD88` or `NOTCH1` the reverse.

Coordinates are therefore required **only for rows that can match an SV** (blank or `sv`); an
`snv` row may leave `chrom`/`start`/`end` empty. The column-level all-or-nothing rule above still
applies, and a blank coordinate on a blank-or-`sv` row is still a hard error.

The column reaches the report through the tool, so a panel that uses it needs
`lrsomatic_report` ≥ 1.6.0 — the tag pinned in
[`modules/local/lrsomaticreport/main.nf`](../../modules/local/lrsomaticreport/main.nf).

## Reference declaration

Panel coordinates are only valid for the reference they were built on, so a coordinate-carrying
panel must declare it, either as a leading comment line:

```
# reference: hg38
gene	chrom	start	end
MYC	chr8	127735434	127742951
```

or as a `reference` column. A declared reference that differs from the one the report is rendered
against is a **hard error**; a panel that declares none loads with a "reference unverified"
footnote. Symbol-only panels need no declaration. Builtin panels ship one file per reference
(`lymphoid.hg38.tsv`, `lymphoid.t2t.tsv`) and are selected by their bare name.

## Supplying a custom panel

A panel does not have to live here: `--report_gene_panel` takes file paths too, staged into the
task alongside these:

```bash
nextflow run . --report_gene_panel /path/to/my_genes.tsv
nextflow run . --report_gene_panel lymphoid,/path/to/my_genes.tsv
```

A TSV in this directory is a _builtin_: its `.hg38`/`.t2t` pair collapses to one entry and it is
embedded in every report. A panel given by path is registered under its filename stem and
embedded only when named. A one-column file of symbols (with or without a `gene` header) is
accepted and gives symbol-only matching. Several panels can be combined, unioned; see
[usage](../../docs/usage.md#applying-several-panels-at-once) for the rules.

## Bundled panels

| File                      | Contents                                                                                                                                  |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `lymphoid.{hg38,t2t}.tsv` | 234 non-Hodgkin lymphoma genes, **scoped per table** — see below                                                                          |
| `sarcoma.hg38.tsv`        | 140 soft-tissue and bone sarcoma genes — tumour suppressors, amplification targets and recurrent fusion partners — GENCODE v46 gene spans |
| `sarcoma.t2t.tsv`         | the same 140 genes, spans from the CHM13v2.0 RefSeq Liftoff v5.1 annotation                                                               |

`sarcoma` carries no `applies_to` column, so every one of its genes filters both tables.

### `lymphoid`

Merged from two curated NHL lists that are deliberately not interchangeable:

- **SV list** — `NHL_Genes_{T2T,hg38}.bed`, 128 genes with per-reference spans, including
  the loci that only make sense as rearrangement partners (`IGH`, `IGK`, `IGL`, `TRA/D`,
  `TRB`, `TRG`, `DUSP22`).
- **Small-variant list** — `twist_genes.tsv`, 149 symbols with no coordinates, the panel's
  coding-mutation content; its `TSG`/`OG` class and remarks are kept in `notes`.

42 genes are on both lists and carry a blank `applies_to`; the 86 bed-only genes are `sv` and the
106 twist-only genes are `snv`, giving 234 rows. The bed coordinates are 1-based inclusive gene
spans like the other builtins (no half-open conversion). `TRA/D` is not an HGNC symbol but is an
`sv` row, so it is matched positionally.

The 106 twist-only genes are `snv`-scoped, so a deletion that removes `MYD88` or `NOTCH1`
entirely does not appear in the panel-filtered SV table; that is a one-cell edit per gene to change.

The two source lists live in the report tool's repository, so this pair is re-synced from its
`assets/gene_lists/` rather than rebuilt here.

#### Corrections applied on top of the upstream lists

This pair is **not** a byte-for-byte copy of the upstream files; re-syncing means re-applying
these, or fixing them upstream first:

| Row                      | Problem                                          | Fix                      |
| ------------------------ | ------------------------------------------------ | ------------------------ |
| `EWSR1`, `KDM6B`, `SIK3` | CHM13 coordinates in the hg38 file               | GENCODE v46 gene spans   |
| `PRKBC`                  | transposition of `PRKCB`; matched nothing        | renamed `PRKCB`          |
| `RCK`                    | obsolete alias; the coordinates are `DDX6`'s     | renamed `DDX6`           |
| `FAM46C`                 | former name of `TENT5C`, which is already listed | dropped (235 rows → 234) |

Regenerating the `sarcoma` spans is mechanical: gene spans keyed on `gene_name`, from `gene`
features (GENCODE) or the min/max of `transcript` features (Liftoff), restricted to
`chr1`–`chr22`, `chrX`, `chrY`:

- hg38: `references/GRCh38.alt-masked-V2/annotation/gencode.v46.basic.annotation.gtf.gz`
- t2t: `references/chm13_v2.0_maskedY.rCRS/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gtf`

Aliases in the source list were mapped by hand (`VEGFR2`→`KDR`, `VEGFR3`→`FLT4`, `MKL2`→`MRTFB`,
`MGEA5`→`OGA`, `HER2`→`ERBB2`, `SYT`→`SS18`, `H3F3A`→`H3-3A`, `H3F3B`→`H3-3B`) and duplicates
merged, leaving 140 genes.

Two `sarcoma` symbols need the T2T annotation handled specially, and both are recorded in
that file's comment header:

- `POU2AF3` — the Liftoff annotation predates the rename and carries it as `COLCA2`.
- `DUX4L10` — the Liftoff annotation has no D4Z4 paralogs at all (only `DUX4` itself), so
  this span comes from `GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz`, whose chromosomes
  are NCBI accessions (`NC_060925.1` = `chr1` … `NC_060948.1` = `chrY`).
