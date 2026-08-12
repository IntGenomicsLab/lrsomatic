# Gene Panel Lists

Each file is a TSV with a required `gene` column (HGNC symbol) and optional metadata columns (`panel`, `notes`).

To supply a custom panel at render time:

```bash
Rscript bin/render_report.R \
  --sample-dir /path/to/sample \
  --sample-id MySample \
  --gene-panel /path/to/my_genes.tsv
```

The minimal format of a custom panel file is one gene symbol per line (no header needed if there is only one column, but a TSV with a `gene` header is preferred).

## Bundled panels

| File | Contents |
|---|---|
| `lymphoid.tsv` | ~70 recurrently mutated genes in B-cell lymphomas (DLBCL, FL, MCL, CLL, BL, MALT) |
