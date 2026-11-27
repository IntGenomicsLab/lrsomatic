---
title: Parameters
render_macros: true
---

# Pipeline parameters

Generated from [`nextflow_schema.json`](https://github.com/IntGenomicsLab/lrsomatic/blob/{{ git_ref }}/nextflow_schema.json).
Pass any of these on the command line (`--param value`) or in a `-params-file`.

{{ render_schema() }}

{{ render_schema(show_hidden=True, section_title="Hidden / advanced parameters") }}
