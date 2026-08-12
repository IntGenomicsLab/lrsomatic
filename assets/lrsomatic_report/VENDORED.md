# Vendored: lrsomatic_report

This directory is a **vendored copy** of the standalone report tool, not a git submodule.
Do not edit it here — fix upstream, tag a release, and re-sync.

| | |
|---|---|
| Upstream | <https://github.com/ljwharbers/lrsomatic_report> |
| Release | `v1.1.0` (`9d660a77d5f23f92e1f7ff34f85da7956f445009`) |
| Vendored commit | `d17a636aeb3f79462b7f58db9102f4030941195b` (`main`) |
| License | MIT (see `LICENSE`) |

The vendored commit is two chore commits ahead of the `v1.1.0` tag. Neither changes
behaviour: `b4cc7620` scrubs real sample identifiers out of the README and the
`--sample-id` help string, `d17a636a` drops a development helper script. Vendoring the
tag itself would publish those identifiers in this repository.

## Why vendored rather than a submodule

`nextflow run IntGenomicsLab/lrsomatic` clones the pipeline repository but does **not**
fetch git submodules, so a gitlink here would be an empty directory for every user who
did not hand-clone with `--recurse-submodules` — and for CI, whose checkout steps do not
pass `submodules: recursive`. Real tracked files work for both.

The tool's *dependencies* (R, Quarto and its R packages) are handled separately, by the
Wave multi-package container declared in `modules/local/lrsomaticreport/main.nf` and
built from that module's `environment.yml`.

## What is included

Only what `bin/render_report.R` needs at run time:

```
bin/  R/  templates/  assets/  LICENSE  README.md
```

Upstream `docs/`, `tests/` and `recipe/` are deliberately excluded — the same set marked
`export-ignore` in the upstream `.gitattributes`.

## Re-syncing on the next upstream release

```bash
TAG=v1.2.0
git clone --depth 1 --branch "$TAG" https://github.com/ljwharbers/lrsomatic_report.git /tmp/lrr
rm -rf assets/lrsomatic_report/{bin,R,templates,assets,LICENSE,README.md}
cp -a /tmp/lrr/{bin,R,templates,assets,LICENSE,README.md} assets/lrsomatic_report/
# then update the table above, and:
#  - modules/local/lrsomaticreport/environment.yml   if upstream recipe/meta.yaml gained a dependency
#  - modules/local/lrsomaticreport/main.nf           container digest + the hard-coded version topic
#  - modules/local/lrsomaticreport/meta.yml          the same version string
```

Rebuild the container after any `environment.yml` change so the image and the file agree:

```bash
wave --conda-file modules/local/lrsomaticreport/environment.yml --freeze --await
```
