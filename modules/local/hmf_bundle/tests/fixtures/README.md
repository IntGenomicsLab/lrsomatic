`hmf_report.test.tar.gz`: the report module's synthetic HMF bundle (`modules/local/lrsomaticreport/tests/hmf`, made-up
coordinates) packed like a hosted bundle: one top-level directory, oncoanalyser layout, built reproducibly
(`tar --sort=name --mtime='2026-12-01 00:00Z' --owner=0 --group=0 --numeric-owner | gzip -n -9`).
