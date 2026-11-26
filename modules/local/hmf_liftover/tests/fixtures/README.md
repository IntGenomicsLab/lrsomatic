Hand-made HMF_LIFTOVER fixtures. source_to_target.chain maps source chr1:0-1000 onto target chr1:100-1100
(+100) and source chr1:2000-2600 onto target chr2:1900-2500 on the reverse strand; source chr1:1000-2000 is a
gap. target.fa is ACGT repeated on chr1 and chr2 (3 kb each). bundle/ is a mini oncoanalyser-style GRCh38 bundle:
GENEA lifts cleanly, GENEB has an exon in the gap (dropped), GENEC sits on the inverted chain (strand flips);
of the hotspots one lifts, one has a REF mismatch on the target, one is in the gap and one is on the inverted chain.
Chains 3 and 4 overlap on the source: a long primary block (source chr1:3000-4000 -> target chr2:0-1000, score 5000)
under 150 one-base blocks of a low-scoring chain (source chr1:3001-3300 -> target chr1). GENED (chr1:3101-3200, not in
the panel) has more than 64 of those blocks before its end, so a fixed look-back would miss the primary block; it must
lift to chr2:101-200.
