// IMPORT MODULES
include { HTSLIB_BGZIPTABIX     } from '../../modules/nf-core/htslib/bgziptabix/main'
include { MODKIT_DMR            } from '../../modules/nf-core/modkit/dmr/main'
include { DMR_HAPLOTYPE_REGIONS } from '../../modules/local/dmr/haplotype_regions/main'
include { DMR_NEAREST_GENE      } from '../../modules/local/dmr/nearest_gene/main'

workflow DMR {

    take:
    modkit_bedgz     // [meta, bed.gz(es)]  -- MODKIT_PILEUP.out.bedgz; only meaningful when modkit_phased is true
    fasta            // [[:], fasta]
    fai              // [[:], fai]
    cpg_islands_bed  // path
    gencode_gene_bed // path

    main:
    ch_versions = channel.empty()

    // Pick the hp1/hp2 pair out of MODKIT_PILEUP's glob-collected output list (which also
    // includes a _combined file when --modkit_phased is set). An unphased run's single
    // unsuffixed file matches neither pattern, so it's dropped by the filter below --
    // defence in depth alongside the modkit_phased gate on the caller's side.
    modkit_bedgz
        .map { meta, files ->
            def flist = files instanceof List ? files : [files]
            def hp1 = flist.find { it.name.endsWith('_hp1.bed.gz') }
            def hp2 = flist.find { it.name.endsWith('_hp2.bed.gz') }
            return [meta, hp1, hp2]
        }
        .filter { _meta, hp1, hp2 -> hp1 && hp2 }
        .set { haplotype_bedmethyl }
    // haplotype_bedmethyl: [meta, hp1_bedgz, hp2_bedgz]

    //
    // MODULE: HTSLIB_BGZIPTABIX (label: process_low)
    // Tabix-index each haplotype's bedMethyl -- modkit dmr pair requires a .tbi alongside each
    // bgzip input. The bedMethyl is already bgzip-compressed (modkit pileup --bgzf), so this
    // only adds the index. hp1/hp2 are tagged onto meta and mixed into one call, then split
    // back apart below.
    //
    haplotype_bedmethyl
        .flatMap { meta, hp1, hp2 ->
            return [
                [meta + [haplotype: 'hp1'], hp1, [], []],
                [meta + [haplotype: 'hp2'], hp2, [], []]
            ]
        }
        .set { bgziptabix_input }
    // bgziptabix_input: [meta+haplotype, bedmethyl, [], []]

    HTSLIB_BGZIPTABIX (
        bgziptabix_input,
        'compress',
        true,
        'bed'
    )

    HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index)
        .map { meta, bedgz, tbi ->
            def haplotype = meta.haplotype
            def sample_meta = meta.findAll { it.key != 'haplotype' }
            return [sample_meta, haplotype, bedgz, tbi]
        }
        .branch { _meta, haplotype, _bedgz, _tbi ->
            hp1: haplotype == 'hp1'
            hp2: haplotype == 'hp2'
        }
        .set { indexed_branched }

    indexed_branched.hp1
        .map { meta, _haplotype, bedgz, tbi -> [meta, bedgz, tbi] }
        .set { hp1_indexed }
    indexed_branched.hp2
        .map { meta, _haplotype, bedgz, tbi -> [meta, bedgz, tbi] }
        .set { hp2_indexed }

    hp1_indexed
        .join(hp2_indexed)
        .set { dmr_haplotype_input }
    // dmr_haplotype_input: [meta, hp1_bedgz, hp1_tbi, hp2_bedgz, hp2_tbi]

    //
    // MODULE: DMR_HAPLOTYPE_REGIONS (label: process_single)
    // Restrict cpg_islands_bed to islands covered in both haplotypes.
    //
    DMR_HAPLOTYPE_REGIONS (
        dmr_haplotype_input,
        cpg_islands_bed,
        fai.first()
    )
    ch_versions = ch_versions.mix(DMR_HAPLOTYPE_REGIONS.out.versions)

    //
    // MODULE: MODKIT_DMR (label: process_medium)
    // Compare methylation between the two haplotypes over the restricted regions.
    //
    dmr_haplotype_input
        .join(DMR_HAPLOTYPE_REGIONS.out.regions_bed)
        .multiMap { meta, hp1_bedgz, hp1_tbi, hp2_bedgz, hp2_tbi, regions_bed ->
            hp1: [meta, hp1_bedgz, hp1_tbi]
            hp2: [meta, hp2_bedgz, hp2_tbi]
            regions: [meta, regions_bed]
        }
        .set { modkit_dmr_input }

    MODKIT_DMR (
        modkit_dmr_input.hp1,
        modkit_dmr_input.hp2,
        modkit_dmr_input.regions,
        fasta.first()
    )

    //
    // MODULE: DMR_NEAREST_GENE (label: process_single)
    // Annotate each DMR with its nearest gene.
    //
    DMR_NEAREST_GENE (
        MODKIT_DMR.out.bed,
        gencode_gene_bed
    )
    ch_versions = ch_versions.mix(DMR_NEAREST_GENE.out.versions)

    emit:
    nearest_gene_bed = DMR_NEAREST_GENE.out.nearest_gene_bed  // [meta, bed]  -- DMRs annotated with nearest gene
    versions         = ch_versions
}
