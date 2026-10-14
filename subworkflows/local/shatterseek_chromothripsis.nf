// IMPORT MODULES
include { SHATTERSEEK as SHATTERSEEK_SEVERUS_ASCAT  } from '../../modules/local/shatterseek/main'
include { SHATTERSEEK as SHATTERSEEK_SEVERUS_WAKHAN } from '../../modules/local/shatterseek/main'
include { SHATTERSEEK as SHATTERSEEK_SAVANA         } from '../../modules/local/shatterseek/main'

//
// ShatterSeek chromothripsis calls with PCAWG confidence tiers for every CN/SV caller pair that produced output for a
// sample: ASCAT + Severus, Wakhan + Severus, and SAVANA on its own. Pass channel.empty() for a caller that did not run.
// The pairs are built as in RECONPLOT_FIGURES, and the callers are read by the same parsers (assets/reconplot/R/parsers).
//
workflow SHATTERSEEK_CHROMOTHRIPSIS {

    take:
    severus_vcf          // [meta, severus_somatic.vcf.gz]
    ascat_segments       // [meta, segments.txt]
    wakhan_solution_dirs // [meta, [solution_rank_N/, ...]]
    savana_cna           // [meta, segmented_absolute_copy_number.tsv]
    savana_vcf           // [meta, classified.somatic.vcf]
    genome               // ShatterSeek genome: hg38, chm13 or hg19

    main:
    ch_versions = channel.empty()

    // A digest of every wrapper file goes into the meta map, so editing either wrapper invalidates -resume
    def wrapper_src = { name ->
        def dir    = file("${projectDir}/assets/${name}", type: 'dir', checkIfExists: true)
        def digest = files("${projectDir}/assets/${name}/**", type: 'file').sort { f -> f.toString() }.collect { f -> f.text }.join('\n').md5()
        return channel.value([[id: name, digest: digest], dir])
    }
    shatterseek_src = wrapper_src('shatterseek')
    reconplot_src   = wrapper_src('reconplot')

    severus_sv_files = severus_vcf.map { meta, vcf -> [meta, [vcf]] }
    // severus_sv_files: [meta, [severus_somatic.vcf.gz]]

    //
    // MODULE: SHATTERSEEK_SEVERUS_ASCAT (label: process_low)
    // Input:  [meta, 'ascat', [segments.txt], 'severus', [vcf]]
    //
    // ASCAT writes an empty segments.txt when it finds no solution: nothing to test for that pair
    ascat_segments
        .filter { _meta, seg -> seg.size() > 0 && seg.countLines() > 1 }
        .map { meta, seg -> [meta, [seg]] }
        .join(severus_sv_files)
        .map { meta, cn, sv -> [meta, 'ascat', cn, 'severus', sv] }
        .set { ascat_input }

    SHATTERSEEK_SEVERUS_ASCAT( ascat_input, shatterseek_src, reconplot_src, genome )
    ch_versions = ch_versions.mix(SHATTERSEEK_SEVERUS_ASCAT.out.versions)

    //
    // MODULE: SHATTERSEEK_SEVERUS_WAKHAN (label: process_low)
    // Input:  [meta, 'wakhan', [integer_profile.bed], 'severus', [vcf]]  -- the top-ranked solution (solution_rank_1/)
    //
    wakhan_solution_dirs
        .map { meta, dirs ->
            def best = [dirs].flatten().find { dir -> dir.name == 'solution_rank_1' }
            def bed = best ? best.resolve('integer_profile.bed') : null
            return [meta, bed?.exists() ? bed : null]
        }
        .filter { _meta, bed -> bed != null }
        .map { meta, bed -> [meta, [bed]] }
        .join(severus_sv_files)
        .map { meta, cn, sv -> [meta, 'wakhan', cn, 'severus', sv] }
        .set { wakhan_input }

    SHATTERSEEK_SEVERUS_WAKHAN( wakhan_input, shatterseek_src, reconplot_src, genome )
    ch_versions = ch_versions.mix(SHATTERSEEK_SEVERUS_WAKHAN.out.versions)

    //
    // MODULE: SHATTERSEEK_SAVANA (label: process_low)
    // Input:  [meta, 'savana', [cna.tsv, somatic.vcf], 'savana', []]  -- single-source mode
    //
    savana_cna
        .join(savana_vcf)
        .map { meta, cna, sv_vcf -> [meta, 'savana', [cna, sv_vcf], 'savana', []] }
        .set { savana_input }

    SHATTERSEEK_SAVANA( savana_input, shatterseek_src, reconplot_src, genome )
    ch_versions = ch_versions.mix(SHATTERSEEK_SAVANA.out.versions)

    emit:
    summary  = SHATTERSEEK_SEVERUS_ASCAT.out.summary.mix(SHATTERSEEK_SEVERUS_WAKHAN.out.summary, SHATTERSEEK_SAVANA.out.summary) // [meta, tsv]
    versions = ch_versions // [versions.yml]
}
