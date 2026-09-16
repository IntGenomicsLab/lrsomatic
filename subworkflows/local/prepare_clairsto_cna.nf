//
// Prepare the ClairS-TO Verdict resource directory (--cna_resource_dir)
//
// Verdict is ClairS-TO's ASCAT-like germline tagger. It ships GRCh38 loci/allele/GC/RT files and has
// no reference-build switch; its only guard is a chr1..chrX contig-name check, which CHM13 passes.
// On CHM13 it therefore runs the GRCh38 resources against T2T coordinates, matches nothing, and tags
// nothing -- silently, because every failure path inside it exits 0. Giving it a matching resource
// set through the (undocumented) --cna_resource_dir flag is what makes it work.
//
// GRCh38 needs nothing from here: the in-image resources are already correct, so this emits [] and
// the module omits the flag.
//

include { UNZIP as UNZIP_CLAIRSTO_RT   } from '../../modules/nf-core/unzip/main'
include { CLAIRSTO_CNA_RESOURCES       } from '../../modules/local/clairsto_cna_resources/main'

workflow PREPARE_CLAIRSTO_CNA {

    take:
        resource_dir   // str:  pre-built Verdict resource directory, or null
        rt_file        // str:  replication timing file (.zip or plain) for this reference, or null
        loci_files     // [path, ...]: per-chromosome ASCAT loci files  (from PREPARE_REFERENCE_FILES)
        allele_files   // [path, ...]: per-chromosome ASCAT allele files
        gc_file        // [path, ...]: ASCAT GC correction file
        genome_name    // str:  reference build name, for the process tag

    main:

        ch_versions = channel.empty()
        cna_resources = channel.value([])

        if (resource_dir) {
            // A directory prepared earlier (for example the one this subworkflow published on a
            // previous run). Validate the layout up front -- ClairS-TO would otherwise just warn
            // and carry on with germline tagging disabled.
            def dir = file(resource_dir, type: 'dir', checkIfExists: true)
            def required = [
                'loci_files/G1000_loci_hg38_chr1.txt',
                'allele_files/G1000_alleles_hg38_chr1.txt',
                'GC_G1000_hg38.txt',
                'RT_G1000_hg38.txt',
            ]
            def missing = required.findAll { f -> !file("${dir}/${f}").exists() }
            if (missing) {
                error("Path provided with --clairsto_cna_resources is not a ClairS-TO Verdict resource directory.\n" +
                      "Missing from ${dir}: ${missing.join(', ')}\n" +
                      "Expected layout (the 'hg38' in these names is ClairS-TO's hard-coded naming, not the build of the contents):\n" +
                      "  loci_files/G1000_loci_hg38_chr{1..22,X}.txt\n" +
                      "  allele_files/G1000_alleles_hg38_chr{1..22,X}.txt\n" +
                      "  GC_G1000_hg38.txt\n" +
                      "  RT_G1000_hg38.txt")
            }
            cna_resources = channel.value(dir)
        }
        else if (rt_file) {
            // Build it from the ASCAT reference files for this genome. The loci, allele and GC files
            // are the same ones the pipeline's own ASCAT step uses; only the replication timing track
            // is specific to this step.
            def ch_rt = channel.empty()
            if (rt_file.endsWith('.zip')) {
                //
                // MODULE: UNZIP_CLAIRSTO_RT (UNZIP alias; label: process_single)
                // Input:  [meta(id=basename), [zip_file]]
                // Output: .unzipped_archive -- [meta, dir]
                //
                UNZIP_CLAIRSTO_RT(channel.fromPath(file(rt_file)).collect().map { it -> [ [ id:it[0].baseName ], it ] })

                ch_rt = UNZIP_CLAIRSTO_RT.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                ch_versions = ch_versions.mix(UNZIP_CLAIRSTO_RT.out.versions)
            } else {
                ch_rt = channel.fromPath(rt_file).collect()
            }

            //
            // MODULE: CLAIRSTO_CNA_RESOURCES (label: process_single)
            // Input:  loci / allele / GC / RT files + genome name
            // Output: .cna_resources -- path to a directory in ClairS-TO's expected layout
            //
            CLAIRSTO_CNA_RESOURCES (
                loci_files,
                allele_files,
                gc_file,
                ch_rt,
                genome_name
            )

            // .first() makes this a value channel: it is one genome-level directory consumed once
            // per sample, and a queue channel would only reach the first sample.
            cna_resources = CLAIRSTO_CNA_RESOURCES.out.cna_resources.first()
            ch_versions = ch_versions.mix(CLAIRSTO_CNA_RESOURCES.out.versions)
        }
        // else: GRCh38 and any other build whose in-image resources are correct -- emit [], and
        // CLAIRSTO leaves --cna_resource_dir off entirely.

    emit:
        cna_resources          // path to the Verdict resource directory, or [] to use the in-image set
        versions = ch_versions
}
