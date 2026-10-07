import java.util.zip.GZIPInputStream

// Record-level checks on the published phased VCFs, shared by the pipeline tests.
class PhasedVcf {

    // Non-header records of a bgzipped VCF, each split into its tab-separated fields
    static List<List<String>> records(String path) {
        new GZIPInputStream(new FileInputStream(path)).withReader { reader ->
            reader.readLines().findAll { !it.startsWith('#') }.collect { it.split('\t') as List }
        }
    }

    static boolean hasFlag(List<String> rec, String flag) {
        rec[7].split(';').contains(flag)
    }

    // Value of a FORMAT field in the first sample column, or null if absent
    static String format(List<String> rec, String key) {
        def i = rec[8].split(':').toList().indexOf(key)
        def values = rec[9].split(':')
        i >= 0 && i < values.size() ? values[i] : null
    }

    static boolean isAlt(List<String> rec) {
        format(rec, 'GT').split(/[\/|]/).any { it != '0' && it != '.' }
    }

    // Default-mode checks: both arms all PASS, somatic records all carry SOMATIC, germline not empty
    static void assertPassOnly(String outdir, List<String> samples, List<String> tumourOnly) {
        samples.each { s ->
            def germ = records("${outdir}/${s}/variants/phased/germline_smallvariants.vcf.gz")
            def som  = records("${outdir}/${s}/variants/phased/somatic_smallvariants.vcf.gz")
            assert !germ.isEmpty() : "${s}: phased germline VCF has no records"
            assert germ.every { it[6] == 'PASS' } : "${s}: non-PASS germline record with the PASS filter on"
            assert som.every { it[6] == 'PASS' } : "${s}: non-PASS somatic record with the PASS filter on"
            assert som.every { hasFlag(it, 'SOMATIC') } : "${s}: somatic record without INFO/SOMATIC"
        }
        // Paired test samples have no PASS somatic calls, so the checks above only bite on tumour-only ones
        tumourOnly.each { s ->
            assert !records("${outdir}/${s}/variants/phased/somatic_smallvariants.vcf.gz").isEmpty() :
                "${s}: tumour-only phased somatic VCF has no records"
        }
    }
}
