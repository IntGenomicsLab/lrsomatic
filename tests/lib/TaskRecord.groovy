// Reads the .ext_record files that tests/smallvar_wiring.config makes tasks write into their work directories.
class TaskRecord {

    // Sample (task tag) -> record (key -> value) of every task whose name ends in ':<process> (<tag>)'
    static Map<String, Map<String, String>> of(String workDir, String process) {
        def records = [:]
        new File(workDir).eachDir { d1 ->
            d1.eachDir { d2 ->
                def run = new File(d2, '.command.run')
                def rec = new File(d2, '.ext_record')
                if (!run.exists() || !rec.exists()) {
                    return
                }
                def header = run.readLines().find { it.startsWith('### name: ') }
                def m = header =~ /^### name: '.*:\Q${process}\E \((.*)\)'$/
                if (m.matches()) {
                    def tag = m.group(1)
                    assert !records.containsKey(tag) : "two ${process} tasks for ${tag}"
                    records[tag] = rec.readLines().collectEntries { line ->
                        def i = line.indexOf('=')
                        [line.substring(0, i), line.substring(i + 1)]
                    }
                }
            }
        }
        return records
    }
}
