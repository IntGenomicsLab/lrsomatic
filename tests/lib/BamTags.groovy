import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.util.zip.GZIPInputStream

// Reads aux tag names directly: htsjdk's SAMRecord cannot show a tag that a record carries twice
class BamTags {

    // [header: String, records: [[name: String, flag: int, tags: [String]]]]
    static Map read(String path) {
        // BGZF is a series of gzip members, which GZIPInputStream reads back to back
        byte[] bytes = new GZIPInputStream(new FileInputStream(path), 65536).withStream { it.bytes }
        def buf = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN)
        if (new String(bytes, 0, 4, 'ISO-8859-1') != 'BAM\u0001') {
            throw new IllegalArgumentException("${path} is not a BAM file")
        }
        buf.position(4)
        int lText = buf.getInt()
        def header = new String(bytes, 8, lText, 'UTF-8').replace('\u0000', '')
        buf.position(8 + lText)
        int nRef = buf.getInt()
        nRef.times {
            int lName = buf.getInt()
            buf.position(buf.position() + lName + 4)    // name, l_ref
        }

        def records = []
        while (buf.remaining() > 0) {
            int end = buf.getInt()
            end += buf.position()
            buf.position(buf.position() + 8)                // refID, pos
            int lReadName = buf.get() & 0xff
            buf.position(buf.position() + 3)                // mapq, bin
            int nCigar = buf.getShort() & 0xffff
            int flag = buf.getShort() & 0xffff
            int lSeq = buf.getInt()
            buf.position(buf.position() + 12)               // next refID, next pos, tlen
            def name = new String(bytes, buf.position(), lReadName - 1, 'ISO-8859-1')
            buf.position(buf.position() + lReadName + 4 * nCigar + (lSeq + 1).intdiv(2) + lSeq)
            def tags = []
            while (buf.position() < end) {
                tags << new String(bytes, buf.position(), 2, 'ISO-8859-1')
                buf.position(buf.position() + 2)
                def type = new String(bytes, buf.position(), 1, 'ISO-8859-1')
                buf.position(buf.position() + 1)
                if (type in ['Z', 'H']) {
                    while (buf.get() != 0) { }
                } else if (type == 'B') {
                    def sub = new String(bytes, buf.position(), 1, 'ISO-8859-1')
                    buf.position(buf.position() + 1)
                    int n = buf.getInt()
                    buf.position(buf.position() + n * size(sub, path))
                } else {
                    buf.position(buf.position() + size(type, path))
                }
            }
            buf.position(end)
            records << [ name: name, flag: flag, tags: tags ]
        }
        [ header: header, records: records ]
    }

    private static int size(String type, String path) {
        switch (type) {
            case ['A', 'c', 'C']: return 1
            case ['s', 'S']: return 2
            case ['i', 'I', 'f']: return 4
            default: throw new IllegalStateException("${path}: unknown aux type '${type}'")
        }
    }
}
