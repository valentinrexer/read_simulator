package readSimulator;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.logging.Logger;

public class IndexedFastaReader {
    private final Path fastaPath;
    private final FastaIndex index;
    private FileChannel channel;
    public static final Logger LOGGER = Logger.getLogger(IndexedFastaReader.class.getName());

    public IndexedFastaReader(Path fastaIndexPath, Path fastaPath) {
        this.index = new FastaIndex(fastaIndexPath);
        this.fastaPath = fastaPath;
        channel = null;
    }

    public void openChannel() throws IOException {
        channel = FileChannel.open(fastaPath, StandardOpenOption.READ);
    }

    public void closeChannel() throws IOException {
        channel.close();
    }

    public boolean isOpen() {
        return channel != null && channel.isOpen();
    }

    public byte[] seekSequence(String chr, long start, long end) throws IOException {
        FastaIndexEntry indexEntry = index.get(chr);
        if (indexEntry == null) return null;

        long basesToRead = end - start + 1;
        int lineBases = indexEntry.lineBases();
        int lineWidth = indexEntry.lineWidth();

        /*
        Calculate the starting byte in the fasta file
        Start at the offset position of the chromosome, skip lines by integer division, then skip remaining positions in
        line with mod operation
         */
        long startingBytePosition = indexEntry.offset()
                + ((start - 1) / lineBases) * lineWidth
                + ((start - 1) % lineBases);

        /*
        estimate how many bytes we have to read from the fasta (including the \n)
        Estimate number of newline characters and add a small buffer (will be discarded later if necessary)
         */
        long bytesToRead = basesToRead + (basesToRead / lineBases) + 2;

        ByteBuffer byteBuffer = ByteBuffer.allocate((int) bytesToRead);
        channel.position(startingBytePosition);
        channel.read(byteBuffer);

        byteBuffer.flip();
        byte[] rawBytes = new byte[byteBuffer.limit()];
        byteBuffer.get(rawBytes);

        ByteBuffer processedBytes = ByteBuffer.allocate((int) bytesToRead);
        for (byte b : rawBytes)
            if (b != '\n' && b != '\r') processedBytes.put(b);

        processedBytes.flip();
        int cleanLength = (int) Math.min(processedBytes.remaining(), basesToRead);
        byte[] result = new byte[cleanLength];
        processedBytes.get(result, 0, cleanLength);
        return result;

    }
}
