package readSimulator;

import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.io.*;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.BlockingQueue;

public class ParallelizedOutputWriter implements Runnable, AutoCloseable {
    private final BlockingQueue<ReadGenerationEventChunk> queue;
    private final BufferedWriter forwardWriter;
    private final BufferedWriter reverseWriter;
    private final BufferedWriter mappingInfoWriter;
    private final AtomicLong id = new AtomicLong(0);
    private volatile boolean running = true;
    private static final HashMap<Integer, String> QUALITY_CACHE = new HashMap<>();
    public static final ReadGenerationEventChunk STOPPING_SIGNAL_CHUNK = new ReadGenerationEventChunk(null);

    public ParallelizedOutputWriter(BlockingQueue<ReadGenerationEventChunk> queue, Path fwFilePath, Path rwFilePath, Path mappingInfoPath) throws IOException {
        this.queue = queue;
        this.forwardWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(fwFilePath.toFile()), StandardCharsets.UTF_8), 1 << 16);
        this.reverseWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(rwFilePath.toFile()), StandardCharsets.UTF_8), 1 << 16);
        this.mappingInfoWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(mappingInfoPath.toFile()), StandardCharsets.UTF_8), 1 << 16);
    }

    @Override
    public void close() throws Exception {
        running = false;
        forwardWriter.close();
        reverseWriter.close();
        mappingInfoWriter.close();
    }

    @Override
    public void run() {
        try {
            while (running) {
                ReadGenerationEventChunk currentChunk = queue.take();
                if (currentChunk == STOPPING_SIGNAL_CHUNK) break;

                for (ReadGenerationEvent event : currentChunk.events()) {
                    long entryId = id.incrementAndGet();
                    writeFastqEntry(forwardWriter, entryId, event.getFwSeqeunce());
                    writeFastqEntry(reverseWriter, entryId, event.getRwSequence());
                    writeMappingInfoEvent(mappingInfoWriter, entryId, event);
                }
            }
        } catch (InterruptedException | IOException e) {
            e.printStackTrace();

        } finally {
            try {
                close();
            } catch (Exception ignore) {
            }
        }

    }

    public void writeFastqEntry(BufferedWriter writer, long id, byte[] seq) throws IOException {
        writer.write("@");
        writer.write(Long.toString(id));
        writer.newLine();
        writer.write(new String(seq));
        writer.newLine();
        writer.write('+');
        writer.write(getQualityString(seq.length));
        writer.newLine();
    }

    public void writeMappingInfoEvent(BufferedWriter writer, long id, ReadGenerationEvent event) throws IOException {
        writer.write(Long.toString(id));
        writer.write("\t");
        writer.write(event.getChromosome());
        writer.write("\t");
        writer.write(event.getGeneId());
        writer.write("\t");
        writer.write(event.getTranscriptId());
        writer.write("\t");
        writer.write(event.fwRegVec);
        writer.write("\t");
        writer.write(event.rwRegVec);
        writer.write("\t");
        writer.write(event.getTranscriptFwRegVec());
        writer.write("\t");
        writer.write(event.getTranscriptRwRegVec());
        writer.write("\t");
        writer.write(event.getFormattedFwMutations());
        writer.write("\t");
        writer.write(event.getFormattedRwMutations());
        writer.newLine();
    }

    public static String getQualityString(int length) {
        return QUALITY_CACHE.computeIfAbsent(length, "I"::repeat);
    }
}
