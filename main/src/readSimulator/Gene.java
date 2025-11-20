package readSimulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.BlockingQueue;

public class Gene {
    private final String geneId;
    private final String chromosome;
    private final HashMap<String, Transcript> transcriptMap;

    public Gene(String geneId, String chromosome) {
        this.geneId = geneId;
        this.chromosome = chromosome;
        this.transcriptMap = new HashMap<>();
    }

    public void makeTranscript(String transcriptId, char strand) {
        transcriptMap.put(transcriptId, new Transcript(transcriptId, chromosome, strand));
    }

    public Transcript getTranscript(String transcriptId) {
        return transcriptMap.get(transcriptId);
    }

    public boolean hasTranscript(String transcriptId) {
        return transcriptMap.containsKey(transcriptId);
    }

    public String getGeneId() {
        return geneId;
    }

    public void buildTranscriptSequences(IndexedFastaReader reader) {
        for (Transcript transcript : transcriptMap.values()) {
            transcript.buildSequence(reader);
        }
    }

    public void generateEventsForAllTranscripts (HashMap<String, Integer> counts,
                                              int initialFragmentLength,
                                              double standardDeviation,
                                              int readLength,
                                              double mutationRate,
                                              BlockingQueue<ReadGenerationEventChunk> queue,
                                              int CHUNK_SIZE) throws InterruptedException {

        for (String transcriptId : transcriptMap.keySet())  {
            int remaining = counts.get(transcriptId);

            while (remaining > 0) {
                int batch = Math.min(remaining, CHUNK_SIZE);

                ReadGenerationEventChunk chunk = generateRandomReadChunkForTranscript(
                        transcriptId,
                        batch,
                        initialFragmentLength,
                        standardDeviation,
                        readLength,
                        mutationRate
                );

                queue.put(chunk);

                remaining -= CHUNK_SIZE;
            }
        }
    }

    public List<ReadGenerationEventChunk> getAllChunks(HashMap<String, Integer> counts,
                                                       int initialFragmentLength,
                                                       double standardDeviation,
                                                       int readLength,
                                                       double mutationRate) {
        List<ReadGenerationEventChunk> chunks = new ArrayList<>();

        for  (String transcriptId : transcriptMap.keySet()) {
            ReadGenerationEventChunk currentChunk = generateRandomReadChunkForTranscript(
                    transcriptId,
                    counts.get(transcriptId),
                    initialFragmentLength,
                    standardDeviation,
                    readLength,
                    mutationRate);
            chunks.add(currentChunk);
        }
        return chunks;
    }

    public ReadGenerationEventChunk generateRandomReadChunkForTranscript(String transcriptId,
                                            int count,
                                            int initialFragmentLength,
                                            double standardDeviation,
                                            int readLength,
                                            double mutationRate) {

        Transcript transcript = getTranscript(transcriptId);
        int[] randomLengths = new int[count];
        int[] randomStartingPositions = new int[count];

        RandomOperationExecutor roe = new RandomOperationExecutor();
        roe.initRandomSamples(count, initialFragmentLength, standardDeviation, transcript.length(), readLength, randomLengths, randomStartingPositions);
        return new ReadGenerationEventChunk(transcript.createEventsForTranscript(
                randomLengths,
                randomStartingPositions,
                readLength,
                mutationRate,
                roe,
                chromosome,
                geneId)
        );
    }
}
