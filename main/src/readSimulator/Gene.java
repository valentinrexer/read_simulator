package readSimulator;

import org.apache.commons.math3.distribution.NormalDistribution;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Gene {
    private final String geneId;
    private final String chromosome;
    private int length;
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

    public void buildTranscriptSequences(IndexedFastaReader reader) {
        for (Transcript transcript : transcriptMap.values()) {
            transcript.buildSequence(reader);
        }
    }

    public String getGeneId() {
        return geneId;
    }

    public String getChromosome() {
        return chromosome;
    }

    public List<Transcript> getTranscripts() {
        return new ArrayList<>(transcriptMap.values());
    }

    public void getRandomReadsForTranscript(String transcriptId, int initialFragmentLength, double standardDeviation, int count) {
        Transcript transcript = getTranscript(transcriptId);
        int[] randomLengths = new int[count];
        int[] randomStartingPositions = new int[count];

        FragmentRandomSampler sampler = new FragmentRandomSampler();
        sampler.initRandomSamples(count, initialFragmentLength, standardDeviation, transcript.length(), randomLengths, randomStartingPositions);

    }
}
