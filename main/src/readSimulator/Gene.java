package readSimulator;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

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
}
