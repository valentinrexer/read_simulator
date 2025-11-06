package readSimulator;


import java.util.HashMap;

public class Gene {
    private String geneId;
    private String chromosome;
    private HashMap<String, Transcript> transcriptMap;

    public Gene(String geneId, String chromosome) {}

    public void addTranscript(Transcript transcript) {
        transcriptMap.put(transcript.getTranscriptId(), transcript);
    }

    public Transcript getTranscript(String transcriptId) {
        return transcriptMap.get(transcriptId);
    }

    public boolean hasTranscript(String transcriptId) {
        return transcriptMap.containsKey(transcriptId);
    }
}
