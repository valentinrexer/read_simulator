package readSimulator;

import java.util.List;

public class ReadGenerationEvent {
    private final String chromosome;
    private final String geneId;
    private final String transcriptId;
    private final byte[] fwSeqeunce;
    private final byte[] rwSequence;
    String fwRegVec;
    String rwRegVec;
    String transcriptFwRegVec;
    String transcriptRwRegVec;
    List<Integer> fwMutations;
    List<Integer> rwMutations;

    public ReadGenerationEvent(String chromosome,
                               String geneId,
                               String transcriptId,
                               byte[] fwSequence,
                               byte[] rwSequence,
                               String fwRegVec,
                               String rwRegVec,
                               String transcriptFwRegVec,
                               String transcriptRwRegVec,
                               List<Integer> fwMutations,
                               List<Integer> rwMutations) {
        this.chromosome = chromosome;
        this.geneId = geneId;
        this.transcriptId = transcriptId;
        this.fwSeqeunce = fwSequence;
        this.rwSequence = rwSequence;
        this.fwRegVec = fwRegVec;
        this.rwRegVec = rwRegVec;
        this.transcriptFwRegVec = transcriptFwRegVec;
        this.transcriptRwRegVec = transcriptRwRegVec;
        this.fwMutations = fwMutations;
        this.rwMutations = rwMutations;
    }

    public byte[] getFwSeqeunce() {
        return fwSeqeunce;
    }

    public byte[] getRwSequence() {
        return rwSequence;
    }

    public String getFwRegVec() {
        return fwRegVec;
    }

    public String getRwRegVec() {
        return rwRegVec;
    }

    public String getTranscriptFwRegVec() {
        return transcriptFwRegVec;
    }

    public String getTranscriptRwRegVec() {
        return transcriptRwRegVec;
    }

    public List<Integer> getFwMutations() {
        return fwMutations;
    }

    public List<Integer> getRwMutations() {
        return rwMutations;
    }

    public String getChromosome() {
        return chromosome;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public static String formattedFwMutations(List<Integer> mutations) {
        StringBuilder s = new StringBuilder();

        for (int i = 0; i < mutations.size(); i++) {
            s.append(mutations.get(i));
            if (i != mutations.size() - 1) {
                s.append(",");
            }
        }
        return s.toString();
    }

    public String getFormattedFwMutations() {
        return formattedFwMutations(fwMutations);
    }

    public String getFormattedRwMutations() {
        return formattedFwMutations(rwMutations);
    }
}
