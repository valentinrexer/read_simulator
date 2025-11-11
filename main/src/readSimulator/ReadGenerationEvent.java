package readSimulator;

import java.util.List;

public class ReadGenerationEvent {
    private byte[] fwSeqeunce;
    private byte[] rwSequence;
    String fwRegVec;
    String rwRegVec;
    String transcriptFwRegVec;
    String transcriptRwRegVec;
    List<Integer> fwMutations;
    List<Integer> rwMutations;

    public ReadGenerationEvent(byte[] fwSequence, byte[] rwSequence) {
        this.fwSeqeunce = fwSequence;
        this.rwSequence = rwSequence;

    }
}
