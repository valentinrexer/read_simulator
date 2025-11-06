package readSimulator;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public class Transcript {
    private final String transcriptId;
    private final List<Coordinates> exonRegions;
    private byte[] sequence;

    public Transcript(String transcriptId) {
        this.transcriptId = transcriptId;
        this.exonRegions = new ArrayList<>();
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void addCoordinates(Coordinates coordinates) {
        this.exonRegions.add(coordinates);
    }

    public void initializeSequence(Path fasta, Path fastaIndex) {
    }
}
