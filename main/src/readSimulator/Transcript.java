package readSimulator;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Transcript {
    private final String transcriptId;
    private final String chromosome;
    private final List<Coordinates> exonRegions;
    private byte[] sequence;

    public Transcript(String transcriptId, String chromosome) {
        this.transcriptId = transcriptId;
        this.chromosome = chromosome;
        this.exonRegions = new ArrayList<>();
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void addCoordinates(Coordinates coordinates) {
        this.exonRegions.add(coordinates);
    }

    public void buildSequence(IndexedFastaReader reader) {
        if (!reader.isOpen()) {
            System.out.println("Cannot build sequence for transcript " + transcriptId + " because reader channel is closed!");
            return;
        }
        exonRegions.sort(Comparator.comparingLong(Coordinates::coordinate1));

        int seqLength = 0;
        List<byte[]> exonVectors = new ArrayList<>();

        try {
            for (Coordinates coordinates : exonRegions) {
                byte[] exonBytes = reader.seekSequence(chromosome, coordinates.coordinate1(), coordinates.coordinate2());
                exonVectors.add(exonBytes);
                seqLength += exonBytes.length;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.sequence = new byte[seqLength];

        int seqPos = 0;
        for (byte[] exonVector : exonVectors) {
            System.arraycopy(exonVector, 0, this.sequence, seqPos, exonVector.length);
            seqPos += exonVector.length;
        }
    }

    public byte[] getSequence() {
        return sequence;
    }
}
