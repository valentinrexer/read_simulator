package readSimulator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Transcript {
    private final String transcriptId;
    private final String chromosome;
    private final List<Coordinates> exonRegions;
    private final char strand;
    private byte[] sequence;

    public Transcript(String transcriptId, String chromosome, char strand) {
        this.transcriptId = transcriptId;
        this.chromosome = chromosome;
        this.exonRegions = new ArrayList<>();
        this.strand = strand;
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

        if (strand == '-') {
            reverseComplementInPlace();
        }
    }

    public byte[] getSequence() {
        return sequence;
    }

    private void reverseComplementInPlace() {
        if (sequence == null) return;

        int left = 0;
        int right = sequence.length - 1;

        while (left <= right) {
            byte leftComplement = complement(sequence[left]);
            byte rightComplement = complement(sequence[right]);

            sequence[left] = rightComplement;
            sequence[right] = leftComplement;

            left++;
            right--;
        }
    }

    private byte complement(byte base) {
        return switch (Character.toUpperCase((char) base)) {
            case 'A' -> (byte) 'T';
            case 'C' -> (byte) 'G';
            case 'G' -> (byte) 'C';
            case 'T' -> (byte) 'A';
            case 'U' -> (byte) 'A';
            default -> (byte) 'N';
        };
    }
}
