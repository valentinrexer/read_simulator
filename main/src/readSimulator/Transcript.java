package readSimulator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class Transcript {
    private final String transcriptId;
    private final String chromosome;
    private final List<Coordinates> exonRegions;
    private final char strand;
    private int[] genomicPositionMapping;
    private byte[] sequence;

    public Transcript(String transcriptId, String chromosome, char strand) {
        this.transcriptId = transcriptId;
        this.chromosome = chromosome;
        this.exonRegions = new ArrayList<>();
        this.strand = strand;
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

        if (strand == '-') reverseComplementInPlace(sequence);

        createGenomicPositionMappingArray();
    }

    public void createGenomicPositionMappingArray() {
        genomicPositionMapping = new int[sequence.length];
        int mappedPos = 0;

        for (Coordinates coordinates : exonRegions) {
            int exonStart = coordinates.coordinate1();
            int exonEnd   = coordinates.coordinate2();

            for (int g = exonStart; g <= exonEnd; g++) {
                genomicPositionMapping[mappedPos++] = g;
            }
        }

    }

    public byte[] getSequence() {
        return sequence;
    }

    public static void reverseComplementInPlace(byte[] sequence) {
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

    public static byte complement(byte base) {
        return switch (Character.toUpperCase((char) base)) {
            case 'A' -> (byte) 'T';
            case 'C' -> (byte) 'G';
            case 'G' -> (byte) 'C';
            case 'T' -> (byte) 'A';
            case 'U' -> (byte) 'A';
            default -> (byte) 'N';
        };
    }

    public int length() {
        return sequence.length;
    }

    public List<ReadGenerationEvent> createEventsForTranscript(int[] fragmentLength,
                                                               int[] startingPosition,
                                                               int readLength,
                                                               double mutationRate,
                                                               RandomOperationExecutor roe,
                                                               String chromosome,
                                                               String geneId) {
        List<ReadGenerationEvent> events = new ArrayList<>();

        for (int i = 0; i < fragmentLength.length; i++) {
            int plusStart = startingPosition[i];
            int plusEnd = startingPosition[i] + readLength;

            int minusStart = startingPosition[i] + fragmentLength[i] - readLength - 1;
            int minusEnd = startingPosition[i] + fragmentLength[i] - 1;

            byte[] fwReadSequence = Arrays.copyOfRange(sequence,
                    plusStart,
                    plusEnd);

            byte[] rwReadSequence = Arrays.copyOfRange(sequence,
                    minusStart,
                    minusEnd);

            reverseComplementInPlace(rwReadSequence);

            String genomicForwardReadVector;
            String genomicReverseReadVector;

            if (strand == '-') {
                genomicForwardReadVector = (genomicPositionMapping[plusStart] + 1) + "-" + (genomicPositionMapping[plusEnd] + 1);
                genomicReverseReadVector = (genomicPositionMapping[minusStart] + 1) + "-" + (genomicPositionMapping[minusEnd] + 1);
            } else {
                genomicForwardReadVector = (genomicPositionMapping[plusStart]) + "-" + (genomicPositionMapping[plusEnd]);
                genomicReverseReadVector = (genomicPositionMapping[minusStart]) + "-" + (genomicPositionMapping[minusEnd]);
            }

            String transcriptForwardReadVector = plusStart + "-" + plusEnd;
            String transcriptReverseReadVector = minusStart + "-" + minusEnd;

            List<Integer> mutatedPositionsForwardRead = roe.mutateInPlace(fwReadSequence, mutationRate);
            List<Integer> mutatedPositionsReverseRead = roe.mutateInPlace(rwReadSequence, mutationRate);

            events.add(new ReadGenerationEvent(
                    chromosome,
                    geneId,
                    transcriptId,
                    fwReadSequence,
                    rwReadSequence,
                    genomicForwardReadVector,
                    genomicReverseReadVector,
                    transcriptForwardReadVector,
                    transcriptReverseReadVector,
                    mutatedPositionsForwardRead,
                    mutatedPositionsReverseRead)
            );
        }
        return events;
    }
}