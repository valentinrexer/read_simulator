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
            case 'T', 'U' -> (byte) 'A';
            default -> (byte) 'N';
        };
    }

    public int length() {
        return sequence.length;
    }

    public List<ReadGenerationEvent> createEventsForTranscript(int[] fragmentLengths,
                                                               int[] startingPositions,
                                                               int readLength,
                                                               double mutationRate,
                                                               RandomOperationExecutor roe,
                                                               String chromosome,
                                                               String geneId) {
        List<ReadGenerationEvent> events = new ArrayList<>();

        for (int i = 0; i < fragmentLengths.length; i++) {
            int plusStrandFirstIndex = startingPositions[i];
            int plusStrandLastIndex = startingPositions[i] + readLength;

            int minusStrandFirstIndex = startingPositions[i] + fragmentLengths[i] - readLength;
            int minusStrandLastIndex = startingPositions[i] + fragmentLengths[i];

            byte[] fwReadSequence = Arrays.copyOfRange(sequence,
                    plusStrandFirstIndex,
                    plusStrandLastIndex);

            byte[] rwReadSequence = Arrays.copyOfRange(sequence,
                    minusStrandFirstIndex,
                    minusStrandLastIndex);

            reverseComplementInPlace(rwReadSequence);

            String transcriptForwardReadVector = plusStrandFirstIndex + "-" + plusStrandLastIndex;
            String transcriptReverseReadVector = minusStrandFirstIndex + "-" + minusStrandLastIndex;

            List<String> genomicCoordinates;
            if (strand == '+') genomicCoordinates = getGenomicCoordinates(plusStrandFirstIndex, plusStrandLastIndex, minusStrandFirstIndex, minusStrandLastIndex);
            else genomicCoordinates = getGenomicCoordinates(minusStrandFirstIndex, minusStrandLastIndex, plusStrandFirstIndex, plusStrandLastIndex);

            String genomicForwardReadVector = genomicCoordinates.get(0);
            String genomicReverseReadVector = genomicCoordinates.get(1);


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

    private List<String> getGenomicCoordinates(int fwFirst,
                                               int fwLast,
                                               int rwFirst,
                                               int rwLast) {

        // Convert to inclusive last indices in transcript coordinates
        int fwLastInclusive = fwLast - 1;
        int rwLastInclusive = rwLast - 1;

        List<Coordinates> fwRegions = new ArrayList<>();
        List<Coordinates> rwRegions = new ArrayList<>();

        if (strand == '+') {
            int transcriptPos = 0; // start of current exon in transcript coordinates

            for (Coordinates exon : exonRegions) {

                int exonStart = exon.coordinate1();
                int exonEnd = exon.coordinate2();
                int exonLen = exonEnd - exonStart + 1;

                int exonT0 = transcriptPos;
                int exonT1 = transcriptPos + exonLen - 1;

                if (fwLastInclusive >= exonT0 && fwFirst <= exonT1) {

                    int overlapStartT = Math.max(fwFirst, exonT0);
                    int overlapEndT = Math.min(fwLastInclusive, exonT1);

                    int genomicStart = exonStart + (overlapStartT - exonT0);
                    int genomicEnd = exonStart + (overlapEndT - exonT0);

                    fwRegions.add(new Coordinates(genomicStart, genomicEnd));
                }

                if (rwLastInclusive >= exonT0 && rwFirst <= exonT1) {

                    int overlapStartT = Math.max(rwFirst, exonT0);
                    int overlapEndT = Math.min(rwLastInclusive, exonT1);

                    int genomicStart = exonStart + (overlapStartT - exonT0);
                    int genomicEnd = exonStart + (overlapEndT - exonT0);

                    rwRegions.add(new Coordinates(genomicStart, genomicEnd));
                }

                transcriptPos += exonLen;
            }
        } else {
            // Transcript sequence is reverse-complemented, so transcript coordinate 0
            // corresponds to the last base of the last exon.
            int transcriptPos = 0;

            for (int i = exonRegions.size() - 1; i >= 0; i--) {
                Coordinates exon = exonRegions.get(i);

                int exonStart = exon.coordinate1();
                int exonEnd = exon.coordinate2();
                int exonLen = exonEnd - exonStart + 1;

                int exonT0 = transcriptPos;
                int exonT1 = transcriptPos + exonLen - 1;

                if (fwLastInclusive >= exonT0 && fwFirst <= exonT1) {
                    int overlapStartT = Math.max(fwFirst, exonT0);
                    int overlapEndT = Math.min(fwLastInclusive, exonT1);

                    int genomicStart = exonEnd - (overlapStartT - exonT0);
                    int genomicEnd = exonEnd - (overlapEndT - exonT0);

                    fwRegions.add(new Coordinates(Math.min(genomicStart, genomicEnd), Math.max(genomicStart, genomicEnd)));
                }

                if (rwLastInclusive >= exonT0 && rwFirst <= exonT1) {
                    int overlapStartT = Math.max(rwFirst, exonT0);
                    int overlapEndT = Math.min(rwLastInclusive, exonT1);

                    int genomicStart = exonEnd - (overlapStartT - exonT0);
                    int genomicEnd = exonEnd - (overlapEndT - exonT0);

                    rwRegions.add(new Coordinates(Math.min(genomicStart, genomicEnd), Math.max(genomicStart, genomicEnd)));
                }

                transcriptPos += exonLen;
            }
        }

        String fwString = coordsToString(fwRegions);
        String rwString = coordsToString(rwRegions);

        return List.of(fwString, rwString);
    }

    private static String coordsToString(List<Coordinates> list) {
        if (list.isEmpty()) return "";
        StringBuilder sb = new StringBuilder();
        for (Coordinates c : list) {
            sb.append(c.coordinate1())
                    .append("-")
                    .append(c.coordinate2())
                    .append("|");
        }
        sb.setLength(sb.length() - 1); // remove final "|"
        return sb.toString();
    }


}