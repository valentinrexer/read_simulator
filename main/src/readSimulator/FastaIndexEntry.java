package readSimulator;

public record FastaIndexEntry(String seqName, long seqLength, long offset, int lineBases, int lineWidth) {

    @Override
    public String toString() {
        return seqName + "\t" + seqLength + "\t" + offset + "\t" + lineBases + "\t" + lineWidth;
    }
}

