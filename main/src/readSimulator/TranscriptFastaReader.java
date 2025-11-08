package readSimulator;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class TranscriptFastaReader {

    /**
     * Reads a gzipped FASTA file and returns a map of transcript IDs to sequences.
     * Key = transcript ID (first token after '>'), Value = full sequence (no newlines)
     */
    public static HashMap<String, String> loadTranscripts(File fastaGzFile) throws IOException {
        HashMap<String, String> transcripts = new HashMap<>();

        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(
                        new GZIPInputStream(new FileInputStream(fastaGzFile)),
                        StandardCharsets.UTF_8
                )
        )) {
            String line;
            String currentId = null;
            StringBuilder seq = new StringBuilder(8000);

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    // save previous transcript
                    if (currentId != null) {
                        transcripts.put(currentId, seq.toString());
                    }

                    // extract transcript ID (first token after '>')
                    currentId = line.substring(1).split("\\s+")[0];
                    seq.setLength(0); // reset for next sequence
                } else {
                    seq.append(line.trim());
                }
            }

            // add last entry
            if (currentId != null) {
                transcripts.put(currentId, seq.toString());
            }
        }

        return transcripts;
    }

    public static HashMap<String, String> getTranscriptMap() {
        // hardcoded path to your FASTA file
        File fasta = new File("/home/valentinrexer/core/uni/bioinformatics/semester7/gobi/data/read_im_Test/Homo_sapiens.GRCh37.75.cdna.all.fa.gz");

        try {
            System.out.println("Reading transcripts from: " + fasta.getAbsolutePath());
            HashMap<String, String> transcripts = loadTranscripts(fasta);
            System.out.println("✅ Loaded " + transcripts.size() + " transcripts.");
            return transcripts;

        } catch (IOException e) {
            System.err.println("❌ Error reading FASTA: " + e.getMessage());
            e.printStackTrace();
        }
        return null;
    }
}
