package readSimulator;

import java.io.*;
import java.nio.file.Path;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

public class ReadCounts {
    private final HashMap<String, HashMap<String, Integer>> countsInfo;
    private static final Logger _LOGGER = Logger.getLogger(ReadCounts.class.getName());

    public ReadCounts(Path filePath) {
        countsInfo = new HashMap<>();

        try (BufferedReader br = new BufferedReader(new FileReader(filePath.toFile()))) {
            String line;
            while ((line = br.readLine()) != null) {
                // Skip header
                if (line.equals("gene\ttranscript\tcount") || line.startsWith("gene")) continue;

                String[] parts = line.split("\t");
                if (parts.length < 3) continue;

                String geneId = parts[0];
                String transcriptId = parts[1];
                int count = Integer.parseInt(parts[2]);

                countsInfo
                        .computeIfAbsent(geneId, k -> new HashMap<>())
                        .put(transcriptId, count);
            }
        } catch (IOException e) {
            _LOGGER.log(Level.SEVERE, "Error reading file: {0}", e.getMessage());
        } catch (NumberFormatException e) {
            _LOGGER.log(Level.WARNING, "Invalid number format in file: {0}", e.getMessage());
        }
    }

    public HashMap<String, HashMap<String, Integer>> getCounts() {
        return countsInfo;
    }

    public List<String> getTranscriptIds() {
        List<String> transcriptIds = new ArrayList<>();
        for (HashMap<String, Integer> transcripts : countsInfo.values())
            transcriptIds.addAll(transcripts.keySet());
        return transcriptIds;
    }

    public List<String> getGeneIds() {
        return new ArrayList<>(countsInfo.keySet());
    }
}
