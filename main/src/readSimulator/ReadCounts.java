package readSimulator;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.nio.file.Path;
import java.util.logging.Logger;
import java.util.logging.Level;

public class ReadCounts {
    private final List<List<String>> countsInfo;
    private final Logger _LOGGER = Logger.getLogger(ReadCounts.class.getName());

    public ReadCounts(Path filePath) {
        countsInfo = new ArrayList<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath.toFile()));
            String line;

            line = br.readLine();
            if (line.equals("gene\ttranscript\tcount")) line = br.readLine();

            while ((line = br.readLine())!= null) {
                List<String> row = new ArrayList<>(Arrays.asList(line.split("\t")));
                countsInfo.add(row);
            }
        } catch (FileNotFoundException e) {
            _LOGGER.log(Level.SEVERE, "File not found: {0}", e.getMessage());
        } catch (IOException e) {
            _LOGGER.log(Level.SEVERE, "Error reading file: {0}", e.getMessage());
        }
    }

    public List<List<String>> getCounts() {
        return countsInfo;
    }

    public List<String> getTranscriptIds() {
        List<String> transcriptIds = new ArrayList<>();

        for  (List<String> row : countsInfo)
            transcriptIds.add(row.get(1));

        return transcriptIds;
    }

    public List<String> getGeneIds(String transcriptId) {
        List<String> geneIds = new ArrayList<>();

        for (List<String> row : countsInfo)
            geneIds.add(row.getFirst());

        return geneIds;
    }

    public void pop(int index) {
        this.countsInfo.remove(index);
    }

    public void pop(String geneId, String transcriptId, String count) {
        List<Integer> toPop = new ArrayList<>();

        for (int i = 0; i < countsInfo.size(); i++) {
            if (!countsInfo.get(i).contains(geneId)) continue;
            if (!countsInfo.get(i).contains(transcriptId)) continue;
            if (!countsInfo.get(i).contains(count)) continue;
            toPop.add(i);
        }

        for (int i = toPop.size() - 1; i >= 0; i--)
            countsInfo.remove(toPop.get(i).intValue());
    }
}
