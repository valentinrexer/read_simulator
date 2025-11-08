package readSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.nio.file.Path;

public class FastaIndex {
    private final HashMap<String, FastaIndexEntry> indexEntries = new HashMap<>();

    public FastaIndex(Path faiPath) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(faiPath.toString()));
            String line;

            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                indexEntries.put(fields[0], new FastaIndexEntry(
                        fields[0],
                        Long.parseLong(fields[1]),
                        Long.parseLong(fields[2]),
                        Integer.parseInt(fields[3]),
                        Integer.parseInt(fields[4])
                        )
                );
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public FastaIndexEntry get(String seqName) {
        return indexEntries.get(seqName);
    }
}
