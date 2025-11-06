package readSimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Gtf {
    private List<Gene> genes;

    public Gtf(Path filePath, List<String> relevantTranscriptIds) {
        HashMap<String, Gene> geneMap = new HashMap<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath.toFile()));
            String line;

            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                if (!fields[2].equals("exon")) continue;

                String transcriptId = Gtf.getAttribute("transcript_id", fields[8]);
                if (!relevantTranscriptIds.contains(transcriptId)) continue;

                String geneId = Gtf.getAttribute("gene_id", fields[8]);

                if (!geneMap.containsKey(geneId)) {
                    Gene gene = new Gene(geneId, fields[0]);
                    Transcript transcript = new Transcript(transcriptId);
                    transcript.addCoordinates(new Coordinates(Integer.parseInt(fields[3]), Integer.parseInt(fields[4])));
                    gene.addTranscript(transcript);
                    geneMap.put(geneId, gene);
                } else {
                    Gene currentGene = geneMap.get(geneId);
                    if (!currentGene.hasTranscript(transcriptId)) {
                        Transcript transcript = new Transcript(transcriptId);
                        transcript.addCoordinates(new Coordinates(Integer.parseInt(fields[3]), Integer.parseInt(fields[4])));
                    } else {
                        Transcript transcript = currentGene.getTranscript(transcriptId);
                        transcript.addCoordinates(new Coordinates(Integer.parseInt(fields[3]), Integer.parseInt(fields[4])));
                    }
                }
            }

            genes = new ArrayList<>(geneMap.values());

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public List<Gene> getGenes() {
        return genes;
    }

    public static String getAttribute(String attribute, String column) {
        Pattern pattern = Pattern.compile(attribute + "\\s+\"([^\"]+)\"");
        Matcher matcher = pattern.matcher(column);
        if (matcher.find()) {
            return matcher.group(1);
        } else {
            System.out.println("Wasn't able to extract '" + attribute + " from string '"+ column +"'");
            return null;
        }
    }
}
