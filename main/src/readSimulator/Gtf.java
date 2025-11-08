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
    private final HashMap<String, Gene> genes;

    public Gtf(Path filePath, List<String> relevantTranscriptIds) {
        genes = new HashMap<>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(filePath.toFile()));
            String line;

            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] fields = line.split("\t");
                if (!fields[2].equals("exon")) continue;

                String transcriptId = Gtf.getAttribute("transcript_id", fields[8]);
                if (!relevantTranscriptIds.contains(transcriptId)) continue;

                String geneId = Gtf.getAttribute("gene_id", fields[8]);

                if (!genes.containsKey(geneId)) {
                    Gene gene = new Gene(geneId, fields[0]);
                    gene.makeTranscript(transcriptId);
                    gene.getTranscript(transcriptId).addCoordinates(new Coordinates(Integer.parseInt(fields[3]), Integer.parseInt(fields[4])));
                    genes.put(geneId, gene);
                } else {
                    Gene currentGene = genes.get(geneId);
                    if (!currentGene.hasTranscript(transcriptId)) {
                        currentGene.makeTranscript(transcriptId);
                        currentGene.getTranscript(transcriptId).addCoordinates(new Coordinates(Integer.parseInt(fields[3]), Integer.parseInt(fields[4])));
                    } else {
                        currentGene.getTranscript(transcriptId).addCoordinates(new Coordinates(Integer.parseInt(fields[3]), Integer.parseInt(fields[4])));
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public Gene getGene(String geneId) {
        return genes.get(geneId);
    }

    public List<Gene> getGenes() {
        return new ArrayList<>(genes.values());
    }

    public void buildGenesTranscriptSequences(IndexedFastaReader reader) {
        for (Gene gene : genes.values()) {
            gene.buildTranscriptSequences(reader);
        }
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
