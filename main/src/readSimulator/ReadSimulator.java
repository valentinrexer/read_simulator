package readSimulator;

import org.apache.commons.cli.*;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;

public class ReadSimulator {

    private final ReadCounts readCounts;
    private final Gtf gtf;
    private final IndexedFastaReader reader;

    public ReadSimulator(Path readCountsPath, Path fastaPath, Path fidxPath, Path gtfPath) throws IOException {
        this.readCounts = new ReadCounts(readCountsPath);
        this.reader = new IndexedFastaReader(fidxPath, fastaPath);
        this.gtf = new Gtf(gtfPath, readCounts.getTranscriptIds());
        this.reader.openChannel();
    }

    public void runSimulation() throws IOException {
        gtf.buildGenesTranscriptSequences(reader);
    }

    public static void main(String[] args) {
        Options options = new Options();

        // Define CLI options
        options.addOption(Option.builder("readcounts")
                .hasArg()
                .argName("path")
                .desc("Path to readcounts.simulation file")
                .required()
                .build());

        options.addOption(Option.builder("fasta")
                .hasArg()
                .argName("path")
                .desc("Path to FASTA file")
                .required()
                .build());

        options.addOption(Option.builder("fidx")
                .hasArg()
                .argName("path")
                .desc("Path to FASTA index (.fai) file")
                .required()
                .build());

        options.addOption(Option.builder("gtf")
                .hasArg()
                .argName("path")
                .desc("Path to GTF annotation file")
                .required()
                .build());

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();

        try {
            CommandLine cmd = parser.parse(options, args);

            Path readCountsPath = Paths.get(cmd.getOptionValue("readcounts"));
            Path fastaPath      = Paths.get(cmd.getOptionValue("fasta"));
            Path fidxPath       = Paths.get(cmd.getOptionValue("fidx"));
            Path gtfPath        = Paths.get(cmd.getOptionValue("gtf"));

            ReadSimulator simulator = new ReadSimulator(readCountsPath, fastaPath, fidxPath, gtfPath);
            simulator.runSimulation();

        } catch (ParseException e) {
            System.err.println("❌ Argument parsing error: " + e.getMessage());
            formatter.printHelp("java readSimulator.ReadSimulator", options, true);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("❌ I/O error: " + e.getMessage());
            e.printStackTrace();
            System.exit(2);
        }
    }
}
