package readSimulator;

import org.apache.commons.cli.*;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

public class ReadSimulator {

    private final ReadCounts readCounts;
    private final Gtf gtf;
    private final IndexedFastaReader reader;

    private final int readLength;
    private final int fragmentLength;
    private final double fragmentSD;
    private final double mutationRate;
    private final Path outputDir;
    private final int CHUNK_SIZE = 30_000;

    public ReadSimulator(
            Path readCountsPath,
            Path fastaPath,
            Path fidxPath,
            Path gtfPath,
            int readLength,
            int fragmentLength,
            double fragmentSD,
            double mutationRate,
            Path outputDir
    ) throws IOException {

        this.readCounts = new ReadCounts(readCountsPath);
        this.reader = new IndexedFastaReader(fidxPath, fastaPath);
        this.gtf = new Gtf(gtfPath, readCounts.getTranscriptIds());
        this.readLength = readLength;
        this.fragmentLength = fragmentLength;
        this.fragmentSD = fragmentSD;
        this.mutationRate = mutationRate;
        this.outputDir = outputDir;

        if (!Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        this.reader.openChannel();
    }

    public void runSimulation() throws IOException {
        System.out.printf(
                "▶ Starting simulation with parameters:%n" +
                        "   Read length: %d%n" +
                        "   Fragment length (mean): %d%n" +
                        "   Fragment SD: %.2f%n" +
                        "   Mutation rate: %.4f%n" +
                        "   Output directory: %s%n",
                readLength, fragmentLength, fragmentSD, mutationRate, outputDir.toAbsolutePath()
        );

        gtf.buildGenesTranscriptSequences(reader);
        int QUEUE_CAPACITY = 500;
        BlockingQueue<ReadGenerationEventChunk> queue  = new ArrayBlockingQueue<>(QUEUE_CAPACITY);

        Path fwPath = outputDir.resolve("fw.fastq");
        Path rwPath = outputDir.resolve("rw.fastq");
        Path mappingInfoPath = outputDir.resolve("read.mappinginfo");

        ParallelizedOutputWriter writer = new ParallelizedOutputWriter(queue, fwPath, rwPath, mappingInfoPath, readLength);
        Thread writerThread = new Thread(writer);
        writerThread.start();

        for (Gene gene : gtf.getGenes()) {
            try {
                gene.generateEventsForAllTranscripts(
                        readCounts.getCounts().get(gene.getGeneId()),
                        fragmentLength,
                        fragmentSD,
                        readLength,
                        mutationRate,
                        queue,
                        CHUNK_SIZE
                );
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }

        /*
        gtf.getGenes().parallelStream().forEach(gene -> {
            try {
                gene.generateEventsForAllTranscripts(
                        readCounts.getCounts().get(gene.getGeneId()),
                        fragmentLength,
                        fragmentSD,
                        readLength,
                        mutationRate,
                        queue,
                        CHUNK_SIZE
                );
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        });
         */

        try {
            queue.put(ParallelizedOutputWriter.STOPPING_SIGNAL_CHUNK);
            writerThread.join();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    public static void main(String[] args) {
        Options options = new Options();

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

        options.addOption(Option.builder("length")
                .hasArg()
                .argName("int")
                .desc("Read length in bases")
                .required()
                .build());

        options.addOption(Option.builder("frlength")
                .hasArg()
                .argName("int")
                .desc("Mean fragment length")
                .required()
                .build());

        options.addOption(Option.builder("SD")
                .hasArg()
                .argName("double")
                .desc("Standard deviation of fragment length")
                .required()
                .build());

        options.addOption(Option.builder("mutationrate")
                .hasArg()
                .argName("double")
                .desc("Mutation rate (0–1)")
                .required()
                .build());

        options.addOption(Option.builder("od")
                .hasArg()
                .argName("path")
                .desc("Output directory")
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
            Path outputDir      = Paths.get(cmd.getOptionValue("od"));

            int readLength      = Integer.parseInt(cmd.getOptionValue("length"));
            int fragmentLength  = Integer.parseInt(cmd.getOptionValue("frlength"));
            double fragmentSD   = Double.parseDouble(cmd.getOptionValue("SD"));
            double mutationRate = Double.parseDouble(cmd.getOptionValue("mutationrate"));

            if (mutationRate < 0 || mutationRate > 1)
                throw new IllegalArgumentException("Mutation rate must be between 0 and 1");

            ReadSimulator simulator = new ReadSimulator(
                    readCountsPath, fastaPath, fidxPath, gtfPath,
                    readLength, fragmentLength, fragmentSD, mutationRate, outputDir
            );
            simulator.runSimulation();

        } catch (ParseException e) {
            System.err.println("❌ Argument parsing error: " + e.getMessage());
            formatter.printHelp("java readSimulator.ReadSimulator", options, true);
            System.exit(1);

        } catch (IllegalArgumentException e) {
            System.err.println("❌ Invalid argument: " + e.getMessage());
            System.exit(2);

        } catch (IOException e) {
            System.err.println("❌ I/O error: " + e.getMessage());
            e.printStackTrace();
            System.exit(3);
        }
    }
}
