package readSimulator;

import java.util.ArrayList;
import java.util.List;
import java.util.SplittableRandom;

public class RandomOperationExecutor {
    private final SplittableRandom rng = new SplittableRandom();
    private static final byte[] BASES = {'A','C','G','T'};

    private double nextGaussian() {
        double u1 = rng.nextDouble();
        double u2 = rng.nextDouble();
        return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    }

    public void initRandomSamples(int n, double meanLength, double standardDeviation, int transcriptLength, int readLength,
                                  int[] fragmentsLengths, int[] startPositions) {
        for (int i = 0; i < n; i++) {
            int fragmentLength;
            do {
                double val = meanLength + standardDeviation * nextGaussian();
                fragmentLength = (int) Math.round(val);
            } while (fragmentLength <= readLength || fragmentLength >= transcriptLength);

            fragmentsLengths[i] = fragmentLength;

            int maxStartPosition = transcriptLength - fragmentLength;
            startPositions[i] = rng.nextInt(maxStartPosition);
        }
    }

    public List<Integer> mutateInPlace(byte[] seq, double mutationRate) {
        if (seq == null || mutationRate <= 0) return null;
        List<Integer> ret = new ArrayList<>();

        if (mutationRate >= 100.0) {
            for (int i = 0; i < seq.length; i++){
                seq[i] = randomDifferentBase(seq[i]);
                ret.add(i);
            }
            return ret ;
        }

        final double log1mP = Math.log(1.0 - mutationRate / 100.0);
        int i = 0;
        final int L = seq.length;

        while (i < L) {
            // geometric skip
            double u = rng.nextDouble();
            int skip = (int) (Math.log(u) / log1mP);
            i += skip;
            if (i >= L) break;
            seq[i] = randomDifferentBase(seq[i]);
            ret.add(i);
            i++;  // move past the mutated site
        }

        return ret;
    }

    private byte randomDifferentBase(byte b) {
        byte old = (byte) Character.toUpperCase((char) b);
        byte newBase;
        do {
            newBase = BASES[rng.nextInt(4)];
        } while (newBase == old);
        return newBase;
    }
}
