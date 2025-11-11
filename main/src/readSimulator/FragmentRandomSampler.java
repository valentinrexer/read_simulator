package readSimulator;

import java.util.SplittableRandom;

public class FragmentRandomSampler {
    private final SplittableRandom random = new SplittableRandom();

    private double nextGaussian() {
        double u1 = random.nextDouble();
        double u2 = random.nextDouble();
        return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    }

    public void initRandomSamples(int n, double meanLength, double standardDeviation, int transcriptLength,
                                  int[] fragmentsLengths, int[] startPositions) {
        for (int i = 0; i < n; i++) {
            int fragmentLength;
            do {
                double val = meanLength + standardDeviation * nextGaussian();
                fragmentLength = (int) Math.round(val);
            } while (fragmentLength < 1 || fragmentLength >= transcriptLength);

            fragmentsLengths[startPositions[i]] = fragmentLength;

            int maxStartPosition = transcriptLength - fragmentLength;
            startPositions[i] = random.nextInt(maxStartPosition);
        }
    }
}
