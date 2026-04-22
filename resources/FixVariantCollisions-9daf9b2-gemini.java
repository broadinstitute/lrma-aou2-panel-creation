import java.io.*;
import java.util.*;

public class FixVariantCollisions {
    private static String WEIGHT_TAG;
    private static boolean WEIGHT_TAG_IN_SAMPLE_COLUMN;
    private static double DEFAULT_WEIGHT;
    private static boolean METHOD;

    enum VariantType { DEL, INV, DUP, INS, SNP, REPLACEMENT }

    enum Genotype {
        PHASED_00, PHASED_01, PHASED_10, PHASED_11,
        UNPHASED_00, UNPHASED_01, UNPHASED_10, UNPHASED_11,
        PHASED_D0, PHASED_0D, PHASED_D1, PHASED_1D, PHASED_DD,
        UNPHASED_D0, UNPHASED_0D, UNPHASED_D1, UNPHASED_1D, UNPHASED_DD;

        public String toStr() {
            switch (this) {
                case PHASED_00: return "0|0"; case PHASED_01: return "0|1";
                case PHASED_10: return "1|0"; case PHASED_11: return "1|1";
                case UNPHASED_00: return "0/0"; case UNPHASED_01: return "0/1";
                case UNPHASED_10: return "1/0"; case UNPHASED_11: return "1/1";
                case PHASED_D0: return ".|0"; case PHASED_0D: return "0|.";
                case PHASED_D1: return ".|1"; case PHASED_1D: return "1|."; case PHASED_DD: return ".|.";
                case UNPHASED_D0: return "./0"; case UNPHASED_0D: return "0/.";
                case UNPHASED_D1: return "./1"; case UNPHASED_1D: return "1/."; case UNPHASED_DD: return "./.";
                default: return "./.";
            }
        }
    }

    private static final int[][] N_GT_COLLISIONS = {
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0},
        {0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0},
        {0,1,1,2, 0,0,0,2, 0,0,1,1,0, 0,0,1,1,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,1,1,2, 0,0,0,2, 0,0,1,1,0, 0,0,1,1,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0},
        {0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0},
        {0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0}
    };

    private static final Genotype[] REMOVE_HAP1_0 = { Genotype.PHASED_00,Genotype.PHASED_01,Genotype.PHASED_00,Genotype.PHASED_01,Genotype.UNPHASED_00,Genotype.UNPHASED_01,Genotype.UNPHASED_10,Genotype.UNPHASED_01,Genotype.PHASED_D0,Genotype.PHASED_0D,Genotype.PHASED_D1,Genotype.PHASED_0D,Genotype.PHASED_DD,Genotype.UNPHASED_D0,Genotype.UNPHASED_0D,Genotype.UNPHASED_D1,Genotype.UNPHASED_1D,Genotype.UNPHASED_DD };
    private static final Genotype[] REMOVE_HAP1_D = { Genotype.PHASED_00,Genotype.PHASED_01,Genotype.PHASED_D0,Genotype.PHASED_D1,Genotype.UNPHASED_00,Genotype.UNPHASED_01,Genotype.UNPHASED_10,Genotype.UNPHASED_D1,Genotype.PHASED_D0,Genotype.PHASED_0D,Genotype.PHASED_D1,Genotype.PHASED_DD,Genotype.PHASED_DD,Genotype.UNPHASED_D0,Genotype.UNPHASED_0D,Genotype.UNPHASED_D1,Genotype.UNPHASED_1D,Genotype.UNPHASED_DD };
    private static final Genotype[] REMOVE_HAP2_0 = { Genotype.PHASED_00,Genotype.PHASED_00,Genotype.PHASED_10,Genotype.PHASED_10,Genotype.UNPHASED_00,Genotype.UNPHASED_01,Genotype.UNPHASED_10,Genotype.UNPHASED_01,Genotype.PHASED_D0,Genotype.PHASED_0D,Genotype.PHASED_D0,Genotype.PHASED_1D,Genotype.PHASED_DD,Genotype.UNPHASED_D0,Genotype.UNPHASED_0D,Genotype.UNPHASED_D1,Genotype.UNPHASED_1D,Genotype.UNPHASED_DD };
    private static final Genotype[] REMOVE_HAP2_D = { Genotype.PHASED_00,Genotype.PHASED_0D,Genotype.PHASED_10,Genotype.PHASED_1D,Genotype.UNPHASED_00,Genotype.UNPHASED_01,Genotype.UNPHASED_10,Genotype.UNPHASED_D1,Genotype.PHASED_D0,Genotype.PHASED_0D,Genotype.PHASED_DD,Genotype.PHASED_1D,Genotype.PHASED_DD,Genotype.UNPHASED_D0,Genotype.UNPHASED_0D,Genotype.UNPHASED_D1,Genotype.UNPHASED_1D,Genotype.UNPHASED_DD };

    private static List<Interval> window = new ArrayList<>();
    private static double maxWeight = 0;

    public static void main(String[] args) throws IOException {
        if (args.length < 4) {
            System.err.println("Usage: java FixVariantCollisions <method> <weight_tag> <weight_loc> <default_w> [<out.windows> <out.hist>] < input.vcf > output.vcf");
            System.exit(1);
        }

        METHOD = args[0].equalsIgnoreCase("1");
        WEIGHT_TAG = args[1];
        WEIGHT_TAG_IN_SAMPLE_COLUMN = Integer.parseInt(args[2]) == 1;
        DEFAULT_WEIGHT = Double.parseDouble(args[3]);
        
        final String OUTPUT_WINDOWS = args.length > 4 ? args[4] : "null";
        final String OUTPUT_HISTOGRAM = args.length > 5 ? args[5] : "null";

        long[] histogram = new long[100];
        int windowLastPos = -1;

        try (
            BufferedReader brVCF = new BufferedReader(new InputStreamReader(System.in));
            BufferedWriter bwVCF = new BufferedWriter(new OutputStreamWriter(System.out))
        ) {
            BufferedWriter bwWindows = OUTPUT_WINDOWS.equalsIgnoreCase("null") ? null : new BufferedWriter(new FileWriter(OUTPUT_WINDOWS));
            
            int nRecords = 0;
            String str;
            
            while ((str = brVCF.readLine()) != null) {
                if (str.startsWith("#")) {
                    bwVCF.write(str); bwVCF.newLine();
                    continue;
                }
                
                Interval tmpInterval = new Interval(str);
                window.add(tmpInterval);
                
                if (!isLastInWindow(windowLastPos)) {
                    Interval nextIv = window.remove(window.size() - 1); // Extract non-overlapping interval
                    
                    for (int i = 0; i < window.size(); i++) window.get(i).inputIndex = i;
                    fixWindow(bwVCF, bwWindows, histogram);
                    
                    nRecords += window.size();
                    if (nRecords % 10000 == 0) System.err.println("Processed " + nRecords + " records");
                    
                    window.clear();
                    window.add(nextIv);
                    windowLastPos = nextIv.last;
                } else {
                    if (tmpInterval.last > windowLastPos) windowLastPos = tmpInterval.last;
                }
            }

            if (!window.isEmpty()) {
                for (int i = 0; i < window.size(); i++) window.get(i).inputIndex = i;
                fixWindow(bwVCF, bwWindows, histogram);
            }

            if (bwWindows != null) bwWindows.close();
            bwVCF.flush();
        }

        if (!OUTPUT_HISTOGRAM.equalsIgnoreCase("null")) {
            try (BufferedWriter bwHistogram = new BufferedWriter(new FileWriter(OUTPUT_HISTOGRAM))) {
                bwHistogram.write("#nCollision \t nHaplotypes\n");
                for (int i = 0; i < histogram.length; i++) {
                    bwHistogram.write(i + "\t" + histogram[i] + "\n");
                }
            }
        }
    }

    private static boolean isLastInWindow(int windowLastPos) {
        if (window.isEmpty() || windowLastPos == -1) return true;
        Interval first = window.get(0);
        Interval last = window.get(window.size() - 1);
        return last.chr == first.chr && last.first <= windowLastPos;
    }

    private static void fixWindow(BufferedWriter bwVCF, BufferedWriter bwWindows, long[] histogram) throws IOException {
        if (window.isEmpty()) return;
        
        window.sort(Comparator.comparingInt(a -> a.last));
        int nSamples = window.get(0).genotypes.length;

        for (int j = 0; j < nSamples; j++) {
            int cols = countCollisions(j);
            int histIdx = Math.min(cols, histogram.length - 1);
            histogram[histIdx]++;
            
            if (cols > 0) {
                if (METHOD) independentSet2(j);
                else independentSet1(j);
            }
        }

        if (bwVCF != null) {
            window.sort(Comparator.comparingInt(a -> a.inputIndex));
            for (Interval iv : window) iv.toVCF(bwVCF);
        }
    }

    private static int countCollisions(int sample) {
        int out = 0;
        for (int i = window.size() - 1; i >= 0; i--) {
            if (!window.get(i).isPresent(sample)) continue;
            for (int j = i - 1; j >= 0; j--) {
                if (window.get(j).last < window.get(i).first) break;
                out += N_GT_COLLISIONS[window.get(j).genotypes[sample].ordinal()][window.get(i).genotypes[sample].ordinal()];
            }
        }
        return out;
    }

    private static void independentSet1(int sample) {
        List<Interval> sampleWindow = new ArrayList<>();
        for (Interval iv : window) {
            if (iv.isPresent(sample)) {
                iv.clearISVariables();
                iv.setWeight(sample);
                sampleWindow.add(iv);
            }
        }
        
        if (sampleWindow.isEmpty()) return;

        maxWeight = 0;
        boolean[] active = new boolean[sampleWindow.size()];
        Arrays.fill(active, true);
        List<Interval> bestIS = new ArrayList<>();

        for (int i = sampleWindow.size() - 1; i >= 0; i--) {
            is1DFS(i, active, new ArrayList<>(), 0.0, sample, sampleWindow, bestIS);
        }

        for (Interval iv : bestIS) iv.inIndependentSet = true;
        markISOverlaps(sampleWindow, sample, true, true);

        for (Interval iv : sampleWindow) {
            if (iv.inIndependentSet) continue;
            Genotype gt = iv.genotypes[sample];
            gt = iv.overlapsIS_hap1 ? REMOVE_HAP1_D[gt.ordinal()] : REMOVE_HAP1_0[gt.ordinal()];
            gt = iv.overlapsIS_hap2 ? REMOVE_HAP2_D[gt.ordinal()] : REMOVE_HAP2_0[gt.ordinal()];
            iv.genotypes[sample] = gt;
        }
    }

    private static void is1DFS(int id, boolean[] active, List<Interval> ancestors, double ancestorsWeight, int sample, List<Interval> sampleWindow, List<Interval> bestIS) {
        Interval current = sampleWindow.get(id);
        double weightPrime = ancestorsWeight + current.weight;
        
        boolean[] activePrime = active.clone();
        activePrime[id] = false;
        
        boolean hasChild = false;
        double upperBound = weightPrime;
        
        for (int i = id - 1; i >= 0; i--) {
            if (activePrime[i] && !sampleWindow.get(i).precedes(current, sample)) activePrime[i] = false;
            if (activePrime[i]) {
                hasChild = true;
                upperBound += sampleWindow.get(i).weight;
            }
        }
        
        if (upperBound < maxWeight) return;
        
        if (hasChild) {
            ancestors.add(current);
            for (int i = id - 1; i >= 0; i--) {
                if (!activePrime[i]) continue;
                if (upperBound < maxWeight) break;
                is1DFS(i, activePrime, ancestors, weightPrime, sample, sampleWindow, bestIS);
                upperBound -= sampleWindow.get(i).weight;
            }
            ancestors.remove(ancestors.size() - 1);
        } else if (weightPrime > maxWeight) {
            maxWeight = weightPrime;
            bestIS.clear();
            bestIS.addAll(ancestors);
            bestIS.add(current);
        }
    }

    private static class EndpointGroup {
        List<Interval> open = new ArrayList<>();
        List<Interval> closed = new ArrayList<>();
    }

    private static void independentSet2(int sample) {
        for (int hap = 1; hap <= 2; hap++) {
            List<Interval> sampleWindow = new ArrayList<>();
            TreeMap<Integer, EndpointGroup> pointsMap = new TreeMap<>();
            
            for (Interval iv : window) {
                boolean onHap = hap == 1 ? iv.onHap1(sample) : iv.onHap2(sample);
                if (onHap) {
                    iv.clearISVariables();
                    iv.setWeight(sample);
                    sampleWindow.add(iv);
                    pointsMap.computeIfAbsent(iv.first, k -> new EndpointGroup()).open.add(iv);
                    pointsMap.computeIfAbsent(iv.last, k -> new EndpointGroup()).closed.add(iv);
                }
            }

            if (sampleWindow.isEmpty()) continue;

            double availableWeight = 0.0;
            Interval previous = null;

            for (EndpointGroup group : pointsMap.values()) {
                for (Interval iv : group.open) {
                    iv.independentSetWeight = availableWeight + iv.weight;
                    iv.independentSetPrevious = previous;
                }
                for (Interval iv : group.closed) {
                    if (iv.independentSetWeight > availableWeight) {
                        availableWeight = iv.independentSetWeight;
                        previous = iv;
                    }
                }
            }

            Interval currTrace = null;
            for (int i = sampleWindow.size() - 1; i >= 0; i--) {
                if (Math.abs(sampleWindow.get(i).independentSetWeight - availableWeight) < 1e-9) {
                    currTrace = sampleWindow.get(i);
                    break;
                }
            }

            while (currTrace != null) {
                currTrace.inIndependentSet = true;
                currTrace = currTrace.independentSetPrevious;
            }

            if (hap == 1) {
                markISOverlaps(sampleWindow, sample, true, false);
                for (Interval iv : sampleWindow) {
                    if (!iv.inIndependentSet) {
                        Genotype gt = iv.genotypes[sample];
                        iv.genotypes[sample] = iv.overlapsIS_hap1 ? REMOVE_HAP1_D[gt.ordinal()] : REMOVE_HAP1_0[gt.ordinal()];
                    }
                }
            } else {
                markISOverlaps(sampleWindow, sample, false, true);
                for (Interval iv : sampleWindow) {
                    if (!iv.inIndependentSet) {
                        Genotype gt = iv.genotypes[sample];
                        iv.genotypes[sample] = iv.overlapsIS_hap2 ? REMOVE_HAP2_D[gt.ordinal()] : REMOVE_HAP2_0[gt.ordinal()];
                    }
                }
            }
        }
    }

    private static void markISOverlaps(List<Interval> sampleWindow, int sample, boolean hap1, boolean hap2) {
        for (Interval iv : sampleWindow) {
            iv.overlapsIS_hap1 = false;
            iv.overlapsIS_hap2 = false;
        }

        for (int i = sampleWindow.size() - 1; i >= 0; i--) {
            Interval ivI = sampleWindow.get(i);
            for (int j = i - 1; j >= 0; j--) {
                Interval ivJ = sampleWindow.get(j);
                if (ivJ.last < ivI.first) break;

                if (ivI.inIndependentSet && !ivJ.inIndependentSet) {
                    if (hap1 && ivI.onHap1(sample) && ivJ.onHap1(sample)) ivJ.overlapsIS_hap1 = true;
                    if (hap2 && ivI.onHap2(sample) && ivJ.onHap2(sample)) ivJ.overlapsIS_hap2 = true;
                } else if (ivJ.inIndependentSet && !ivI.inIndependentSet) {
                    if (hap1 && ivJ.onHap1(sample) && ivI.onHap1(sample)) ivI.overlapsIS_hap1 = true;
                    if (hap2 && ivJ.onHap2(sample) && ivI.onHap2(sample)) ivI.overlapsIS_hap2 = true;
                }
            }
        }
    }

    private static class Interval {
        VariantType variantType;
        int chr, first, last, inputIndex;
        double weight;
        String[] vcfRecord;
        Genotype[] genotypes;
        
        boolean inIndependentSet;
        double independentSetWeight;
        Interval independentSetPrevious;
        boolean overlapsIS_hap1, overlapsIS_hap2;

        public Interval(String record) {
            this.vcfRecord = record.split("\t", -1); 
            
            String svTypeStr = getInfoField(vcfRecord[7], "SVTYPE");
            if (svTypeStr != null) {
                variantType = parseSVType(svTypeStr);
            } else {
                variantType = vcfRecord[3].length() == 1 ? (vcfRecord[4].length() > 1 ? VariantType.INS : VariantType.SNP) 
                                                         : (vcfRecord[4].length() == 1 ? VariantType.DEL : VariantType.REPLACEMENT);
            }

            chr = parseChr(vcfRecord[0]);
            int pos = Integer.parseInt(vcfRecord[1]);
            
            String svLenStr = getInfoField(vcfRecord[7], "SVLEN");
            int length = svLenStr != null ? Math.abs(Integer.parseInt(svLenStr)) 
                                          : (variantType == VariantType.REPLACEMENT ? vcfRecord[3].length() - 1 
                                                                                    : Math.max(vcfRecord[3].length(), vcfRecord[4].length()) - 1);

            first = pos; last = pos;
            if (variantType == VariantType.DEL || variantType == VariantType.INV || variantType == VariantType.DUP || variantType == VariantType.REPLACEMENT) {
                first = pos + 1; last = pos + length;
            } else if (variantType == VariantType.INS) {
                last = pos + 1;
            }

            int nSamples = vcfRecord.length - 9;
            genotypes = new Genotype[nSamples];
            for (int i = 0; i < nSamples; i++) genotypes[i] = parseGT(vcfRecord[9 + i]);
        }

        public void setWeight(int sample) {
            String value = null;
            if (WEIGHT_TAG_IN_SAMPLE_COLUMN) {
                int p = vcfRecord[8].indexOf(WEIGHT_TAG);
                if (p >= 0) {
                    int j = 0;
                    for (int i = 0; i < p; i++) {
                        if (vcfRecord[8].charAt(i) == ':') j++;
                    }
                    String gt = vcfRecord[9 + sample];
                    int gtLength = gt.length();
                    for (int i = 0; i < gtLength; i++) {
                        if (gt.charAt(i) != ':') continue;
                        j--;
                        if (j > 0) continue;
                        p = gt.indexOf(':', i + 1);
                        value = p >= 0 ? gt.substring(i + 1, p) : gt.substring(i + 1);
                        break;
                    }
                }
            } else {
                value = getInfoField(vcfRecord[7], WEIGHT_TAG);
            }
            this.weight = value != null ? Double.parseDouble(value) : DEFAULT_WEIGHT;
        }

        public void clearISVariables() {
            inIndependentSet = false;
            independentSetWeight = 0.0;
            independentSetPrevious = null;
            overlapsIS_hap1 = false;
            overlapsIS_hap2 = false;
        }

        public boolean isPresent(int sample) {
            Genotype gt = genotypes[sample];
            return gt != Genotype.PHASED_00 && gt != Genotype.UNPHASED_00 && gt != Genotype.PHASED_DD && 
                   gt != Genotype.UNPHASED_DD && gt != Genotype.PHASED_D0 && gt != Genotype.UNPHASED_D0 && 
                   gt != Genotype.PHASED_0D && gt != Genotype.UNPHASED_0D;
        }

        public boolean onHap1(int sample) {
            Genotype gt = genotypes[sample];
            return gt == Genotype.PHASED_10 || gt == Genotype.PHASED_1D || gt == Genotype.PHASED_11 || gt == Genotype.UNPHASED_11;
        }

        public boolean onHap2(int sample) {
            Genotype gt = genotypes[sample];
            return gt == Genotype.PHASED_01 || gt == Genotype.PHASED_D1 || gt == Genotype.PHASED_11 || gt == Genotype.UNPHASED_11;
        }

        public boolean precedes(Interval next, int sample) {
            return N_GT_COLLISIONS[genotypes[sample].ordinal()][next.genotypes[sample].ordinal()] == 0 || last < next.first;
        }

        public void toVCF(BufferedWriter bw) throws IOException {
            bw.write(String.join("\t", Arrays.copyOfRange(vcfRecord, 0, 9)));
            
            String formatStr = vcfRecord[8];
            int p = formatStr.indexOf("GT");
            int gtIndex = 0;
            for (int i = 0; i < p; i++) {
                if (formatStr.charAt(i) == ':') gtIndex++;
            }
            
            for (int i = 0; i < genotypes.length; i++) {
                bw.write("\t");
                String gt = vcfRecord[9 + i];
                int gtLength = gt.length();
                int k = 0, p2 = 0;
                for (int j = 0; j < gtLength; j++) {
                    if (gt.charAt(j) == ':') {
                        k++;
                        if (k == gtIndex) { p2 = j + 1; break; }
                    }
                }
                if (p2 > 0) bw.write(gt.substring(0, p2));
                bw.write(genotypes[i].toStr());
                int q = gt.indexOf(':', p2);
                if (q >= 0) bw.write(gt.substring(q));
            }
            bw.newLine();
        }

        private static Genotype parseGT(String gt) {
            if (gt.length() < 3) return Genotype.UNPHASED_DD;
            char a = gt.charAt(0), b = gt.charAt(1), c = gt.charAt(2);
            if (b == '/') {
                if (a == '.') return c == '.' ? Genotype.UNPHASED_DD : (c == '0' ? Genotype.UNPHASED_D0 : Genotype.UNPHASED_D1);
                if (a == '0') return c == '.' ? Genotype.UNPHASED_0D : (c == '0' ? Genotype.UNPHASED_00 : Genotype.UNPHASED_01);
                return c == '.' ? Genotype.UNPHASED_1D : (c == '0' ? Genotype.UNPHASED_10 : Genotype.UNPHASED_11);
            } else {
                if (a == '.') return c == '.' ? Genotype.PHASED_DD : (c == '0' ? Genotype.PHASED_D0 : Genotype.PHASED_D1);
                if (a == '0') return c == '.' ? Genotype.PHASED_0D : (c == '0' ? Genotype.PHASED_00 : Genotype.PHASED_01);
                return c == '.' ? Genotype.PHASED_1D : (c == '0' ? Genotype.PHASED_10 : Genotype.PHASED_11);
            }
        }

        private static String getInfoField(String info, String field) {
            String target = field + "=";
            int p = info.indexOf(target);
            if (p < 0) return null;
            int valStart = p + target.length();
            int valEnd = info.indexOf(";", valStart);
            return valEnd < 0 ? info.substring(valStart) : info.substring(valStart, valEnd);
        }

        private static int parseChr(String chr) {
            String s = chr.toUpperCase().replace("CHR", "");
            switch(s) {
                case "X": return 23; case "Y": return 24; case "M": case "MT": return 25;
                default: try { return Integer.parseInt(s); } catch(NumberFormatException e) { return -1; }
            }
        }

        private static VariantType parseSVType(String type) {
            String t = type.toUpperCase();
            if (t.equals("DEL") || t.equals("DEL:ME")) return VariantType.DEL;
            if (t.equals("INV")) return VariantType.INV;
            if (t.equals("DUP") || t.equals("DUP:TANDEM") || t.equals("DUP:INT") || t.equals("CNV")) return VariantType.DUP;
            if (t.equals("INS") || t.equals("INS:ME") || t.equals("INS:NOVEL")) return VariantType.INS;
            return VariantType.REPLACEMENT;
        }
    }
}
