
package genotypecalling;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.TreeSet;

/**
 *
 * @author dashazhernakova
 */
public class SNVMixToGenConverter {
    String[] files;
    String[] samples;
    private HashMap<String, Integer> snpMap;
    private byte[][] alleles;
    private TreeSet<String> snpSet;
    private int numInds;
    private int numSNPs;
    byte[][] matrix;
    private float p_threshold = 0.8f;

    static final String[] chromosomes = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"};

    public void setPthreshold(float p){
        if ((p >= 0) && (p <= 1)){
            p_threshold = p;
        }

    }

    /**
     * Parses files whose names match snvmixFnamePattern (regex), creates output directory
     *
     * @param dir
     * @param outPrefix
     * @param snvmixFnamePattern regexp pattern for matching genotype files
     * @throws IOException
     */
    public void parse(String dir, String outPrefix, String snvmixFnamePattern, String conversionVCFFile) throws IOException {

        System.out.println("\nInput dir: " + dir);
        System.out.println("Output prefix: " + outPrefix);
        System.out.println("Pattern: " + snvmixFnamePattern + "\n");

        if (!Gpio.exists(dir)) {
            throw new IOException("ERROR: could not find dir: " + dir);
        }

        getFnames(dir, snvmixFnamePattern);

        System.out.println("Found " + files.length + " SNVMix genotype files:");

        if (files.length == 0) {
            System.out.println("ERROR: No samples detected!");
            System.exit(-1);
        }
        for (String f : files)
            System.out.println(f);

        numInds = files.length;

        try {
            writeIndividuals(outPrefix);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: Failed to write individuals");
        }

        convert(outPrefix);

        System.out.println("\nFinished successfully");

    }

    public void parse(String fileListPath, String outPrefix, String conversionVCFFile) {
        try {
            getFnames(fileListPath);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("\nERROR: Failed to get file names. Exiting!");
            System.exit(-1);
        }

        System.out.println("\nFile with samples and file paths: " + fileListPath);
        System.out.println("Output prefix: " + outPrefix);

        System.out.println("Found " + files.length + " SNVMix genotype files:");

        if (files.length == 0) {
            System.out.println("ERROR: No samples detected!");
            System.exit(-1);
        }
        for (String f : files)
            System.out.println(f);


        numInds = files.length;

        try {
            writeIndividuals(outPrefix);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: Failed to write individuals");
        }

        convert(outPrefix);

        System.out.println("\nFinished successfully");

    }
    public void convert(String outPrefix){
        System.out.println("\n\nStarted conversion");
        long st = System.currentTimeMillis();
        for (String chr : chromosomes) {
            System.out.println("\n---------------- chr" + chr + " ----------------");
            try {
                getAllSNPs(chr);
                makeMatrix(chr);
                //makeMatrix(chr);
            } catch (IOException e) {
                e.printStackTrace();
                System.out.println("\nERROR: Failed to make genotype matrix. Exiting!");
                System.exit(-1);
            }

            try {
                writeGenotypes(outPrefix, chr);
            } catch (IOException e) {
                e.printStackTrace();
                System.out.println("\nERROR: Failed to write the genotypes. Exiting!");
                System.exit(-1);
            }
        }
        long end = System.currentTimeMillis();
        System.out.println("\n--------------------------------");
        System.out.println("Finished the conversion in: " + (int)((end-st) / 1000) + " sec");
    }

    private void getAllSNPs(String chr) throws IOException {
        System.out.println("\nGetting all SNPs for chr" + chr);

        long st = System.currentTimeMillis();
        snpSet = new TreeSet<String>(COMPARE_BY_POS);
        snpMap = new HashMap<String, Integer>();

        String[] els;
        int gen, snpCtr = 0;
        float probability;
        String chrom, snpId;
        boolean chrRead;
        TextFile tf;

        for (String f : files){
            //System.out.println(f);
            tf = new TextFile(f, TextFile.R);
            chrRead = false;

            while ((els = tf.readLineElems(TextFile.tab)) != null){
                snpId = els[0].replaceFirst("chr", "");
                chrom = snpId.split(":")[0];

                if (chrom.equals(chr)){
                    chrRead = true;
                    gen = Integer.parseInt(els[3].split(",")[5]);
                    probability = Float.parseFloat(els[3].split(",")[1 + gen]);

                    if (probability >= p_threshold){
                        if (!snpMap.containsKey(snpId)){
                            snpMap.put(snpId, snpCtr);
                            snpCtr++;
                        }
                    }

                }
                else if (chrRead){
                    break;
                }
            }
        }

        snpSet.addAll(snpMap.keySet());
        numSNPs = snpSet.size();

        System.out.println("Loaded " + numSNPs + " SNPs");
        long end = System.currentTimeMillis();
        System.out.println("Finished in: " + (int)((end-st) / 1000) + " sec");
    }

    private void makeMatrix(String chr) throws IOException{
        System.out.println("\nGenerating genotype matrix for chr" + chr);

        long st = System.currentTimeMillis();

        //initialize genotype matrix
        matrix = new byte[numSNPs][numInds];
        //alleles = new HashMap<Integer, byte[]>();
        alleles = new byte[numSNPs][2];
        int indCnt = 0;

        SNVMix snvmix = new SNVMix();

        String [] els;
        String snpId, chrom;
        boolean chrRead;
        int gen;
        float probability;
        float[] genProbs;
        byte ref, alt;
        int snpIntId;

        for (String f : files){
            //System.out.println(f);
            TextFile tf = new TextFile(f, TextFile.R);
            chrRead = false;

            while ((els = tf.readLineElems(TextFile.tab)) != null){
                snpId = els[0].replaceFirst("chr", "");
                chrom = snpId.split(":")[0];

                if (chrom.equals(chr)){
                    chrRead = true;
                    gen = Integer.parseInt(els[3].split(",")[5]);
                    probability = Float.parseFloat(els[3].split(",")[1 + gen]);

                    if (probability >= p_threshold){
                        snpIntId = snpMap.get(snpId);

                        //ref = els[1].getBytes()[0];
                        //alt = els[2].getBytes()[0];
                        //alleles.put(snpIntId, new byte[] {ref, alt});

                        //genProbs = new float[] {Float.parseFloat(els[3].split(",")[2]), Float.parseFloat(els[3].split(",")[3]), Float.parseFloat(els[3].split(",")[4])};;
                        byte[] genotype = snvmix.getByteGenotype(els);
                        //matrix[snpId][0][indCnt] = genotype[0];
                        //matrix[snpId][1][indCnt] = genotype[1];

                        //code the genotype as 0 if ref-ref, 1 if ref-alt, 2 if alt-alt; put ref allele and alt allele before that
                        alleles [snpIntId][0] = els[1].getBytes()[0];
                        alleles [snpIntId][1] = els[2].getBytes()[0];
                        matrix [snpIntId][indCnt] = byteGenotype(els);
                    }
                }
                else if (chrRead){
                    break;
                }
            }
            indCnt++;
            tf.close();
        }

        long end = System.currentTimeMillis();
        System.out.println("Finished in: " + (int)((end-st) / 1000) + " sec");

    }

    /**
     * Gets paths to files if a list of samples and corresponding file paths is given
     * @param fileList - path to the file with sample names in the first column and paths to snvmix files in the second column (tab as a separator)
     * @throws IOException
     */
    private void getFnames(String fileList) throws IOException {
        TextFile flist = new TextFile(fileList, false);
        files = flist.readAsArray(1, TextFile.tab);
        samples = new String[files.length];
        flist.close();
        flist = new TextFile(fileList, false);
        samples = flist.readAsArray(0, TextFile.tab);
        flist.close();
    }


    /**
     * Gets files from folder dirName whose names match snvmixFnamePattern (regex)
     *
     * @param dirName
     * @param snvmixFnamePattern regexp pattern for matching genotype files
     *
     */
    private void getFnames(String dirName, String snvmixFnamePattern){
        File dir = new File(dirName);
        ArrayList<String> fnames = new ArrayList<String>();
        ArrayList<String> sampleNames = new ArrayList<String>();
        for (File ch : dir.listFiles())
            if (ch.isDirectory())
                for (File child : ch.listFiles())
                    if ((child.getName().matches(snvmixFnamePattern)) && (child.length() != 0)){
                        fnames.add(child.getPath());
                        sampleNames.add(ch.getName());
                    }

        files = fnames.toArray(new String[0]);
        samples = sampleNames.toArray(new String[0]);
    }

    private byte byteGenotype(String[] els){
        int type = Integer.parseInt(els[3].split(",")[5]);
        return (byte) type;
    }

    private void writeIndividuals(String outPrefix) throws IOException {
        System.out.println("\nTotal number of detected individuals: " + samples.length);
        System.out.println("Writing individuals to output file...");

        TextFile sample = new TextFile(outPrefix + ".sample", true);
        sample.writeln("ID_1 ID_2 missing");
        sample.writeln("0 0 0");
        for (String ind : samples) {
            sample.writeln(ind + " " + ind + " 0");
        }
        sample.close();
    }
    /**
     * Creates a genotype matrix containing genotype probabilities: matrix[SNPid][sampleId][type] (type: 0 - homozygote by ref, 1 - heterozygote, 2 - homozygote by alt)
     * @throws IOException
     */

    private void writeGenotypes (String outPrefix, String chr) throws IOException {
        System.out.println("\nWriting genotypes...");

        TextFile gen = new TextFile(outPrefix + chr + ".gen", true);
        for (String snp : snpSet) {
            int snpId = snpMap.get(snp);
            gen.writeln(printGenotypeLine(snp, snpId));
        }
        gen.close();

        System.out.println("Finished writing.");
    }

    /**
     * Writes the genotype line for a snp in the gen format
     * @param snp - snp name
     * @return String with the genotype line
     */
    private String printGenotypeLine(String snp, int snpId){
        String snpName = snp,
                chr = snp.split(":")[0],
                pos = snp.split(":")[1];

        String out = chr + " " + snpName + " " + pos + " " + (char) alleles[snpId][0] + " " + (char) alleles[snpId][1];
        byte type;
        for (int indId = 0; indId < numInds; indId++) {
            type = matrix[snpId][indId];
            if (type == 0)
                out += " 0 0 0";
            else if (type == 1)
                out += " 1 0 0";
            else if (type == 2)
                out += " 0 1 0";
            else
                out += " 0 0 1";
        }
        return out;
    }


    public static Comparator<String> COMPARE_BY_POS = new Comparator<String>(){
        @Override
        public int compare(String snp1, String snp2){
            String chr1 = snp1.split(":")[0];
            String chr2 = snp2.split(":")[0];
            String pos1 = snp1.split(":")[1];
            String pos2 = snp2.split(":")[1];
            int chrCompare = chr1.compareTo(chr2);
            if (chrCompare != 0)
                return chrCompare;
            return pos1.compareTo(pos2);
        }
    };

    public static void main(String[] args) throws IOException {
        SNVMixToGenConverter c = new SNVMixToGenConverter();
        c.parse("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked",
            "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/test", 
            "reads_unique_hits.*cov5.snvmix",
            null);
        /*c.parse("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/test",
                null);      */
    }
}            


