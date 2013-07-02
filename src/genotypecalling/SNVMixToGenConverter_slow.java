
package genotypecalling;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author dashazhernakova
 */
public class SNVMixToGenConverter_slow {
    String[] files;
    String[] samples;
    private HashMap<String, Integer> individualMap;
    private ArrayList<String> individuals;
    private HashMap<String, Integer> snpMap;
    private HashMap<Integer, byte[]> alleles;
    private TreeSet<String> snpSet;
    private int numInds;
    private int numSNPs;
    float[][][] matrix;
    private TreeSet<String> chromosomes;
    private HashMap<String, String> SNPIdConverter = null;
    //private static final float P_THRESHOLD = 0.95f;
    private float p_threshold = 0.8f;

    
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
                throw new IOException("Error: could not find dir: " + dir);
            }

            getFnames(dir, snvmixFnamePattern);

            System.out.println("Found " + files.length + " SNVMix genotype files");
            if (files.length == 0) {
                    System.exit(0);
            }

            numInds = files.length;

            for (String f : files)
                System.out.println(f);

            
            makeMatrix();

            writeGenotypesByChr(outPrefix, conversionVCFFile);
        }
        
    public void parse(String fileListPath, String outPrefix, String conversionVCFFile) {
        try {
            getFnames(fileListPath);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("\nFailed to get file names. Exiting!");
            System.exit(-1);
        }

        System.out.println("\nFile with samples and file paths: " + fileListPath);
        System.out.println("Output prefix: " + outPrefix);

        System.out.println("Found " + files.length + " SNVMix genotype files");
        if (files.length == 0) {
                System.exit(0);
        }

        numInds = files.length;

        for (String f : files)
            System.out.println(f);

        try {
            makeMatrix();
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("\nFailed to make genotype matrix. Exiting!");
            System.exit(-1);
        }

        try {
            writeGenotypesByChr(outPrefix, conversionVCFFile);
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("\nFailed to write the genotypes. Exiting!");
            System.exit(-1);
        }
        System.out.println("\nFinished successfully");

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
                       sampleNames.add(dirName);
                   }
                
        files = fnames.toArray(new String[0]);
        samples = sampleNames.toArray(new String[0]);
     }

    private void getAllSNPs() throws IOException{
        System.out.println("Getting all SNPs");
        //get a set of all SNPs
        //HashSet<String> allSNPs = new HashSet<String>();
        snpMap = new HashMap<String, Integer>();
        long st = System.currentTimeMillis();
        int snpCtr = 0;
        for (String f : files){
            System.out.println(f);
            TextFile tf = new TextFile(f, TextFile.R);
            String [] elems;
            String snpId;
            int gen;
            float probability;
            while ((elems = tf.readLineElems(TextFile.tab)) != null){
                snpId = elems[0].replaceFirst("chr", "");
                gen = Integer.parseInt(elems[3].split(",")[5]);
                probability = Float.parseFloat(elems[3].split(",")[1 + gen]);

                if (probability >= p_threshold){
                    if (! snpMap.containsKey(snpId)){
                        //allSNPs.add(snpId);
                        snpMap.put(snpId, snpCtr);
                        snpCtr++;
                    }

                }
            }
            tf.close();
        }


        //make SNP map
        System.out.println("\nMaking a SNP set");
        snpSet = new TreeSet<String>(COMPARE_BY_POS);
        snpSet.addAll(snpMap.keySet());
        //String[] uniqueSNPsArray = allSNPs.toArray(new String[0]);
        //allSNPs.clear();

        /*snpMap = new HashMap<String, Integer>();
        //snpChrMap = new HashMap<String, Byte>();
        //snpChrPosMap = new HashMap<String, Integer>();
        snpSet = new TreeSet<String>(COMPARE_BY_POS);
        snpSet.addAll(allSNPs);
        allSNPs.clear();
        //chromosomes = new

        for (String snp : snpSet){
            //Byte chr = ChrAnnotation.parseChr(snp.split(":")[0]);
            //Integer pos = Integer.parseInt(snp.split(":")[1]);

            //chromosomes.add(snp.split(":")[0]);
            //snpSet.add(snp);
            snpMap.put(snp, snpCtr);
            //snpChrPosMap.put(snp, pos);
            //snpChrMap.put(snp, chr);
            snpCtr++;
        } */
        numSNPs = snpSet.size();
        long end = System.currentTimeMillis();
        System.out.println("Time spent: " + (int)((end-st) / 1000) + " sec");
        System.out.println("Number of SNPs: " + numSNPs);
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
     * Gets all SNPs that pass the p-value threshold and are diallelic
     * @throws IOException
     */
    private void getAllSNPsExcludeTriallelic() throws IOException{
        System.out.println("Getting all SNPs, excluding SNPs with more than 2 alleles");
        //get a map of all SNPs mapped to their alleles
        HashMap<String, String> allSNPs = new HashMap<String, String>();
        HashSet<String> triallelic = new HashSet<String>();
        for (String f : files){
            TextFile tf = new TextFile(f, TextFile.R);
            String [] els;
            int gen;
            float probability;
            String alleles;
            while ((els = tf.readLineElems(TextFile.tab)) != null){
                if (els[0].contains("chr"))
                    els[0] = els[0].replaceFirst("chr", "");
                gen = Integer.parseInt(els[3].split(",")[5]);
                probability = Float.parseFloat(els[3].split(",")[1 + gen]);
                alleles = els[1] + els[2];
                if (probability >= p_threshold) {
                    if (triallelic.contains(els[0])) {
                        continue;
                    }
                    if (allSNPs.containsKey(els[0])){
                        if (!allSNPs.get(els[0]).equals(alleles)) {//alleles not equal to previous instance of that SNP => not diallelic
                            allSNPs.remove(els[0]);
                            triallelic.add(els[0]);
                        }
                    }
                    else
                        allSNPs.put(els[0], alleles);
                }
            }
            tf.close();
        }


        //make SNP map and position maps
        String[] uniqueSNPsArray = allSNPs.keySet().toArray(new String[0]);
        allSNPs.clear();

        snpMap = new HashMap<String, Integer>();
        //snpChrMap = new HashMap<String, Byte>();
        //snpChrPosMap = new HashMap<String, Integer>();
        snpSet = new TreeSet<String>(COMPARE_BY_POS);
        int snpCtr = 0;
        for (String snp : uniqueSNPsArray){
            Byte chr = ChrAnnotation.parseChr(snp.split(":")[0]);
            Integer pos = Integer.parseInt(snp.split(":")[1]);

            snpSet.add(snp);
            snpMap.put(snp, snpCtr);
            //snpChrPosMap.put(snp, pos);
            //snpChrMap.put(snp, chr);
            snpCtr++;
        }
        numSNPs = snpSet.size();
    }

    /**
     * Creates a genotype matrix containing genotype probabilities: matrix[SNPid][sampleId][type] (type: 0 - homozygote by ref, 1 - heterozygote, 2 - homozygote by alt)
     * @throws IOException
     */
    private void makeMatrix() throws IOException{
        System.out.println("Making genotype matrix");
        //get all SNPs passing p-value threshold
        getAllSNPs();

        matrix = new float[numSNPs][numInds][2];
        individuals = new ArrayList<String>();
        individualMap = new HashMap<String, Integer>();
        alleles = new HashMap<Integer, byte[]>();

        int indCnt = 0;

        for (String sampleFile : files){
            TextFile sFile = new TextFile(sampleFile, TextFile.R);

            String[] els;
            int gen;
            float[] prob;
            byte ref, alt;

            while ( (els = sFile.readLineElems(TextFile.tab)) != null ){

                if (els[0].contains("chr"))
                    els[0] = els[0].replaceFirst("chr", "");

                ref = els[1].getBytes()[0];
                alt = els[2].getBytes()[0];
                gen = Integer.parseInt(els[3].split(",")[5]);
                prob = new float[] {Float.parseFloat(els[3].split(",")[2]), Float.parseFloat(els[3].split(",")[3]), Float.parseFloat(els[3].split(",")[4])};;

                if (prob[gen - 1] >= p_threshold){
                    //if (! snpMap.containsKey(els[0])) //skip if triallelic
                    //continue;

                    Integer snpId = snpMap.get(els[0]);
                    alleles.put(snpId, new byte[] {ref, alt});

                    matrix [snpId][indCnt][0] = prob[0];
                    matrix [snpId][indCnt][1] = prob[1];
                    //matrix [snpId][indCnt][2] = prob[2];
                }

             }

            sFile.close();

            String sampleId = samples[indCnt];
            individuals.add(sampleId);
            individualMap.put(sampleId, indCnt);
            indCnt++;
        }

    }
    /*
    private void writeGenotypes (String outPrefix, String conversionVCFFile) throws IOException {

        System.out.println("\nTotal number of detected SNPs: " + snpSet.size());
        System.out.println("Writing SNPs to output file...");

        TextFile gen = new TextFile(outPrefix + ".gen", true);

        System.out.println("\nWriting genotypes...");

        if (conversionVCFFile != null){ //convert to rs ids
            System.out.println("Will be using " + conversionVCFFile + " to convert SNP ids");
            makeSNPIdConverter(conversionVCFFile);
        }

        for (String snp : snpSet){
            gen.writeln(printGenotypes(snp));
        }

        gen.close();

        System.out.println("\nTotal number of detected individuals: " + individuals.size());
        System.out.println("Writing individuals to output file...");

        TextFile sample = new TextFile(outPrefix + ".sample", true);
        sample.writeln("ID_1 ID_2 missing");
        sample.writeln("0 0 0");
        for (String ind : individuals) {
                sample.writeln(ind + " " + ind + " 0");
        }

        sample.close();
        System.out.println("\nFinished.");
    }
    */
    private void writeGenotypesByChr (String outPrefix, String conversionVCFFile) throws IOException {

        //Write samples
        System.out.println("\nTotal number of detected individuals: " + individuals.size());
        System.out.println("Writing individuals to output file...");

        TextFile sample = new TextFile(outPrefix + ".sample", true);
        sample.writeln("ID_1 ID_2 missing");
        sample.writeln("0 0 0");
        for (String ind : individuals) {
                sample.writeln(ind + " " + ind + " 0");
        }
        sample.close();


        //Write genotypes
        System.out.println("\nTotal number of detected SNPs: " + snpSet.size());
        System.out.println("\nWriting genotypes by chromosomes...");

        if (conversionVCFFile != null){ //convert to rs ids
            System.out.println("Will be using " + conversionVCFFile + " to convert SNP ids");
            makeSNPIdConverter(conversionVCFFile);
        }
        else
            System.out.println("Will not convert SNP ids");

        String curChr = "", chr;
        int nSNPs = 0;
        TextFile gen = null;
        for (String snp : snpSet){
            chr = snp.split(":")[0];
            if (curChr.equals(chr)){
                gen.writeln(printGenotypes(snp));
                nSNPs++;
            }
            else{
                if (! curChr.isEmpty()){
                    gen.close();
                    System.out.println(" (" + nSNPs + " SNPs)");
                }
                //System.out.println(" (" + nSNPs + " SNPs)");
                nSNPs = 0;
                //next chr
                curChr = chr;
                System.out.print("Writing genotypes on chr" + curChr);
                gen = new TextFile(outPrefix + curChr + ".gen", true);
                gen.writeln(printGenotypes(snp));
                nSNPs++;
            }
        }
        gen.close();
        System.out.println(" (" + nSNPs + " SNPs)");
        System.out.println("\nFinished.");
    }

    /*
     * Makes a HashMap to use when converting snvmix position ids (chr:pos) to rs ids
     * conversionFile - path to the VCF file
     */
    private void makeSNPIdConverter(String conversionVCFFile) throws IOException {
        TextFile tf = new TextFile(conversionVCFFile, false);
        String[] els;
        SNPIdConverter = new HashMap<String, String>();
                        
        while ((els = tf.readLineElems(TextFile.tab))[0].startsWith("#"))
            continue;
        while ((els = tf.readLineElems(TextFile.tab)) != null){
            if (snpSet.contains(els[0] + ":" + els[1]))
                SNPIdConverter.put(els[0] + ":" + els[1], els[2]);
        }
        tf.close();
    }
    private String getSampleNameFromFname(String fname){
        String[] spl = fname.split("/");
        return spl[spl.length - 2];
    }
    
    /**
     * Writes the genotype line for a snp in the gen format
     * @param snp - snp name
     * @return String with the genotype line
     */
    private String printGenotypes(String snp){
        int snpId = snpMap.get(snp);
        String snpName = snp,
        chr = snp.split(":")[0],
        pos = snp.split(":")[1];

        if (SNPIdConverter != null){
            snpName = SNPIdConverter.get(snp);
        }
                
        //String out = snpChrMap.get(snp) + " " + snpName + " " + snpChrPosMap.get(snp) + " " + alleles.get(snpId).charAt(0)+ " " + alleles.get(snpId).charAt(1);
        String out = chr + " " + snpName + " " + pos + " " + (char) alleles.get(snpId)[0] + " " + (char) alleles.get(snpId)[1];
        for (int indId = 0; indId < numInds; indId++) {
            out += " " + matrix[snpId][indId][0] + matrix[snpId][indId][1] + (float)(1 - matrix[snpId][indId][0] - matrix[snpId][indId][1]);
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
        SNVMixToGenConverter_slow c = new SNVMixToGenConverter_slow();
        /*c.parse("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked",
            "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/test", 
            "reads_unique_hits.*cov5.filtered.snvmix", 
            null); */
        c.parse("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/test",
                null);
    }
}            


