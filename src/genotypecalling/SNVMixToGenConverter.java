
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import org.apache.commons.lang.StringUtils;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author dashazhernakova
 */
public class SNVMixToGenConverter {
    private HashMap<String, Integer> individualMap;
    private ArrayList<String> individuals;
    private HashMap<String, Integer> snpMap;
    private HashMap<String, Byte> snpChrMap;
    private HashMap<String, Integer> snpChrPosMap;
    private HashMap<Integer, String> alleles;
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

            String[] files = getFnames(dir, snvmixFnamePattern);

            System.out.println("Found " + files.length + " SNVMix genotype files");
            if (files.length == 0) {
                    System.exit(0);
            }

            numInds = files.length;

            for (String f : files)
                System.out.println(f);

            
            makeMatrix(files);

            writeGenotypesByChr(outPrefix, conversionVCFFile);
        }
        /**
	 * Gets files whose names match snvmixFnamePattern (regex)
	 * 
	 * @param dir
	 * @param snvmixFnamePattern regexp pattern for matching genotype files
	 * 
	 */
        private String[] getFnames(String dirName, String snvmixFnamePattern){
        File dir = new File(dirName);
        ArrayList<String> fnames = new ArrayList<String>();
        for (File ch : dir.listFiles())
            if (ch.isDirectory())
                for (File child : ch.listFiles())
                   if ((child.getName().matches(snvmixFnamePattern)) && (child.length() != 0))
                       fnames.add(child.getPath());
                
        return fnames.toArray(new String[0]);           
        }
        
        private void getAllSNPs(String[] files) throws IOException{
            
            //get a set of all SNPs
            HashSet<String> allSNPs = new HashSet<String>();
            for (String f : files){
                TextFile tf = new TextFile(f, TextFile.R);
                String [] elems;
                
                while ((elems = tf.readLineElems(TextFile.tab)) != null){
                    if (elems[0].contains("chr"))
                        elems[0] = elems[0].replaceFirst("chr", "");
                    int gen = Integer.parseInt(elems[3].split(",")[5]);
                    float probability = Float.parseFloat(elems[3].split(",")[1 + gen]);
                    
                    if (probability >= p_threshold)
                        allSNPs.add(elems[0]);
                }
                tf.close();
            }
            
            //make SNP map
            String[] uniqueSNPsArray = allSNPs.toArray(new String[0]);
            allSNPs.clear();
            snpMap = new HashMap<String, Integer>();
            snpChrMap = new HashMap<String, Byte>();
            snpChrPosMap = new HashMap<String, Integer>();
            snpSet = new TreeSet<String>(COMPARE_BY_POS);
            //chromosomes = new 
            int snpCtr = 0;
            for (String snp : uniqueSNPsArray){
                Byte chr = ChrAnnotation.parseChr(snp.split(":")[0]);
		Integer pos = Integer.parseInt(snp.split(":")[1]);
                
                //chromosomes.add(snp.split(":")[0]);
                snpSet.add(snp);
                snpMap.put(snp, snpCtr);
                snpChrPosMap.put(snp, pos);
                snpChrMap.put(snp, chr);
                snpCtr++;
            }
            numSNPs = snpSet.size();
        }
        
        private void makeMatrix(String[] files) throws IOException{
            
            getAllSNPs(files);
            
            matrix = new float[numSNPs][numInds][3];
            
            individuals = new ArrayList<String>();
            individualMap = new HashMap<String, Integer>();
            alleles = new HashMap<Integer, String>();
            
            SNVMix snvmix = new SNVMix();
            int indCnt = 0;
            
            for (String sampleFile : files){
                TextFile sFile = new TextFile(sampleFile, TextFile.R);
                
                String[] els;
                int gen;
                float[] prob;
                String ref, alt;
                
                while ( (els = sFile.readLineElems(TextFile.tab)) != null ){
                    
                    if (els[0].contains("chr"))
                        els[0] = els[0].replaceFirst("chr", "");
                    
                    ref = els[1];
                    alt = els[2];
                    
                    gen = Integer.parseInt(els[3].split(",")[5]);
                    prob = new float[] {Float.parseFloat(els[3].split(",")[2]), Float.parseFloat(els[3].split(",")[3]), Float.parseFloat(els[3].split(",")[4])};;
                    
                    if (prob[gen - 1] >= p_threshold){
                        Integer snpId = snpMap.get(els[0]);
                        
                        alleles.put(snpId, ref + alt);
                        
                        matrix [snpId][indCnt][0] = prob[0];
                        matrix [snpId][indCnt][1] = prob[1];
                        matrix [snpId][indCnt][2] = prob[2];
                    }
                    
                 }
                
                sFile.close();
                
                String sampleId = getSampleNameFromFname(sampleFile);
                individuals.add(sampleId);
                individualMap.put(sampleId, indCnt);
                indCnt++;
            }
            
        }
        
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
        
        private void writeGenotypesByChr (String outPrefix, String conversionVCFFile) throws IOException {
            
            System.out.println("\nTotal number of detected individuals: " + individuals.size());
            System.out.println("Writing individuals to output file...");
            
            TextFile sample = new TextFile(outPrefix + ".sample", true);
            sample.writeln("ID_1 ID_2 missing");
            sample.writeln("0 0 0");
            for (String ind : individuals) {
                    sample.writeln(ind + " " + ind + " 0");
            }
                        
            sample.close();
            
            
            System.out.println("\nTotal number of detected SNPs: " + snpSet.size());
            System.out.println("Writing SNPs to output file...");
            
                      
            System.out.println("\nWriting genotypes by chomosomes...");
                        
            if (conversionVCFFile != null){ //convert to rs ids
                System.out.println("Will be using " + conversionVCFFile + " to convert SNP ids");
                makeSNPIdConverter(conversionVCFFile);
            }
            
            String curChr = "", chr;
            TextFile gen = null;
            for (String snp : snpSet){
                chr = snp.split(":")[0];
                if (curChr.equals(chr)){
                    gen.writeln(printGenotypes(snp));
                }
                else{
                    if (! curChr.isEmpty()){
                        gen.close();
                    }
                    
                    //next chr
                    curChr = chr;
                    System.out.println("Writing genotypes on chr" + curChr);
                    gen = new TextFile(outPrefix + "_" + curChr + ".gen", true);
                    gen.writeln(printGenotypes(snp));
                }
                
            }
            
            gen.close();   
            
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
        String snpName = snp;
        if (SNPIdConverter != null){
            snpName = SNPIdConverter.get(snp);
        }
                
        String out = snpChrMap.get(snp) + " " + snpName + " " + snpChrPosMap.get(snp) + " " + alleles.get(snpId).charAt(0)+ " " + alleles.get(snpId).charAt(1);
        
        for (float[] f : matrix [snpId]){
            for (float ff : f){
                out += " " + ff;
            }
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
            "reads_unique_hits.*cov5.filtered.snvmix", 
            null);
    }
}            


