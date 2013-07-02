
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;
import umcg.genetica.io.trityper.util.ChrAnnotation;

/**
 *
 * @author dashazhernakova
 */

public class SNVMixToTriTyperConverter {

	private HashMap<String, Integer> individualMap;
	private ArrayList<String> individuals;
	private HashMap<String, Integer> snpMap;
	private HashMap<String, Byte> snpChrMap;
	private HashMap<String, Integer> snpChrPosMap;
	private ArrayList<String> snpList;
	private int numInds;
        private int numSNPs;
	byte[][][] matrix;
        private HashMap<String, String> SNPIdConverter = null;
	//private static final float P_THRESHOLD = 0.95f;
        private float p_threshold = 0.95f;
        
        /**
	 * Gets files whose names match snvmixFnamePattern (regex)
	 * 
	 * @param dirName
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
       
        public void setPthreshold(float p){
            if ((p >= 0) && (p <= 1)){
                p_threshold = p;
            }
            
        }
	/**
	 * Parses files whose names match snvmixFnamePattern (regex), creates output directory
	 * 
	 * @param dir
	 * @param outdir
	 * @param snvmixFnamePattern regexp pattern for matching genotype files
         * @param dosage boolean: whether to write dosage values
	 * @throws IOException 
	 */
	public void parse(String dir, String outdir, String snvmixFnamePattern, String conversionVCFFile, boolean dosage) throws IOException {

		System.out.println("\nInput dir: " + dir);
                System.out.println("Output dir: " + outdir);
                System.out.println("Pattern: " + snvmixFnamePattern + "\n");
		System.out.println("Writing dosage values? " + dosage);
                
                  if (!Gpio.exists(dir)) {
			throw new IOException("Error: could not find dir: " + dir);
		}

                if (!outdir.endsWith("/")) {
			outdir += "/";
		}

		Gpio.createDir(outdir);


		String[] files = getFnames(dir, snvmixFnamePattern);

		System.out.println("Found " + files.length + " SNVMix genotype files");
		if (files.length == 0) {
			System.exit(0);
		}
                
        numInds = files.length;

        for (String f : files)
            System.out.println(f);

        if (! dosage){
            makeMatrix(files);
            writeGenotypes(outdir, conversionVCFFile);
        }
        else{
            makeMatrixWithDosage(files);
            writeGenotypesDosage(outdir, conversionVCFFile);
        }

	}
        
        /**
         * Gets a set of all SNPs that pass the probability threshold
         * @param files
         * @throws IOException 
         */
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
            snpList = new ArrayList<String>();
            int snpCtr = 0;
            for (String snp : uniqueSNPsArray){
                Byte chr = ChrAnnotation.parseChr(snp.split(":")[0]);
		        Integer pos = Integer.parseInt(snp.split(":")[1]);
                
                snpList.add(snp);
                snpMap.put(snp, snpCtr);
                snpChrPosMap.put(snp, pos);
                snpChrMap.put(snp, chr);
                snpCtr++;
            }
            numSNPs = snpList.size();
        }
        
        /**
         * Makes a byte matrix of genotypes with [SNPid][alleleNumber][individualId]
         * @param files
         * @throws IOException 
         */
        private void makeMatrix(String[] files) throws IOException{
            
            getAllSNPs(files);
            
            matrix = new byte[numSNPs][2][numInds];
            
            individuals = new ArrayList<String>();
            individualMap = new HashMap<String, Integer>();
            
            SNVMix snvmix = new SNVMix();
            int indCnt = 0;
            
            for (String sampleFile : files){
                TextFile sFile = new TextFile(sampleFile, TextFile.R);
                
                String[] els;
                int gen;
                float probability;
                
                while ( (els = sFile.readLineElems(TextFile.tab)) != null ){
                    if (els[0].contains("chr"))
                        els[0] = els[0].replaceFirst("chr", "");
                    
                    Integer snpId = snpMap.get(els[0]); //null if prob < threshold
                    
                    gen = Integer.parseInt(els[3].split(",")[5]);
                    probability = Float.parseFloat(els[3].split(",")[1 + gen]);
                    
                    if ((snpId != null) && (probability >= p_threshold) ){
                        byte[] genotype = snvmix.getByteGenotype(els);
                        matrix[snpId][0][indCnt] = genotype[0];
                        matrix[snpId][1][indCnt] = genotype[1];
                    }
                }
                
                sFile.close();
                
                String sampleId = getSampleNameFromFname(sampleFile);
                individuals.add(sampleId);
                individualMap.put(sampleId, indCnt);
                indCnt++;
            }
            
        }
        
        /**
         * Makes a byte matrix of dosages with [SNPid][0,1][individualId] containing genotype bytes, and [SNPid][2][individualId] containing the dosage value
         * @param files
         * @throws IOException 
         */
        private void makeMatrixWithDosage(String[] files) throws IOException {
            getAllSNPs(files);
            
            matrix = new byte[numSNPs][3][numInds];
            
            individuals = new ArrayList<String>();
            individualMap = new HashMap<String, Integer>();
            
            SNVMix snvmix = new SNVMix();
            int indCnt = 0;
            
            for (String sampleFile : files){
                TextFile sFile = new TextFile(sampleFile, TextFile.R);
                
                String[] els;
                int gen;
                double[] probabilities = new double[3];
                
                while ( (els = sFile.readLineElems(TextFile.tab)) != null ){
                    if (els[0].contains("chr"))
                        els[0] = els[0].replaceFirst("chr", "");
                    
                    Integer snpId = snpMap.get(els[0]); //null if prob < threshold
                    
                    gen = Integer.parseInt(els[3].split(",")[5]);
                    //probability = Float.parseFloat(els[3].split(",")[1 + gen]);
                    probabilities[0] = Float.parseFloat(els[3].split(",")[2]);
                    probabilities[1] = Float.parseFloat(els[3].split(",")[3]);
                    probabilities[2] = Float.parseFloat(els[3].split(",")[4]);

                    if ((snpId != null) && (probabilities[gen - 1] >= p_threshold) ){ 
                        byte[] genotype = snvmix.getByteGenotype(els);
                        matrix[snpId][0][indCnt] = genotype[0];
                        matrix[snpId][1][indCnt] = genotype[1];
                        double dosageValue = probabilities[1] + 2 * probabilities[2];
                        int dosageInt = (int) Math.round(dosageValue * 100d);
                        byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                        if (dosageInt < 0 || dosageInt > 200) {
                            System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpId + "\t" + probabilities[0] + "\t" + probabilities[1] + "\t" + probabilities[2]);
                        } else {
                            matrix[snpId][2][indCnt] = dosageByte;
                        }
                    }
                }
                
                sFile.close();
                
                String sampleId = getSampleNameFromFname(sampleFile);
                individuals.add(sampleId);
                individualMap.put(sampleId, indCnt);
                indCnt++;
            }
        }
        private String getSampleNameFromFname(String fname){
            String[] spl = fname.split("/");
            return spl[spl.length - 2];
        }
        
        private void writeGenotypes (String outdir, String conversionVCFFile) throws IOException {
            
            System.out.println("\nTotal number of detected SNPs: " + snpList.size());
            System.out.println("Writing SNPs to output directory...");
            
            TextFile snpFile = new TextFile(outdir + "SNPs.txt", TextFile.W);
            TextFile snpMappings = new TextFile(outdir + "SNPMappings.txt", TextFile.W);
            
            if (conversionVCFFile == null){ // do not convert to rs ids
                for (String snp : snpList){
                    snpFile.writeln(snp);
                    snpMappings.writeln(snp.split(":")[0] + "\t" + snp.split(":")[1] + "\t" + snp);
                }
            }
            else{ //convert to rs ids
                System.out.println("Will be using " + conversionVCFFile + " to convert SNP ids");
                makeSNPIdConverter(conversionVCFFile);
                for (String snp : snpList){
                    String rsId = SNPIdConverter.get(snp);
                    if (rsId == null)
                        rsId = snp;
                    snpFile.writeln(rsId);
                    snpMappings.writeln(snp.split(":")[0] + "\t" + snp.split(":")[1] + "\t" + rsId);
                }
            }
            snpFile.close();
            snpMappings.close();
            
            
            System.out.println("\nTotal number of detected individuals: " + individuals.size());
            System.out.println("Writing individuals to output directory...");
            
            TextFile indOut = new TextFile(outdir + "Individuals.txt", TextFile.W);
            TextFile indPhenoOut = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);
            
            for (String ind : individuals) {
                    indOut.writeln(ind);
                    indPhenoOut.writeln(ind + "\tunknown\tinclude\tunknown");
            }
            indOut.close();
            indPhenoOut.close();
            
            
            System.out.println("\nWriting genotypes...");
            
            String outfilename = outdir + "GenotypeMatrix.dat";
            WGAFileMatrixGenotype genotypefile = new WGAFileMatrixGenotype(numSNPs, individuals.size(), new File(outfilename), false);
            
            for (int snpId = 0; snpId < numSNPs; snpId++){
                genotypefile.setAlleles(snpId, matrix[snpId][0], matrix[snpId][1]);
            }
            
            genotypefile.close();
            
            System.out.println("\nFinished.");
        }
        
        private void writeGenotypesDosage (String outdir, String conversionVCFFile) throws IOException {
            
            System.out.println("\nTotal number of detected SNPs: " + snpList.size());
            System.out.println("Writing SNPs to output directory...");
            
            TextFile snpFile = new TextFile(outdir + "SNPs.txt", TextFile.W);
            TextFile snpMappings = new TextFile(outdir + "SNPMappings.txt", TextFile.W);
            
            if (conversionVCFFile == null){ // do not convert to rs ids
                for (String snp : snpList){
                    snpFile.writeln(snp);
                    snpMappings.writeln(snp.split(":")[0] + "\t" + snp.split(":")[1] + "\t" + snp);
                }
            }
            else{ //convert to rs ids
                System.out.println("Will be using " + conversionVCFFile + " to convert SNP ids");
                makeSNPIdConverter(conversionVCFFile);
                for (String snp : snpList){
                    String rsId = SNPIdConverter.get(snp);
                    if (rsId == null)
                        rsId = snp;
                    snpFile.writeln(rsId);
                    snpMappings.writeln(snp.split(":")[0] + "\t" + snp.split(":")[1] + "\t" + rsId);
                }
            }
            snpFile.close();
            snpMappings.close();
            
            
            System.out.println("\nTotal number of detected individuals: " + individuals.size());
            System.out.println("Writing individuals to output directory...");
            
            TextFile indOut = new TextFile(outdir + "Individuals.txt", TextFile.W);
            TextFile indPhenoOut = new TextFile(outdir + "PhenotypeInformation.txt", TextFile.W);
            
            for (String ind : individuals) {
                    indOut.writeln(ind);
                    indPhenoOut.writeln(ind + "\tunknown\tinclude\tunknown");
            }
            indOut.close();
            indPhenoOut.close();
            
            
            System.out.println("\nWriting genotypes...");
            
            String outfilename = outdir + "GenotypeMatrix.dat";
            WGAFileMatrixGenotype genotypefile = new WGAFileMatrixGenotype(numSNPs, individuals.size(), new File(outfilename), false);
            WGAFileMatrixImputedDosage matrixDosage = new WGAFileMatrixImputedDosage(numSNPs, individuals.size(), new File(outdir + "/ImputedDosageMatrix.dat"), false);

            for (int snpId = 0; snpId < numSNPs; snpId++){
                genotypefile.setAlleles(snpId, matrix[snpId][0], matrix[snpId][1]);
                for (int sample = 0; sample < individuals.size(); sample++) {
                    matrixDosage.setDosage(snpId, sample, matrix[snpId][2][sample]);
                }
            }
            
            genotypefile.close();
            matrixDosage.close();
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
        HashSet<String> snpSet = new HashSet<String>(snpList);
                
        while ((els = tf.readLineElems(TextFile.tab))[0].startsWith("#"))
            continue;
        while ((els = tf.readLineElems(TextFile.tab)) != null){
            if (snpSet.contains(els[0] + ":" + els[1]))
                SNPIdConverter.put(els[0] + ":" + els[1], els[2]);
        }
        tf.close();
    }
    
    
    
public static void main(String[] args) throws IOException {
    SNVMixToTriTyperConverter c = new SNVMixToTriTyperConverter();
    c.parse("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked", 
            "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp2/", 
            "reads_unique_hits.*cov5.snvmix", 
            null,
            true);
    
        
}

    

    
}


