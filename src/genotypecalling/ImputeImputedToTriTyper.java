/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;
import umcg.genetica.io.trityper.WGAFileMatrixImputedDosage;

/**
 *
 * @author harmjan
 */
public class ImputeImputedToTriTyper {
    ArrayList<String> vectorSNP;
    ArrayList<String> vectorSNPMappings;
    ArrayList<String> vectorInd;
    HashSet<String> SNPsToInclude;
    HashSet<String> duplicates;
    /**
     * Runs everything
     * @param inputDir - folder with 1 .gen file per chromosome, 1 _info file per chromosome and one or more .sample files
     * @param outputDir - TriTyper folder
     * @throws IOException 
     */
    
    public void convert(String inputDir, String outputDir) throws IOException {
        
        Gpio.createDir(outputDir);
        vectorSNP = new ArrayList<String>();
        vectorSNPMappings = new ArrayList<String>();
        vectorInd = getIndividuals(inputDir);
        writeIndividuals(outputDir);
        rmdup(inputDir);
        getSNPsToInclude(inputDir, 0.3f);
        
        writeSNPs(outputDir, vectorSNP, vectorSNPMappings);

        writeGenotypeMatrix(inputDir, outputDir);
        
    }

    private void writeGenotypeMatrix(String inputDir, String outputDir) throws IOException{
        int nrSNPs = SNPsToInclude.size();
        int nrSamples = vectorInd.size();
        int nSNPs = 0;
        int numAllSNPs = 0;
        System.out.println(vectorSNP.size() + "\t" + SNPsToInclude.size());
        
        File genMatrix = new File(outputDir + "GenotypeMatrix.dat");
        File dosageMatrix = new File(outputDir + "/ImputedDosageMatrix.dat");
        genMatrix.delete();
        dosageMatrix.delete();
        
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrSamples, genMatrix, false);
        WGAFileMatrixImputedDosage matrixImputedDosage = new WGAFileMatrixImputedDosage(nrSNPs, nrSamples, dosageMatrix, false);

        int snpIndex = 0;
        for (int chr = 1; chr <= 22; chr++) {
            // make a file list of batches for this chr....
            String fileName = getGenFile(inputDir, chr);
            System.out.println("Processing file:\t" + fileName);
            long startTime = System.currentTimeMillis();
            try {
                TextFile in = new TextFile(fileName, TextFile.R);
                String str;
                while ((str = in.readLine()) != null) {

                    String data[] = str.split(" ");
                    String snpId = chr + ":" + data[2];
                    
                    //if SNP doesn't pass info threshold: skip it
                    if (SNPsToInclude.contains(snpId)){
                        byte[] allele1 = new byte[nrSamples];
                        byte[] allele2 = new byte[nrSamples];

                        byte[] alleles = new byte[2];
                        alleles[0] = data[3].getBytes()[0];
                        alleles[1] = data[4].getBytes()[0];

                        byte[] dosage = new byte[nrSamples];

                        for (int sample = 0; sample < nrSamples; sample++) {
                            double dosageValue = Double.parseDouble(data[sample * 3 + 5 + 1]) + 2 * Double.parseDouble(data[sample * 3 + 5 + 2]);
                            int dosageInt = (int) Math.round(dosageValue * 100d);
                            byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                            if (dosageInt < 0 || dosageInt > 200) {
                                System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 5] + "\t" + data[sample * 3 + 5 + 1] + "\t" + data[sample * 3 + 5 + 2]);
                            } else {
                                //dosage[sample] = (byte) dosageByte;
                                matrixImputedDosage.setDosage(snpIndex, sample, (byte) dosageByte);
                            }
                            if (dosageValue < 0.5) {
                                allele1[sample] = alleles[0];
                                allele2[sample] = alleles[0];
                            } else {
                                if (dosageValue > 1.5) {
                                    allele1[sample] = alleles[1];
                                    allele2[sample] = alleles[1];
                                } else {
                                    allele1[sample] = alleles[0];
                                    allele2[sample] = alleles[1];
                                }
                            }
                        }
                        fileMatrixGenotype.setAlleles(snpIndex, allele1, allele2);



                        nSNPs ++;
                        snpIndex++;
                    }
                }
                    in.close();
                } catch (IOException e) {
                    System.out.println("Error parsing file:\t" + e.getMessage());
                    e.printStackTrace();
                }
            System.out.println("Num SNPs: " + nSNPs);
            numAllSNPs += nSNPs;
            long endTime = System.currentTimeMillis();
            long duration = endTime - startTime;
            SimpleDateFormat dateFormat = new SimpleDateFormat("HH:mm:ss");
            System.out.println("Finished in " + dateFormat.format(new Date(duration)));
        }
        System.out.println("Num " + numAllSNPs);
        System.out.println("");
        fileMatrixGenotype.close();
        matrixImputedDosage.close();
    }
    /**
     * Gets imputed gen file names for the chromosome
     * @param inputDir
     * @param chr
     * @return 
     */
    private String getGenFile(String inputDir, int chr) {
        File dir = new File(inputDir);
        String[] files = dir.list();

        //ArrayList< String> filelist = new ArrayList<String>();
        String fname = null;
        for (File f : dir.listFiles()) {
            if (f.getPath().matches(".*/chr" + chr + ".*gen.gz")) {
                fname = f.getPath();
                System.out.println(fname);
                if (new File(fname.replace(".gen", ".rmdup.gen")).exists()){
                    fname = fname.replace(".gen", ".rmdup.gen");
                    System.out.println("Using file without duplicates!");
                }
                break;
            }
        }

        return fname;
    }
    
    private void writeSNPs(String outputDir, ArrayList<String> vectorSNP, ArrayList<String> vectorSNPMappings) throws IOException{
        System.out.println("\nWriting SNP mappings to file:");

        TextFile mappings = new TextFile(outputDir + "SNPMappings.txt", true);
        for (String mapp : vectorSNPMappings) {
            mappings.writeln(mapp);
        }
        System.out.println("");
        mappings.close();

        System.out.println("\nWriting SNPs to file:");

        TextFile outSNPFile = new TextFile(outputDir + "SNPs.txt", true);
        for (String snp : vectorSNP) {
            outSNPFile.writeln(snp);
        }
        System.out.println("");
        outSNPFile.close();
    }
    
    
    
    /**
     * gets SNPs that pass info threshold
     * @param genFname
     * @param infoThresh
     * @return
     * @throws IOException 
     */
    private void getSNPsToInclude(String inputDir, float infoThresh) throws IOException{
        SNPsToInclude = new HashSet<String>();
        duplicates = new HashSet<String>();
        int numAllSNPs = 0;
        for (int chr = 1; chr <= 22; chr++) {
            // make a file list of batches for this chr....
            String fileName = getGenFile(inputDir, chr);
            System.out.println("Processing file:\t" + fileName);
            if (fileName != null){
                try {
                    TextFile info = new TextFile(fileName.replace(".gen", "_info"), false);
                    String[] els = info.readLineElems(TextFile.space);
                    float inf;
                    int numSNPs = 0;
                    System.out.println(vectorSNP.size() + "\t" + SNPsToInclude.size() + "\t" + vectorSNPMappings.size());
                    while( (els = info.readLineElems(TextFile.space)) != null ){
                        inf = Float.parseFloat(els[4]);
                        if (inf > infoThresh){
                            vectorSNPMappings.add(chr + "\t" + els[2] + "\t" + chr + ":" + els[2]);
                            vectorSNP.add(chr + ":" + els[2]);
                            SNPsToInclude.add(chr + ":" + els[2]);
                            numSNPs++;
                        }
                        
                    }
                    numAllSNPs += numSNPs;
                    info.close();
                    System.out.println("On chr " + chr + " detected " + numSNPs + " SNPs");
                } catch(IOException e) {
                        System.out.println("Error parsing file:\t" + e.getMessage());
                        e.printStackTrace();
                    }
                }
        }
        System.out.println(vectorSNP.size() + "\t" + SNPsToInclude.size());
        System.out.println("Overall number of SNPs: " + numAllSNPs);
    }
    
    /**
     * Gets individuals from .sample file
     * @param inputDir
     * @return
     * @throws IOException 
     */
    private ArrayList<String> getIndividuals(String inputDir) throws IOException{
        File dir = new File(inputDir);
        //String[] files = dir.listFiles();
        String fName = "";
        for (File f : dir.listFiles()){
            if (f.getPath().endsWith(".sample")){
                fName = f.getPath();
                break;
            }
        }
        if (fName.isEmpty())
            fName = "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/impute/CR0.5_pr0.8/chr1.sample";
        TextFile sample = new TextFile(fName, false);
        ArrayList<String> inds = new ArrayList<String>();
        String[] els = sample.readLineElems(TextFile.space);
        els = sample.readLineElems(TextFile.space);
        while ((els = sample.readLineElems(TextFile.space)) != null){
            inds.add(els[0]);
        }
        sample.close();
        return inds;
    }
    
    private void writeIndividuals(String outputDir) throws IOException{
        TextFile inds = new TextFile(outputDir + "/Individuals.txt", true);
        TextFile phen = new TextFile(outputDir + "/PhenotypeInformation.txt", true);
        for (String ind : vectorInd){
            inds.writeln(ind);
            phen.writeln(ind + "\tunknown\tinclude\tunknown");
        }
        inds.close();
        phen.close();
    }
    
    /**
     * If there are duplicate SNPs creates rmdup files without them
     * @param inputDir 
     */
    private void rmdup(String inputDir){
        for (int chr = 1; chr <= 22; chr++) {
            String fname = getGenFile(inputDir, chr);
            if (fname.contains("rmdup.gen")){
                System.out.println("Duplicates already removed for chr" + chr);
            }
            else if (containsDuplicates(fname, chr)){
                try {
                    removeDuplicates(fname, chr);
                } catch (IOException ex) {
                    Logger.getLogger(ImputeImputedToTriTyper.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }
    /**
     * Goes through the _info file and checks if there are duplicate SNPs
     * @param inputDir
     * @param chr
     * @return 
     */
    private boolean containsDuplicates(String fileName, int chr){
        System.out.println("Processing file:\t" + fileName);
        HashSet<String> snps = new HashSet<String>();
        HashSet<String> dups = new HashSet<String>();
        TextFile info;
        if (fileName != null){
            try {
                info = new TextFile(fileName.replace(".gen", "_info"), false);
                String[] els = info.readLineElems(TextFile.space);
                while( (els = info.readLineElems(TextFile.space)) != null ){
                    if (! snps.contains(els[2])){
                        snps.add(els[2]);
                    }
                    else{
                        System.out.println("Duplicates detected for chr" + chr);
                        return true;
                    }
                }
                info.close();
            } catch (IOException e) {
                System.out.println("Error parsing file:\t" + e.getMessage());
                e.printStackTrace();
            }


        }
        else{
            System.out.println("No file for chr" + chr);
        }
        return false;
    }
    
    /**
     * Creates rmdup files without duplicates (taking the first occurrence and removing all other)
     * @param inputDir
     * @param chr
     * @throws IOException 
     */
    private void removeDuplicates(String fileName, int chr) throws IOException{
        
        System.out.println("Creating rmdup file:\t" + fileName);
        
        TextFile info = new TextFile(fileName.replace(".gen", "_info"), false);
        TextFile info_rmdup = new TextFile(fileName.replace(".gen", ".rmdup_info"), true);
        TextFile gen = new TextFile(fileName, false);
        TextFile gen_rmdup = new TextFile(fileName.replace(".gen", ".rmdup.gen"), true);
        
        HashSet<String> snps = new HashSet<String>();
        
        String info_line = info.readLine();
        info_rmdup.writeln(info_line); //write header
        info_line = info.readLine();
        String gen_line = gen.readLine();
        String id;
        while(info_line != null){
            id = info_line.split(" ")[2];
            if (! gen_line.split(" ")[2].equals(id)){
                System.out.println("Something wrong with gen and info alignment!");
                return;
            }
            if (! snps.contains(id)){
                snps.add(id);
                info_rmdup.writeln(info_line);
                gen_rmdup.writeln(gen_line);
            }
            else{
                System.out.println("Skipping duplicate: " + id);
            }
            info_line = info.readLine();
            gen_line = gen.readLine();
        }
        info.close();
        gen.close();
        info_rmdup.close();
        gen_rmdup.close();
            
    }
    
    public static void main(String[] args) throws IOException {
        ImputeImputedToTriTyper i = new ImputeImputedToTriTyper();
        i.convert("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/impute/CR0.5_pr0.8/", 
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/impute/CR0.5_pr0.8/TriTyper_info0.3/");
    }
}
