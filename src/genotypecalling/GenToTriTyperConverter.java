
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package genotypecalling;

import umcg.genetica.io.Gpio;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.WGAFileMatrixGenotype;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Converts gen genotypes to TriTyper
 *
 */
public class GenToTriTyperConverter {
    ArrayList<String> vectorSNP;
    ArrayList<String> vectorInd;
    String[] fileNames;
    int nrSNPs;
    int nrInds;
    static final String[] chromosomes = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"};

    /**
     * @param fileName - path to the only genotype file to convert or to the path to file list
     * @param sampleFile
     */
    public GenToTriTyperConverter(String fileName, String sampleFile){
        this(fileName, sampleFile, false);
    }

    /**
     * @param fileName - path to the only genotype file to convert or to the path to file list
     * @param sampleFile
     * @param fileList - boolean; true if fileName is a path to a list of files rather than a gen file
     */
    public GenToTriTyperConverter(String fileName, String sampleFile, boolean fileList) {
        if (! fileList)
            fileNames = new String[] {fileName};
        else{
            FileGetter filegetter = new FileGetter();
            try {
                filegetter.loadFileNames(fileName);
                fileNames = filegetter.fileNames;
            } catch (IOException e) {
                System.out.println("Failed to load files from file list: " + fileName);
                e.printStackTrace();
                System.exit(-1);
            }
        }
        getIndividuals(sampleFile);
    }

    /**
     *
     * @param inDir - folder with gen files
     * @param pattern - gen file name pattern
     * @param sampleFile
     */
    public GenToTriTyperConverter(String inDir, String pattern, String sampleFile){
        FileGetter filegetter = new FileGetter();
        filegetter.loadFileNamesByChr(inDir, pattern);
        fileNames = filegetter.fileNames;

        getIndividuals(sampleFile);
    }

    /**
     * Runs everything
     * @param outputDir - TriTyper folder
     * @throws IOException
     */

    public void convert(String outputDir) throws IOException {
        if (fileNames.length == 1){
            convertFile(fileNames[0], outputDir);
        }
		else{
        	convertPerChr(outputDir);
		}
    }
    /*public void convertPerFile(String outputDir) throws IOException {
        for (int chr = 0; chr < 24; chr++){
            String fileName = fileNames[chr];
            String[] spl = fileName.split("/");
            String sampleOutDir = outputDir + "/" + spl[spl.length - 1].replace(".gen","") + "/";

            System.out.println("\n------------------------------------------");
            System.out.println("Processing input file " + fileName);
            System.out.println("Sample output folder: " + sampleOutDir);

            Gpio.createDir(sampleOutDir);


            writeIndividuals(sampleOutDir);


            getSNPs(fileName);
            writeSNPs(sampleOutDir);

            writeGenotypeMatrix(fileName, sampleOutDir);
        }
    }*/

    private void convertFile(String fileName, String outputDir) throws IOException {
        System.out.println("Processing input file " + fileName);
        System.out.println("Sample output folder: " + outputDir);

        Gpio.createDir(outputDir);


        writeIndividuals(outputDir);


        getSNPs(fileName);
        writeSNPs(outputDir);

        writeGenotypeMatrix(fileName, outputDir);
    }
    /**
     * Converts gen files to TriTyper per chromosome
     * @param outputDir - output folder to contain TriTyper folders per chromosome
     * @throws IOException
     */
    public void convertPerChr(String outputDir) throws IOException {

        String chrOutDir, chrString;
        //for (String chr : chromosomes){
        for (int chr = 1; chr < 24; chr++){
            if (chr == 23)
                chrString = "X";
            else
                chrString = String.valueOf(chr);

            System.out.println("\n--------------------- chr" + chrString + "---------------------");
            chrOutDir = outputDir + "/chr" + chrString + "/";
            Gpio.createDir(chrOutDir);

            writeIndividuals(chrOutDir);

            String fileName = fileNames[chr];
            getSNPs(fileName);
            writeSNPs(chrOutDir);

            writeGenotypeMatrix(fileName, chrOutDir);
        }
    }

    /**
     * Gets a List of SNPs from the gen file
     * @param fileName
     * @throws IOException
     */
    private void getSNPs(String fileName) throws IOException {
        TextFile info = new TextFile(fileName, false);
        String[] els;
        vectorSNP = new ArrayList<String>();
        while( (els = info.readLineElems(TextFile.space)) != null ){
            //vectorSNPMappings.add(chr + "\t" + els[2] + "\t" + chr + ":" + els[2]);
            vectorSNP.add(els[1]);
        }
        nrSNPs = vectorSNP.size();
        System.out.println("Detected " + nrSNPs + " SNPs");
    }

    /**
     * Makes a byte genotype matrix from gen file and writes it to GenotypeMatrix.dat
     * @param fileName
     * @param outputDir
     * @throws IOException
     */
    private void writeGenotypeMatrix(String fileName, String outputDir) throws IOException{
        File genMatrix = new File(outputDir + "/GenotypeMatrix.dat");
        genMatrix.delete();
        WGAFileMatrixGenotype fileMatrixGenotype = new WGAFileMatrixGenotype(nrSNPs, nrInds, genMatrix, false);

        int snpIndex = 0;

        System.out.println("Processing file:\t" + fileName);
        long startTime = System.currentTimeMillis();
        try {
            TextFile in = new TextFile(fileName, TextFile.R);
            String str;

            while ((str = in.readLine()) != null) {
                String data[] = str.split(" ");

                byte[] allele1 = new byte[nrInds];
                byte[] allele2 = new byte[nrInds];

                byte[] alleles = new byte[2];
                alleles[0] = data[3].getBytes()[0];
                alleles[1] = data[4].getBytes()[0];

                for (int sample = 0; sample < nrInds; sample++) {
                    double dosageValue = Double.parseDouble(data[sample * 3 + 5 + 1]) + 2 * Double.parseDouble(data[sample * 3 + 5 + 2]);
                    int dosageInt = (int) Math.round(dosageValue * 100d);
                    byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
                    if (dosageInt < 0 || dosageInt > 200) {
                        System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snpIndex + "\t" + data[sample * 3 + 5] + "\t" + data[sample * 3 + 5 + 1] + "\t" + data[sample * 3 + 5 + 2]);
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

                snpIndex++;
            }
            in.close();
        } catch (IOException e) {
            System.out.println("Error parsing file:\t" + e.getMessage());
            e.printStackTrace();
        }
        System.out.println("Num SNPs: " + snpIndex);
        long endTime = System.currentTimeMillis();
        long duration = endTime - startTime;
        System.out.println("Finished in " + (int)(duration / 1000) + " sec");

        fileMatrixGenotype.close();

    }
    /**
     * Gets imputed gen file names for the chromosome
     * @param inputDir
     * @param chr
     * @return
     */
    private String getGenFile(String inputDir, String chr, String pattern) {
        File dir = new File(inputDir);
        String[] files = dir.list();

        String fname = null;
        String fpattern = pattern.replace("(chr)", chr).replace(".", "\\.");
        System.out.println("Pattern: " + fpattern);
        for (File f : dir.listFiles()) {
            if (f.getName().matches(fpattern)) {
                fname = f.getPath();
                System.out.println("Processing file: " + fname);

                break;
            }
        }

        return fname;
    }

    /**
     * Writes SNPs.txt and SNPMappings.txt (SNP positions determined from their names!)
     * @param outputDir
     * @throws IOException
     */
    private void writeSNPs(String outputDir) throws IOException{
        System.out.println("Writing SNPs to folder: " + outputDir + "\n");
        TextFile mappings = new TextFile(outputDir + "/SNPMappings.txt", true);
        TextFile outSNPFile = new TextFile(outputDir + "/SNPs.txt", true);
        for (String snp : vectorSNP) {
            String[] spl = snp.split(":");
            outSNPFile.writeln(snp);
            mappings.writeln(spl[0] + "\t" + spl[1] + "\t" + snp);
        }
        outSNPFile.close();
        mappings.close();
    }


    /**
     * Gets individuals from .sample file
     * @param sampleFile - path to the sample file
     * @return
     * @throws IOException
     */
    private void getIndividuals(String sampleFile){
        vectorInd = new ArrayList<String>();

        System.out.println("Getting individuals from " + sampleFile);
        TextFile sample = null;
        try {
            sample = new TextFile(sampleFile, false);

            String[] els = sample.readLineElems(TextFile.space);
            els = sample.readLineElems(TextFile.space);
            while ((els = sample.readLineElems(TextFile.space)) != null){
                vectorInd.add(els[0]);
            }
            sample.close();
        } catch (IOException e) {
            System.out.println("\nFailed to find a .sample file!");
            System.out.println("sampleFile=" + sampleFile);
            e.printStackTrace();
            System.exit(-1);
        }
        nrInds = vectorInd.size();
        System.out.println("Found " + nrInds + " individuals");
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

    public static void main(String[] args) throws IOException {
        GenToTriTyperConverter i = new GenToTriTyperConverter("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp/",
                "test[1-9X]+.gen",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp/test.sample");
        i.convert("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/test_trityper/");
    }
}
