
package genotypecalling;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Merges 2 TriTyper genotypes. The first dir that is input into constructor has higher priority
 * @author dashazhernakova
 */
public class AddToTriTyper {
    TriTyperGenotypeData genotypeData2;
    TriTyperGenotypeData genotypeData1;
    SNPLoader loader1;
    SNPLoader loader2;
    
    private int numInds;
    private int numSNPs;
    private int numSNPs1;
    private ArrayList<String> absentIn1List;

    public AddToTriTyper(String dir1, String dir2) throws IOException{
        genotypeData1 = new TriTyperGenotypeData();
        genotypeData1.load(dir1);
        loader1 = genotypeData1.createSNPLoader();
        
        genotypeData2 = new TriTyperGenotypeData();
        genotypeData2.load(dir2);
        loader2 = genotypeData2.createSNPLoader();
        
        int numInds1 = genotypeData1.getIndividuals().length;
        int numInds2 = genotypeData1.getIndividuals().length;
        
        if (numInds1 == numInds2){
            numInds = numInds1;
        }
        
    }
    /**
     * Makes a new TriTyper folder with dir1 genotypes + genotypes from dir2
     * @param outDir
     */
    public void add(String outDir) throws IOException{
        HashMap<String, Integer> snps1 = genotypeData1.getSnpToSNPId();
        HashMap<String, Integer> snps2 = genotypeData2.getSnpToSNPId();
        numSNPs1 = snps1.size();
        absentIn1List = new ArrayList<String>();
        HashSet<String> absentIn1 = new HashSet<String>();
        
        //find SNPs present in this TriTyper, but absent in the second TriTyper
        for (String snp : snps2.keySet()){
            if (! snps1.containsKey(snp)){
                //absentIn1.add(snp);
                absentIn1List.add(snp);
            }
        }
        
        System.out.println("Number of SNPs from the first dataset absent in the second: " + absentIn1List.size());
        
        numSNPs = numSNPs1 + absentIn1List.size();
        
        writeIndividuals(outDir);
        
        writeSNPs(outDir);
        
        addGenotypes(outDir);
    }
    
    private void writeIndividuals(String outDir) throws IOException{
        TextFile outFile = new TextFile(outDir + "/Individuals.txt", true);
        String[] inds1 = genotypeData1.getIndividuals();
        String[] inds2 = genotypeData2.getIndividuals();
        for (int i=0; i < numInds; i++){
            if (inds1[i].equals(inds2[i])){
                outFile.writeln(inds1[i]);
            }
            else{
                System.out.println("Individuals are not aligned!!!");
            }
        }
        outFile.close();
    }
    
    /**
     * Writes new SNPs.txt and SNPMappings.txt files
     * @param outDir
     * @throws IOException 
     */
    private void writeSNPs(String outDir) throws IOException {
        TextFile SNPs = new TextFile(outDir + "/SNPs.txt", true);
        TextFile SNPMappings = new TextFile(outDir + "/SNPMappings.txt", true);
        
        //write dataset1 SNPs
        SNP snp;
        for (int snpid = 0; snpid < numSNPs1; snpid++){
            snp = genotypeData1.getSNPObject(snpid);
            SNPs.writeln(snp.getName());
            SNPMappings.writeln(snp.getChr() + "\t" + snp.getChrPos() + "\t" + snp.getName());
        }
        
        //write absent SNPs
        for (String snpName : absentIn1List){ 
            int snpid = genotypeData2.getSnpToSNPId().get(snpName);
            snp = genotypeData2.getSNPObject(snpid);
            SNPs.writeln(snpName);
            SNPMappings.writeln(snp.getChr() + "\t" + snp.getChrPos() + "\t" + snpName);
        }
        
        SNPMappings.close();
        SNPs.close();
    }
    
    /**
     * Writes new genotypes adding those from dataset2
     * @param outdir
     * @throws IOException
     */
    private void addGenotypes(String outdir) throws IOException {
        WGAFileMatrixGenotype genotypeMatrix = new WGAFileMatrixGenotype(numSNPs, numInds, new File(outdir + "GenotypeMatrix.dat"), false);
        WGAFileMatrixImputedDosage dosageMatrix = new WGAFileMatrixImputedDosage(numSNPs, numInds, new File(outdir + "/ImputedDosageMatrix.dat"), false);
        
        SNP snp;
        for (int snpid = 0; snpid < numSNPs1; snpid++){
            snp = genotypeData1.getSNPObject(snpid);
            loader1.loadDosage(snp);
            
            genotypeMatrix.setAlleles(snpid, snp.getAllele1(), snp.getAllele2());
            setDosage(dosageMatrix, snp, snpid);
        }
        
        //add absent in 1
        for (String snpName : absentIn1List){ 
            int snpid = genotypeData2.getSnpToSNPId().get(snpName);
            snp = genotypeData2.getSNPObject(snpid);
            loader1.loadDosage(snp);
            
            genotypeMatrix.setAlleles(snpid, snp.getAllele1(), snp.getAllele2());
            setDosage(dosageMatrix, snp, snpid);
        }  
        genotypeMatrix.close();
        dosageMatrix.close();
    
    } 
    
    private void setDosage(WGAFileMatrixImputedDosage dosageMatrix, SNP snp, int snpid) throws IOException{
        double[] dosageValues = snp.getDosageValues();
        for (int ind = 0; ind < numInds; ind++){
            int dosageInt = (int) Math.round(dosageValues[ind] * 100d);
            byte dosageByte = (byte) (Byte.MIN_VALUE + dosageInt);
            if (dosageInt < 0 || dosageInt > 200) {
                System.out.println("Warning, incorrect dosage!:\t" + dosageInt + "\t" + snp.getName());
            } else {
                dosageMatrix.setDosage(snpid, ind, dosageByte);
            }
        }
    }

    public static void main(String[] args) throws IOException {
        AddToTriTyper a = new AddToTriTyper("",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/impute/CR0.8_pr0.8/TriTyper/");
    }
    
}
