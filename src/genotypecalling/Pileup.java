package genotypecalling;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 * 
 * Genotypes called by samtools mpileup tool and then converted by Marijke's script computeCountsFromPileup.pl
 */
public class Pileup {
    HashMap<String, float[]> calledSNPsMap = new HashMap<String, float[]>();
    TreeMap<String, String> SNP2genotype = new TreeMap<String, String>();
    /*
     * reads the output of Marijke's script
     */
    public void readFromPileup(String fName) throws IOException{
        TextFile pileup = new TextFile(fName, false);
        //HashMap<String, String> pos2rs = pos2SNPid(SNPmappingsF);
        String[] calledSNPs = pileup.readLineElems(TextFile.tab);
        while ( (calledSNPs = pileup.readLineElems(TextFile.tab)) != null){
            float[] alleles = {Float.valueOf(calledSNPs[7])/Float.valueOf(calledSNPs[3]), 
                Float.valueOf(calledSNPs[8])/Float.valueOf(calledSNPs[3]), 
                Float.valueOf(calledSNPs[9])/Float.valueOf(calledSNPs[3]), 
                Float.valueOf(calledSNPs[10])/Float.valueOf(calledSNPs[3])};
            calledSNPsMap.put(calledSNPs[0] + ":" + calledSNPs[1], alleles);    
        }
        int het = 0, all = 0;
        for (String snp : calledSNPsMap.keySet()){
            float[] alleles = calledSNPsMap.get(snp);
            /*if (isHeterozygous(alleles, 0.3f)){
                //System.out.println(snp + " : " + alleles[0] + ", " + alleles[1] + ", " + alleles[2] + ", " + alleles[3]);
                het++;
                System.out.println(snp + " : " + getGenotype(alleles, 0.3f));
            }
            all++;
             * 
             */
        }
        //System.out.println(het + " heterozygous out of " + all + " (" + 100*het/all + "%)");
             
        pileup.close();
        
    }
    
    public void readGenotypes(String fName, float threshold) throws IOException{
        TextFile pileup = new TextFile(fName, false);
        //HashMap<String, String> pos2rs = pos2SNPid(SNPmappingsF);
        String[] calledSNPs = pileup.readLineElems(TextFile.tab);
        String genotype;
        while ( (calledSNPs = pileup.readLineElems(TextFile.tab)) != null){
            float[] alleles = {Float.valueOf(calledSNPs[7])/Float.valueOf(calledSNPs[3]), 
                Float.valueOf(calledSNPs[8])/Float.valueOf(calledSNPs[3]), 
                Float.valueOf(calledSNPs[9])/Float.valueOf(calledSNPs[3]), 
                Float.valueOf(calledSNPs[10])/Float.valueOf(calledSNPs[3])};
            genotype = getGenotype(alleles, threshold);
            if (! genotype.isEmpty())
                SNP2genotype.put(calledSNPs[0] + ":" + calledSNPs[1], genotype);
        }
        pileup.close();
        
    }
    
    public boolean isHeterozygous(float[] alleles, float threshold){
        boolean hetero = false;
        for (float al : alleles){
            if ((al > threshold) && (hetero))
                return true;
            else if (al > threshold)
                hetero = true;
        }
        return false;
    }
    
    public String indexToNucleotide(int ind){
        if (ind == 0) return "A";
        if (ind == 1) return "C";
        if (ind == 2) return "G";
        if (ind == 3) return "T";
        return null;
    }
    
    /*
     * gets genotype as a list of alleles from the fractions of number of reads having each SNP variant
     * alleles - 4 coverage numbers
     * threshold - threshold to consider a caled variant
     * 
     */
    public ArrayList<String> getArrayGenotype(float[] alleles, float threshold){
        ArrayList<String> genotype = new ArrayList<String>(2);
        for (int i = 0; i < alleles.length; i++){
            if (alleles[i] > threshold)
                genotype.add(indexToNucleotide(i)); 
        }
        if (genotype.size() == 1){ //if homozygous
            genotype.add(genotype.get(0));
        }
        //System.out.println(genotype);
        return genotype;
    }
    
    /*
     * gets genotype as a String AB from the fractions of number of reads having each SNP variant
     * alleles - 4 coverage numbers
     * threshold - threshold to consider a caled variant
     * 
     */
    public String getGenotype(float[] alleles, float threshold){
        String genotype = "";
        for (int i = 0; i < alleles.length; i++){
            if (alleles[i] > threshold)
                genotype += indexToNucleotide(i); 
        }
        if (genotype.length() == 1){ //if homozygous
            genotype += genotype.charAt(0);
        }
        //System.out.println(genotype);
        return genotype;
    }
    
    /*
     * Runs Marijke's script to convert raw mpileup output to counts
     */
    public void computeCountsFromPileup(String pileup, String out, String minCoverage){
        try{
            Process p = Runtime.getRuntime().exec("perl " + System.getProperty("user.dir") + "/src/genotypecalling/computeCountsFromPileup.pl " + pileup + " " + out + " " + minCoverage);
        }catch(IOException e){
            System.out.println(e);
        }
    }
    public static void main(String[] args) {
        Pileup p = new Pileup();
        
    }
}
