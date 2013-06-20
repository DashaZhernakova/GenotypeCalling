package genotypecalling;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.util.BaseAnnot;



/**
 *
 * @author dashazhernakova
 */
public class SNVMix {
    TreeMap<String, String> SNP2genotype;
    //HashMap<String, String> SNP2line;
    HashSet<SNP> snps;
    
    public SNVMix(){
        SNP2genotype = new TreeMap<String, String>();
        //SNP2line = new HashMap<String, String>();
        snps = new HashSet<SNP>();
        
    }
    
    public Genotypes readGenotypes2(String fName) throws IOException{
        TextFile genotypes = new TextFile(fName, false);
        Genotypes gen = new Genotypes();
        String[] els;
        String line;
        while ( (line = genotypes.readLine()) != null ){
            els = line.split("\t");
            if (els[0].contains("chr"))
                els[0] = els[0].replaceFirst("chr", "");
            gen.SNP2genotype.put(els[0], getGenotype(els));
            //SNP2line.put(els[0], line);
            
        }
        
        genotypes.close();
        return gen;
    }
    
    public void readGenotypes(String fName) throws IOException{
        TextFile genotypes = new TextFile(fName, false);
        
        String[] els;
        String line;
        while ( (line = genotypes.readLine()) != null ){
            els = line.split("\t");
            if (els[0].contains("chr"))
                els[0] = els[0].replaceFirst("chr", "");
            SNP2genotype.put(els[0], getGenotype(els));
            //SNP2line.put(els[0], line);
            
        }
        
        genotypes.close();
    }
    public void readGenotypesFilterCoverage(String fName, int N) throws IOException{
        TextFile genotypes = new TextFile(fName, false);
        
        String[] els;
        String line;
        int numReads = 0;
        int i = 0, j = 0;
        while ( (line = genotypes.readLine()) != null ){
            els = line.split("\t");
            if (els[0].contains("chr"))
                els[0] = els[0].replaceFirst("chr", "");
            numReads = Integer.parseInt(els[3].split(",|\\:")[1]) + Integer.parseInt(els[3].split(",|\\:")[3]);
            if (numReads >= N){
                SNP2genotype.put(els[0], getGenotype(els));
                //SNP2line.put(els[0], line);
                i++;
            }
            else
                j++;
            
        }
        System.out.println("pass: " + i);
        System.out.println("fail: " + j);
        genotypes.close();
    }
    public void readGenotypesFilterProbability(String fName, float P) throws IOException{
        TextFile genotypes = new TextFile(fName, false);
        
        String[] els;
        String line;
        int gen;
        float pr;
        while ( (line = genotypes.readLine()) != null ){
            els = line.split("\t");
            if (els[0].contains("chr"))
                els[0] = els[0].replaceFirst("chr", "");
            gen = Integer.parseInt(els[3].split(",")[5]);
            pr = Float.parseFloat(els[3].split(",")[1 + gen]);
            if (pr >= P){
                SNP2genotype.put(els[0], getGenotype(els));
                //SNP2line.put(els[0], line);
            }
        }
        
        genotypes.close();
    }
    
    /**
     * Gets the genotype probability
     * @param snvmixEls - SNVMix line split by tabs
     * @return genotype probability
     */
    public float getProbability (String[] snvmixEls){
        int gen = Integer.parseInt(snvmixEls[3].split(",")[5]);
        return Float.parseFloat(snvmixEls[3].split(",")[1 + gen]);
    }
    
    /**
     * Gets coverage for reference, alternative alleles and the overall coverage
     * @param snvmixEls
     * @return [ref allele coverage, alt allele coverage, overall coverage]
     */
    public int[] getCoverage(String[] snvmixEls){
        String[] spl = snvmixEls[3].split(",");
        int[] coverage = new int[3];
        coverage[0] = Integer.valueOf(spl[0].split(":")[1]);
        coverage[1] = Integer.valueOf(spl[1].split(":")[1]);
        coverage[2] = coverage[0] + coverage[1];
        
        return coverage;
    }
    
    public byte[] getByteGenotype(String[] genotypeLine){
        String type = genotypeLine[3].split(",")[5];
        byte refb = BaseAnnot.toByte(genotypeLine[1]);
        byte altb = BaseAnnot.toByte(genotypeLine[2]);
        if (type.equals("1"))
            return new byte[] {refb, refb};
        if (type.equals("2"))
            return new byte[] {refb, altb};
        if (type.equals("3"))
            return new byte[] {altb, altb};
        return null;
    }
    
    public String getGenotype(String[] genotypeLine){
        String type = genotypeLine[3].split(",")[5];
        if (type.equals("1"))
            return genotypeLine[1] + genotypeLine[1];
        if (type.equals("2"))
            return genotypeLine[1] + genotypeLine[2];
        if (type.equals("3"))
            return genotypeLine[2] + genotypeLine[2];
        return null;
    }
    /*
     * gets genotype as a String, but sorted in lexicological order
     */
    public String getSortedGenotype(String[] genotypeLine){
        String type = genotypeLine[3].split(",")[5];
        if (type.equals("1"))
            return genotypeLine[1] + genotypeLine[1];
        if (type.equals("2")){
            char[] alleles = (genotypeLine[1] + genotypeLine[2]).toCharArray();
            Arrays.sort(alleles);            
            return new String(alleles);
        }
        if (type.equals("3"))
            return genotypeLine[2] + genotypeLine[2];
        return null;
    }
    
}
