
package genotypecalling;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import umcg.genetica.containers.Pair;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author dashazhernakova
 */
public class LD {
        TriTyperGenotypeData genotypeData;
        SNPLoader loader;
        HashMap<String, String> pos2id;
    /*public LD() throws IOException{
        genotypeData = new TriTyperGenotypeData();
        //genotypeData.load("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/");
        genotypeData.load("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/iCHIP_genotypes/TriTyper/");
        loader = genotypeData.createSNPLoader();
        }
    */
     public LD(String TriTyperDir) throws IOException{
        genotypeData = new TriTyperGenotypeData();
        //genotypeData.load("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/");
        genotypeData.load(TriTyperDir);
        loader = genotypeData.createSNPLoader();
        makePos2id(TriTyperDir + "/SNPMappings.txt");
        }
     
    public boolean equalGenotypes(String g1, String g2){
        HashSet<String> gt1 = new HashSet (Arrays.asList( String.valueOf(g1.charAt(0)), String.valueOf(g1.charAt(1)) ));
        HashSet<String> gt2 = new HashSet (Arrays.asList( String.valueOf(g2.charAt(0)), String.valueOf(g2.charAt(1)) ));
        return gt1.equals(gt2);
    }
    public boolean checkLD(String pos1, String genotype1, String pos2, String genotype2) throws IOException{
        //TriTyperGenotypeData d = new TriTyperGenotypeData();
        //d.load("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/");
        
        try{
        String id1 = pos2id.get(pos1);
        String id2 = pos2id.get(pos2);
        Integer snpid1 = genotypeData.getSnpToSNPId().get(id1);
        Integer snpid2 = genotypeData.getSnpToSNPId().get(id2);
        SNP snp1 = genotypeData.getSNPObject(snpid1);
        SNP snp2 = genotypeData.getSNPObject(snpid2);
        
        loader.loadGenotypes(snp1);
        loader.loadGenotypes(snp2);
        
        DetermineLD ldcalc = new DetermineLD();
        Pair<Double, Double> pair = ldcalc.getLD(snp1, snp2, genotypeData, DetermineLD.INCLUDE_CASES_AND_CONTROLS, false);
        double dp = pair.getLeft();
        double r2 = pair.getRight();
        String predicted = "";
        if (r2 >= 0.8){
            predicted = predictGenotype(snp1, snp2, genotype1);
            if (equalGenotypes(predicted, genotype2)){
                System.out.println("CORRECT LD STRUCTURE" + "\t" + id1 + "\t" + snp1.getChr() + ":" + snp1.getChrPos() + "\t" + genotype1 + "\t" + id2 + "\t" + snp2.getChr() + ":" + snp2.getChrPos() + "\t" + genotype2);
                return true;
            }
            else{
                System.out.println("WRONG LD STRUCTURE" + "\t" + id1 + "\t" + snp1.getChr() + ":" + snp1.getChrPos() + "\t" + genotype1 + "\t" + id2 + "\t" + snp2.getChr() + ":" + snp2.getChrPos() + "\t" + genotype2 + "\tPredicted genotype: " + predicted);
                return false;
            }
                
        }
        
        snp1.clearGenotypes();
        snp2.clearGenotypes();
        //System.out.println("Not in LD");
        } catch (Exception e){
            //System.out.println("No such SNP: " + id1 +" or " + id2);
        }
        return true;
    }
    public String predictGenotype(SNP snp1, SNP snp2, String genotype1){
        byte[] gen = genotype1.getBytes();
        byte[] res_gen = new byte[2];
        byte minor_al1 = snp1.getMinorAllele(),
                minor_al2 = snp2.getMinorAllele();
        
        if (gen[0] == gen[1]){
            if (gen[0] == minor_al1){
                res_gen[0] = res_gen[1] = minor_al2;
            }
            else{
                res_gen[0] = res_gen[1] = getAlternativeAllele(snp2.getAlleles(), minor_al2);
            }
        }
        else{
            res_gen[0] = snp2.getAlleles()[0];
            res_gen[1] = snp2.getAlleles()[1];
        }
        return new String(res_gen);
        
    }
    private byte getAlternativeAllele(byte[] alleles, byte al){
        for (byte b : alleles){
            if (b != al)
                return b;
        }
        return 0;
    }
    public double getMAF(String pos){
        try{
            String id = pos2id.get(pos);
            Integer snpid = genotypeData.getSnpToSNPId().get(id);
            
            SNP snp = genotypeData.getSNPObject(snpid);
            loader.loadGenotypes(snp);
            return snp.getMAF();
        } catch(Exception e){}
        return 4;
        
    }
    private void makePos2id(String snp_mappings_file) throws IOException{
        TextFile snp_mappings = new TextFile(snp_mappings_file, false);
        String [] els;
        String s;
        SNP snp;
        
        // generate a map of SNP ids chr:pos to rs ids
        pos2id = new HashMap<String, String>();
        while ( (els = snp_mappings.readLineElems(TextFile.tab)) != null ){
            pos2id.put(els[0] + ":" + els[1], els[2]);
        }
    }
    public static void main(String[] args) throws IOException {
        //LD ld = new LD();
        //ld.checkLD("rs10276096", "TT", "rs7792674", "TC");
        
        /*
        char[] ch = new char[] {'g', 'c', 'e', 'c', 'i'};
        byte[] b = new byte[] {64, 65, 66};
        String s = new String(b);
        System.out.println(s);
        b = s.getBytes();
        for (byte bb : b){
            System.out.println(bb);
        }
        */
    }
    
}
