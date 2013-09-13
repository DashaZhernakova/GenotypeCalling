
package genotypecalling;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.TreeMap;

/**
 * Generic genotype class (TreeMap<String, String> SNP2genotype)
 * @author dashazhernakova
 */
public class Genotypes {
    TreeMap<String, String> SNP2genotype;
    HashSet<SNP> snps;

    public Genotypes(){
        SNP2genotype = new TreeMap<String, String>();
    }
    public Genotypes(TreeMap<String, String> genotypes){
        SNP2genotype = genotypes;
        System.out.println("Loaded " + genotypes.size() + " SNPs");
    }
    public void doLDcheck(String TriTyperDir) throws IOException{
        LD ld = new LD(TriTyperDir);
        
        String genotype1, genotype2;
        int wrong = 0, correct = 0;
        String cur_start = SNP2genotype.firstKey();
        int pos1, pos2;
        boolean LDstructureOk;
        
        for (String snp1 : SNP2genotype.keySet()){
            wrong = 0; correct = 0;
            genotype1 = SNP2genotype.get(snp1);
            for (String snp2 : SNP2genotype.tailMap(cur_start, false).keySet()){ 
                if (snp1.equals(snp2))
                    continue;
                
                genotype2 = SNP2genotype.get(snp2);
                pos1 = Integer.parseInt(snp1.split(":")[1]);
                pos2 = Integer.parseInt(snp2.split(":")[1]);
                if (getDistance(snp1, snp2) < 1000000){ //if within 1mb, check LD structure
                    LDstructureOk = ld.checkLD(snp1, genotype1, snp2, genotype2);
                    
                    /*if (! LDstructureOk)
                        System.out.println("WRONG LD STRUCTURE: " + snp1 + "\t" + genotype1 + "\t" + snp2 + "\t" + genotype2);
                    else
                        System.out.println("GOOD LD STRUCTURE: " + snp1 + "\t" + genotype1 + "\t" + snp2 + "\t" + genotype2);
                    */
                }
                else if (pos2 < pos1) //get to the closest snp to the left of the 1mb window
                    cur_start = snp2;
                else if (pos2 > pos1) //surpassed the window, go to the next snp
                    break;
                
            }
            if (wrong + correct > 0)
                System.out.println(snp1 + "\twrong: " + wrong + "\tcorrect: " + correct);
        }
    }

    private int getDistance(String snp1, String snp2){
        String chr1 = snp1.split(":")[0],
                chr2 = snp2.split(":")[0],
                pos1 = snp1.split(":")[1],
                pos2 = snp2.split(":")[1];
        
        if (! chr1.equals(chr2))
            return Integer.MAX_VALUE;
        return Math.abs(Integer.parseInt(pos1) - Integer.parseInt(pos2));
    }
    /*
     * checks if the 2 genotypes are equal (Strings AA, BB, AB)
     */
    public boolean equalGenotypes(String g1, String g2){
        HashSet<String> gt1 = new HashSet (Arrays.asList( String.valueOf(g1.charAt(0)), String.valueOf(g1.charAt(1)) ));
        HashSet<String> gt2 = new HashSet (Arrays.asList( String.valueOf(g2.charAt(0)), String.valueOf(g2.charAt(1)) ));
        return gt1.equals(gt2);
    }
    /*
     * checks if the 2 genotypes are equal (Strings AA, BB, AB). + Checks the complementary genotype if not equal
     */
    public boolean equalGenotypesComplementary(String g1, String g2){
        HashSet<String> gt1 = new HashSet (Arrays.asList( String.valueOf(g1.charAt(0)), String.valueOf(g1.charAt(1)) ));
        HashSet<String> gt2 = new HashSet (Arrays.asList( String.valueOf(g2.charAt(0)), String.valueOf(g2.charAt(1)) ));
        if (gt1.equals(gt2))
            return true;
        else
            return gt1.equals(complement(gt2));
        
    }
    
    //Complements the genotype
    public HashSet<String> complement(HashSet<String> genotype){
        HashSet<String> new_g = new HashSet<String>();
        for (String al : genotype){
            if (al.equals("A"))
                new_g.add("T");
            else if (al.equals("T"))
                new_g.add("A");
            else if (al.equals("G"))
                new_g.add("C");
            else if (al.equals("C"))
                new_g.add("G");
        }
        return new_g;
    }

    public String compare(Genotypes genotypes2){
        int numShared = 0, numSame = 0;
        String gen1, gen2;
        String res = "";
        for (String snp : genotypes2.SNP2genotype.keySet()){
            
            if (SNP2genotype.containsKey(snp)){
                gen1 = SNP2genotype.get(snp);
                gen2 = genotypes2.SNP2genotype.get(snp);
                numShared++;
                if (equalGenotypesComplementary(gen1, gen2))
                    numSame++;
                else
                    System.out.println(snp + "\tFirst: " + gen1 + "\tSecond: " + gen2);
            }
               
        }
        res = numShared + "\t" + numSame + "\t" + 100*numSame/numShared + "%";
        System.out.println("Number of shared SNPs: " + numShared);
        System.out.println("Number of same SNPs: " + numSame + " (" + 100*numSame/numShared + "%)" );
        return res;
    }
    
    public void compare(Genotypes genotypes2, String outFile, String id) throws IOException{
        int numShared = 0, numSame = 0;
        String gen1, gen2;
        File out = new File(outFile);
        if (!out.exists())
            out.createNewFile();
        FileWriter writer = new FileWriter(out, true);
        BufferedWriter bufferedWritter = new BufferedWriter(writer);
        for (String snp : genotypes2.SNP2genotype.keySet()){
            if ( (snp.equals("18:678947")) || (snp.equals("18:675903")) || (snp.equals("18:685797")) ){
                System.out.println(snp);
            }


            if (SNP2genotype.containsKey(snp)){
                gen1 = SNP2genotype.get(snp);
                gen2 = genotypes2.SNP2genotype.get(snp);
                numShared++;
                if (equalGenotypesComplementary(gen1, gen2))
                    numSame++;
                else
                    System.out.println(snp + "\tFirst: " + gen1 + "\tSecond: " + gen2);
            }
        }
        if (numShared == 0)
            System.out.println(id + "\tNo shared genotypes!");
        else
            bufferedWritter.write(id + "\t" +numShared + "\t" + numSame + "\t" + 100*numSame/numShared + "%" + "\n");
        //System.out.println("Number of shared SNPs: " + numShared);
        //System.out.println("Number of same SNPs: " + numSame + " (" + 100*numSame/numShared + "%)" );
        bufferedWritter.close();
    }


    public static void main(String[] args) throws IOException {
        
        String trityper = args[1], snvmix_dir = args[2], gte_file = args[3];

        //TriTyper d = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/dna-seq/TriTyper/all_chr/");
        TriTyper d = new TriTyper(trityper);

        TextFile gte = new TextFile(gte_file, false);
        SNVMix r = new SNVMix();

        Genotypes dg = d.readGenotypes("NA12812");
        Genotypes rg = r.readGenotypes2("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked/ERR188022/reads_unique_hits.sorted.mpileup.cov5.snvmix");
        rg.compare(dg);

        
    }
}
