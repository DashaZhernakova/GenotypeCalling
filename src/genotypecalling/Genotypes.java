
package genotypecalling;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;

/**
 *
 * @author dashazhernakova
 */
public class Genotypes {
    TreeMap<String, String> SNP2genotype;
    //HashMap<String, String> SNP2line;
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
            //if (! snp1.startsWith("4:"))
            //    continue;
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
        //TextFile out = new TextFile("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/comparisonTriTyperToSNVMix.txt", true);
        File out = new File(outFile);
        if (!out.exists())
            out.createNewFile();
        FileWriter writer = new FileWriter(out, true);
        BufferedWriter bufferedWritter = new BufferedWriter(writer);
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
        bufferedWritter.write(id + "\t" +numShared + "\t" + numSame + "\t" + 100*numSame/numShared + "%" + "\n");
        //System.out.println("Number of shared SNPs: " + numShared);
        //System.out.println("Number of same SNPs: " + numSame + " (" + 100*numSame/numShared + "%)" );
        bufferedWritter.close();
    }
    
    public static void main(String[] args) throws IOException {
        
        /*
        SNVMix s = new SNVMix();
        //s.readGenotypesFilterProbability("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked/ERR188032/reads_unique_hits.sorted.mpileup.cov5.filtered.snvmix", 0.95f);
        s.readGenotypesFilterProbability("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked/ERR188032/reads_unique_hits.sorted.mpileup.cov5.snvmix", 0.95f);
        Genotypes snvmix = new Genotypes(s.SNP2genotype);
         
        
        TriTyper t = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp/");
        Genotypes trityper = t.readGenotypes("ERR188032");
        //snvmix.compare(trityper);
        
        t = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/TriTyper/Chr10/");
        Genotypes trityper1 = t.readGenotypes("NA12275");
        snvmix.compare(trityper1);
        trityper.compare(trityper1);
        
        
        */
        
        TriTyper t_rna_dir = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/rna-seq/SNVMix/");
        TriTyper t_dna_dir = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/dna-seq/TriTyper/all_chr/");
        Genotypes t_rna;
        Genotypes t_dna;
        
        TextFile gte = new TextFile("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/ERRtoNA.txt", false);
        Map<String, String> conv = gte.readAsHashMap(0, 1);
        gte.close();
        
        TextFile out = new TextFile("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/SNVMixToTriTyperComparison_QCfilter2_withCovFilter.txt", true);
        for (Entry<String, String> entry : conv.entrySet()){
            try{
                t_rna = t_rna_dir.readGenotypes(entry.getKey(), true);
                t_dna = t_dna_dir.readGenotypes(entry.getValue(), true);
                out.writeln(entry.getKey() + "\t" + t_dna.compare(t_rna));
            }catch (Exception e) {}
        }
        
        out.close();
        
        
        
        
        //snvmix.doLDcheck("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/TriTyper/Chr6/");
        
        
        //Process the whole folder:
        /*
        File dir = new File("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/TriTyper/");
        TriTyper tri;
        for (File f : dir.listFiles()){
            tri = new TriTyper(f.getPath());
            Genotypes trityper = tri.readGenotypes("NA12275");
            snvmix.compare(trityper, "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/comparisonNumbers.txt", f.getName());
            
        }
        */
        
    }
}
