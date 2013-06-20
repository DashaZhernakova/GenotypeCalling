package genotypecalling;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;



/**
 *
 * @author dashazhernakova
 */
public class OLD_PedMap {
    HashMap<String, Integer> snp2line = new HashMap<String, Integer>(); //map of SNP position in map file to its chrNum_position
    
    /*
     * makes a map of SNP position in map file to its rs number
     */
    /*public HashMap<String, String> pos2SNPid(String SNPmappingsF) throws IOException{
        TextFile SNPmappings = new TextFile(SNPmappingsF, false);
        HashMap<String, String> pos2rs = new HashMap<String, String>();
        String[] els;
        while ( (els = SNPmappings.readLineElems(TextFile.tab)) != null){
            pos2rs.put(els[0] + "_" + els[1], els[2]);
        }
        SNPmappings.close();
        return pos2rs;
    }*/
    
    /*
     * chooses shared SNPs 
     */
    public HashMap<Integer, String> getSNPsToCompare(HashMap<String, Integer> snp2line, Set<String> RNAseqSNP){
        HashMap<Integer, String> pos = new HashMap<Integer, String>();
        int posInMap = 0;
        for (String snp : RNAseqSNP){
            if (snp2line.containsKey(snp))
                pos.put(snp2line.get(snp), snp);
            posInMap++;
        }
        return pos;
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
    /*
     * Read genotypes from ped and map files and compare them to mpileup file processed by Marijke's script
     * fName - path and prefix of the ped and map files
     * sampleId - sample id for which to compare genotypes
     * SNP2genotype - genotype to compare with 
     */
    public void readAndCompare(String fName, String sampleId, TreeMap<String, String> SNP2genotype) throws IOException{
        TextFile ped = new TextFile(fName + ".ped", false);
        TextFile map = new TextFile(fName + ".map", false);
        
        //LD ld = new LD("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/RNA-Seq_Sebo/iCHIP_individual212-0169/pedMap/iCHIP_SNP_mappings_hg19.txt");
        
        //get line numbers for SNPs from map file
        String[] els;
        int lineNum = 0;
        
        //check column separator (may be space or tab)
        Pattern sep = TextFile.space;
        els = map.readLineElems(sep);
        if (els.length < 4){
            sep = TextFile.tab;
            els = els[0].split("\t");
            if (els.length < 4)
                System.out.println("Something wrong with field separators: " + els[0]);
        }
        while ( (els) != null){
            
            snp2line.put(els[0] + ":" + els[3], 2*lineNum + 6 );
            lineNum++;
            els = map.readLineElems(sep);
        }
        
        //compare
        int nEq = 0;
        TreeMap<String, String> wrongGenotypes = new TreeMap<String, String>();
        TreeMap<String, String> correctGenotypes = new TreeMap<String, String>();
        HashMap<Integer, String> pos = getSNPsToCompare(snp2line, SNP2genotype.keySet()); //maps lineNumbers to chr:position
        String RNAseq_genotype = "",
                array_genotype = "";
        while ( (els = ped.readLineElems(TextFile.space)) != null){
            if (els[1].equals(sampleId)){
                for (int p : pos.keySet()){
                    RNAseq_genotype = SNP2genotype.get(pos.get(p));
                    array_genotype = els[p] + els[p + 1];
                    if (RNAseq_genotype.isEmpty())
                        System.out.println("");
                    //if (! equalGenotypes(RNAseq_genotype, array_genotype)){
                    if (! equalGenotypesComplementary(RNAseq_genotype, array_genotype)){
                        //System.out.println(pos.get(p) + "\t Array genotype: " + array_genotype + "\t\tRNA-seq genotype: " + RNAseq_genotype +"\tMAF: " + ld.getMAF(pos.get(p)));
                        System.out.println(pos.get(p) + "\t Array genotype: " + array_genotype + "\t\tRNA-seq genotype: " + RNAseq_genotype);
                        nEq ++;
                        wrongGenotypes.put(pos.get(p), RNAseq_genotype);
                    }
                    else
                        correctGenotypes.put(pos.get(p), RNAseq_genotype);
                }
                break;
            }
        }
        
        Genotypes g = new Genotypes();
        g.SNP2genotype = wrongGenotypes;
        g.doLDcheck("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/RNA-Seq_Sebo/iCHIP_individual212-0169/pedMap/iCHIP_SNP_mappings_hg19.txt");
        System.out.println("Overall number of detected SNPs in RNA-seq: " + SNP2genotype.keySet().size());
        System.out.println("Out of " + pos.keySet().size() + " shared SNPs " + nEq + " have a different genotype (" + 100*nEq/pos.keySet().size() + "%)");
        System.out.println("\n\n\n");
        g.SNP2genotype = correctGenotypes;
        //g.doLDcheck("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/RNA-Seq_Sebo/iCHIP_individual212-0169/pedMap/iCHIP_SNP_mappings_hg19.txt");
        
        Mpileup bq = new Mpileup();
        bq.makeMpileupMap("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/L2/accepted_hits.filtered.iCHIP.mpileup.cov5");
        float[] out;
        /*
        for (String posit : correctGenotypes.keySet()){
            out = bq.getQualFrequencies(posit);
            if (posit.equals("10:99080585"))
                System.out.println("10:99080585");
            if (out != null){
                System.out.println(posit + "\t" + correctGenotypes.get(posit) + "\t" + out[0] + "\t" + out[1] + "\t" + out[2]);
                //bq.getQualFrequencies(posit);
            }
        }
        */
        ped.close();
        map.close();
    }
    public static void main(String[] args) throws IOException {
        OLD_PedMap c = new OLD_PedMap();
        //c.readAndCompare("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/genotypeCalls_chr7/output", "NA18486", pileup.SNP2genotype);
        /*System.out.println("\n\n########## SNVMix2 ############\n");
        GetSharedSNPs g = new GetSharedSNPs();
        g.getSharedFromSNVMix("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/", 
                "accepted_hits.filtered.iCHIP.mpileup.cov5.snvmix.gz");
        c.readAndCompare("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/RNA-Seq_Sebo/iCHIP_individual212-0169/pedMap/212-0169_hg19",
                "212-0169",
                g.getSNP2genotype());
        */
        
        
        VCF vcf = new VCF();
        vcf.readGenotypes("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/L2/VarScan/7_accepted_hits.filtered.mpileup.cov5.varscan.vcf");
        
        SNVMix snvmix = new SNVMix();
        
        snvmix.readGenotypesFilterProbability("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/L2/7_accepted_hits.filtered.mpileup.cov5.snvmix", 0.97f);
        //snvmix.readGenotypes("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/L2/iCHIP_masked/accepted_hits.filtered.masked.mpileup.cov5.snvmix");
        System.out.println(snvmix.SNP2genotype.size());
        c.readAndCompare("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/RNA-Seq_Sebo/iCHIP_individual212-0169/pedMap/212-0169_hg19",
                "212-0169",
                snvmix.SNP2genotype);
                
        /*Pileup pileup = new Pileup();
        pileup.readGenotypes("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out/NA18486_yale/7_accepted_hits.filtered.pileup.cov10.out", 0.3f);
        System.out.println("\n\n########## Naive Pileup ############\n");
        PedMap c = new PedMap();
        //c.readAndCompare("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/genotypeCalls_chr7/output", "NA18486", pileup.SNP2genotype);
        System.out.println("\n\n########## SNVMix2 ############\n");
        SNVMix snvmix = new SNVMix();
        snvmix.readGenotypes("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out/NA18486_yale/7_accepted_hits.filtered.pileup.cov3.snvmix2.snvmix1");
        c.readAndCompare("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/genotypeCalls_chr7/output", "NA18486", snvmix.SNP2genotype);
        System.out.println("\n\n########## BCFtools ############\n");
        VCF vcf = new VCF();
        //vcf.readGenotypes("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out/NA18486_yale/7_accepted_hits.filtered.var.flt.cov1.vcf");
        System.out.println(pileup.SNP2genotype.keySet().size() + " SNPs;");
        //c.readAndCompare("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/genotypeCalls_chr7/output", "NA18486", vcf.SNP2genotype);
       */
    }
}
