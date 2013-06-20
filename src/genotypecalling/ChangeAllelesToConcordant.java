
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class ChangeAllelesToConcordant {
    static final String[] CHROM = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"};
    
    public void change(String inDir, String legendDir, String genPat, String legendPat) throws IOException{
        HashMap<String, byte[]> pos2alleles;
        TextFile gen, out;
        
        for (String chr : CHROM) {
            String genPattern = ".*[/._(chr)]" + chr + genPat;
            String legendPattern = ".*[/._(chr)]" + chr + legendPat;
            
            String genFname = getFname(inDir, genPattern);
            if (genFname == null){
                System.out.println("No .gen file for chromosome " + chr);
                continue;
            }
                
            gen = new TextFile(getFname(inDir, genPattern), false);
            out = new TextFile(genFname.replace(".gen", ".corr.gen"), true);
            pos2alleles = getPos2Alleles(getFname(legendDir, legendPattern));
            
            String[] els;
            int refN = 0, altN = 0, allN = 0;
            while((els = gen.readLineElems(TextFile.space)) != null){
                if (pos2alleles.containsKey(els[2])){
                    byte[] alleles = new byte[2];
                    alleles[0] = els[3].getBytes()[0];
                    alleles[1] = els[4].getBytes()[0];

                    byte[] refAlleles = pos2alleles.get(els[2]);
                    
                    if ( (alleles[0] == refAlleles[0]) && ((alleles[1] == refAlleles[1]))){
                        out.writelnDelimited(els, TextFile.space);
                    }
                    else if ( (alleles[0] == refAlleles[0]) && ((alleles[1] != refAlleles[1]))){
                        els[4] = String.valueOf((char) refAlleles[1]);
                        out.writelnDelimited(els, TextFile.space);
                        altN++;
                    }
                    else if ( (alleles[0] != refAlleles[0]) && ((alleles[1] == refAlleles[1]))){
                        els[3] = String.valueOf((char) refAlleles[0]);
                        out.writelnDelimited(els, TextFile.space);
                        refN++;
                    }
                    else{
                        System.out.println("All alleles not equal for pos " + els[2] + 
                                "\tref: " + String.valueOf((char) refAlleles[0]) + "," + String.valueOf((char) refAlleles[1]) + 
                                        "\tsnvmix: " + String.valueOf((char) alleles[0]) +","+ String.valueOf((char) alleles[1]));
                        allN++;
                    }
                }
                else{
                    out.writelnDelimited(els, TextFile.space);
                }
                
            }
            System.out.println("Chr" + chr + ":\tNum dif alt: " + altN + "; Num dif ref " + refN + "; Num dif all: " + allN);
            gen.close();
            out.close();
        }
    }
    /**
     * Creates a map of positions mapped to alleles according to the reference panel
     * @param legendFname - file name of the legend file
     * @return
     * @throws IOException 
     */
    private HashMap<String, byte[]> getPos2Alleles(String legendFname) throws IOException{
        HashMap<String, byte[]> pos2alleles = new HashMap<String, byte[]>();
        TextFile legend = new TextFile(legendFname, false);
        String[] els = legend.readLineElems(TextFile.space);
        
        while((els = legend.readLineElems(TextFile.space)) != null){
            String pos = els[1];
            byte[] alleles = new byte[2];
            alleles[0] = els[2].getBytes()[0];
            alleles[1] = els[3].getBytes()[0];
            
            pos2alleles.put(pos, alleles);
        }
        
        legend.close();
        return pos2alleles;
    }
    private String getFname(String dirName, String pattern){
        File dir = new File(dirName);
        for (File child : dir.listFiles())
            if (child.getName().matches(pattern))
                return child.getPath();
        return null;
    }
    
    public static void main(String[] args) throws IOException {
        ChangeAllelesToConcordant c = new ChangeAllelesToConcordant();
        c.change("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/impute/CR0.5_pr0.8/tmp_chr10/", 
                "/Users/dashazhernakova/Documents/UMCG/hg19/GoNL/",
                ".filtered_maf0.01_hwe0.0001_cr0.5.gen",
                ".legend.*");
    }
}
