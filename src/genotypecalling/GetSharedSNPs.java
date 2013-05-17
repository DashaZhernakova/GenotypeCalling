
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class GetSharedSNPs {
    //HashMap<String, String> snps = new HashMap<String, String>();
    HashSet<String> snps = new HashSet<String>(); //chr:position:sorted genotype
    ArrayList<String> fileNames = new ArrayList<String>();
    
    /*
     * Scans through output files rom SNVMix located in dirName, with file name fnamePattern, retains only shared SNPs (same genotypes) in the HashSet snps
     */
    public void getSharedFromSNVMix(String dirName, String fnamePattern) throws IOException{
        getFilenames(dirName, fnamePattern);
        HashSet<String> cur_snps;
        String[] els;
        SNVMix snvmix = new SNVMix();
        //String p = "/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/", f = "accepted_hits.filtered.iCHIP.mpileup.cov5.snvmix.gz";
        //fileNames.add(p+"L5/"+f);
        //fileNames.add(p+"L4/"+f);
        for (String fName : fileNames){
            System.out.println("Processing " + fName.split("/")[fName.split("/").length - 2]);
            TextFile in = new TextFile(fName, false);
            cur_snps = new HashSet<String>();
            while ( (els = in.readLineElems(TextFile.tab)) != null ){
                    cur_snps.add(els[0] + ":" + snvmix.getSortedGenotype(els));
            }
            //if first file
            if (snps.isEmpty()){
                snps.addAll(cur_snps);
                System.out.println("Number of SNPs: " + snps.size());
            }
            else{ //retain only shared
                System.out.println("Number of SNPs: " + snps.size());
                for (Iterator<String> it = snps.iterator(); it.hasNext(); ){
                    String snp = it.next();
                    if (! cur_snps.contains(snp))
                        it.remove();
                }
                System.out.println("Number of shared SNPs with all previous samples: " + snps.size());
            }
            in.close();
        }
    }
    private void getFilenames(String dirName, String fnamePattern) throws IOException{
        System.out.println("Started");
        File dir = new File(dirName);
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory()) {
                for (File child : ch.listFiles()){
                   if (child.getName().equals(fnamePattern) ){
                       fileNames.add(child.getPath());
                  }
                }
            }
        }
    }
    public TreeMap<String, String> getSNP2genotype(){
        TreeMap<String, String> SNP2genotype = new TreeMap<String, String>();
        for (String snp : snps){
            SNP2genotype.put(snp.split(":")[0] + ":" + snp.split(":")[1], snp.split(":")[2]);
        }
        return SNP2genotype;
    }
    
    public static void main(String[] args) throws IOException {
        GetSharedSNPs g = new GetSharedSNPs();
        g.getSharedFromSNVMix("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/", 
                "accepted_hits.filtered.mpileup.cov5.snvmix.gz");
    }
}
