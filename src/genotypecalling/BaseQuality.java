
package genotypecalling;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import jsc.independentsamples.TwoSampleTtest;
import org.apache.commons.lang.ArrayUtils;
import umcg.genetica.io.text.TextFile;
/**
 *
 * @author dashazhernakova
 * 
 * For a SNP gets from a mpileup file the base calling qualities (Phred)
 * Returns the average quality for bases supporting reference allele and supporting alternative allele
 * Also returns the p-value for T-test comparing those two
 */
public class BaseQuality {
    HashMap<String, String> pos2mpileup;
    
    /*
     * Parses the 5th column of mpileup file line
     * Returns the String with insertions, deletions and $ and ^. symbols removed
     */
    private String parseMpileupString(String line){
        String out = "";
        String nucleotides = "acgtACGT";
        
        for (int i = 0; i < line.length(); i++){
            char c = line.charAt(i);
            if ((c == '.') || (c == ',') || (c == '>') || (c == '<')) //reference allle or reference skip
                out+=c;
            else if (nucleotides.contains(String.valueOf(c))) //alternative allele
                out+=c;
            else if (c == '$') //end of read seqment
                continue;
            else if (c == '^') //start of read segment
                i++; //skip the mapping quality following ^
            else if ((c == '+') || (c == '-')){ // skip insertion/deletion
                i++;
                int n = Integer.parseInt(String.valueOf(line.charAt(i))); //number of bases inserted/deleted
                i = i + n;
            }
        }
        return out;
        
    }
    /*
     * Checks whether all elements in an array are equal
     */
    private boolean allElementsEqual(double[] array){
        HashSet<Double> set = new HashSet<Double>();
        for (double d : array){
            set.add(d);
        }
        if (set.size() == 1)
            return true;
        return false;
    }
    /*
     * Gets the average quality for reference allele bases, for alternative allele and the p-value for T-test testin those two
     * Returns float[] with those values 
     */
    public float[] getQualitiesMpileup(String line){
        
        char[] call = parseMpileupString(line.split("\t")[4]).toCharArray();
        byte[] qual = line.split("\t")[5].getBytes();
        
        ArrayList<Double> ref = new ArrayList<Double>();
        ArrayList<Double> alt = new ArrayList<Double>();
        String nucleotides = "acgtACGT";
        int sum_ref = 0;
        int sum_alt = 0;
        float[] freq;
        
        for (int i = 0; i < call.length; i++){
            char c = call[i];
            if ((c == '.') || (c == ',')){
                ref.add((double) qual[i]);
            }
            else if (nucleotides.contains(String.valueOf(c))){
                alt.add((double) qual[i]);
            }
            
        }
        for (double d : ref){
            sum_ref += d;
        }
        for (double d : alt){
            sum_alt += d;
        }
        
        //Do T-test for difference in qualities of reference vs alternative base calls
        double[] ref_ar = ArrayUtils.toPrimitive(ref.toArray(new Double[0])); 
        double[] alt_ar = ArrayUtils.toPrimitive(alt.toArray(new Double[0])); 
        TwoSampleTtest t = null;
        
        //t-test doesn't work when all elements in array are equal. So test that:
        if ((allElementsEqual(alt_ar)) || (allElementsEqual(ref_ar)))
            return new float[] {(float)sum_ref/ref.size(), (float)sum_alt/alt.size(), 0};
        try{
            t = new TwoSampleTtest(ref_ar, alt_ar);
            freq = new float[] {(float)sum_ref/ref.size(), (float)sum_alt/alt.size(), (float) t.getSP()};
        }
        catch(Exception e){
            freq = new float[] {(float)sum_ref/ref.size(), (float)sum_alt/alt.size(), 0};
        }
        return freq;
    }
   
    /*
     * Constructs a HashMap with SNP chr:position mapped to mpileup lines
     * fname - mpileup filename for a sample
     */
    public void makeMpileupMap(String fname) throws IOException{
        TextFile mpileup = new TextFile(fname, false);
        pos2mpileup = new HashMap<String, String>();
        String line, pos;
        while ((line = mpileup.readLine()) != null){
            pos = line.split("\t")[0] + ":" + line.split("\t")[1];
            pos2mpileup.put(pos, line);
        }
        mpileup.close();
    }
    /*
     * Runs getQualitiesMpileup
     * 
     * TODO:remove it
     */
    
    public float[] getQualFrequencies(String pos){
        float[] out = getQualitiesMpileup(pos2mpileup.get(pos));
        if (out[2] == 0.0)
            return null;
        //System.out.println(out[0] + "\t" + out[1] + "\t" + out[2]);
        return out;
    }
    public static void main(String[] args) throws IOException {
         BaseQuality bq = new BaseQuality();
        bq.makeMpileupMap("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/L2/accepted_hits.filtered.iCHIP.mpileup.cov5");
        //bq.getQualitiesMpileup("3	49315589	G	89	>><><><<><<<,..TTT.T.TT.T,,T.,,,,..T.T..T,,..,.T,,..T.TT,,,,,,,,,T,T,,,.,.,,,,,,.,....^~.^~.^~.	gc`gcecieiiccOTB_cEXLBcEE`fBBiiiiFNFFENNQihRFfNFhiNZNNBQiihihicfiU^Fgg^FeFdheccbcZcfce_ba	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        System.out.println(bq.getQualFrequencies("3:49315589"));
        
         
    }
}
