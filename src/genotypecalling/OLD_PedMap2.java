
package genotypecalling;

import java.io.IOException;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Pattern;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class OLD_PedMap2 {
    HashMap<String, Integer> snp2line; //map of SNP id (chr:position) to its column number in ped file 
    
    public OLD_PedMap2(){
        snp2line = new HashMap<String, Integer>();
    }
    /*
     * Checks column separator (may be space or tab)
     */
    private Pattern getSeparator(String mapFname) throws IOException{
        TextFile map = new TextFile(mapFname + ".map", false);
        String[] els;
        
        Pattern sep = TextFile.space;
        els = map.readLineElems(sep);
        if (els.length < 4){
            sep = TextFile.tab;
            els = els[0].split("\t");
            if (els.length < 4){
                System.out.println("Something wrong with field separators: " + els[0]);
                return null;
            }
        }
        map.close();
        return sep;
    }
    /*
     * Read genotypes from ped and map files and compare them to mpileup file processed by Marijke's script
     * fName - path and prefix of the ped and map files
     * sampleId - sample id for which to compare genotypes
     * SNP2genotype - genotype to compare with 
     */
    public void readGenotypes(String fName, String sampleId) throws IOException{
        TextFile map = new TextFile(fName + ".map", false);
        
        //get line numbers for SNPs from map file
        String[] els;
        int lineNum = 0;
        
        //check column separator (may be space or tab)
        Pattern sep = getSeparator(fName + ".map");
        
        while ( (els = map.readLineElems(sep)) != null){
            snp2line.put(els[0] + ":" + els[3], 2*lineNum + 6 );
            lineNum++;
        }
        map.close();
        
        // process the genotypes
        TextFile ped = new TextFile(fName + ".ped", false);
        Genotypes gen = new Genotypes();
        
        while ( (els = ped.readLineElems(TextFile.space)) != null){
            if (els[1].equals(sampleId)){
                for (String snp : snp2line.keySet()){
                    int col = snp2line.get(snp);
                    gen.SNP2genotype.put(snp, els[col] + els[col + 1]);
                }
                break;
            }
        }
        
        ped.close();
        
    }
    
    
}
