package genotypecalling;

import java.io.IOException;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.converters.TriTyperToVCF;
import umcg.genetica.io.trityper.converters.VCFToTriTyper;

/**
 *
 * @author dashazhernakova
 */
public class VCF {
    TreeMap<String, String> SNP2genotype = new TreeMap<String, String>();
    
    private int getSampleColNum(String[] els, String sampleId){
        for (int i = 0; i < els.length; i++){
            String colName = els[i];
            if (colName.equals(sampleId)){
                return i;
            }
        }
        return 0;
    }
    public Genotypes readGenotypes2(String fName, String sampleId) throws IOException{
        Genotypes gen = new Genotypes();
        TextFile file = new TextFile(fName, false);
        String[] els;
        int sampleColNum = 0;
        while ( (els = file.readLineElems(TextFile.tab))[0].startsWith("#") ){
            //skip the header
            if (els[0].startsWith("#CHROM")){
                //get the column number for this sample
                sampleColNum = getSampleColNum(els, sampleId);
                break;
            }
        }
        if (sampleColNum == 0){
            System.out.println("No sample with sample id: " + sampleId + ", exiting!");
            return null;
        }
        
        //process the genotype lines
        int numSNPs = 0;
        while ( (els = file.readLineElems(TextFile.tab)) != null ){
            if (els[4].length() == 1){ //to skip indels
            gen.SNP2genotype.put(els[0] + ":" + els[1], getGenotype(els, sampleColNum));
            numSNPs++;
            }
        }
        System.out.println("Loaded " + numSNPs + " SNPs");
        file.close();
        return gen;
    }
    
    public void readGenotypes(String fName) throws IOException{
        TextFile genotypes = new TextFile(fName, false);
        
        String[] els;
        while ( (els = genotypes.readLineElems(TextFile.tab))[0].startsWith("#") ){}
        while ( (els = genotypes.readLineElems(TextFile.tab)) != null ){
            SNP2genotype.put(els[0] + ":" + els[1], getGenotype(els));
        }
        genotypes.close();
    }
    
    /*
     * reads the genotype as a String of 2 nucleotides if there is one sample per vcf file
     */
    private String getGenotype(String[] spl_line){
        String ref = spl_line[3], 
                alt = spl_line[4], 
                gen = spl_line[9].split(":")[0];
        String al1 = gen.split("/|[|]")[0],
                al2 = gen.split("/|[|]")[1];
        return numToNucleotide(ref, alt, al1) + numToNucleotide(ref, alt, al2);
    }
    
    /*
     * reads the genotype as a String of 2 nucleotides if there are multiple samples per vcf file
     */
    private String getGenotype(String[] spl_line, int sampleColNum){
        String ref = spl_line[3], 
                alt = spl_line[4], 
                gen = spl_line[sampleColNum].split(":")[0];
        String al1 = gen.split("/|[|]")[0],
                al2 = gen.split("/|[|]")[1];
        return numToNucleotide(ref, alt, al1) + numToNucleotide(ref, alt, al2);
    }
    
    /*
     * Converts the vcf-style genotype to a String of alleles
     */
    private String numToNucleotide(String ref, String alt, String num){
        if (num.equals("0"))
            return ref;
        String[] alts = alt.split(",");
        if (alts[Integer.parseInt(num) - 1] != null)
            return alts[Integer.parseInt(num) - 1];
        else
            System.out.println("Bad allele number! " + num + " (alt = " + alt);
        return null;
    }
    
    public void convertFromTriTyper(String pathToTriTyper, String outPath) throws IOException{
        TriTyperToVCF conv = new TriTyperToVCF();
        conv.convert(pathToTriTyper, outPath, "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/iCHIP_genotypes/TriTyper/SNPs.txt");
    }
    
    public void convertToTriTyper(String pathToTriTyper, String vcfPath, String snppattern) throws IOException{
        VCFToTriTyper t = new VCFToTriTyper();
        t.parse(vcfPath, pathToTriTyper, snppattern);
    }
    
    public static void main(String[] args) throws IOException {
       VCF vcf = new VCF();
       /*vcf.readGenotypes("/Users/dashazhernakova/Documents/UMCG/data/lincRNA_Sebo/mappedData/L2/7_var.flt.vcf");
       int i = 0;
       for (String pos : vcf.SNP2genotype.keySet()){
           System.out.println(pos + "\t" + vcf.SNP2genotype.get(pos));
           i++;
           //if (i > 10) break;
       }
        System.out.println(i);
    }
    * */
       vcf.convertFromTriTyper("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/iCHIP_genotypes/TriTyper/", 
               "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/iCHIP_genotypes/TriTyper/SNPs.vcf");
}
}
