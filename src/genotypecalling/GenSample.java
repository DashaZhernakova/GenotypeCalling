
package genotypecalling;

import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

/**
 * WRONG CHROMOSOME NAMES
 * Impute2 genotype format
 * @author dashazhernakova
 */
public class GenSample {
    HashMap<String, String> SNPIdConverter;
    HashSet<String> SNPsToInclude;
    HashSet<String> SNPsBeforeImputation;

    public Genotypes readGenotypes(String genFname, String sampleFname, String sample) throws IOException {
        Genotypes gen = new Genotypes();

        // Get sample column number
        TextFile sampleFile = new TextFile(sampleFname, false);
        String[] els = sampleFile.readLineElems(TextFile.space);
        els = sampleFile.readLineElems(TextFile.space);
        int sampleCol = -1, i = 0;

        while ((els = sampleFile.readLineElems(TextFile.space)) != null){
            if (els[0].equals(sample)){
                sampleCol = i;
                break;
            }
            i++;
        }
        sampleFile.close();
        if (sampleCol == -1){
            System.out.println("No such sample! " + sample);
            return null;
        }

        // Read genotypes for the sample
        TextFile genFile = new TextFile(genFname, false);
        String genString = null;
        while ((els = genFile.readLineElems(TextFile.space)) != null){
            genString = getGenotypeString(els, sampleCol);
            if (genString != null)
                gen.SNP2genotype.put(els[1], genString);
        }
        genFile.close();

        return gen;
    }

    public void compareWithSNVMix(String filePrefix, String sampleId, TreeMap<String, String> dnaGenotype, String chr) throws IOException{
        TextFile gen = new TextFile(filePrefix + ".gen", false);
        TextFile sample = new TextFile(filePrefix + ".sample", false);

        //Get sample column name
        String[] els;
        int sampleCol = 0;
        //two header lines
        sample.readLineElems(TextFile.space);
        sample.readLineElems(TextFile.space);
        
        while((els = sample.readLineElems(TextFile.space)) != null){
            if (els[0].equals(sampleId))
                break;
            sampleCol++;
        }
        sample.close();
        
        //make the conversion file
        //makeSNPIdConverter("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/dna-seq/VCF/ALL.chr10.phase1_release_v3.20101123.snps_indels_svs.genotypes_filtered_common_snps_FIN_samp.recode.vcf.gz", getSNPSet(filePrefix + ".gen"));
        
        Genotypes g = new Genotypes();
        //read genotypes
        String snpId = "", dnaGen, rnaGen;
        int eq = 0, neq = 0, eq_b = 0, neq_b = 0;
        while((els = gen.readLineElems(TextFile.space)) != null){
            snpId = els[1];
            if (snpId.startsWith("rs")){
                snpId = chr + ":" + els[2];
            }
            else{
                snpId = snpId.replace("chr", "");
                /*try{
                    snpId = SNPIdConverter.get(snpId.replace("chr", ""));
                } catch (Exception e) {
                    //System.out.println("No rs SNP id for SNP: " + snpId);
                }*/
            }
            if (! SNPsToInclude.contains(snpId))
                continue;
            if (dnaGenotype.containsKey(snpId)){
                dnaGen = dnaGenotype.get(snpId);
                rnaGen = getGenotypeString(els, sampleCol);
                if (rnaGen == null)
                    continue;
                if (g.equalGenotypes(dnaGen, rnaGen)){
                    eq++;
                    //if (SNPsBeforeImputation.contains(snpId))
                      //  eq_b++;
                    //System.out.println("EQUAL " + snpId);
                }
                else{
                    neq++;
                    //if (SNPsBeforeImputation.contains(snpId))
                      //  neq_b++;
                    //System.out.println(snpId + "\t" + dnaGen + "\t" + rnaGen + "\t" + els[3*sampleCol + 5] + ", " + els[3*sampleCol + 6] + "," + els[3*sampleCol + 7]);
                }
            }
            
            
        }
        gen.close();
        int sum = eq + neq;
        System.out.println("Overall: " + sum + " Equal: " + eq + " (" + 100*eq/(eq+neq) + "%)\nNot equal " + neq + " (" + 100*neq/(eq+neq) + "%)");
        System.out.println("Out of those: equal in initial: " + eq_b + " not equal: " + neq_b);
    }
    /**
     * Gets a String with individual genotype according to the max probability
     * @param els - SNP line of gen file
     * @param sampleCol - sample position number
     */
    private String getGenotypeString(String[] els, int sampleCol){
        String gen = null;
        float[] genPr = new float [3];

        genPr[0] = Float.parseFloat(els[3*sampleCol + 5]);
        genPr[1] = Float.parseFloat(els[3*sampleCol + 6]);
        genPr[2] = Float.parseFloat(els[3*sampleCol + 7]);
        float maxPr = Math.max(Math.max(genPr[0], genPr[1]), genPr[2]);
        
        if (maxPr == 0){
            return null;
        }
        if (maxPr == genPr[0])
            gen = els[3] + els[3];
        else if (maxPr == genPr[1])
            gen = els[3] + els[4];
        else if (maxPr == genPr[2])
            gen = els[4] + els[4];
        
        return gen;
    }
    
    private HashSet<String> getSNPSet(String genFname) throws IOException{
        TextFile gen = new TextFile(genFname, false);
        HashSet<String> snpSet = new HashSet<String>();
        String[] els;
        
        while ((els = gen.readLineElems(TextFile.space)) != null){
            if (! els[1].startsWith("rs")){
                snpSet.add(els[1].replace("chr", ""));
            }
        }
        
        gen.close();
        return snpSet;
    }
    
    public void getSNPsPassQC(String dirName, float infoThresh) throws IOException{
        File dir = new File(dirName);
        SNPsToInclude = new HashSet<String>();
        
        for (File child : dir.listFiles()){
            if ((child.getName().endsWith("_info"))){
                TextFile info = new TextFile(child.getPath(), false);
                String[] els = info.readLineElems(TextFile.space);
                float inf;
                while( (els = info.readLineElems(TextFile.space)) != null ){
                    inf = Float.parseFloat(els[4]);
                    if (inf > infoThresh)
                        SNPsToInclude.add("10:" + els[2]);
                    
                }
                info.close();
            }
        }
            
    }
    
    public void getSNPsBeforeImputation(String path, int sampleCol) throws IOException{
        SNPsBeforeImputation = new HashSet<String>();
        TextFile notImputed = new TextFile(path, false);
        
        String[] els;
        while ((els = notImputed.readLineElems(TextFile.space)) != null){
            if (!( (els[5 + sampleCol].equals("0")) && (els[6 + sampleCol].equals("0")) && (els[7 + sampleCol].equals("0")))){
                SNPsBeforeImputation.add(els[1]);
            }
            
        }
        
        notImputed.close();
    }
    /*
     * Makes a HashMap to use when converting snvmix position ids (chr:pos) to rs ids
     * conversionFile - path to the VCF file
     */
    private void makeSNPIdConverter(String conversionVCFFile, HashSet<String> snpSet) throws IOException {
        TextFile tf = new TextFile(conversionVCFFile, false);
        String[] els;
        SNPIdConverter = new HashMap<String, String>();
                
        while ((els = tf.readLineElems(TextFile.tab))[0].startsWith("#"))
            continue;
        while ((els = tf.readLineElems(TextFile.tab)) != null){
            if (snpSet.contains(els[0] + ":" + els[1]))
                SNPIdConverter.put(els[0] + ":" + els[1], els[2]);
        }
        tf.close();
    }
    
    public static void main(String[] args) throws IOException {
        
        GenSample g = new GenSample();
        Genotypes gen = new Genotypes();
        gen = g.readGenotypes("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp/test1.gen",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp/test.sample", "ERR188032");
        System.out.println(gen.SNP2genotype.get("1:100317394"));
       }
}
