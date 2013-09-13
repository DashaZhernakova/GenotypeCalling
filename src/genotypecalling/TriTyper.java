
package genotypecalling;

import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author dashazhernakova
 */
public class TriTyper {
    TriTyperGenotypeData genotypeData;
    SNPLoader loader;
    HashMap<String, String> pos2id;
    
    public TriTyper(String dirName) throws IOException{
        genotypeData = new TriTyperGenotypeData();
        //genotypeData.load("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/");
        genotypeData.load(dirName);
        loader = genotypeData.createSNPLoader();
        }
    
    public Genotypes readGenotypes(String sampleId) throws IOException{
        return readGenotypes(sampleId, false);
    }
    
    public Genotypes readGenotypes(String sampleId, boolean filterQC) throws IOException{
        Genotypes gen = new Genotypes();
        Integer snpid;
        SNP snp;
        int samplePos, i = 0;
        if (genotypeData.getIndividualToId().containsKey(sampleId))
              samplePos  = genotypeData.getIndividualToId().get(sampleId);
        else{
            System.out.println("No such sample!");
            return null;
        }
        for (String id : genotypeData.getSNPs()){
            snpid = genotypeData.getSnpToSNPId().get(id);
            snp = genotypeData.getSNPObject(snpid);
            loader.loadGenotypes(snp);
            
            if ((snp.getAllele1()[samplePos] == 0) ||(snp.getAllele2()[samplePos] == 0))
                continue;
            if (filterQC){
                if (passesQC(snp)){ //CR > 0.5, MAF > 0.05, HWEP > 0.001
                    i++;
                    if (snp.getChr() != 23)
                        gen.SNP2genotype.put(snp.getChr() + ":" + snp.getChrPos(), new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
                    else
                        gen.SNP2genotype.put("X:" + snp.getChrPos(), new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
            
                }
            }
            else{
                i++;

                //process X chromosome correctly
                if (snp.getChr() != 23)
                    gen.SNP2genotype.put(snp.getChr() + ":" + snp.getChrPos(), new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
                else
                    gen.SNP2genotype.put("X:" + snp.getChrPos(), new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
                //System.out.println(snp.getChr() + ":" + snp.getChrPos() + "\t" + (char)snp.getAllele1()[samplePos] + (char)snp.getAllele2()[samplePos]);
            }
        }
        System.out.println("Loaded " + i + " SNPs");
        return gen;
        
    }
    
    public Genotypes readGenotypesAsRs(String sampleId, boolean filterQC) throws IOException{
        Genotypes gen = new Genotypes();
        Integer snpid;
        SNP snp;
        int samplePos, i = 0;
        if (genotypeData.getIndividualToId().containsKey(sampleId))
              samplePos  = genotypeData.getIndividualToId().get(sampleId);
        else{
            System.out.println("No such sample!");
            return null;
        }
        for (String id : genotypeData.getSNPs()){
            snpid = genotypeData.getSnpToSNPId().get(id);
            snp = genotypeData.getSNPObject(snpid);
            loader.loadGenotypes(snp);
            
            if ((snp.getAllele1()[samplePos] == 0) ||(snp.getAllele2()[samplePos] == 0))
                continue;
            if (filterQC){
                if (passesQC(snp)){ //CR > 0.5, MAF > 0.05, HWEP > 0.001
                    i++;
                    if (snp.getChr() != 23)
                        gen.SNP2genotype.put(id, new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
                    else
                        gen.SNP2genotype.put(id, new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
            
                }
            }
            else{
                i++;

                //process X chromosome correctly
                if (snp.getChr() != 23)
                    gen.SNP2genotype.put(snp.getChr() + ":" + snp.getChrPos(), new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
                else
                    gen.SNP2genotype.put("X:" + snp.getChrPos(), new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]}));
                //System.out.println(snp.getChr() + ":" + snp.getChrPos() + "\t" + (char)snp.getAllele1()[samplePos] + (char)snp.getAllele2()[samplePos]);
            }
        }
        System.out.println("Loaded " + i + " SNPs");
        return gen;
        
    }
    
    public SNP getSNPByPosId(String snpPos) throws IOException{
        String id = pos2id.get(snpPos);
        
        Integer snpid = genotypeData.getSnpToSNPId().get(id);
        SNP snp = genotypeData.getSNPObject(snpid);
        loader.loadGenotypes(snp);
        
        return snp;
    }
    private Map<String, String> converter(String gteFile) throws IOException{
        TextFile tf = new TextFile(gteFile, false);
        Map<String, String> conv = tf.readAsHashMap(1, 0);
        tf.close();
        return conv;
    }
    public void compare(TriTyper t2, String gteFile, String outPath) throws IOException{
        Integer snpid1, snpid2;
        SNP snp2, snp1;
        int samplePos1, samplePos2, i = 0;
        String sampleId2;
        Map<String, String> convert = converter(gteFile);
        TextFile out = new TextFile(outPath, true)        ;
        for (String sampleId1 : genotypeData.getIndividuals()){
            int numEqQC = 0, numDifQC = 0, numShared = 0, numDif = 0, numEq = 0;
            sampleId2 = convert.get(sampleId1);
            if (t2.genotypeData.getIndividualToId().containsKey(sampleId2)){
                  samplePos1  = genotypeData.getIndividualToId().get(sampleId1);
                  samplePos2  = t2.genotypeData.getIndividualToId().get(sampleId2);
            }
            else{
                System.out.println(sampleId1 + " : No such sample!");
                break;
            }
            System.out.println("samples: " + sampleId1 + "; " + sampleId2);
            List<String> snps2 = Arrays.asList(t2.genotypeData.getSNPs());
            for (String id : genotypeData.getSNPs()){
                snpid1 = genotypeData.getSnpToSNPId().get(id);
                snp1 = genotypeData.getSNPObject(snpid1);
                
                String id2 = snp1.getChr() + ":" + snp1.getChrPos();
                
                if (snps2.contains(id2)){
                    snpid2 = t2.genotypeData.getSnpToSNPId().get(id2);
                    snp2 = t2.genotypeData.getSNPObject(snpid2);
                    t2.loader.loadGenotypes(snp2);
                    loader.loadGenotypes(snp1);
                   
                    numShared++;
                    if (numShared % 1000 == 0) System.out.println(numShared);
                    if ((snp1.getAllele1()[samplePos1] != 0) && (snp1.getAllele2()[samplePos1] != 0) && (snp2.getAllele1()[samplePos2] != 0) && (snp2.getAllele2()[samplePos2] != 0)){
                    if ((passesQC(snp1)) && (passesQC(snp2))){
                        if (compareAlleles(snp1, samplePos1, snp2, samplePos2)){
                            //System.out.println("EQUAL\t" + snp1.getChr() + ":" + snp1.getChrPos() + "\t" + new String(new char[] {(char)snp1.getAllele1()[samplePos1], (char)snp1.getAllele2()[samplePos1]})
                              //      + "\t" + new String(new char[] {(char)snp2.getAllele1()[samplePos2], (char)snp2.getAllele2()[samplePos2]}));
                            numEqQC++;
                        }
                        else{
                            //System.out.println("NOT EQUAL\t" + snp1.getChr() + ":" + snp1.getChrPos() + "\t" + new String(new char[] {(char)snp1.getAllele1()[samplePos1], (char)snp1.getAllele2()[samplePos1]})
                              //      + "\t" + new String(new char[] {(char)snp2.getAllele1()[samplePos2], (char)snp2.getAllele2()[samplePos2]}));
                            numDifQC++;
                        }
                    }
                    else{
                        if (compareAlleles(snp1, samplePos1, snp2, samplePos2)){
                            //System.out.println("EQUAL\t" + snp1.getChr() + ":" + snp1.getChrPos() + "\t" + new String(new char[] {(char)snp1.getAllele1()[samplePos1], (char)snp1.getAllele2()[samplePos1]})
                              //      + "\t" + new String(new char[] {(char)snp2.getAllele1()[samplePos2], (char)snp2.getAllele2()[samplePos2]}));
                            numEq++;
                        }
                        else{
                            //System.out.println("NOT EQUAL\t" + snp1.getChr() + ":" + snp1.getChrPos() + "\t" + new String(new char[] {(char)snp1.getAllele1()[samplePos1], (char)snp1.getAllele2()[samplePos1]})
                              //      + "\t" + new String(new char[] {(char)snp2.getAllele1()[samplePos2], (char)snp2.getAllele2()[samplePos2]}));
                            numDif++;
                        }
                        }
                    }
                    
                }
            }
            out.writeln(sampleId1 + "\t" + numShared + "\t" + (int) (numEq + numDif) + "\t" + numEq + "\t" + numDif + "\t" + (int)(numEqQC + numDifQC) + "\t" + numEqQC + "\t" + numDifQC);
            //break;
        }
        out.close();
    }
    
    public boolean compareAlleles(SNP snp1, int samplePos1, SNP snp2, int samplePos2){
        if (((snp1.getAllele1()[samplePos1] == snp2.getAllele1()[samplePos2]) && (snp1.getAllele2()[samplePos1] == snp2.getAllele2()[samplePos2])) || 
                                ((snp1.getAllele1()[samplePos1] == snp2.getAllele2()[samplePos2]) && (snp1.getAllele2()[samplePos1] == snp2.getAllele1()[samplePos2]))){
            
            return true;
        }
        return false;
    }
    public boolean passesQC (SNP snp){
        if ( (snp.getCR() > 0.5) && (snp.getMAF() > 0.05) && (snp.getHWEP() > 0.001))
            return true;
        return false;
    }
    public void test() throws IOException{
        Genotypes gen = new Genotypes();
        Integer snpid;
        SNP snp;
        for (String id : genotypeData.getSNPs()){
            snpid = genotypeData.getSnpToSNPId().get(id);
            snp = genotypeData.getSNPObject(snpid);
            loader.loadGenotypes(snp);
            System.out.println(snp);
        }
    } 
    
    public void makePosIdToRsMap(){
        pos2id = new HashMap<String, String>();
        Integer snpid;
        SNP snp;
                
        for (String id : genotypeData.getSNPs()){
            snpid = genotypeData.getSnpToSNPId().get(id);
            snp = genotypeData.getSNPObject(snpid);
            if (snp.getChr() == 23){
                pos2id.put("X:" + snp.getChrPos(), id);
            }
            else
                pos2id.put(snp.getChr() + ":" + snp.getChrPos(), id);
        }
          
    }
    
    public int getSamplePos(String sampleId){
        int samplePos = -1;
        if (genotypeData.getIndividualToId().containsKey(sampleId))
              samplePos  = genotypeData.getIndividualToId().get(sampleId);
        else{
            System.out.println("No such sample!");
            
        }
        return samplePos;
    }
    
    public String getStringGenotype(SNP snp, String sampleId){
        int samplePos = getSamplePos(sampleId);
        return new String(new char[] {(char)snp.getAllele1()[samplePos], (char)snp.getAllele2()[samplePos]});
    }
    
    
    public static void main(String[] args) throws IOException {
        TriTyper r = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/SNVMix-TriTyper/");
        TriTyper d = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/dna-seq/TriTyper_pos");

        r.compare(d,
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/FIN/gte_HGtoERR_FIN.txt",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/rna-seq/SNVMix-TriTyper/comp2dnaseq.txt");
        //TriTyper t2 = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/dna-seq/TriTyper");
        //t2.readGenotypes("HG00360");
        //TriTyper t = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/RNA-seq/SNVMix/");
        //TriTyper t2 = new TriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/RNA-seq/SNVMix/");
        //t.test();
        //t2.compare(t);
        
    }
        
}
