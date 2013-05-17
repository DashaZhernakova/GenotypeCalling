
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang.StringUtils;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;

/**
 *
 * @author dashazhernakova
 */
public class GetAllStatisticsRNAseq {
    HashSet<String> exonicSNPs;
    TriTyper rnaTrityperReader;
    TriTyper dnaTrityperReader;
    public GetAllStatisticsRNAseq(String rnaTrityperDir, String dnaTrityperDir) throws IOException{
        rnaTrityperReader = new TriTyper(rnaTrityperDir);
        rnaTrityperReader.makePosIdToRsMap();
        
        dnaTrityperReader = new TriTyper(dnaTrityperDir);
        dnaTrityperReader.makePosIdToRsMap();
        
        
    }
    /**
     * Launches SNP statistics gatherer (getSNPStatisticsPerSample) function for each snvmix file in the sample folders from dirName
     * @param dirName - folder with sample subfolders 
     * @param fnamePattern - regex to match the snvmix output
     */
    public void runForEachSample(String dirName, String fnamePattern){
        System.out.println("Getting all files with reads per transcripts counts. Files names finish with " + fnamePattern);
        File dir = new File(dirName);
        for (File ch : dir.listFiles()) {
            if (ch.isDirectory())
                for (File child : ch.listFiles()){
                   if (child.getName().matches(fnamePattern)){
                        try {
                            getSNPStatisticsPerSample(child.getPath());
                        } catch (IOException ex) {
                            Logger.getLogger(GetAllStatisticsRNAseq.class.getName()).log(Level.SEVERE, null, ex);
                        }
                  }
            }
        }
    }
    
    /**
     * Gathers information about all SNPs called: coverage, exon/intron, average base quality, genotype probability from SNVMix, CR, MAF, HWEP, etc
     * @param fName - path to snvmix file (mpileup file name is the same but without "snvmix" extention)
     */
    public void getSNPStatisticsPerSample(String fName) throws IOException{
        File mpileupFile = new File(fName.replace(".snvmix", "").replace(".cov5", ""));
        if (! mpileupFile.exists()){
            mpileupFile = new File(fName.replace(".snvmix", ".gz"));
            if (! mpileupFile.exists()){
                System.out.println("Could not find mpileup file for " + fName + ". \nExiting!");
                return;
            }
        }
        
        String sampleId = fName.split("/")[fName.split("/").length - 2];
        
        TextFile snvmix = new TextFile(fName, false);
        TextFile mpileup = new TextFile(mpileupFile.getPath(), false);
        TextFile out = new TextFile(fName + ".stats", true);
        writeHeader(out);
        
        String[] snvmixEls;
        String[] mpileupEls = mpileup.readLineElems(TextFile.tab);
        
        SNVMix snvmixReader = new SNVMix();
        
        String snpId, rnaGenotype, dnaGenotype, coverage;
        boolean isExonic, isConcordant;
        float prob;
        float[] baseQuals;
        SNP rnaSNP, dnaSNP;
        
        BaseQuality bq = new BaseQuality();
        
        while ((snvmixEls = snvmix.readLineElems(TextFile.tab)) != null){
            snpId = snvmixEls[0].replace("chr", "");
            rnaGenotype = snvmixReader.getSortedGenotype(snvmixEls);
            prob = snvmixReader.getProbability(snvmixEls);
            isExonic = exonicSNPs.contains(snpId);
            
            rnaSNP = rnaTrityperReader.getSNPByPosId(snpId);
            dnaSNP = dnaTrityperReader.getSNPByPosId(snpId);
            isConcordant = rnaTrityperReader.compareAlleles(dnaSNP, dnaTrityperReader.getSamplePos(sampleId), rnaSNP, rnaTrityperReader.getSamplePos(sampleId));
            
            //find the SNP in mpileup file
            while ((mpileupEls[0].replace("chr", "") + ":" + mpileupEls[1]) != snpId){
                mpileupEls = mpileup.readLineElems(TextFile.tab);
            }
            
            coverage = mpileupEls[3];
            baseQuals = bq.getQualitiesMpileup(StringUtils.join(mpileupEls, "\t"));
            
            out.writelnTabDelimited(new String[]{
                snpId, 
                rnaGenotype, 
                dnaTrityperReader.getStringGenotype(dnaSNP, sampleId), 
                String.valueOf(isConcordant), 
                String.valueOf(prob),
                coverage,
                String.valueOf(rnaSNP.getCR()),
                String.valueOf(rnaSNP.getMAF()),
                String.valueOf(rnaSNP.getHWEP()),
                String.valueOf(isExonic), 
                String.valueOf(baseQuals[0]),
                String.valueOf(baseQuals[1])});
        }
        
        
        out.close();
        mpileup.close();
        snvmix.close();
    }
    
    /**
     * Makes a HashSet of SNPs that are located in exons
     * @param bedFileName - path to a bed file with only exonic SNPs
     */
    private void makeSetOfExonicSNPs(String bedFileName) throws IOException{
        TextFile bed = new TextFile(bedFileName, false);
        exonicSNPs = new HashSet<String>();
        
        String[] els;
        
        while ((els = bed.readLineElems(TextFile.tab)) != null){
            exonicSNPs.add(els[0] + ":" + els[2]);
        }
        
        bed.close();
    }
    
    private void writeHeader(TextFile out) throws IOException{
        out.writelnTabDelimited(new String[]{
                "SNP_id", 
                "RNA-seq_genotype", 
                "DNA-seq_genotype", 
                "Genotypes_concordant", 
                "SNVMix_probability",
                "Coverage",
                "CR",
                "MAF",
                "HWEP",
                "SNP_is_exonic", 
                "Average_base_quality_ref",
                "Average_base_quality_alt"});
    }
    
    public static void main(String[] args) throws IOException {
        GetAllStatisticsRNAseq stat = new GetAllStatisticsRNAseq("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/rna-seq/SNVMix-Trityper-noCovThres", 
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/dna-seq/TriTyper/all_chr/");
        
        stat.getSNPStatisticsPerSample("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked/ERR188032/reads_unique_hits.sorted.mpileup.cov5.snvmix");
    }
    
}
