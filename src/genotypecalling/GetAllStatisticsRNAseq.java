
package genotypecalling;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
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
    Map<String, String> IDconverter;
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
    
    /*
     * checks if the 2 genotypes are equal (Strings AA, BB, AB)
     */
    private boolean equalGenotypes(String g1, String g2){
        HashSet<String> gt1 = new HashSet (Arrays.asList( String.valueOf(g1.charAt(0)), String.valueOf(g1.charAt(1)) ));
        HashSet<String> gt2 = new HashSet (Arrays.asList( String.valueOf(g2.charAt(0)), String.valueOf(g2.charAt(1)) ));
        return gt1.equals(gt2);
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
        String dnaSampleId = sampleId;
        try{
            dnaSampleId = IDconverter.get(sampleId);
        } catch (Exception e){}
        
        TextFile snvmix = new TextFile(fName, false);
        TextFile mpileup = new TextFile(mpileupFile.getPath(), false);
        TextFile out = new TextFile(fName + ".stats", true);
        writeHeader(out);
        
        String[] snvmixEls;
        String[] mpileupEls = mpileup.readLineElems(TextFile.tab);
        
        SNVMix snvmixReader = new SNVMix();
        
        String snpId, rnaGenotype, dnaGenotype, allCov, refCov, altCov;
        boolean isExonic, isConcordant;
        float prob;
        float[] baseQuals;
        int[] coverage;
        SNP rnaSNP, dnaSNP;
        String CR = "NA", MAF = "NA", HWEP = "NA";
        
        Mpileup mp = new Mpileup();
        
        while ((snvmixEls = snvmix.readLineElems(TextFile.tab)) != null){
            snpId = snvmixEls[0].replace("chr", "");
            rnaGenotype = snvmixReader.getSortedGenotype(snvmixEls);
            prob = snvmixReader.getProbability(snvmixEls);
            if (snpId.equals("X:200860")){
                System.out.println(snpId);
            }
            coverage = snvmixReader.getCoverage(snvmixEls);
            allCov = String.valueOf(coverage[2]).replace("NaN", "NA");
            refCov = String.valueOf(coverage[0]).replace("NaN", "NA");
            altCov = String.valueOf(coverage[1]).replace("NaN", "NA");
            
            isExonic = exonicSNPs.contains(snpId);
            
            try{
                rnaSNP = rnaTrityperReader.getSNPByPosId(snpId);
                CR = String.valueOf(rnaSNP.getCR());
                MAF = String.valueOf(rnaSNP.getMAF());
                HWEP = String.valueOf(rnaSNP.getHWEP());
            } catch (Exception e){
                //System.out.println("The SNP " + snpId + " couldn't be found in the TriTyper genotype data");
            }
            dnaSNP = dnaTrityperReader.getSNPByPosId(snpId);
            //isConcordant = rnaTrityperReader.compareAlleles(dnaSNP, dnaTrityperReader.getSamplePos(sampleId), rnaSNP, rnaTrityperReader.getSamplePos(sampleId));
            isConcordant = equalGenotypes(rnaGenotype, dnaTrityperReader.getStringGenotype(dnaSNP, dnaSampleId));
            //find the SNP in mpileup file
            while (!(mpileupEls[0].replace("chr", "") + ":" + mpileupEls[1]).equals(snpId)){
                mpileupEls = mpileup.readLineElems(TextFile.tab);
            }
            
            //coverage = String.valueOf(mp.getRealCoverage(StringUtils.join(mpileupEls, "\t"))); coverage from mpileup file
            
            baseQuals = mp.getQualitiesMpileup(StringUtils.join(mpileupEls, "\t"));
            
            out.writelnTabDelimited(new String[]{
                snpId, 
                rnaGenotype, 
                dnaTrityperReader.getStringGenotype(dnaSNP, dnaSampleId), 
                getNumFromBoolean(isConcordant), 
                String.valueOf(prob),
                allCov,
                refCov,
                altCov,
                CR,
                MAF,
                HWEP,
                getNumFromBoolean(isExonic), 
                String.valueOf(baseQuals[0]),
                String.valueOf(baseQuals[1])});
        }
        
        
        out.close();
        mpileup.close();
        snvmix.close();
    }
    
    /**
     * Makes a HashSet of SNPs that are located in exons
     * in the DNA-seq genotype folder run: awk 'BEGIN {FS="\t"}; {OFS="\t"}; {print $1, $2 - 1, $2, $3}' SNPMappings.txt | intersectBed -a stdin -b exons_v70_withGeneTypes.bed -u > exonicSNPs.bed
     * @param bedFileName - path to a bed file with only exonic SNPs
     * 
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
    
    /**
     * Makes a map for conversion of sample IDs (used if they are different in expression (and RNA-seq genotype) and genotype)
     * @param fName - path to genotype to expression coupling file 
     * @throws IOException 
     */
    private void makeIDconverter(String fName) throws IOException{
        TextFile gte = new TextFile(fName, false);
        IDconverter = gte.readAsHashMap(1, 0);
        gte.close();
    }
    
    private void writeHeader(TextFile out) throws IOException{
        out.writelnTabDelimited(new String[]{
                "SNP_id", 
                "RNA-seq_genotype", 
                "DNA-seq_genotype", 
                "Genotypes_concordant", 
                "SNVMix_probability",
                "Overall_coverage",
                "Ref_coverage",
                "Alt_coverage",
                "CR",
                "MAF",
                "HWEP",
                "SNP_is_exonic", 
                "Average_base_quality_ref",
                "Average_base_quality_alt"});
    }
    
    private String getNumFromBoolean(boolean b){
        if (b == true)
            return "1";
        return "0";
    }
    public static void main(String[] args) throws IOException {
        GetAllStatisticsRNAseq stat = new GetAllStatisticsRNAseq("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/rna-seq/SNVMix-Trityper-noCovThres", 
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/dna-seq/TriTyper/all_chr/");
        stat.makeSetOfExonicSNPs("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/CEU/dna-seq/TriTyper/all_chr/exonicSNPs.bed");
        stat.makeIDconverter("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/CEU/gte_coupling_CEU.txt");
        stat.getSNPStatisticsPerSample("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/mappedData_masked/ERR188032/reads_unique_hits.sorted.mpileup.cov5.snvmix");
    }
    
}
