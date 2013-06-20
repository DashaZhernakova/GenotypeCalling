/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package genotypecalling;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import umcg.genetica.io.trityper.converters.TriTyperToVCF;

/**
 *
 * @author dashazhernakova
 */
public class GenotypeCalling {

    /**
     * @param args the command line arguments
     */
    public static void usage(){
        System.out.println("Input arguments:\n--mode");
        System.out.println("\tSNVMixToTriTyper - converter from SNVMix files per sample to joint TriTyper format");
        System.out.println("\t\t--in - path to the folder containing folders with snvmix genotype files");
        System.out.println("\t\t--out - output folder");
        System.out.println("\t\t--pattern - SNVMix file name pattern");
        System.out.println("\t\t--toRsIds - path to a VCF file to convert from snvmix chr:pos ids to rs ids");
        System.out.println("\t\t--p-value - SNVMix p-value threshold to use. If genotype probability is lower than this value, this SNP is skipped. (Default: 0.95)");
        
        System.out.println("\tSNVMixToGen - converter from SNVMix files per sample to joint Gen + sample format for Impute2 and ShapeIt");
        System.out.println("\t\t--in - path to the folder containing folders with snvmix genotype files");
        System.out.println("\t\t--out - output prefix");
        System.out.println("\t\t--pattern - SNVMix file name pattern");
        System.out.println("\t\t--toRsIds - path to a VCF file to convert from snvmix chr:pos ids to rs ids");
        System.out.println("\t\t--p-value - SNVMix p-value threshold to use. If genotype probability is lower than this value, this SNP is skipped. (Default: 0.8)");
        
        System.out.println("\tcorrectAlleles - tool to change the wrong alleles called from RNA-seq to reference panel alleles. The sample alleles are NOT checked, only the alleles are changed to reference panel alleles");
        System.out.println("\t\t--genDir - path to the folder with .gen genotype files");
        System.out.println("\t\t--legendDir - path to the folder with .legend files");
        System.out.println("\t\t--genPattern - .gen file name pattern (after chr name)");
        System.out.println("\t\t--legendPattern - .legend file name pattern (after chr name)");
        
    }
    public static void main(String[] args) {
        /*TriTyperToVCF c = new TriTyperToVCF();
        try {
            c.convert("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/",
                    "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/vcf/",
                    "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/newGenotypes/SNPlist_chr7.txt");
        } catch (IOException ex) {
            Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
        }
        */
        SNVMixToTriTyperConverter converter;
        SNVMixToGenConverter convToGen;
        String mode = "";
        int i = 0;
        for (i = 0; i < args.length; i++) {
	    String arg = args[i];
	    String val = null;

	    if (i + 1 < args.length) {
		val = args[i + 1];
	    }

	    if (arg.equals("--mode")) {
		mode = val;
                break;
	    }
            
	}
        
        if (mode == null){ 
	    System.out.println("ERROR: Please supply --mode");
            usage();
        }
        
        else if(mode.equals("SNVMixToTriTyper")){
            converter = new SNVMixToTriTyperConverter();
            String in = null, out = null, pattern = null, vcfForConversion = null, p_val = null;
            boolean dosage = false;
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;
                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) {
                    in = val;
                }
                if (arg.equals("--out")) {
                    out = val;
                }
                if (arg.equals("--pattern")) {
                    pattern = val;
                }
                if (arg.equals("--p-value")) {
                    p_val = val;
                }
                if (arg.equals("--dosage")) {
                    dosage = Boolean.parseBoolean(val);
                }
                if (arg.equals("--toRsIds")) {
                    vcfForConversion = val;
                    if ((! vcfForConversion.endsWith(".vcf")) && (! vcfForConversion.endsWith(".vcf.gz"))){
                        System.out.println("To convert the ids to rs ids a VCF file is needed. Check you file format and extension (needs to be .vcf)!\nExiting.");
                        return;
                    }
                }
            }
            try {
                if (p_val != null)
                    converter.setPthreshold(Float.valueOf(p_val));
                converter.parse(in, out, pattern, vcfForConversion, dosage);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("Input arguments :\nin = " + in +"\nout = " + out + "\npattern = " + pattern + "\ntoRsIds = " + vcfForConversion);
                usage();
            }
            
        }
        else if(mode.equals("SNVMixToGen")){
            convToGen = new SNVMixToGenConverter();
            String in = null, out = null, pattern = null, vcfForConversion = null, p_val = null;
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;
                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) {
                    in = val;
                }
                if (arg.equals("--out")) {
                    out = val;
                }
                if (arg.equals("--pattern")) {
                    pattern = val;
                }
                if (arg.equals("--p-value")) {
                    p_val = val;
                }
                if (arg.equals("--toRsIds")) {
                    vcfForConversion = val;
                    if (! vcfForConversion.endsWith(".vcf")){
                        System.out.println("To convert the ids to rs ids a VCF file is needed. Check you file format and extension (needs to be .vcf)!\nExiting.");
                        return;
                    }
                }
            }
            try {
                if (p_val != null)
                    convToGen.setPthreshold(Float.valueOf(p_val));
                convToGen.parse(in, out, pattern, vcfForConversion);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("Input arguments :\nin = " + in +"\nout = " + out + "\npattern = " + pattern + "\ntoRsIds = " + vcfForConversion);
                usage();
            }
            
        }
        else if(mode.equals("correctAlleles")){
            ChangeAllelesToConcordant c = new ChangeAllelesToConcordant();
            String genDir = null, legendDir = null, genPattern = null, legendPattern = null, p_val = null;
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;
                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--genDir")) {
                    genDir = val;
                }
                if (arg.equals("--legendDir")) {
                    legendDir = val;
                }
                if (arg.equals("--genPattern")) {
                    genPattern = val;
                }
                if (arg.equals("--legendPattern")) {
                    legendPattern = val;
                }
                
            }
            try {
                System.out.println("Input arguments :\ngenDir = " + genDir +"\nlegendDir = " + legendDir + "\ngenPattern = " + genPattern + "\nlegendPattern = " + legendPattern);
                c.change(genDir, legendDir, genPattern, legendPattern);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
                System.out.println("Input arguments :\ngenDir = " + genDir +"\nlegendDir = " + legendDir + "\ngenPattern = " + genPattern + "\nlegendPattern = " + legendPattern);
                usage();
            }
            
        }
        else
            usage();
        /*
        else if (mode.equals("comparePileup2pedMap")) {
        
            String snps = null, mappings = null, pedMapDir = null, sample = null;
            float thres = 0.3f;
            PedMap pedmap = new PedMap();
            Pileup pileup = new Pileup();
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--genotypes")) {
                    snps = val;
                }
                if (arg.equals("--SNPmappings")) {
                    mappings = val;
                }
                if (arg.equals("--pedMap")) {
                    pedMapDir = val;
                }
                if (arg.equals("--sample")) {
                    sample = val;
                }
                if (arg.equals("--threshold")) {
                    sample = val;
                }
            }
            System.out.println("Reading genotypes from pileup file " + snps);
            System.out.println("SNP mappings file: " + mappings);
            try {
                pileup.readGenotypes(snps, thres);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
            }
            System.out.println("reading genotypes from ped and map files " + pedMapDir);
            System.out.println("Comparing genotype calls for sample " + sample);
            try {
                pedmap.readAndCompare(pedMapDir, sample, pileup.SNP2genotype);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        else if (mode.equals("comparePileup2SNVMix")) {
        
            String snps = null, mappings = null, pedMapDir = null, sample = null;
            PedMap pedmap = new PedMap();
            SNVMix snvmix = new SNVMix();
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--genotypes")) {
                    snps = val;
                }
                if (arg.equals("--SNPmappings")) {
                    mappings = val;
                }
                if (arg.equals("--pedMap")) {
                    pedMapDir = val;
                }
                if (arg.equals("--sample")) {
                    sample = val;
                }
                
            }
            System.out.println("Reading genotypes from pileup file " + snps);
            System.out.println("SNP mappings file: " + mappings);
            try {
                snvmix.readGenotypes(snps);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
            }
            System.out.println("reading genotypes from ped and map files " + pedMapDir);
            System.out.println("Comparing genotype calls for sample " + sample);
            try {
                pedmap.readAndCompare(pedMapDir, sample, snvmix.SNP2genotype);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        else if (mode.equals("TriTyperToVCF")){
            TriTyperToVCF converter = new TriTyperToVCF();
            String in = null, out = null, list = null;
            for (int j = i; j < args.length; j++) {
                String arg = args[j];
                String val = null;

                if (j + 1 < args.length) {
                    val = args[j + 1];
                }
                if (arg.equals("--in")) {
                    in = val;
                }
                if (arg.equals("--out")) {
                    out = val;
                }
                if (arg.equals("--list")) {
                    list = val;
                }
            }
            try {
                converter.convert(in, out, list);
            } catch (IOException ex) {
                Logger.getLogger(GenotypeCalling.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        * */
    }
}
