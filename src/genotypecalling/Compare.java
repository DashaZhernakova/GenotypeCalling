package genotypecalling;


import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.HashMap;

public class Compare {
    public void compareSNVMixAndTriTyper(String trityper, String snvmix_dir, String snvmix_fname, Float prob, String gte_file, String out) throws IOException {

        TriTyper d = new TriTyper(trityper);
        SNVMix r = new SNVMix();

        TextFile gte = new TextFile(gte_file, false);
        HashMap<String, String> conv = new HashMap<String, String>(gte.readAsHashMap(0,1));
        Genotypes dg, rg;
        for (String dnasample : conv.keySet()){
            try{
                dg = d.readGenotypes(dnasample, true);
                rg = r.readGenotypesFilterProbability(snvmix_dir + "/" + conv.get(dnasample) + "/" + snvmix_fname, prob);
                rg.compare(dg, out, dnasample);
            } catch (IOException e){
                e.printStackTrace();
                System.out.println("No genotype: " + dnasample + " : " + conv.get(dnasample));
            }

        }

    }

    public void compareGenAndTriTyper(String trityper, String fList, String sampleFile, String gte_file, String out) throws IOException {
        GenSample r = new GenSample();
        TriTyper d = new TriTyper(trityper);

        TextFile gte = new TextFile(gte_file, false);
        HashMap<String, String> conv = new HashMap<String, String>(gte.readAsHashMap(0,1));
        Genotypes dg = null, rg = null;

        TextFile fileList = new TextFile(fList, false);
        String[] genFiles = fileList.readAsArray(0, TextFile.tab);
        fileList.close();

        for (String dnasample : conv.keySet()){
            System.out.println(dnasample);
            for (String genFile : genFiles){
                String name = genFile.split("/")[genFile.split("/").length - 1];
                try{
                      rg = r.readGenotypes(genFile, sampleFile, conv.get(dnasample));
                    if ( (rg.SNP2genotype.containsKey("18:678947")) || (rg.SNP2genotype.containsKey("18:675903")) || (rg.SNP2genotype.containsKey("18:685797")) ){
                        System.out.println("c");
                    }
                      if (rg != null)
                          dg = d.readGenotypes(dnasample, true);
                    if ( (dg.SNP2genotype.containsKey("18:678947")) || (dg.SNP2genotype.containsKey("18:675903")) || (dg.SNP2genotype.containsKey("18:685797")) ){
                        System.out.println("c");
                    }
                      if ((rg != null) && (dg != null))
                        dg.compare(rg, out, dnasample + " (" + name + ")");
                  } catch (IOException e){
                      e.printStackTrace();
                      System.out.println("No genotype: " + dnasample + " : " + conv.get(dnasample));
                  }
              }
        }
    }
    public static void main(String[] args) throws IOException {
        Compare c = new Compare();
        c.compareGenAndTriTyper("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/YRI/dna-seq/TriTyper_pos/",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp.txt",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/YRI/rna-seq/gen/yri_0.8.sample",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/YRI/gte_YRI.txt",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/tmp.comp.txt");
    }
}
