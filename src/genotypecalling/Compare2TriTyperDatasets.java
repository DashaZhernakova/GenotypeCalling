package genotypecalling;

import org.apache.commons.lang.ArrayUtils;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.RandomAccessGenotypeDataReaderFormats;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.io.text.TextFile;

import java.io.IOException;
import java.util.*;


/**
 *
 * @author Patrick Deelen
 */
public class Compare2TriTyperDatasets {
	RandomAccessGenotypeData dataset1;
	RandomAccessGenotypeData dataset2;
	HashMap<String, Integer> dataset1SharedSampleMap;
	HashMap<String, Integer> dataset2SharedSampleMap;
	ArrayList<SharedSample> sharedSamplesList;
	Map<String, String> gte;

	public Compare2TriTyperDatasets(){}

	public Compare2TriTyperDatasets(String f1, String f2){

		try {
			dataset1 = new TriTyperGenotypeData(f1);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 1: " + f1);

		}

		try {
			System.out.println("Second");
			dataset2 = new TriTyperGenotypeData(f2);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 2: " + f2);

		}

	}

	public Compare2TriTyperDatasets(String f1, String f2, HashSet<String> includeFirstDataset){

		try {
			VariantFilter variantFilter = new VariantIdIncludeFilter(includeFirstDataset);
			//dataset1 = new TriTyperGenotypeData(f1, 1000, variantFilter);
			//RandomAccessGenotypedDataReaderFormats
			dataset1 = RandomAccessGenotypeDataReaderFormats.TRITYPER.createFilteredGenotypeData(f1, 1000, variantFilter, null);
			//RandomAccessGenotypeData genotypeData1 = new TriTyperGenotypeData()
		} catch (IOException ex) {
			System.err.println("Error reading dataset 1: " + f1);

		}

		try {
			System.out.println("Second");
			dataset2 = new TriTyperGenotypeData(f2);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 2: " + f2);

		}

	}

	private void getSharedSamples(String path) throws IOException {
		TextFile gte_coupling = new TextFile(path, false);
		int col;
		gte = gte_coupling.readAsHashMap(0, 1);

		ArrayList<String> samples1 = new ArrayList<String>();
		ArrayList<String> samples2 = new ArrayList<String>();
		for (Sample s : dataset1.getSamples())
			samples1.add(s.getId());
		for (Sample s : dataset2.getSamples())
			samples2.add(s.getId());


		//making a correct gte map
		for(String sample : samples1){
			if (gte.containsKey(sample)){
				break;
			}
			if (gte.containsValue(sample)){
				gte = gte_coupling.readAsHashMap(1, 0);
				break;
			}
		}
		gte_coupling.close();


		//getting shared samples
		HashSet<String> sharedSamples = new HashSet<String>();

		for (String sample : samples1){
			String translation = gte.get(sample);
			if (translation != null){
				if (samples2.contains(translation)){
					sharedSamples.add(sample);
					sharedSamples.add(translation);
				}
			}
		}

		System.out.println("\nNumber of shared samples: " + sharedSamples.size()/2 + "\n");

		dataset1SharedSampleMap = new HashMap<String, Integer>();
		dataset2SharedSampleMap = new HashMap<String, Integer>();

		int i = 0;
		for(Sample sample : dataset1.getSamples()){
			String id = sample.getId();
			if (sharedSamples.contains(id))
				dataset1SharedSampleMap.put(id, i);
			++i;
		}

		i = 0;
		for(Sample sample : dataset2.getSamples()){
			String id = sample.getId();
			if (sharedSamples.contains(id))
				dataset2SharedSampleMap.put(id, i);
			++i;
		}

		sharedSamplesList = new ArrayList<SharedSample>();
		for (String id : dataset1SharedSampleMap.keySet()){
			sharedSamplesList.add(new SharedSample(id));
		}

	}

	private boolean isHeterogygous(Alleles alleles){
		char[] alls = alleles.getAllelesAsChars();

		if (alls[0] == alls[1])
			return true;
		return false;
	}

	private boolean passQC(GeneticVariant snp){
		if ((snp.getCallRate() > 0.3) && (snp.getMinorAlleleFrequency() > 0.05))
			return true;
		return false;
	}


	public void compare(){
		//System.out.println("variantId\tcountIdentical\tcountMismatch");
		//ArrayList<String> sharedSamplesList = new ArrayList<String>(dataset1SharedSampleMap.keySet());
		int numShared = sharedSamplesList.size();
		int[] sharedPerSample = new int[numShared];
		int[] samePerSample = new int[numShared];

		for(GeneticVariant dataset1Variant : dataset1){
			GeneticVariant dataset2Variant = dataset2.getSnpVariantByPos(dataset1Variant.getSequenceName(), dataset1Variant.getStartPos());

			if (dataset2Variant == null){
				//skipping non shared variants
				continue;
			}


			/*if ((! passQC(dataset1Variant)) || (! passQC(dataset2Variant))){
				//skipping SNPs not passing QC'
				continue;
			}
*/
			List<Alleles> dataset1VariantAlleles = dataset1Variant.getSampleVariants();
			List<Alleles> dataset2VariantAlleles = dataset2Variant.getSampleVariants();


			int countIdentical = 0;
			int countMismatch = 0;
			int sample = 0;
			for(SharedSample sharedSample : sharedSamplesList){

				//Get alleles for this shared sample
				Alleles x = dataset1VariantAlleles.get(dataset1SharedSampleMap.get(sharedSample.id));
				Alleles y = dataset2VariantAlleles.get(dataset2SharedSampleMap.get(gte.get(sharedSample.id)));
				if ((! isCalled(x) || (! isCalled(y)))){
					continue;
				}
				//Compare if same alleles. In this case AT == TA
				sharedSample.numShared++;
				//sharedPerSample[sample]++;
				if(x.sameAlleles(y)){
					sharedSample.numConcordant++;
					//samePerSample[sample]++;
				}
				sample++;

			}

			//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + countIdentical + "\t" + countMismatch);

		}
		int num = 0;
		for (SharedSample sample : sharedSamplesList){
			System.out.println(sample);
			num++;
		}
	}

	private int getIntGenotype(char[] alleles){
		int genotype = -1;
		if ((alleles[0] == '0') || (alleles[1] == '0'))
			return -1;
		if (alleles[0] == alleles[1])
			return 1;
		return 2;


	}
	private double calculateCorrelationSharedSamples(List<Alleles> dataset1VariantAlleles, List<Alleles> dataset2VariantAlleles){
		int[] alleles1 = new int[sharedSamplesList.size()];
		int[] alleles2 = new int[sharedSamplesList.size()];
		ArrayList<Integer> allelesList1 = new ArrayList<Integer>();
		ArrayList<Integer> allelesList2= new ArrayList<Integer>();
		int cnt = 0;
		boolean empty = true;
		for(SharedSample sharedSample : sharedSamplesList){

			//Get alleles for this shared sample
			Alleles x = dataset1VariantAlleles.get(dataset1SharedSampleMap.get(sharedSample.id));
			Alleles y = dataset2VariantAlleles.get(dataset2SharedSampleMap.get(gte.get(sharedSample.id)));
			if ((! isCalled(x) || (! isCalled(y)))){
				continue;
			}
			empty = false;

			allelesList1.add(getIntGenotype(x.getAllelesAsChars()));
			allelesList2.add(getIntGenotype(y.getAllelesAsChars()));

			alleles1[cnt] = getIntGenotype(x.getAllelesAsChars());
			alleles2[cnt] = getIntGenotype(y.getAllelesAsChars());
			cnt++;
		}

		if (! empty){
			System.out.print(JSci.maths.ArrayMath.correlation(alleles1, alleles2));
			alleles1 = ArrayUtils.toPrimitive(allelesList1.toArray(new Integer[cnt]));
			alleles2 = ArrayUtils.toPrimitive(allelesList2.toArray(new Integer[cnt]));
			System.out.println("\t" + JSci.maths.ArrayMath.correlation(alleles1, alleles2));
			return JSci.maths.ArrayMath.correlation(alleles1, alleles2);
		}
		else
			return -1;

		//return corr;
	}

	private boolean isCalled(Alleles al){
		if ((al.getAllelesAsChars()[0] == '0') || (al.getAllelesAsChars()[1] == '0'))
			return false;
		return true;
	}

	public void compareHighlyCorrelatedSNPs(){
		//System.out.println("variantId\tcountIdentical\tcountMismatch");
		//ArrayList<String> sharedSamplesList = new ArrayList<String>(dataset1SharedSampleMap.keySet());
		int numShared = sharedSamplesList.size();
		int[] sharedPerSample = new int[numShared];
		int[] samePerSample = new int[numShared];
		int numPass = 0;

		for(GeneticVariant dataset1Variant : dataset1){
			GeneticVariant dataset2Variant = dataset2.getSnpVariantByPos(dataset1Variant.getSequenceName(), dataset1Variant.getStartPos());

			if (dataset2Variant == null){
				//skipping non shared variants
				continue;
			}


			if ((! passQC(dataset1Variant)) || (! passQC(dataset2Variant))){
				//skipping SNPs not passing QC'
				continue;
			}
			numPass++;
			List<Alleles> dataset1VariantAlleles = dataset1Variant.getSampleVariants();
			List<Alleles> dataset2VariantAlleles = dataset2Variant.getSampleVariants();

			int countIdentical = 0;
			int countMismatch = 0;
			int sample = 0;
			int hetero = 0;
			double corr = calculateCorrelationSharedSamples(dataset1VariantAlleles, dataset2VariantAlleles);
			if (corr > 0.9){
				for(SharedSample sharedSample : sharedSamplesList){

					//Get alleles for this shared sample
					Alleles x = dataset1VariantAlleles.get(dataset1SharedSampleMap.get(sharedSample.id));
					Alleles y = dataset2VariantAlleles.get(dataset2SharedSampleMap.get(gte.get(sharedSample.id)));

					if ((! isCalled(x) || (! isCalled(y)))){
						continue;
					}
					//Compare if same alleles. In this case AT == TA
					sharedSample.numShared++;
					//sharedPerSample[sample]++;
					if(x.sameAlleles(y)){
						sharedSample.numConcordant++;
						countIdentical++;
					}
					else{
						countMismatch++;
					}
					if (isHeterogygous(x)){
						sharedSample.numHeterozygous++;
						hetero++;
					}
					sample++;

				}
				//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + corr + "\t" + countIdentical + "\t" + countMismatch + "\t" + hetero);
				//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + countIdentical + "\t" + countMismatch);
			}
		}
		int num = 0;
		for (SharedSample sample : sharedSamplesList){
			System.out.println(sample.toFullString());
			num++;
		}
		System.out.println("number of SNPs " + numPass);
	}

	private void getAllSNPPassQC(){
		for(GeneticVariant var : dataset2){
			if (passQC(var)){
				//System.out.println(var.getPrimaryVariantId());
				GeneticVariant var2 = dataset1.getSnpVariantByPos(var.getSequenceName(), var.getStartPos());
				if (var2 != null){
					System.out.println(var.getSequenceName() +" " + var.getStartPos());
				}
			}
		}
	}
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		/*String f2 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/LL/genotypes/SNVMix-TriTyper/",
				f1="/Users/dashazhernakova/Documents/UMCG/data/BBMRI/LL/genotypes/TriTyper_pos/";

		Compare2TriTyperDatasets compare = new Compare2TriTyperDatasets(f1, f2);
		compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/LL/expression_table/gte_coupling.txt");
		//String f1 = args[1], f2 = args[2];
		*/

		/*String f1 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/CODAM-imputed-20130828-trityper/",
				f2 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/SNVMix-TriTyper/";

		HashSet<String> includedSNPs = new HashSet<String>();
		includedSNPs.add("rs6517457");
		includedSNPs.add("rs1494558");

		Compare2TriTyperDatasets compare = new Compare2TriTyperDatasets(f1, f2, includedSNPs);
		compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/gte_coupling_codam.txt");
		*/
		/*
		String f2 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/SNVMix-TriTyper/",
				f1="/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/CODAM-imputed-20130828-trityper/";

		Compare2TriTyperDatasets compare = new Compare2TriTyperDatasets(f1, f2);
		compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM//gte_coupling_codam.txt");
		 */

		Compare2TriTyperDatasets compare;
		if (args.length < 4){
			String f2 = "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/GBR/rna-seq/SNVMix-TriTyper/",
					f1="/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/GBR/dna-seq/TriTyper/";

			compare = new Compare2TriTyperDatasets(f1, f2);
			compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/GBR/gte_HGtoERR_GBR.txt");
		}
		else{
			compare = new Compare2TriTyperDatasets(args[1], args[2]);
			compare.getSharedSamples(args[3]);
		}
		//compare.getAllSNPPassQC();
		compare.compareHighlyCorrelatedSNPs();



	}
}

