package genotypecalling;

/**
 * Created with IntelliJ IDEA.
 * User: dashazhernakova
 * Date: 09.09.13
 * Time: 11:58
 * To change this template use File | Settings | File Templates.
 */
public class SharedSample {
	public String id;
	public int numShared;
	public int numConcordant;
	public int numHeterozygous;

	public SharedSample(String sampleId){
		id = sampleId;
		numShared = 0;
		numConcordant = 0;
		numHeterozygous = 0;

	}

	public String toString(){

		return id + " (" + numShared + " " + String.format("%.2f", 100*numConcordant/numShared) + "% " + String.format("%.2f", 100*numHeterozygous/numShared) + "%)";
	}
	public String toFullString(){

		return id + "\t" + numShared + "\t" + numConcordant + "\t" + String.format("%.5f", numConcordant/numShared) + "\t" + numHeterozygous + "\t" + String.format("%.5f", numHeterozygous/numShared);
	}

}
