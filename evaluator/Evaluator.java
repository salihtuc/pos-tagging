package evaluator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import utils.FileUtils;
import utils.StringCoder;
import utils.StringUtils;


public class Evaluator {
	
	private int[] clusters;
	private int numClusters;
	
	private int[] goldTags;
	private int numGoldTags;
	
	int[][] coocCounts;
	int[] clusterCounts;
	int[] goldTagCounts;
	
	public void setData(int[] clusters, int[] goldTags){
		this.clusters = clusters;
		this.goldTags = goldTags;
		
		//Induce the number of gold tags
		//TODO: put this in the constructor
		numGoldTags = 0;
		List<Integer> seenTags = new ArrayList<Integer>();
		for (int goldTag:goldTags){
			if (!seenTags.contains(goldTag)){
				numGoldTags++;
				seenTags.add(goldTag);
			}
		}
		
		//Use the max cluster index instead
		numClusters = -1;
		for (int cluster:clusters){
			if (cluster > numClusters) numClusters = cluster;
		}
		
		numClusters++;
		
		clusterCounts = new int[numClusters];
		goldTagCounts = new int[numGoldTags];
		coocCounts = new int[numClusters][numGoldTags];
		
		//Make sure that both array have the same number of words
		assert(clusters.length == goldTags.length);
		
		//Initialise arrays
		Arrays.fill(clusterCounts, 0);
		Arrays.fill(goldTagCounts, 0);
		for (int[] cluster:coocCounts) {
			Arrays.fill(cluster, 0);
		}
		
		//Count the co-occurrences of cluster i with goldTag j
		for (int word = 0; word < clusters.length; word++){
			clusterCounts[clusters[word]]++;
			goldTagCounts[goldTags[word]]++;
			coocCounts[clusters[word]][goldTags[word]]++;
		}
		
	}
	
	public double manyToOne(){
		Map<Integer, Integer> manyToOneMap = new HashMap<Integer, Integer>();
		
		//Many-to-one mapping
		for (int cluster = 0; cluster < coocCounts.length; cluster++){
			int mostFreqTag = 0; int mostFreqTagCount = 0;
			//Find the most frequent tag
			for (int goldTag = 0; goldTag < coocCounts[cluster].length; goldTag++){
				if (coocCounts[cluster][goldTag] > mostFreqTagCount){
					mostFreqTag = goldTag;
					mostFreqTagCount = coocCounts[cluster][goldTag];
				}
			}
			manyToOneMap.put(cluster, mostFreqTag);
		}
		
		int correctTags = 0;
		for (int word = 0; word < goldTags.length; word++) {
			if (goldTags[word] == manyToOneMap.get(clusters[word])) correctTags++;
		}
		
		return ((double)correctTags/(double)goldTags.length)*100;
	}
	
	public double oneToOne(){
		Map<Integer, Integer> oneToOneMap = new HashMap<Integer, Integer>();
		
		int bestClust, mostFreqTag;
		int mostFreqTagCount;
		List<Integer> usedClust = new ArrayList<Integer>();
		List<Integer> usedTags = new ArrayList<Integer>();
		//Greedy 1-to-1: Until we run out of either clusters or gold tags 
		for (int i=0; i<numClusters && i<numGoldTags; i++){
			bestClust=0; mostFreqTag=0;
			mostFreqTagCount=-1;
			for (int cluster = 0; cluster < coocCounts.length; cluster++){
				//Do not allow cluster re-use
				if (usedClust.contains(cluster)) continue;
				//Find the most frequent tag
				for (int goldTag = 0; goldTag < coocCounts[cluster].length; goldTag++){
					//Do not allow gold-tag re-use
					if (usedTags.contains(goldTag)) continue;
					if (coocCounts[cluster][goldTag] > mostFreqTagCount){
						mostFreqTagCount = coocCounts[cluster][goldTag];
						mostFreqTag=goldTag;
						bestClust=cluster;
					}
				}
			}
			oneToOneMap.put(bestClust, mostFreqTag);
			usedTags.add(mostFreqTag);
			usedClust.add(bestClust);
		}
		
		int correctTags = 0;
		for (int word = 0; word < goldTags.length; word++) {
			//If the map doesn't contain the key, then word has no tag 
			//(therefore the tagging is wrong)
			if (!oneToOneMap.containsKey(clusters[word])) continue;
			if (goldTags[word] == oneToOneMap.get(clusters[word])) correctTags++;
		}
		
		return ((double)correctTags/(double)goldTags.length)*100;
	}
	
	public double VMeasure(){
		int totalWords = clusters.length;
		double H_T=0, H_CL=0, I=0;
		
		//H(CL): Cluster entropy
		for (int cluster = 0; cluster < numClusters; cluster++){
			double prob = (double)clusterCounts[cluster]/totalWords;
			if(prob!=0){
				H_CL -= prob*Math.log(prob)/Math.log(2);
			}
		}
		//H(T): Tag entropy
		for (int goldTag=0; goldTag < numGoldTags; goldTag++){
			double prob = (double)goldTagCounts[goldTag]/totalWords;
			if(prob!=0){
				H_T -= prob*Math.log(prob)/Math.log(2);
			}
		}
		//I(CL,T): Mutual information
		for (int cluster = 0; cluster < numClusters; cluster++){
			double clProb = (double)clusterCounts[cluster]/totalWords;			
			for (int goldTag=0; goldTag < numGoldTags; goldTag++){
				double tagProb = (double)goldTagCounts[goldTag]/totalWords;
				double prob = 0;
				prob = (double)coocCounts[cluster][goldTag]/totalWords;
				if(prob!=0){
					I += prob*Math.log(prob/(tagProb*clProb))/Math.log(2);			
				}
			}
		}
		
		//H(CL|T): Conditional cluster entropy
		double H_CL_T=H_CL-I;
		//H(T|CL): Conditional tag entropy
		double H_T_CL=H_T-I;
		//h=1-H(CL|T)/H(CL)
		double c=1-(H_CL_T/H_CL);
		//c=1-H(T|CL)/H(T)
		double h=1-(H_T_CL/H_T);
		//V-Measure = (2*h*c)/(h+c)
		double VM=(2*h*c)/(h+c);
		
		return VM*100;
	}
	
	public double VI(){
		int totalWords = clusters.length;
		double H_T=0, H_CL=0, I=0;
		
		//H(CL): Cluster entropy
		for (int cluster = 0; cluster < numClusters; cluster++){
			double prob = (double)clusterCounts[cluster]/totalWords;
			if(prob!=0){
				H_CL -= prob*Math.log(prob)/Math.log(2);
			}
		}
		//H(T): Tag entropy
		for (int goldTag=0; goldTag < numGoldTags; goldTag++){
			double prob = (double)goldTagCounts[goldTag]/totalWords;
			if(prob!=0){
				H_T -= prob*Math.log(prob)/Math.log(2);
			}
		}
		//I(CL,T): Mutual information
		for (int cluster = 0; cluster < numClusters; cluster++){
			double clProb = (double)clusterCounts[cluster]/totalWords;			
			for (int goldTag=0; goldTag < numGoldTags; goldTag++){
				double tagProb = (double)goldTagCounts[goldTag]/totalWords;
				double prob = 0;
				prob = (double)coocCounts[cluster][goldTag]/totalWords;
				if(prob!=0){
					I += prob*Math.log(prob/(tagProb*clProb))/Math.log(2);			
				}
			}
		}
		
		//H(CL|T): Conditional cluster entropy
		double H_CL_T=H_CL-I;
		//H(T|CL): Conditional tag entropy
		double H_T_CL=H_T-I;
		//VI(CL,T) = H(CL|T)+H(T|CL)
		double VI=H_CL_T+H_T_CL;
		
		return VI;
	}
	
	//****** TEST METHODS ******
	public static void main(String[] args){
		//new Evaluator(args[0], args[1]);
		new Evaluator("udEnglish12kBrown.txt","udEnglishTagged12k.txt");
	}
	
	public Evaluator(String fileGold, String fileInduced){
		//Read the gold-tags
		int[] gT = readCorpus(new File(fileGold));
		//Read the clusters
		int[] cl = readCorpus(new File(fileInduced));
		setData(cl, gT);
		
		System.out.println("M-1:\t"+manyToOne());
		System.out.println("1-1:\t"+oneToOne());
		System.out.println("VM:\t"+VMeasure());
		System.out.println("VI:\t"+VI());
	}
	
	private int[] readCorpus(File file){
		FileUtils f = new FileUtils();
		StringUtils s = new StringUtils();
		StringCoder tagsCoder = new StringCoder();
		StringCoder wordsCoder = new StringCoder();
		
		String line;
		int sentenceInd = 0, totalWords = 0;

		//Get the corpus statistics
		try {
			BufferedReader in = f.createIn(file);
			while ((line = in.readLine())!=null){
				sentenceInd++;
				totalWords+=line.split("\\s+").length;
			}
			in.close();
		} 
		catch (IOException e) {
			System.err.println("Error reading the corpus file.");
			System.exit(-1);
		}
		
		int[][] corpus = new int[sentenceInd][];
		int[][] corpusGoldTags = new int[sentenceInd][];
		int[] tags = new int[totalWords];

		sentenceInd = 0;
		//Read the corpus
		try {
			BufferedReader in = f.createIn(file);
			while ((line = in.readLine())!=null){
				//Map word types to integers
				int[] lineWords = new int[line.split("\\s+").length];
				//Map tags to integers
				int[] lineTags = new int[line.split("\\s+").length];
				int wordInd = 0;
				for (String word:line.split("\\s+")){
					String tag = null;
					tag = s.extractTag(word);
					word = s.extractWord(word);
					//... and then get the tag index
					lineTags[wordInd] = tagsCoder.encode(tag);

					//Get the word index
					lineWords[wordInd] = wordsCoder.encode(word);
					wordInd++;
					totalWords++;
				}
				corpus[sentenceInd] = lineWords;
				corpusGoldTags[sentenceInd] = lineTags;
				sentenceInd++;
			}
			in.close();
		}
		catch (IOException e) {
			System.err.println("Error reading the corpus file.");
			System.exit(-1);
		}
		
		//Populate tags
		int wordInd = 0;
		for (int sentence = 0; sentence < corpus.length; sentence++){
			for (int word = 0; word < corpus[sentence].length; word++){
				tags[wordInd] = corpusGoldTags[sentence][word];
				wordInd++;
			}
		}
		return tags;
	}
}
