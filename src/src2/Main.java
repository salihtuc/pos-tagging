
/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import lbfgsb.LBFGSBException;
//import lbfgsb.Minimizer;
import lbfgsb.Result;

public class Main {

	public static PrintWriter pw = null;
	
	public static int tagSize = 13;
	public static ArrayList<String> sentences = new ArrayList<>();
	public static ArrayList<String> words = new ArrayList<>();
	
	public static ArrayList<String> allWords = new ArrayList<>();
	public static HashSet<String> uniqueValues;
	
	public static HashMap<Integer, HashMap<Integer, Node>> globalMap;	// Holds all lattices
	
	public static ArrayList<HashMap<Integer, HashMap<Integer, Node>>> sentenceGlobals = new ArrayList<>();
	
	
	public static HashMap<String, Double> initialProbabilities = new HashMap<>();
	
	public static HashMap<String, Double> transitionProbabilities = new HashMap<>();
	public static HashMap<String, Double> emissionProbabilities = new HashMap<>();
	
	public static ArrayList<String> tagList = new ArrayList<>();	// Holds tags
	public static ArrayList<String> featureList = new ArrayList<>();	// Holds features
	
	public static HashMap<String, Integer> tagFeature2Index = new HashMap<>();
	public static HashMap<String, Integer> gradFeature2Index = new HashMap<>();
	
	public static double[] tagFeatureWeights = null;	// Feature sized weight array that we use in LBFGS-B
	public static double[] tagFeatureGradients = null;

	public static void main(String[] args) {
		
		// Time operations. Just using for information.
		long startTime = System.nanoTime();
		try {
			pw = new PrintWriter(new File("outp.txt"));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		
		//try (BufferedReader br = new BufferedReader(new FileReader("PennCorpus12.txt"))) {
		try (BufferedReader br = new BufferedReader(new FileReader("input.txt"))) {
			String start = "<start> ";
			String line;

			while ((line = br.readLine()) != null) {
				
				if(line.split(" ").length > 2) {
					sentences.add(start + line.toLowerCase());
					
					allWords.addAll(Arrays.asList((start + line.toLowerCase()).split(" ")));
				}
				
			}

			allWords.add("<end>");
			uniqueValues = new HashSet<>(allWords);
			allWords.clear();
			allWords.addAll(uniqueValues);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//String sentence = "<s> Natural language is a delicate thing";
		
		
		// Initialization of features and tags.
		
		fillTagList();	// Creating tags (filling tagList)
		fillTransitionMap();	// Creating transitionProbabilities map
		
		for(String sentence : sentences)
			fillEmissionMap(sentence);	// Creating emissionProbabilities map
		
		fillEmissionMap("<end>");
		fillInitials("initialProb.txt");
		
//		System.out.println(transitionProbabilities);
		
		tagFeatureWeights = createWeightsArray(transitionProbabilities, emissionProbabilities, gradFeature2Index);
//		System.out.println(transitionProbabilities);
//		System.out.println(emissionProbabilities);
//		System.out.println(transitionProbabilities.size());
//		System.out.println(emissionProbabilities.size());
//		System.out.println(tagFeatureWeights.length);
//		tagFeatureGradients = createZeroArray(tagFeatureWeights.length);
//		fillFeatures(transitionProbabilities, emissionProbabilities);
		
		
//		System.out.println("Transitions:\n" + transitionProbabilities);
//		System.out.println("------------\nEmissions:\n" + emissionProbabilities + "\n-----------");
		
		for(String sentence : sentences) {
			
			globalMap = new HashMap<>();
			// Create all lattices and put them into the globalMap
			HashMap<Integer, Node> latticeMap = new HashMap<>();	// Original sentence
			fillLattice(sentence.split(" "), latticeMap, 0);
			globalMap.put(0, latticeMap);
			
			HashMap<Integer, Node> latticeMap1 = new HashMap<>();	// Negative sample 1 (DEL1WORD)
//			fillLattice(sentence.split(" "), latticeMap1, 1);
//			globalMap.put(1, latticeMap1);
			
			HashMap<Integer, Node> latticeMap2 = new HashMap<>();	// Negative sample 2 (TRANS)
			fillLattice(sentence.split(" "), latticeMap2, 2);
			globalMap.put(1, latticeMap2);
			
//			HashMap<Integer, Node> latticeMap3 = new HashMap<>();	// Negative sample 3 (DEL1SUBSEQ)
//			fillLattice(sentence.split(" "), latticeMap3, 3);
//			globalMap.put(3, latticeMap3);
		
		
			sentenceGlobals.add(globalMap);
		
//		double[] a = createUniformArray(tagSize);
//		
//		for(int i = 0; i < tagSize; i++){
//			System.out.println(a[i]);
//		}
		
			
			for(int j = 0; j < globalMap.size(); j++) {
				HashMap<Integer, Node> targetMap = globalMap.get(j);
//				System.out.println("Lattice: " + j + "\n" + targetMap);
//				System.out.println("----------------------");
				
				// Iterate over targetMap
				iterateFromStartToEnd(targetMap);
				iterateFromEndToStart(targetMap);
				
			}
		}
		
		/* LBFGS-B part */
		Minimizer alg = new Minimizer();
        alg.getStopConditions().setMaxIterations(10000);
        alg.setDebugLevel(1);
        
		Result ret;
		try {
			double dSum = 0;
			for(double d : tagFeatureWeights) {
				dSum += d;
			}
			System.out.println("First: " + dSum);
			pw.println("First " + dSum);
			
			ret = alg.run(new Function(), tagFeatureWeights);
			
			double finalValue = ret.functionValue;
	        double [] finalGradient = ret.gradient;
	        
	        System.out.println("Final Value: " + finalValue);
	        Main.pw.println("Final Value: " + finalValue);
	        System.out.println("Gradients:");
	        printDoubleArray(finalGradient);
	        double dSum2 = 0;
			for(double d : finalGradient) {
				dSum2 += d;
			}
			System.out.println("Last: " + dSum2);
			pw.println("Last: " + dSum2);
	        
	        updateProbabilities2(finalGradient, gradFeature2Index);
	        tagFeatureWeights = finalGradient;
	        
			for(String s : Main.transitionProbabilities.keySet()) {
				Main.pw.println(s + " " + Main.transitionProbabilities.get(s)); 
			}
			
			for(String s : Main.emissionProbabilities.keySet()) {
				Main.pw.println(s + " " + Main.emissionProbabilities.get(s)); 
			}
	        
		} catch (LBFGSBException e) {
			e.printStackTrace();
		}
		
//		List<Node> list = new ArrayList<Node>(targetMap.values());
//		ListIterator itr = list.listIterator(list.size());
//		
//		while(itr.hasPrevious()){
//			System.out.println(itr.previous());
//		}
		
		// Printing the scores
//		for(Node n : targetMap.values()){
//			System.out.println(n.tagScores);
//		}
		
		// Time operations. Just using for information.
		long endTime = System.nanoTime();
		long duration = (endTime - startTime);  //divide by 1000000 to get milliseconds.
		System.out.println("\nRunning time: " + duration + " nanoseconds ~ " + duration/1000000 + " milliseconds");
		pw.close();
	}
	
	protected static void printDoubleArray(double[] array) {
		for(double d : array) {
			System.out.print(d + " ");
			//Main.pw.print(d + " ");
		}
		System.out.println();
		//Main.pw.println();
	}
	
	private static void fillInitials(String fileName){
		try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
			String line;
			
			while ((line = br.readLine()) != null) {
				String[] values = line.split(" ");
				String tag = values[0].replaceAll("\\|", "-");
				double prob = Double.parseDouble(values[1]);
				
				if(tag.endsWith("_t")) {
					tag = tag.substring(0, tag.length()-2);
//					System.out.println(tag);
					
					if(transitionProbabilities.containsKey(tag)) {
						transitionProbabilities.put(tag, prob);
					}
				}
				else {
					if(!tag.endsWith("<start>") && !tag.endsWith("<end>")) {
						tag = tag.substring(0, tag.lastIndexOf("/"));
//						System.out.println(tag);
					}
					if(emissionProbabilities.containsKey(tag)) {
						emissionProbabilities.put(tag, prob);
					}
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void fillTagList(){
		tagList.add("t1");
		tagList.add("t2");
		tagList.add("t3");
		tagList.add("t4");
		tagList.add("t5");
		tagList.add("t6");
		tagList.add("t7");
		tagList.add("t8");
		tagList.add("t9");
		tagList.add("t10");
		tagList.add("t11");
		tagList.add("t12");
		tagList.add("<s>");
		tagList.add("</s>");
		
		tagSize = tagList.size();
	}
	
	private static void fillTransitionMap(){
		
//		double value = 1.0 / (tagSize);	// For uniform values
		Random r = new Random();
		double value = 0.00000001;	// For zero values
		
		for(int i = 0; i < (tagSize); i++){
			for(int j = 0; j < (tagSize); j++){
//				double value = (r.nextInt(100) / 10000.0);
				
				String key = (tagList.get(i) + "-" + tagList.get(j));

				if(!key.endsWith("<s>") && !key.startsWith("</s>")) {
					if(key.equals("<s>-</s>")) {
						transitionProbabilities.put(key, 0.0);
					}
					else {
						transitionProbabilities.put(key, value);
					}
				}
				else if(key.endsWith("</s>")) {
					if(key.startsWith("</s>") || key.startsWith("<s>")) {
						transitionProbabilities.put(key, 0.0);
					}
					else {
						transitionProbabilities.put(key, value);
					}
				}
				else {
					transitionProbabilities.put(key, 0.0);
				}
				
			}
		}
	}
	
	private static void fillEmissionMap(String sentence){
		words = new ArrayList<String>(Arrays.asList(sentence.split(" ")));
		//double value = 1.0 / (allWords.size());	// For uniform values
		double value = 0.00000001;	// For zero values
		Random r = new Random();
		
		
		for(int i = 0; i < tagSize; i++){
			for(String word : words){
//				double value = (r.nextInt(100) / 10000.0);
				String key = tagList.get(i) + "-" + word;
				if(!emissionProbabilities.containsKey(key)) {
					
					if(key.startsWith("<s>") || key.endsWith("<start>")) {
						if(key.startsWith("<s>") && key.endsWith("<start>")) {
							emissionProbabilities.put(key, 1.0);
						}
						else {
							emissionProbabilities.put(key, 0.0);
						}
					}
					else if(key.contains("</s>")) {
						if(key.equals("</s>-<end>")) {
							emissionProbabilities.put(key, 1.0);
						}
						else {
							emissionProbabilities.put(key, 0.0);
						}
					}
					else {
						emissionProbabilities.put(key, value);
					}
				}
			}
		}
	}
	
//	private static void fillFeatures(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap){
//		int i = 0;
//		for(String s : transitionMap.keySet()) {
//			featureList.add(s);
//			tagFeature2Index.put(s, i);
//			i++;
//		}
//		
//		for(String s : emissionMap.keySet()) {
//			featureList.add(s);
//			tagFeature2Index.put(s, i);
//			i++;
//		}
//	}
//	
//	private static void fillFeatures2(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap, HashMap<String, Integer> gradFeature2Index){
//		int i = 0;
//		for(String s : transitionMap.keySet()) {
//			featureList.add(s);
//			tagFeature2Index.put(s, i);
//			i++;
//		}
//		
//		for(String s : emissionMap.keySet()) {
//			featureList.add(s);
//			tagFeature2Index.put(s, i);
//			i++;
//		}
//	}
//	
	protected static double[] createWeightsArray(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap, HashMap<String, Integer> gradFeature2Index) {
		int size = transitionMap.size() + emissionMap.size();
		
		double[] weights = new double[size];
		
		int i = 0;
		for (String s : transitionMap.keySet()) {
			weights[i] = transitionMap.get(s);
			
//			if(s.equals("<s>-<start>") || s.equals("</s>-<end>")) {
//				s = s + "t"; //TODO
//			}
			gradFeature2Index.put(s, i);
			i++;	
		}
		for (String s : emissionMap.keySet()) {
			weights[i] = emissionMap.get(s);
			gradFeature2Index.put(s, i);
			i++;
		}
		
		return weights;
	}
	
	// Old version
//	protected static double[] createWeightsArray(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap) {
//		int size = transitionMap.size() + emissionMap.size();
//		
//		double[] weights = new double[size];
//		
//		int i = 0;
//		for(double d : transitionMap.values()) {
//			weights[i] = d;
//			i++;
//		}
//		
//		for(double d : emissionMap.values()) {
//			weights[i] = d;
//			i++;
//		}
//		
//		return weights;
//	}
	
//	protected static double[] createWeightsArray2(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap, HashMap<String, Integer> gradFeature2Index) {
//		int size = (tagSize * (tagSize-1)) + ((allWords.size()-1) * 12 + 1); // transition + emission
//		
//		double[] weights = new double[size];
//		
//		int i = 0;
//		
//		for(String s : transitionMap.keySet()) {
//			if(s.endsWith("<s>")) {
//				
//			}
//			else {
//				weights[i] = transitionMap.get(s);
//				gradFeature2Index.put(s, i);
//				i++;
//			}
//			
//		}
//		for(String s : emissionMap.keySet()) {
//			if(s.endsWith("<s>") || s.startsWith("<s>")) { // tag || word
//				if(s.endsWith("<s>") && s.startsWith("<s>")) {
//					weights[i] = 1.0;
//					gradFeature2Index.put(s, i);
//					i++;
//				}
//				else {
//				}
//				
//			}
//			else {
//				weights[i] = emissionMap.get(s);
//				gradFeature2Index.put(s, i);
//				i++;
//			}
//			
//		}
//		
//		return weights;
//	}
//	
	protected static double[] createGradArray(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap, HashMap<String, Integer> gradFeature2Index) {
		//int size = (tagSize * (tagSize-1)) + ((allWords.size()-1) * 12 + 1); // transition + emission
		int size = transitionMap.size() + emissionMap.size();
		
		double[] weights = new double[size];
		int index = 0;
		
		for(String s : transitionMap.keySet()) {
//			System.out.println(s);
			if(s.equals("</s>-</s>") || s.equals("<s>-<s>")) {
				index = gradFeature2Index.get(s);
				weights[index] = 0.0;
			}
			else if(s.startsWith("</s>") || s.endsWith("<s>")) {
				index = gradFeature2Index.get(s);
				weights[index] = 0.0;
			}
			else {
				index = gradFeature2Index.get(s);
				weights[index] = transitionMap.get(s);
			}
		}
		
		for(String s : emissionMap.keySet()) {
			if(s.startsWith("<s>") || s.endsWith("<start>")) {
				if(s.startsWith("<s>") && s.endsWith("<start>")) {
					index = gradFeature2Index.get(s);
					weights[index] = 1.0;
				}
				else {
					index = gradFeature2Index.get(s);
					weights[index] = 0.0;
				}
			}
			else if(s.startsWith("</s>") || s.endsWith("<end>")) {
				if(s.equals("</s>-<end>")) {
					index = gradFeature2Index.get(s);
					weights[index] = 1.0;
				}
				else {
					index = gradFeature2Index.get(s);
					weights[index] = 0.0;
				}
			}
			else {
				index = gradFeature2Index.get(s);
				weights[index] = emissionMap.get(s);
			}
		}
		
		return weights;
	}
	
	
	// This function is using for creating lattices.
	static int sentenceCount = 0;
	public static void fillLattice(String[] words, HashMap<Integer, Node> lattice, int latticeType) {
		int N = words.length - 1;
		//if (N > 2) {
			sentenceCount++;
			if (latticeType == 0) { // Original Lattice
				for (int i = 0; i < N + 1; i++) {
					Node n = new Node(i, words[i]);

					if ((i + 1) <= N) {
						n.next.add(i + 1);
					}
					if ((i - 1) >= 0) {
						n.prev.add(i - 1);
					}

					if (i == N) {
						n.isEndState = true;
					}
					lattice.put(i, n);

				}
			} else if (latticeType == 1) { // Negative Lattice 1 : Del1Word
				for (int i = 0; i < N + 1; i++) {
					if ((i + N - 1) <= N) {
						Node n = new Node(i, words[i]);

						n.next.add(i + 1);
						n.next.add(i + N + 1);

						if (i != 0) {
							n.prev.add(i - 1);
						}

						if (i == N || i == N - 1) {
							n.isEndState = true;
						}

						lattice.put(i, n);
					} 
					else {
						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(i + N - 1, words[i]);

						n1.prev.add(i - 1);

						if (i + 1 <= N) {
							n1.next.add(i + 1);
						}
						if (i + N + 1 < 2 * N) {
							n1.next.add(i + N + 1);
						}

						if (i == N) {
							n1.isEndState = true;
							n2.isEndState = true;
						} else if (i == N - 1) {
							n1.isEndState = true;
						}

						int j = i + N - 1;

						if (j - 1 == N) {
							n2.prev.add(j - N - 1);
						} else {
							n2.prev.add(j - 1);
							n2.prev.add(j - N - 1);
						}

						if (j + 1 < 2 * N) {
							n2.next.add(j + 1);
						}

						lattice.put(i, n1);
						lattice.put(j, n2);
					}
				}
			} else if (latticeType == 2) { // Negative Lattice 2: Trans1
				for (int i = 0; i < N + 1; i++) {
					if (i == 0) {
						Node n = new Node(i, words[i]);
						n.next.add(i + 1);
						n.next.add(i + N + 1);

						lattice.put(i, n);
					} else if (i == 1) {
						int j = 2 * N;
						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(j, words[i]);

						n1.next.add(i + 1);
						n1.next.add(i + N + 1);
						n1.prev.add(i - 1);

						n2.next.add(j + N - 1);
						n2.prev.add(j - N + 1);

						lattice.put(i, n1);
						lattice.put(j, n2);

					} else if (i == 2) {
						int j = i + N - 1;
						int k = i + (2 * N) - 1;
						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(j, words[i]);
						Node n3 = new Node(k, words[i]);

						n1.next.add(i + 1);
						if(N > 3)
							n1.next.add(i + N + 1);
						n1.prev.add(i - 1);

						n2.next.add(j + N - 1);
						n2.prev.add(j - N - 1);

						if(N > 3)
							n3.next.add(k + N - 1);
						n3.prev.add(k - N + 1);

						lattice.put(i, n1);
						lattice.put(j, n2);
						lattice.put(k, n3);
					} else if (i == N) {
						int j = i + N - 1;
						int k = i + (3 * N) - 4;
						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(j, words[i]);
						Node n3 = new Node(k, words[i]);

						n1.prev.add(i - 1);

						n2.next.add(j + N - 1);
						n2.prev.add(j - N - 1);

						n3.prev.add(k - N + 1);
						if(N > 3)
							n3.prev.add(k - 1);

						n1.isEndState = true;
						n3.isEndState = true;

						lattice.put(i, n1);
						lattice.put(j, n2);
						lattice.put(k, n3);

					} else {
						int j = i + N - 1;
						int k = i + (2 * N) - 1;
						int m = i + (3 * N) - 4;
						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(j, words[i]);
						Node n3 = new Node(k, words[i]);
						Node n4 = new Node(m, words[i]);

						if (i == N - 1) {
							n1.next.add(i + 1);
							n1.prev.add(i - 1);

							n2.next.add(j + N - 1);
							n2.prev.add(j - N - 1);

							n3.prev.add(k - N + 1);
							n3.isEndState = true;

							n4.next.add(m + 1);
							n4.prev.add(m - N + 1);
							n4.prev.add(m - 1);
						} else {
							n1.next.add(i + 1);
							n1.next.add(i + N + 1);
							n1.prev.add(i - 1);

							n2.next.add(j + N - 1);
							n2.prev.add(j - N - 1);

							n3.next.add(k + N - 1);
							n3.prev.add(k - N + 1);

							n4.next.add(m + 1);
							n4.prev.add(m - N + 1);

							if(i != 3) {
								n4.prev.add(m - 1);
							}
						}

						lattice.put(i, n1);
						lattice.put(j, n2);
						lattice.put(k, n3);
						lattice.put(m, n4);
					}
				}
			} else if (latticeType == 3) {
				for (int i = 0; i < N + 1; i++) {
					if (i == 0 || i == 1) {
						Node n1 = new Node(i, words[i]);

						if (i == 1) {
							n1.prev.add(i - 1);
							n1.isEndState = true;
						}
						n1.next.add(i+1);
						for (int j = i + N + 1; j < (2 * N); j++) {
							n1.next.add(j);
						}

						lattice.put(i, n1);
					} else if (i == N || i == N - 1) {
						int j = i + N - 1;

						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(j, words[i]);

						if (i == N - 1) {
							n1.next.add(i + 1);
							n2.next.add(j + 1);
							n1.isEndState = true;
						}

						if (i == N) {
							n1.isEndState = true;
							n2.isEndState = true;
						}

						n1.prev.add(i - 1);
						
						if(N > 3)
						n2.prev.add(j - 1);

						for (int k = 0; k < (j - N); k++) {
							n2.prev.add(k);
						}

						lattice.put(i, n1);
						lattice.put(j, n2);
					} else {
						int j = i + N - 1;

						Node n1 = new Node(i, words[i]);
						Node n2 = new Node(j, words[i]);

						n1.next.add(i + 1);
						n2.next.add(j + 1);

						n1.isEndState = true;

						n1.prev.add(i - 1);
						if (j - 1 != N)
							n2.prev.add(j - 1);

						for (int m = i + N + 1; m < (2 * N); m++) {
							n1.next.add(m);
						}

						for (int k = 0; k < (j - N); k++) {
							n2.prev.add(k);
						}

						lattice.put(i, n1);
						lattice.put(j, n2);
					}
				}
			}
		//}
			int size = lattice.size();
			Node finalNode = new Node(size, "<end>");
			for(Node n : returnEndStates(lattice)) {
				n.next.add(finalNode.stateNum);
				finalNode.prev.add(n.stateNum);
				
				lattice.put(n.stateNum, n);
			}
			lattice.put(finalNode.stateNum, finalNode);
			
	}

//	private static double[] createZeroArray(int size){
//		double[] array = new double[size];
//		
//		for(int i = 0; i < size; i++){
//			array[i] = 0.0;
//		}
//		
//		return array;
//	}

	// For creating list values uniformly
	private static ArrayList<Double> createUniformList(int size, boolean isStartState){
		ArrayList<Double> list = new ArrayList<>();

		if(isStartState) {
			for(int i = 0; i < size; i++){
				if(i == size-2) {
					list.add((1.0));
				}
				else {
					list.add((0.0));
				}
			}
		}
		else {
			for(int i = 0; i < size; i++){
				if(i == size-1) {
					list.add((1.0));
				}
				else {
					list.add((0.0));
				}
			}
		}
		
		
		return list;
	}
	
	// Iteration from start to end. Using for alpha values.
	public static void iterateFromStartToEnd(HashMap<Integer, Node> latticeMap){
		for(Node node : latticeMap.values()){
//			System.out.println(node.word);
			if(node.prev.isEmpty()){	// Start node
				node.alpha = createUniformList(tagSize, node.isStartState);
			}
			else if(node.next.isEmpty()){	// End node
				node.beta = createUniformList(tagSize, node.isStartState);
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha");
				node.tagScores = multiply(node.alpha, node.beta);
			}
			else{	// Others
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha");
			}
			
			latticeMap.put(node.stateNum, node);
		}
		//System.out.println("Lattice: " + latticeMap);
	}
	
	// Iteration from end to start. Using for beta values
	public static void iterateFromEndToStart(HashMap<Integer, Node> latticeMap){
		List<Node> list = new ArrayList<Node>(latticeMap.values());
		ListIterator<Node> iterator = list.listIterator(list.size());
		
		while(iterator.hasPrevious()){
			Node node = (Node) iterator.previous();
			if(!node.next.isEmpty()){
				node.beta = calculateValue(node.next, node, latticeMap, "beta");
				node.tagScores = multiply(node.alpha, node.beta);
			}
			
			latticeMap.put(node.stateNum, node);
		}
	}
	
	// Calculates values for alpha or beta.
	public static ArrayList<Double> calculateValue(List<Integer> neighbors, Node n, HashMap<Integer, Node> latticeMap, String decide){
		if(decide.equals("alpha")){
			n.alpha.clear();
			
//			if(n.next.isEmpty()) {
//				for(int i = 0; i < tagSize-1; i++){
//					n.alpha.add(i, 0.0);
//				}
//				double finalResult = 0;
//				for(int counter : neighbors){
//					Node n2 = latticeMap.get(counter);
//					finalResult += calculate(n, n2, tagSize-1, decide);
//					
//				}
//				n.alpha.add(tagSize-1, finalResult);
//			}
//			else {
				for(int i = 0; i < tagSize; i++){
					double finalResult = 0;
					for(int counter : neighbors){
						Node n2 = latticeMap.get(counter);
						finalResult += calculate(n, n2, i, decide);
						
					}
					n.alpha.add(i, finalResult);
					
				}
//				n.alpha.add(tagSize-2, (0.0));
//				n.alpha.add(tagSize-1, (0.0));
//			}
			
			return n.alpha;
		}
		else{
			n.beta.clear();
//			if(n.next.isEmpty()) {
//				for(int i = 0; i < tagSize-1; i++){
//					n.beta.add(i, 0.0);
//				}
//				n.beta.add(tagSize-1, 1.0);
//			}
//			else {
				for(int i = 0; i < tagSize; i++){
					double finalResult = 0;
					for(int counter : neighbors){
						Node n2 = latticeMap.get(counter);
						finalResult += calculate(n, n2, i, decide);
						
					}
					n.beta.add(i, finalResult);
				}
//				n.beta.add(tagSize-2, (0.0));
//				n.beta.add(tagSize-1, (0.0));
//			}
			
			return n.beta;
		}
	}
	
/*	public static ArrayList<Double> calculateValue2(List<Integer> neighbors, Node n, HashMap<Integer, Node> latticeMap, String decide){
		if(decide.equals("alpha")){
			for(int i = 0; i < tagSize-1; i++){
				double finalResult = 0;
				for(int counter : neighbors){
					Node n2 = latticeMap.get(counter);
					finalResult += Math.exp(calculate(n, n2, i, decide));
					
				}
				n.alpha.add(i, finalResult);
			}
			n.alpha.add(tagSize-1, 0.0);
			
			return n.alpha;
		}
		else{
			for(int i = 0; i < tagSize-1; i++){
				double finalResult = 0;
				for(int counter : neighbors){
					Node n2 = latticeMap.get(counter);
					finalResult += calculate(n, n2, i, decide);
					
				}
				n.beta.add(i, finalResult);
			}
			n.beta.add(tagSize-1, 0.0);
			
			return n.beta;
		}
	}*/
	
	// Using for calculating values of a node
	private static double calculate(Node n, Node n2, int tagNumber, String decide){
		double sum = 0;
		for (int j = 0; j < tagSize; j++) {
			
			if(decide.equals("alpha"))
				sum += (n2.alpha.get(j)) * (transitionProbabilities.get(tagList.get(j) + "-" + tagList.get(tagNumber)));
			else {
				//System.out.println(tagList.get(tagNumber) + "-" + tagList.get(j));
				sum += (n2.beta.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j))) * (emissionProbabilities.get(tagList.get(j)  + "-" + n2.word));
			}
		}
		
		if(decide.equals("alpha")) {
			return sum* emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word);
		}
		else {
			return sum;
		}
		
//		return Math.exp(sum* emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word));
	}
	
/*	private static double calculate2(Node n, Node n2, int tagNumber, String decide){
		double sum = 0;
		ArrayList<Double> list = new ArrayList<>();
		for (int j = 0; j < tagSize; j++) {
			
			if(decide.equals("alpha")) {
				sum += Math.log(n2.alpha.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j)) * emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word));
				list.add(Math.log(n2.alpha.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j)) * emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word)));
			}
			else {
				//System.out.println(tagList.get(tagNumber) + "-" + tagList.get(j));
				sum += Math.log(n2.beta.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j)) * emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word));
				list.add(Math.log(n2.beta.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j)) * emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word)));
			}
		}
		
		double max = Collections.max(list);

		return Math.exp(Math.log(Math.exp(sum-max)) + max);
	}*/
	
	
	// Using for calculating scores. No need!
	private static ArrayList<Double> multiply(ArrayList<Double> alpha, ArrayList<Double> beta){
		ArrayList<Double> scores = new ArrayList<>();
		
		for(int i = 0; i < tagSize; i++){
			double score = ((alpha.get(i)) * (beta.get(i)));
			
			scores.add(i, score);
		}
		
		return scores;
	}
	
	public static Node returnEndState(HashMap<Integer, Node> lattice){
		Node returnNode = null;
		for(Node n : lattice.values()) {
			if(n.next.isEmpty()) {
				returnNode = n;
			}
		}
		return returnNode;
	}
	
	public static ArrayList<Node> returnEndStates(HashMap<Integer, Node> lattice){
		
		ArrayList<Node> endStates = new ArrayList<>();
		
		for(Node n : lattice.values()) {
			if(n.isEndState) {
				endStates.add(n);
			}
		}
		
		return endStates;
	}
	
//	public static void updateProbabilities(double[] weights) {
//		for(int i = 0; i < weights.length; i++) {
//			if(i < transitionProbabilities.size()) {
//				if(featureList.get(i).endsWith("<s>")) {
//					transitionProbabilities.put(featureList.get(i), 0.0);
//				}
//				else {
//					transitionProbabilities.put(featureList.get(i), weights[i]);
//				}
//			}
//			else {
//				if(featureList.get(i).startsWith("<s>") || featureList.get(i).endsWith("<s>")) {
//					if(featureList.get(i).startsWith("<s>") && featureList.get(i).endsWith("<s>")) {	// <s>-<s> emission
//						emissionProbabilities.put(featureList.get(i), 10.0);
//					}
//					else {
//						emissionProbabilities.put(featureList.get(i), 0.0);
//					}
//				}
//				else {
//					emissionProbabilities.put(featureList.get(i), weights[i]);
//				}
//			}
//		}
//	}
	public static void updateProbabilities2(double[] weights, HashMap<String, Integer> gradFeature2Index) {
		for(String s : gradFeature2Index.keySet()) {
			int index = gradFeature2Index.get(s);
			
//			if(s.equals("<s>-<s>t") || s.equals("</s>-</s>t")) {
//				s = s.substring(0, s.length()-1);
//			}
			
			if(transitionProbabilities.containsKey(s)) {
				transitionProbabilities.put(s, weights[index]);
			}
			else if(emissionProbabilities.containsKey(s)) {
				emissionProbabilities.put(s, weights[index]);
//				emissionProbabilities.put("<s>-<s>", 10.0);
			}
		}
	}
	
}