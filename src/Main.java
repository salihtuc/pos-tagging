
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
	
	public static HashMap<String, Double> transitionProbabilities = new HashMap<>();
	public static HashMap<String, Double> emissionProbabilities = new HashMap<>();
	
	public static ArrayList<String> tagList = new ArrayList<>();	// Holds tags
	public static ArrayList<String> featureList = new ArrayList<>();	// Holds features
	
	public static HashMap<String, Integer> tagFeature2Index = new HashMap<>();
	
	public static double[] tagFeatureWeights = null;	// Feature sized weight array that we use in LBFGS-B
	public static double[] tagFeatureGradients = null;

	public static void main(String[] args) {
		
		// Time operations. Just using for information.
		long startTime = System.nanoTime();
		try {
			pw = new PrintWriter(new File("outp.txt"));
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		//try (BufferedReader br = new BufferedReader(new FileReader("PennCorpus12.txt"))) {
		try (BufferedReader br = new BufferedReader(new FileReader("input.txt"))) {
			String start = "<s> ";
			String line;

			while ((line = br.readLine()) != null) {
				sentences.add(start + line.toLowerCase());
				
				allWords.addAll(Arrays.asList((start + line.toLowerCase()).split(" ")));
				
			}

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
		
		tagFeatureWeights = createWeightsArray(transitionProbabilities, emissionProbabilities);
		tagFeatureGradients = createZeroArray(tagFeatureWeights.length);
		fillFeatures(transitionProbabilities, emissionProbabilities);
		
		
//		System.out.println("Transitions:\n" + transitionProbabilities);
//		System.out.println("------------\nEmissions:\n" + emissionProbabilities + "\n-----------");
		
		
		for(String sentence : sentences) {
			
			globalMap = new HashMap<>();
			// Create all lattices and put them into the globalMap
			HashMap<Integer, Node> latticeMap = new HashMap<>();	// Original sentence
			fillLattice(sentence.split(" "), latticeMap, 0);
			globalMap.put(0, latticeMap);
			
			HashMap<Integer, Node> latticeMap1 = new HashMap<>();	// Negative sample 1 (DEL1WORD)
			fillLattice(sentence.split(" "), latticeMap1, 1);
			globalMap.put(1, latticeMap1);
			
			HashMap<Integer, Node> latticeMap2 = new HashMap<>();	// Negative sample 2 (TRANS)
			fillLattice(sentence.split(" "), latticeMap2, 2);
			globalMap.put(2, latticeMap2);
			
			HashMap<Integer, Node> latticeMap3 = new HashMap<>();	// Negative sample 3 (DEL1SUBSEQ)
			fillLattice(sentence.split(" "), latticeMap3, 3);
			globalMap.put(3, latticeMap3);
		
		
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
	        
	        updateProbabilities(finalGradient);
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
		
		tagSize = tagList.size();
	}
	
	private static void fillTransitionMap(){
		
		double value = 1.0 / (tagSize);	// For uniform values
		
		//double value = 0.0;	// For zero values
		
		for(int i = 0; i < (tagSize); i++){
			for(int j = 0; j < (tagSize); j++){
				transitionProbabilities.put((tagList.get(i) + "-" + tagList.get(j)), value);
			}
		}
	}
	
	private static void fillEmissionMap(String sentence){
		String[] word = sentence.split(" ");
		words= new ArrayList<String>(Arrays.asList(word));
		double value = 1.0 / (allWords.size() * tagSize);	// For uniform values
		//double value = 0.0;	// For zero values
		
		for(int i = 0; i < tagSize; i++){
			for(String wordd : words){
				if(!emissionProbabilities.containsKey(tagList.get(i) + "-" + wordd))
					emissionProbabilities.put(tagList.get(i) + "-" + wordd, value);
			}
		}
	}
	
	private static void fillFeatures(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap){
		int i = 0;
		for(String s : transitionMap.keySet()) {
			featureList.add(s);
			tagFeature2Index.put(s, i);
			i++;
		}
		
		for(String s : emissionMap.keySet()) {
			featureList.add(s);
			tagFeature2Index.put(s, i);
			i++;
		}
	}
	
	protected static double[] createWeightsArray(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap) {
		int size = transitionMap.size() + emissionMap.size();
		
		double[] weights = new double[size];
		
		int i = 0;
		for(double d : transitionMap.values()) {
			weights[i] = d;
			i++;
		}
		
		for(double d : emissionMap.values()) {
			weights[i] = d;
			i++;
		}
		
		return weights;
	}
	
	// This function is using for creating lattices.
	public static void fillLattice(String[] words, HashMap<Integer, Node> lattice, int latticeType){
		int N = words.length - 1;
		if(latticeType == 0){	// Original Lattice
			for(int i = 0; i < N+1; i++){
				Node n = new Node(i, words[i]);
				
				if((i+1) <= N){
					n.next.add(i+1);
				}
				if((i-1) >= 0){
					n.prev.add(i-1);
				}
				
				if(i == N) {
					n.isEndState = true;
				}
				lattice.put(i, n);
				
			}
		}
		else if(latticeType == 1){	// Negative Lattice 1 : Del1Word
			for(int i = 0; i < N+1; i++){
				if((i+N - 1) <= N){
					Node n = new Node(i, words[i]);
					
					n.next.add(i+1);
					n.next.add(i+N+1);
					
					if(i != 0){
						n.prev.add(i-1);
					}
					
					if(i == N || i == N-1) {
						n.isEndState = true;
					}
					
					lattice.put(i, n);
				}
				else{
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(i+N-1, words[i]);
					
					n1.prev.add(i-1);
					
					if(i+1 <= N){
						n1.next.add(i+1);
					}
					if(i+N+1 < 2*N){
						n1.next.add(i+N+1);
					}
					
					if(i == N) {
						n1.isEndState = true;
						n2.isEndState = true;
					}
					else if(i == N-1) {
						n2.isEndState = true;
					}
					
					int j = i+N-1;
					
					if(j-1 == N){
						n2.prev.add(j-N-1);
					}
					else{
						n2.prev.add(j-1);
						n2.prev.add(j-N-1);
					}
					
					if(j+1 < 2*N){
						n2.next.add(j+1);
					}
					
					lattice.put(i, n1);
					lattice.put(j, n2);
				}
			}
		}
		else if(latticeType == 2){	// Negative Lattice 2: Trans1
			for(int i = 0; i < N+1; i++){
				if (i == 0) {
					Node n = new Node(i, words[i]);
					n.next.add(i+1);
					n.next.add(i+N+1);
					
					lattice.put(i, n);
				} 
				else if (i == 1) {
					int j = 2*N;
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(j, words[i]);
					
					n1.next.add(i+1);
					n1.next.add(i+N+1);
					n1.prev.add(i-1);
					
					n2.next.add(j+N-1);
					n2.prev.add(j-N+1);
					
					lattice.put(i, n1);
					lattice.put(j, n2);
					
				} 
				else if (i == 2) {
					int j = i+N-1;
					int k = i+(2*N)-1;
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(j, words[i]);
					Node n3 = new Node(k, words[i]);
					
					n1.next.add(i+1);
					n1.next.add(i+N+1);
					n1.prev.add(i-1);
					
					n2.next.add(j+N-1);
					n2.prev.add(j-N-1);
					
					n3.next.add(k+N-1);
					n3.prev.add(k-N+1);
					
					lattice.put(i, n1);
					lattice.put(j, n2);
					lattice.put(k, n3);
				} 
				else if (i == N) {
					int j = i+N-1;
					int k = i+(3*N)-4;
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(j, words[i]);
					Node n3 = new Node(k, words[i]);
					
					n1.prev.add(i-1);
					
					n2.next.add(j+N-1);
					n2.prev.add(j-N-1);
					
					n3.prev.add(k-N+1);
					n3.prev.add(k-1);
					
					n1.isEndState = true;
					n3.isEndState = true;
					
					lattice.put(i, n1);
					lattice.put(j, n2);
					lattice.put(k, n3);
					
				} 
				else {
					int j = i+N-1;
					int k = i+(2*N)-1;
					int m = i+(3*N)-4;
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(j, words[i]);
					Node n3 = new Node(k, words[i]);
					Node n4 = new Node(m, words[i]);
					
					if(i == N-1){
						n1.next.add(i+1);
						n1.prev.add(i-1);
						
						n2.next.add(j+N-1);
						n2.prev.add(j-N-1);
						
						n3.prev.add(k-N+1);
						n3.isEndState = true;
						
						n4.next.add(m+1);
						n4.prev.add(m-N+1);
						n4.prev.add(m-1);
					}
					else{
						n1.next.add(i+1);
						n1.next.add(i+N+1);
						n1.prev.add(i-1);
						
						n2.next.add(j+N-1);
						n2.prev.add(j-N-1);
						
						n3.next.add(k+N-1);
						n3.prev.add(k-N+1);
						
						n4.next.add(m+1);
						n4.prev.add(m-N+1);
						
						// Look at it!!!
						n4.prev.add(m-1);
					}
					
					lattice.put(i, n1);
					lattice.put(j, n2);
					lattice.put(k, n3);
					lattice.put(m, n4);
				}
			}
		}
		else if(latticeType == 3){
			for(int i = 0; i < N+1; i++){
				if(i == 0 || i == 1){
					Node n1 = new Node(i, words[i]);
					
					if(i == 1){
						n1.prev.add(i-1);
						n1.isEndState = true;
					}
					
					for(int j = i+N+1; j < (2*N); j++){
						n1.next.add(j);
					}
					
					lattice.put(i, n1);
				}
				else if(i == N || i == N-1){
					int j = i+N-1;
					
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(j, words[i]);
					
					if(i == N-1){
						n1.next.add(i+1);
						n2.next.add(j+1);
						n1.isEndState = true;
					}
					
					if(i == N) {
						n1.isEndState = true;
						n2.isEndState = true;
					}
					
					n1.prev.add(i-1);
					n2.prev.add(j-1);
					
					for(int k = 0; k < (j-N); k++){
						n2.prev.add(k);
					}
					
					lattice.put(i, n1);
					lattice.put(j, n2);
				}
				else{
					int j = i+N-1;
					
					Node n1 = new Node(i, words[i]);
					Node n2 = new Node(j, words[i]);
					
					n1.next.add(i+1);
					n2.next.add(j+1);
					
					n1.isEndState = true;
					
					n1.prev.add(i-1);
					if(j-1 != N)
						n2.prev.add(j-1);
					
					for(int m = i+N+1; m < (2*N); m++){
						n1.next.add(m);
					}
					
					for(int k = 0; k < (j-N); k++){
						n2.prev.add(k);
					}
					
					lattice.put(i, n1);
					lattice.put(j, n2);
				}
			}
		}
	}

	private static double[] createZeroArray(int size){
		double[] array = new double[size];
		
		for(int i = 0; i < size; i++){
			array[i] = 0.0;
		}
		
		return array;
	}

	// For creating list values uniformly
	private static ArrayList<Double> createUniformList(int size, boolean isStartState, String decide){
		ArrayList<Double> list = new ArrayList<>();
		double initialValue;
		
		if(decide.equals("alpha")) {
			initialValue = Math.PI;
		}
		else {
			initialValue = 1.0;
		}
		
		if(isStartState) {
			for(int i = 0; i < size; i++){
				if(i == size-1) {
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
					list.add((0.0));
				}
				else {
					list.add((initialValue/(size-1)));
				}
			}
		}
		
		
		return list;
	}
	
	// Iteration from start to end. Using for alpha values.
	public static void iterateFromStartToEnd(HashMap<Integer, Node> latticeMap){
		for(Node node : latticeMap.values()){
			if(node.prev.isEmpty()){	// Start node
				node.alpha = createUniformList(tagSize, node.isStartState, "alpha");
			}
			else if(node.next.isEmpty()){	// End node
				node.beta = createUniformList(tagSize, node.isStartState, "beta");
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha");
				node.tagScores = multiply(node.alpha, node.beta);
			}
			else{	// Others
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha");
			}
			
			latticeMap.put(node.stateNum, node);
		}
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
			for(int i = 0; i < tagSize-1; i++){
				double finalResult = 0;
				for(int counter : neighbors){
					Node n2 = latticeMap.get(counter);
					finalResult += calculate(n, n2, i, decide);
					
				}
				n.alpha.add(i, finalResult);
			}
			n.alpha.add(tagSize-1, (0.0));
			
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
			n.beta.add(tagSize-1, (0.0));
			
			return n.beta;
		}
	}
	
	public static ArrayList<Double> calculateValue2(List<Integer> neighbors, Node n, HashMap<Integer, Node> latticeMap, String decide){
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
	}
	
	// Using for calculating values of a node
	private static double calculate(Node n, Node n2, int tagNumber, String decide){
		double sum = 0;
		for (int j = 0; j < tagSize; j++) {
			
			if(decide.equals("alpha"))
				sum += (n2.alpha.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j)));
			else {
				//System.out.println(tagList.get(tagNumber) + "-" + tagList.get(j));
				sum += (n2.beta.get(j)) * (transitionProbabilities.get(tagList.get(tagNumber) + "-" + tagList.get(j)));
			}
		}

		return Math.exp(sum* emissionProbabilities.get(tagList.get(tagNumber) + "-" + n.word));
	}
	
	private static double calculate2(Node n, Node n2, int tagNumber, String decide){
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
	}
	
	
	// Using for calculating scores. No need!
	private static ArrayList<Double> multiply(ArrayList<Double> alpha, ArrayList<Double> beta){
		ArrayList<Double> scores = new ArrayList<>();
		
		for(int i = 0; i < tagSize; i++){
			double score = (Math.log(alpha.get(i)) * (beta.get(i)));
			
			scores.add(i, score);
		}
		
		return scores;
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
	
	public static void updateProbabilities(double[] weights) {
		for(int i = 0; i < weights.length; i++) {
			if(i < transitionProbabilities.size()) {
				transitionProbabilities.put(featureList.get(i), weights[i]);
			}
			else {
				emissionProbabilities.put(featureList.get(i), weights[i]);
			}
		}
	}
	
	public static ArrayList<String> returnUniqueWords(HashMap<Integer, Node> lattice) {
		ArrayList<String> uniqueWords = new ArrayList<>();
		
		for(Node n : lattice.values()) {
			if(!uniqueWords.contains(n.word)) {
				uniqueWords.add(n.word);
			}
		}
		
		return uniqueWords;
	}
}
