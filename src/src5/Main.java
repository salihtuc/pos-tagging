
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
import lbfgsb.Minimizer;
//import lbfgsb.Minimizer;
import lbfgsb.Result;

public class Main {

	public static PrintWriter pw = null;

	public static int tagSize = 12;
	public static ArrayList<String> sentences = new ArrayList<>();
	public static ArrayList<String> words = new ArrayList<>();

	public static ArrayList<String> allWords = new ArrayList<>();
	public static HashSet<String> uniqueValues;

	public static HashMap<Integer, HashMap<Integer, Node>> globalMap; // Holds all lattices

	public static ArrayList<HashMap<Integer, HashMap<Integer, Node>>> sentenceGlobals = new ArrayList<>();

	public static HashMap<String, Double> initialProbabilitiesFromFile = new HashMap<>();

	public static HashMap<String, Double> transitionProbabilities = new HashMap<>();
	public static HashMap<String, Double> emissionProbabilities = new HashMap<>();
	public static HashMap<String, Double> initialProbabilities = new HashMap<>();

	public static ArrayList<String> tagList = new ArrayList<>(); // Holds tags
	public static ArrayList<String> featureList = new ArrayList<>(); // Holds features

	public static HashMap<String, Integer> tagFeature2Index = new HashMap<>();
	public static HashMap<String, Integer> gradFeature2Index = new HashMap<>();

	public static double[] tagFeatureWeights = null; // Feature sized weight array that we use in LBFGS-B
	public static double[] tagFeatureGradients = null;

	public static void main(String[] args) {

		// Time operations. Just using for information.
		long startTime = System.nanoTime();
		try {
			pw = new PrintWriter(new File("outp.txt"));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}

		// try (BufferedReader br = new BufferedReader(new
		// FileReader("PennCorpus12.txt"))) {
		try (BufferedReader br = new BufferedReader(new FileReader("input.txt"))) {
			String start = "<start> ";
			String line;

			while ((line = br.readLine()) != null) {

				if (line.split(" ").length > 2) {
					sentences.add(start + line.toLowerCase());
//					sentences.add(line.toLowerCase());

					allWords.addAll(Arrays.asList((line.toLowerCase()).trim().split(" ")));
					// allWords.addAll(Arrays.asList((line.toLowerCase()).split("
					// ")));
				}

			}
//			// System.out.println(allWords);
//			for (int i = 0; i < allWords.size(); i++) {
//				if (Collections.frequency(allWords, allWords.get(i)) == 1) {
//					allWords.set(i, "unknown");
//				}
//			}
//			int sum = 0;
//			for (int j = 0; j < sentences.size(); j++) {
//				String sentence = sentences.get(j);
//				int count = sentence.split(" ").length;
//				String modifiedSentence = "";
//				for (int k = sum; k < sum + count; k++) {
//					modifiedSentence += allWords.get(k) + " ";
//
//				}
//				sum += count;
//				sentences.set(j, "<start> " + modifiedSentence);
//				// System.out.println("<start> "+modifiedSentence);
//			}
//			 allWords.add("<end>");
			uniqueValues = new HashSet<>(allWords);
			allWords.clear();
			allWords.addAll(uniqueValues);

		} catch (IOException e) {
			e.printStackTrace();
		}

		// String sentence = "<s> Natural language is a delicate thing";

		// Initialization of features and tags.

		fillTagList(); // Creating tags (filling tagList)
		fillTransitionMap(); // Creating transitionProbabilities map
		fillInitialMap();

		for (String sentence : sentences)
			fillEmissionMap(sentence); // Creating emissionProbabilities map

		// fillEmissionMap("<end>");
		fillInitialFromFile("initialProb.txt");

		tagFeatureWeights = createWeightsArray(transitionProbabilities, initialProbabilities, emissionProbabilities,
				gradFeature2Index);

		for (String sentence : sentences) {

			globalMap = new HashMap<>();
			// Create all lattices and put them into the globalMap
			HashMap<Integer, Node> latticeMap = new HashMap<>(); // Original sentence
			fillLattice(sentence.split(" "), latticeMap, 0);
			globalMap.put(0, latticeMap);

			HashMap<Integer, Node> latticeMap1 = new HashMap<>(); // Negative sample 1 (DEL1WORD)
			// fillLattice(sentence.split(" "), latticeMap1, 1);
			// globalMap.put(1, latticeMap1);

			HashMap<Integer, Node> latticeMap2 = new HashMap<>(); // Negative sample 2 (TRANS)
			fillLattice(sentence.split(" "), latticeMap2, 2);
			globalMap.put(1, latticeMap2);

			// HashMap<Integer, Node> latticeMap3 = new HashMap<>(); // Negative sample 3
			// (DEL1SUBSEQ)
			// fillLattice(sentence.split(" "), latticeMap3, 3);
			// globalMap.put(3, latticeMap3);

			sentenceGlobals.add(globalMap);

			for (int j = 0; j < globalMap.size(); j++) {
				HashMap<Integer, Node> targetMap = globalMap.get(j);

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
			for (double d : tagFeatureWeights) {
				dSum += d;
			}
			System.out.println("First: " + dSum);
			pw.println("First " + dSum);

			ret = alg.run(new Function(), tagFeatureWeights);

			double finalValue = ret.functionValue;
			double[] finalGradient = ret.gradient;

			System.out.println("Final Value: " + finalValue);
			Main.pw.println("Final Value: " + finalValue);
			System.out.println("Gradients:");
			printDoubleArray(finalGradient);
			double dSum2 = 0;
			for (double d : finalGradient) {
				dSum2 += d;
			}
			System.out.println("Last: " + dSum2);
			pw.println("Last: " + dSum2);

			updateProbabilities(finalGradient, gradFeature2Index);
			tagFeatureWeights = finalGradient;

			for (String s : Main.transitionProbabilities.keySet()) {
				Main.pw.println(s + " " + Main.transitionProbabilities.get(s));
			}
			for (String s : Main.initialProbabilities.keySet()) {
				Main.pw.println(s + " " + Main.initialProbabilities.get(s));
			}
			for (String s : Main.emissionProbabilities.keySet()) {
				Main.pw.println(s + " " + Main.emissionProbabilities.get(s));
			}

		} catch (LBFGSBException e) {
			e.printStackTrace();
		}

		// Time operations. Just using for information.
		long endTime = System.nanoTime();
		long duration = (endTime - startTime); // divide by 1000000 to get milliseconds.
		System.out.println("\nRunning time: " + duration + " nanoseconds ~ " + duration / 1000000 + " milliseconds");
		pw.close();
	}

	protected static void printDoubleArray(double[] array) {
		for (double d : array) {
			System.out.print(d + " ");
			// Main.pw.print(d + " ");
		}
		System.out.println();
		// Main.pw.println();
	}

	private static void fillInitialFromFile(String fileName) {
		try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
			String line;
			int i = 1;
			while ((line = br.readLine()) != null) {
				String[] values = line.split(" ");
				String key = values[0];
				double prob = Double.parseDouble(values[1]);

				if (key.endsWith("_t")) { // Transition or Initial

					key = key.substring(0, key.length() - 2);

					if (key.contains("<s>") || key.contains("</s>")) {
						initialProbabilities.put(key, prob);
					} else {
						transitionProbabilities.put(key, prob);
					}
				} else {
					if (emissionProbabilities.containsKey(key)) {
						emissionProbabilities.put(key, prob);
					}
				}

				i++;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void fillInitialMap() {
		double value = 0.1;

		for (int i = 0; i < tagList.size(); i++) {
			String key1 = "<s>|" + tagList.get(i);
			String key2 = tagList.get(i) + "|</s>";

			initialProbabilities.put(key1, value);
			initialProbabilities.put(key2, value);
		}
	}

	private static void fillTagList() {
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
		// tagList.add("<s>");
		// tagList.add("</s>");

		tagSize = tagList.size();
	}

	private static void fillTransitionMap() {

		// double value = 1.0 / (tagSize); // For uniform values
		Random r = new Random();
		double value = 0.1; // For zero values

		for (int i = 0; i < (tagSize); i++) {
			for (int j = 0; j < (tagSize); j++) {
				// double value = (r.nextInt(100) / 10000.0);

				String key = (tagList.get(i) + "|" + tagList.get(j));

				transitionProbabilities.put(key, value);

			}
		}
	}

	private static void fillEmissionMap(String sentence) {
		words = new ArrayList<String>(Arrays.asList(sentence.split(" ")));
		// double value = 1.0 / (allWords.size()); // For uniform values
		double value = 0.1; // For zero values
		Random r = new Random();

		for (int i = 0; i < tagSize; i++) {
			for (String word : words) {
				// double value = (r.nextInt(100) / 10000.0);
				String key = tagList.get(i) + "|" + word;
				if (!emissionProbabilities.containsKey(key)) {
					emissionProbabilities.put(key, value);
				}
			}
		}
	}

	// protected static double[] createWeightsArray(HashMap<String, Double>
	// transitionMap, HashMap<String, Double> initialMap, HashMap<String, Integer>
	// gradFeature2Index) {
	// int size = transitionMap.size() + initialMap.size();
	//
	// double[] weights = new double[size];
	//
	// int i = 0;
	// for (String s : transitionMap.keySet()) {
	// weights[i] = transitionMap.get(s);
	//
	// gradFeature2Index.put(s, i);
	// i++;
	// }
	// for (String s : initialMap.keySet()) {
	// weights[i] = initialMap.get(s);
	//
	// gradFeature2Index.put(s, i);
	// i++;
	// }
	//
	// return weights;
	// }

	// SAME FUNCTION WITH EMISSION PROBS
	protected static double[] createWeightsArray(HashMap<String, Double> transitionMap,
			HashMap<String, Double> initialMap, HashMap<String, Double> emissionMap,
			HashMap<String, Integer> gradFeature2Index) {
		int size = transitionMap.size() + initialMap.size() + emissionMap.size();

		double[] weights = new double[size];

		int i = 0;
		for (String s : transitionMap.keySet()) {
			weights[i] = transitionMap.get(s);

			gradFeature2Index.put(s, i);
			i++;
		}
		for (String s : initialMap.keySet()) {
			weights[i] = initialMap.get(s);

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

	protected static double[] createGradArray(HashMap<String, Double> transitionMap, HashMap<String, Double> initialMap,
			HashMap<String, Double> emissionMap, HashMap<String, Integer> gradFeature2Index) {
		int size = gradFeature2Index.size();

		double[] weights = new double[size];
		int index = 0;

		for (String s : transitionMap.keySet()) {
			index = gradFeature2Index.get(s);
			weights[index] = transitionMap.get(s);
		}

		for (String s : initialMap.keySet()) {
			index = gradFeature2Index.get(s);
			weights[index] = initialMap.get(s);
		}

		for (String s : emissionMap.keySet()) {
			index = gradFeature2Index.get(s);
			weights[index] = emissionMap.get(s);
		}

		return weights;
	}

	// This function is using for creating lattices.
	static int sentenceCount = 0;

	public static void fillLattice(String[] words, HashMap<Integer, Node> lattice, int latticeType) {
		int N = words.length - 1;
		if (N > 2) {
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
					} else {
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
						if (N > 3)
							n1.next.add(i + N + 1);
						n1.prev.add(i - 1);

						n2.next.add(j + N - 1);
						n2.prev.add(j - N - 1);

						if (N > 3)
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
						if (N > 3)
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

							if (i != 3) {
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
						n1.next.add(i + 1);
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

						if (N > 3)
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
		}
		// int size = lattice.size();
		// Node finalNode = new Node(size, "<end>");
		// for(Node n : returnEndStates(lattice)) {
		// n.next.add(finalNode.stateNum);
		// finalNode.prev.add(n.stateNum);
		//
		// lattice.put(n.stateNum, n);
		// }
		// lattice.put(finalNode.stateNum, finalNode);

	}

	// For creating list values uniformly
	private static ArrayList<Double> createInitialList(String word, String decide) {
		int size = tagSize;
		ArrayList<Double> list = new ArrayList<>();

		if (decide.equals("alpha")) {
			for (int i = 0; i < size; i++) {
				double value = Math.exp(initialProbabilities.get("<s>|" + tagList.get(i))
						+ (emissionProbabilities.get(tagList.get(i) + "|" + word)));

				// we put only sum of initial value + emission
				list.add(value);
			}
		} else {
			for (int i = 0; i < size; i++) {
				double value = Math.exp(initialProbabilities.get(tagList.get(i) + "|</s>"));// we do not add emission </s> |
																					// <end> +
																					// (emissionProbabilities.get(tagList.get(i)+
																					// "|" + word));

				// we put only sum of initial value
				list.add(value);
			}
		}

		return list;
	}

	// Iteration from start to end. Using for alpha values.
	public static void iterateFromStartToEnd(HashMap<Integer, Node> latticeMap) {
		ArrayList<Integer> startStates = returnStartStates(latticeMap);
		for (Node node : latticeMap.values()) {

			// System.out.println(node.word);
			if (startStates.contains(node.stateNum)) { // First nodes for calculation
				node.alpha = createInitialList(node.word, "alpha");
			} else if (node.next.isEmpty()) { // End node
				node.beta = createInitialList(node.word, "beta");
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha");
				node.tagScores = multiply(node.alpha, node.beta);
			} else if (node.prev.isEmpty()) { // Start node
				// TODO
			} else { // Others
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha");
			}

			latticeMap.put(node.stateNum, node);
		}
	}

	// Iteration from end to start. Using for beta values
	public static void iterateFromEndToStart(HashMap<Integer, Node> latticeMap) {
		List<Node> list = new ArrayList<Node>(latticeMap.values());
		ListIterator<Node> iterator = list.listIterator(list.size());

		while (iterator.hasPrevious()) {
			Node node = (Node) iterator.previous();
			if (!node.next.isEmpty() && !node.prev.isEmpty()) {
				node.beta = calculateValue(node.next, node, latticeMap, "beta");
				node.tagScores = multiply(node.alpha, node.beta);
			}

			latticeMap.put(node.stateNum, node);
		}
	}

	// Calculates values for alpha or beta.
	public static ArrayList<Double> calculateValue(List<Integer> neighbors, Node n, HashMap<Integer, Node> latticeMap,
			String decide) {

		if (decide.equals("alpha")) {
			n.alpha.clear();

			for (int i = 0; i < tagSize; i++) {
				double finalResult = 0;
				for (int counter : neighbors) {
					Node n2 = latticeMap.get(counter);
					finalResult += Math.exp(calculate(n, n2, i, decide));

				}
				n.alpha.add((finalResult));

			}
			return n.alpha;
		}

		else {
			n.beta.clear();
			for (int i = 0; i < tagSize; i++) {
				double finalResult = 0;
				for (int counter : neighbors) {
					Node n2 = latticeMap.get(counter);
					finalResult += Math.exp(calculate(n, n2, i, decide));

				}
				n.beta.add((finalResult));
			}
			return n.beta;
		}
	}

	// Using for calculating values of a node
	private static double calculate(Node n, Node n2, int tagNumber, String decide) {
		double sum = 0;
		ArrayList<Double> list = new ArrayList<>();
		for (int j = 0; j < tagSize; j++) {

			if (decide.equals("alpha")) {

				sum = transitionProbabilities.get(tagList.get(j) + "|" + tagList.get(tagNumber))
						+ (emissionProbabilities.get(tagList.get(tagNumber) + "|" + n.word)) + Math.log(n2.alpha.get(j));
				list.add(sum);

			} 
			else {

				sum =  transitionProbabilities.get(tagList.get(tagNumber) + "|" + tagList.get(j))
						+ (emissionProbabilities.get(tagList.get(j) + "|" + n2.word)) + Math.log(n2.beta.get(j));
				list.add(sum);

			}

		}
		sum = logSumOfExponentials(list);
		list.clear();
		return sum;

	}

	// Using for calculating scores.
	private static ArrayList<Double> multiply(ArrayList<Double> alpha, ArrayList<Double> beta) {
		ArrayList<Double> scores = new ArrayList<>();

		for (int i = 0; i < tagSize; i++) {
			double score = (Math.log(alpha.get(i)) + Math.log(beta.get(i)));

			scores.add(i, score);
		}

		return scores;
	}

	public static Node returnEndState(HashMap<Integer, Node> lattice) {
		Node returnNode = null;
		for (Node n : lattice.values()) {
			if (n.next.isEmpty()) {
				returnNode = n;
			}
		}
		return returnNode;
	}

	public static ArrayList<Integer> returnStartStates(HashMap<Integer, Node> lattice) {

		ArrayList<Integer> startStates = new ArrayList<>();

		for (Node n : lattice.values()) {
			if (n.prev.contains(0)) {
				startStates.add(n.stateNum);
			}
		}

		return startStates;
	}

	public static ArrayList<Node> returnEndStates(HashMap<Integer, Node> lattice) {

		ArrayList<Node> endStates = new ArrayList<>();

		for (Node n : lattice.values()) {
			if (n.isEndState) {
				endStates.add(n);
			}
		}

		return endStates;
	}

	public static void updateProbabilities(double[] weights, HashMap<String, Integer> gradFeature2Index) {
		for (String s : gradFeature2Index.keySet()) {
			int index = gradFeature2Index.get(s);

			if (transitionProbabilities.containsKey(s)) {
				transitionProbabilities.put(s, weights[index]);
			} else if (initialProbabilities.containsKey(s)) {
				initialProbabilities.put(s, weights[index]);
			} else if (emissionProbabilities.containsKey(s)) {
				emissionProbabilities.put(s, weights[index]);
			}
		}
	}

	public static double max(double[] values) {
		double max = -Double.MAX_VALUE;
		for (double value : values) {
			if (value > max)
				max = value;
		}
		return max;
	}

	public static double logSumOfExponentials(double[] xs) {
		if (xs.length == 1)
			return xs[0];
		double max = max(xs);
		double sum = 0.0;
		for (double x : xs)
			if (x != Double.NEGATIVE_INFINITY)
				sum += Math.exp(x - max);
		return max + java.lang.Math.log(sum);
	}

	public static double logSumOfExponentials(ArrayList<Double> x) {
		double[] xs = new double[x.size()];
		for (int i = 0; i < x.size(); i++)
			xs[i] = x.get(i);
		return logSumOfExponentials(xs);
	}

}
