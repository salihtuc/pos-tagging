import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.MutablePair;
import org.apache.commons.lang3.tuple.Pair;

public class JointModel {

	// Variables from params.properties
	static String wordVectorFile = "wordVectorsEng.txt";
	public static double generalInitialValue = 0.0;
	static double STOP_FACTOR = 1;
	static int FREQ_THRESHOLD = 0;

	// enum constants
	static int STOP = 1;
	static int REPEAT = 2;
	static int MODIFY = 3;
	static int SUFFIX = 4;
	static int PREFIX = 5;
	static int COMPOUND = 6;
	static int DELETE = 7;

	// conditions for features/run types
	static final boolean REPEAT_ON = true;
	static final boolean MODIFY_ON = false;
	static final boolean DELETE_ON = true;
	static final boolean INDUCTIVE = true;
	static final boolean DOT = true;
	static boolean AFFIX_NEIGHBORS = true;
	static boolean AFFIX_CONTEXT = false; // NOT for neighbors
	static boolean TEST = false; // turn on while testing to avoid creating new features - don't change manually

	static int AFFIX_FREQ_THRESHOLD;
	static int MIN_SEG_LENGTH;
	static int TOP_AFFIX_SELECT;
	static int MAX_AFFIX_LENGTH = 4;
	static int TOP_AFFIX_NEIGHBORS = 1;

	static int HEURISTIC_FREQ_THRESHOLD;
	static int VECTOR_SIZE = 200;
	static int tagSize;

	// language specific
	static String ALPHABET = "qwertyuiopasdfghjklzxcvbnm";

	// General lists about words & sentences
	public static ArrayList<String> sentences = new ArrayList<>();
	public static ArrayList<String> words = new ArrayList<>();
	public static ArrayList<String> allWords = new ArrayList<>();
	public static HashSet<String> uniqueValues;

	// Tag and Feature lists
	public static ArrayList<String> tagList = new ArrayList<>(); // Holds tags
	public static ArrayList<String> featureList = new ArrayList<>(); // Holds features

	// Structures related with Lattices
	public static HashMap<Integer, HashMap<Integer, Node>> globalMap; // Holds all lattices
	public static ArrayList<HashMap<Integer, HashMap<Integer, Node>>> sentenceGlobals = new ArrayList<>();

	// Word vectors
	static HashMap<String, ArrayList<Double>> wordVec = new HashMap<String, ArrayList<Double>>();

	// "Word2..." structures
	public static HashMap<String, Integer> word2Cnt = new HashMap<>();
	static HashMap<String, ArrayList<String>> word2Neighbors = new HashMap<String, ArrayList<String>>();
	static HashMap<String, Double> word2MaxDot = new HashMap<String, Double>();

	// Feature weights / Probabilities
	public static HashMap<String, Double> transitionProbabilities = new HashMap<>();
	public static HashMap<String, Double> coarseProbabilities = new HashMap<>();
	public static HashMap<String, Double> initialProbabilities = new HashMap<>();

	public static HashMap<String, Double> generalEmissionProbabilities = new HashMap<>();
	public static HashMap<String, Double> generalEmissionProbabilitiesNegative = new HashMap<>();

	// index <-> feature
	static HashMap<String, Integer> feature2Index = new HashMap<String, Integer>();
	static ArrayList<String> index2Feature = new ArrayList<String>();
	public static HashMap<String, Integer> tagFeature2Index = new HashMap<>();
	public static HashMap<String, Integer> gradFeature2Index = new HashMap<>();
	public static HashMap<String, ArrayList<Integer>> feature2CoarseIndex = new HashMap<>();
	
	// affixes
	static LinkedHashSet<String> prefixes = new LinkedHashSet<String>();
	static LinkedHashSet<String> suffixes = new LinkedHashSet<String>();
	static HashMap<String, Map<String, Double>> suffixNeighbor = new HashMap<String, Map<String, Double>>();
	static HashMap<String, Map<String, Double>> prefixNeighbor = new HashMap<String, Map<String, Double>>();

	// weights
	public static double[] tagFeatureWeights = null; // Feature sized weight array that we use in LBFGS-B
	static ArrayList<Double> weights = new ArrayList<Double>();

	// caching of features
	static HashMap<String, HashMap<String, HashMap<Integer, HashMap<Integer, Double>>>> w2P2TypeFeatures = new HashMap<String, HashMap<String, HashMap<Integer, HashMap<Integer, Double>>>>();

	public static PrintWriter pw = null;

	/* Methods for initializing & Pre-computation */

	static void initialize() {
		try {
			readWordVectors();
			readInputFile();
			fillTagList(); // Creating tags (filling tagList)
			fillTransitionMap(); // Creating transitionProbabilities map
			fillInitialMap();
			fillEmissionMap();

			fillInitialFromFile("initialProb.txt");
			
			selectMostFrequentAffixes();

			int i = 0;

			System.err.println("Initializing features....");
			for (String word : word2Cnt.keySet()) {
				if (word2Cnt.get(word) < FREQ_THRESHOLD)
					continue;

				for (String neighbor : getNeighbors(word)) {
					for (Pair<String, Integer> parentAndType : getCandidates(neighbor)) {
						getFeatures(neighbor, parentAndType.getKey(), parentAndType.getValue());
					}
				}
				System.err.print("\r" + (i++));
			}
			System.err.println();
			
			fillCoarseMap();

			initLattice();
			
			tagFeatureWeights = createWeightsArray(transitionProbabilities, initialProbabilities, coarseProbabilities, gradFeature2Index);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	static void readInputFile() {
		try (BufferedReader br = new BufferedReader(new FileReader("input.txt"))) {
			String start = "<start> ";
			String line;

			while ((line = br.readLine()) != null) {

				if (line.split(" ").length > 2) {
					sentences.add(start + line.toLowerCase());
					fillGeneralEmissions(line.toLowerCase());

					allWords.addAll(Arrays.asList((line.toLowerCase()).trim().split(" ")));
				}

			}

			uniqueValues = new HashSet<>(allWords);

			for (String s : uniqueValues) {
				word2Cnt.put(s, Collections.frequency(allWords, s));
			}
//			System.out.println(word2Cnt);

			allWords.clear();
			allWords.addAll(uniqueValues);

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	static void readWordVectors() throws IOException, InterruptedException {
		BufferedReader br = new BufferedReader(new FileReader(wordVectorFile));
		try {
			StringBuilder sb;
			String line = br.readLine();
			while (line != null) {
				sb = new StringBuilder();
				sb.append(line);
				String[] parts = sb.toString().split(" ");
				ArrayList<Double> vector = new ArrayList<Double>();
				String word = parts[0];
				for (int i = 1; i < Math.min(VECTOR_SIZE + 1, parts.length); i++) {
					vector.add(Double.parseDouble(parts[i]));
				}
				wordVec.put(word, vector);
				line = br.readLine();
			}
		} finally {
			br.close();
		}
		System.err.println("Read in " + Integer.toString(wordVec.size()) + " vectors");
	}

	private static void fillTagList() {

		tagList.add("t0");
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
		tagList.add("t13");
		tagList.add("t14");
		tagList.add("t15");
		tagList.add("t16");

		// tagList.add("ADJ");
		// tagList.add("ADV");
		// tagList.add("CONJ");
		// tagList.add("DET");
		// tagList.add("ENDPUNC");
		// tagList.add("INPUNC");
		// tagList.add("LPUNC");
		// tagList.add("RPUNC");
		// tagList.add("N");
		// tagList.add("POS");
		// tagList.add("PREP");
		// tagList.add("PRT");
		// tagList.add("TO");
		// tagList.add("W");
		// tagList.add("V");
		// tagList.add("VBN");
		// tagList.add("VBG");

		tagSize = tagList.size();
	}

	private static void fillTransitionMap() {

		// double value = 1.0 / (tagSize); // For uniform values
		Random r = new Random();
		double value = 0.000000001; // For zero values
		value = generalInitialValue;

		for (int i = 0; i < (tagSize); i++) {
			for (int j = 0; j < (tagSize); j++) {
				// value = (r.nextInt(100) / 10000.0);
				// value = r.nextGaussian();
				String key = (tagList.get(i) + "|" + tagList.get(j));

				transitionProbabilities.put(key, value);

			}
		}
	}

	private static void fillInitialMap() {
		double value = 0.000000001;
		Random r = new Random();
		value = generalInitialValue;

		for (int i = 0; i < tagList.size(); i++) {
			// value = r.nextGaussian();
			// value = (r.nextInt(100) / 10000.0);
			String key1 = "<s>|" + tagList.get(i);
			String key2 = tagList.get(i) + "|</s>";

			initialProbabilities.put(key1, value);
			initialProbabilities.put(key2, value);
		}
	}
	
	private static void fillEmissionMap() {
		Random r = new Random();
		double value = 0.000000001; // For zero values
		value = generalInitialValue;

		for (int i = 0; i < (tagSize); i++) {
			for (String word : allWords) {
				String key;
				
				for(String neighbor : getNeighbors(word)) {
					key = (tagList.get(i) + "|" + neighbor);
					generalEmissionProbabilitiesNegative.put(key, value);
				}
				
				key = (tagList.get(i) + "|" + word);
				generalEmissionProbabilities.put(key, value);
			}
		}
	}
	
	private static void fillCoarseMap() {
		for(String feature : feature2Index.keySet()) {
			coarseProbabilities.put(feature, 0.0);
		}
	}
	
	private static void fillInitialFromFile(String fileName) {
		try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
			String line;
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
					if (generalEmissionProbabilities.containsKey(key)) {
						generalEmissionProbabilities.put(key, prob);
					}
				}
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	static void initLattice() {
		for (String sentence : sentences) {

			globalMap = new HashMap<>();
			// Create all lattices and put them into the globalMap
			HashMap<Integer, Node> latticeMap = new HashMap<>(); // Original sentence
			fillLattice(sentence.split(" "), latticeMap, 0);
			globalMap.put(0, latticeMap);

			// HashMap<Integer, Node> latticeMap1 = new HashMap<>(); // Negative sample 1
			// (DEL1WORD)
			// fillLattice(sentence.split(" "), latticeMap1, 1);
			// globalMap.put(1, latticeMap1);

			HashMap<Integer, Node> latticeMap2 = new HashMap<>(); // Negative sample 2 (TRANS)
			fillLattice(sentence.split(" "), latticeMap2, 2);
			globalMap.put(1, latticeMap2);

			// HashMap<Integer, Node> latticeMap3 = new HashMap<>(); // Negative sample 3
			// (DEL1SUBSEQ)
			// fillLattice(sentence.split(" "), latticeMap3, 3);
			// globalMap.put(1, latticeMap3);

			sentenceGlobals.add(globalMap);

			for (int j = 0; j < globalMap.size(); j++) {
				HashMap<Integer, Node> targetMap = globalMap.get(j);

				String decideOriginal;

				if (j == 0) {
					decideOriginal = "original";
				} else {
					decideOriginal = "negative";
				}

				// Iterate over targetMap
				iterateFromStartToEnd(targetMap, decideOriginal);
				iterateFromEndToStart(targetMap, decideOriginal);

			}
		}
	}

	public static void fillLattice(String[] words, HashMap<Integer, Node> lattice, int latticeType) {
		int N = words.length - 1;
		if (N > 2) {
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
	}

	// Iteration from start to end. Using for alpha values.
	public static void iterateFromStartToEnd(HashMap<Integer, Node> latticeMap, String decideOriginal) {
		ArrayList<Integer> startStates = returnStartStates(latticeMap);
		for (Node node : latticeMap.values()) {

			if (startStates.contains(node.stateNum)) { // First nodes for calculation
				node.alpha = createInitialList(node.word, "alpha");
			} else if (node.next.isEmpty()) { // End node
				node.beta = createInitialList(node.word, "beta");
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha", decideOriginal);
				node.tagScores = multiply(node.alpha, node.beta);
			} else if (node.prev.isEmpty()) { // Start node

			} else { // Others
				node.alpha = calculateValue(node.prev, node, latticeMap, "alpha", decideOriginal);
			}

			latticeMap.put(node.stateNum, node);
		}
	}

	// Iteration from end to start. Using for beta values
	public static void iterateFromEndToStart(HashMap<Integer, Node> latticeMap, String decideOriginal) {
		List<Node> list = new ArrayList<Node>(latticeMap.values());
		ListIterator<Node> iterator = list.listIterator(list.size());

		while (iterator.hasPrevious()) {
			Node node = (Node) iterator.previous();
			if (!node.next.isEmpty() && !node.prev.isEmpty()) {
				node.beta = calculateValue(node.next, node, latticeMap, "beta", decideOriginal);
				node.tagScores = multiply(node.alpha, node.beta);
			}

			latticeMap.put(node.stateNum, node);
		}
	}

	// Calculates values for alpha or beta.
	public static ArrayList<Double> calculateValue(List<Integer> neighbors, Node n, HashMap<Integer, Node> latticeMap,
			String decide, String decideOriginal) {

		if (decide.equals("alpha")) {
			n.alpha.clear();

			for (int i = 0; i < tagSize; i++) {
				double finalResult = 0;
				for (int counter : neighbors) {
					Node n2 = latticeMap.get(counter);
					finalResult += Math.exp(calculate(n, n2, i, decide, decideOriginal));

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
					finalResult += Math.exp(calculate(n, n2, i, decide, decideOriginal));

				}
				n.beta.add(finalResult);
			}
			return n.beta;
		}
	}

	// Using for calculating values of a node
	private static double calculate(Node n, Node n2, int tagNumber, String decide, String decideOriginal) {
		double sum = 0;
		ArrayList<Double> list = new ArrayList<>();
		for (int j = 0; j < tagSize; j++) {

			double emission = 0, transition = 0;

			if (decide.equals("alpha")) {

				transition = transitionProbabilities.get(tagList.get(j) + "|" + tagList.get(tagNumber));

				if (decideOriginal.equals("original"))
					emission = generalEmissionProbabilities.get(tagList.get(tagNumber) + "|" + n.word);
				else
					emission = generalEmissionProbabilitiesNegative.get(tagList.get(tagNumber) + "|" + n.word);

				sum = transition + emission + Math.log(n2.alpha.get(j));
				list.add(sum);

			} else {

				transition = transitionProbabilities.get(tagList.get(tagNumber) + "|" + tagList.get(j));

				if (decideOriginal.equals("original"))
					emission = generalEmissionProbabilities.get(tagList.get(j) + "|" + n2.word);
				else
					emission = generalEmissionProbabilitiesNegative.get(tagList.get(j) + "|" + n2.word);

				sum = transition + emission + Math.log(n2.beta.get(j));
				list.add(sum);

			}

		}
		sum = Tools.logSumOfExponentials(list);
		list.clear();
		return sum;

	}

	// Using for calculating scores.
	private static ArrayList<Double> multiply(ArrayList<Double> alpha, ArrayList<Double> beta) {
		ArrayList<Double> scores = new ArrayList<>();

		for (int i = 0; i < tagSize; i++) {
			double score = Math.log(alpha.get(i)) + Math.log(beta.get(i));

			scores.add(i, score);
		}

		return scores;
	}

	// For creating list values uniformly
	private static ArrayList<Double> createInitialList(String word, String decide) {
		int size = tagSize;
		ArrayList<Double> list = new ArrayList<>();

		if (decide.equals("alpha")) {
			for (int i = 0; i < size; i++) {

				double emission = 0;

				if (generalEmissionProbabilities.containsKey(tagList.get(i) + "|" + word)) {
					emission = generalEmissionProbabilities.get(tagList.get(i) + "|" + word);
				}

				double value = Math.exp(initialProbabilities.get("<s>|" + tagList.get(i)) + (emission));

				// we put only sum of initial value + emission
				list.add(value);
			}
		} else {
			for (int i = 0; i < size; i++) {
				double value = Math.exp(initialProbabilities.get(tagList.get(i) + "|</s>"));// we do not add emission
																							// </s> |
				// <end> +
				// (emissionProbabilities.get(tagList.get(i)+
				// "|" + word));

				// we put only sum of initial value
				list.add(value);
			}
		}

		return list;
	}

	private static void fillGeneralEmissions(String sentence) {
		words = new ArrayList<String>(Arrays.asList(sentence.split(" ")));
		// double value = 1.0 / (allWords.size()); // For uniform values
		double value = 0.000000001; // For zero values

		value = generalInitialValue;
		Random r = new Random();

		for (int i = 0; i < tagSize; i++) {
			for (String word : words) {
				// value = r.nextGaussian();
				// value = (r.nextInt(100) / 10000.0);

				String key = tagList.get(i) + "|" + word;
				if (!generalEmissionProbabilities.containsKey(key)) {
					generalEmissionProbabilities.put(key, value);
				}
			}
		}
	}

	static void updateGeneralEmissions() {

		for (int i = 0; i < tagSize; i++) {
			for (String word : allWords) {
				ArrayList<Double> wordProbList = new ArrayList<>();
				for (Pair<String, Integer> candidate : getCandidates(word, false)) {
					HashMap<Integer, Double> features = getFeatures(word, candidate.getKey(),
							candidate.getValue());

					double sum = 0;
					
					for (int featureIndex : features.keySet()) {
						String feature = index2Feature.get(featureIndex);
						
						if (feature.startsWith(tagList.get(i)) || !feature.contains("|")) {
							sum += coarseProbabilities.get(feature);
						}
					}
					
					wordProbList.add(sum);
				}
				generalEmissionProbabilities.put(tagList.get(i) + "|" + word, Math.exp(Tools.logSumOfExponentials(wordProbList)));
				
				ArrayList<String> neighbors = JointModel.getNeighbors(word);
				ArrayList<Double> wordProbListNeg = new ArrayList<>();
				for (String neighbor : neighbors) {

					for (Pair<String, Integer> candidate : getCandidates(neighbor, false)) {
						HashMap<Integer, Double> features = getFeatures(neighbor, candidate.getKey(),
								candidate.getValue());

						double sum = 0;
						
						for (int featureIndex : features.keySet()) {
							String feature = index2Feature.get(featureIndex);

							if (feature.startsWith(tagList.get(i)) || !feature.contains("|")) {
								sum += coarseProbabilities.get(feature);
							}
						}
						wordProbListNeg.add(sum);
					}

				}
				generalEmissionProbabilitiesNegative.put(tagList.get(i) + "|" + word, Math.exp(Tools.logSumOfExponentials(wordProbListNeg)));
			}
		}
	}

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

	public static ArrayList<Integer> returnStartStates(HashMap<Integer, Node> lattice) {

		ArrayList<Integer> startStates = new ArrayList<>();

		for (Node n : lattice.values()) {
			if (n.prev.contains(0)) {
				startStates.add(n.stateNum);
			}
		}

		return startStates;
	}

	/* Get Methods */

	static ArrayList<String> getNeighbors(String word) {

		if (word2Neighbors.containsKey(word))
			return word2Neighbors.get(word);

		// get a random subset of neighbors (maybe by transposing characters?)
		int NUM_NEIGHBORS = 5;
		ArrayList<String> neighbors = new ArrayList<String>();
		neighbors.add(word); // word is in its own neighbors set

		// generate neighbors
		HashSet<Integer> positions = new HashSet<Integer>();

		// prefix positions
		for (int i = 0; i < NUM_NEIGHBORS; i++) {
			if (word.length() > i + 1)
				positions.add(i);
			if (word.length() - i - 2 >= 0)
				positions.add(word.length() - i - 2);
		}

		for (int pos : positions) {
			// permute and add
			String newWord = word.substring(0, pos) + word.charAt(pos + 1) + word.charAt(pos);
			if (pos + 2 < word.length())
				newWord += word.substring(pos + 2);
			if (!newWord.equals(word))
				neighbors.add(newWord);
		}

		if (word.length() >= 4) {
			int n = word.length();
			if (word.length() >= 4)
				neighbors.add("" + word.charAt(1) + word.charAt(0) + word.substring(2, n - 2) + word.charAt(n - 1)
						+ word.charAt(n - 2));
			if (word.length() >= 5)
				neighbors.add("" + word.charAt(0) + word.charAt(2) + word.charAt(1) + word.substring(3, n - 2)
						+ word.charAt(n - 1) + word.charAt(n - 2));
			if (word.length() >= 6)
				neighbors.add("" + word.charAt(0) + word.charAt(2) + word.charAt(1) + word.substring(3, n - 3)
						+ word.charAt(n - 2) + word.charAt(n - 3) + word.charAt(n - 1));
			if (word.length() >= 5)
				neighbors.add("" + word.charAt(1) + word.charAt(0) + word.substring(2, n - 3) + word.charAt(n - 2)
						+ word.charAt(n - 3) + word.charAt(n - 1));
		}

		word2Neighbors.put(word, neighbors);

		return neighbors;
	}

	// generate all possible parent candidates - if heuristic is true, will use the
	// heuristic to prune
	static ArrayList<Pair<String, Integer>> getCandidates(String word, boolean heuristic) {
		ArrayList<Pair<String, Integer>> candidates = new ArrayList<Pair<String, Integer>>();

		for (int i = 1; i < word.length(); i++) {
			// suffix case
			String parent = word.substring(0, i);
			if (2 * parent.length() >= word.length()) { // be careful here
				if (!heuristic || checkHeuristic(word, parent, SUFFIX))
					candidates.add(new MutablePair<String, Integer>(parent, SUFFIX));
			}

			// REPEAT and MODIFY case
			if (2 * parent.length() >= word.length()) { // be careful here
				if (REPEAT_ON)
					if (word.charAt(i - 1) == word.charAt(i)) {
						if (!heuristic || checkHeuristic(word, parent, REPEAT))
							candidates.add(new MutablePair<String, Integer>(parent, REPEAT));
					}
				int n = parent.length();
				if (MODIFY_ON && ALPHABET.contains(parent.substring(n - 1))) {

					for (int q = 0; q < ALPHABET.length(); q++) {
						if (ALPHABET.charAt(q) == parent.charAt(n - 1))
							continue;
						String newParent = parent.substring(0, n - 1) + ALPHABET.charAt(q);

						if (word2Cnt.containsKey(newParent) && word2Cnt.get(newParent) > HEURISTIC_FREQ_THRESHOLD
								&& Tools.dot(word, newParent) > 0.2)
							candidates.add(new MutablePair<String, Integer>(newParent, MODIFY));

					}
				}
				if (DELETE_ON && i < word.length() - 1 && suffixes.contains(word.substring(i))) {

					for (int q = 0; q < ALPHABET.length(); q++) {
						String newParent = parent + ALPHABET.charAt(q);
						if (word.contains(newParent))
							continue;
						if (word2Cnt.containsKey(newParent) && word2Cnt.get(newParent) > HEURISTIC_FREQ_THRESHOLD
								&& Tools.dot(word, newParent) > 0)
							candidates.add(new MutablePair<String, Integer>(newParent, DELETE));
					}
				}
			}

			// prefix case
			parent = word.substring(i);
			if (2 * parent.length() >= word.length())
				if (!heuristic || checkHeuristic(word, parent, PREFIX))
					candidates.add(new MutablePair<String, Integer>(parent, PREFIX));
		}

		// stop case
		if (!heuristic || candidates.size() == 0)
			candidates.add(new MutablePair<String, Integer>(word, STOP));
		return candidates;
	}

	static ArrayList<Pair<String, Integer>> getCandidates(String word) {
		return getCandidates(word, false); // no heuristics used
	}

	static HashMap<Integer, Double> getFeatures(String word, String parent, int type) {
		// check feature cache
		if (checkCacheExists(word, parent, type))
			return w2P2TypeFeatures.get(word).get(parent).get(type);

		// stop features
		if (type == STOP || parent.equals(word)) {
			HashMap<Integer, Double> stopFeatures = getStopFeatures(word);
			cacheFeatures(word, parent, type, stopFeatures);
			return stopFeatures;
		}

		HashMap<Integer, Double> features = new HashMap<Integer, Double>();

		// DOT
		double cosine = Tools.dot(word, parent);

		// maxDot caching
		if (!word2MaxDot.containsKey(word))
			word2MaxDot.put(word, cosine);
		else {
			double maxDotOld = word2MaxDot.get(word);
			if (cosine > maxDotOld)
				word2MaxDot.put(word, cosine);
		}

		if (DOT)
			Tools.addFeature(features, "DOT", cosine, "other?");

		// affix
		String affix = "";
		String inVocab = "";

		if (word2Cnt.containsKey(parent) && word2Cnt.get(parent) > HEURISTIC_FREQ_THRESHOLD) {
			Tools.addFeature(features, "_IV_", Math.log(word2Cnt.get(parent)), "other?");
		} else
			Tools.addFeature(features, "_OOV_", 1., "other");

		if (type == SUFFIX) {
			// suffix case
			affix = word.substring(parent.length());
			if (affix.length() > MAX_AFFIX_LENGTH || !suffixes.contains(affix))
				affix = "UNK";
			if (!affix.equals("UNK")) {
				Tools.addFeature(features, inVocab + "SUFFIX_" + affix, 1., "tagDependent");
			}

			if (AFFIX_NEIGHBORS && suffixNeighbor.containsKey(affix)) {
				int i = 0;
				for (String neighbor : suffixNeighbor.get(affix).keySet()) {
					if (word2Cnt.containsKey(parent + neighbor)) {
						Tools.addFeature(features, "COR_S_" + affix, 1., "tagDependent");
						break;
					}
					i++;
					if (i == TOP_AFFIX_NEIGHBORS)
						break;
				}
			}

			// context features
			if (AFFIX_CONTEXT) {
				Tools.addFeature(features,
						inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 1), 1.,
						"tagDependent");

				if (parent.length() >= 2) {
					Tools.addFeature(features,
							inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 2), 1.,
							"tagDependent");
				}
			}
		} else if (type == REPEAT) {
			// assuming affix is only on the right side
			affix = word.substring(parent.length() + 1);
			if (!suffixes.contains(affix))
				affix = "UNK";
			if (!affix.equals("UNK")) {
				Tools.addFeature(features, inVocab + "SUFFIX_" + affix, 1., "tagDependent");
			}

			if (AFFIX_NEIGHBORS && suffixNeighbor.containsKey(affix)) {
				int i = 0;
				for (String neighbor : suffixNeighbor.get(affix).keySet()) {
					if (word2Cnt.containsKey(parent + neighbor)) {
						Tools.addFeature(features, "COR_S_" + affix, 1., "tagDependent");
						break;
					}
					i++;
					if (i == TOP_AFFIX_NEIGHBORS)
						break;
				}
			}

			// context features
			if (AFFIX_CONTEXT) {
				Tools.addFeature(features,
						inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 1), 1.,
						"tagDependent");
				if (parent.length() >= 2)
					Tools.addFeature(features,
							inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 2), 1.,
							"tagDependent");
			}
			// REPEAT specific features
			int parentLen = parent.length();
			Tools.addFeature(features, inVocab + "REPEAT_" + word.charAt(parentLen), 1., "other");

		} else if (type == MODIFY) { // change last letter of parent
			// assuming affix is only on the right side
			affix = word.substring(parent.length());
			if (!suffixes.contains(affix))
				affix = "UNK";
			if (!affix.equals("UNK"))
				Tools.addFeature(features, inVocab + "SUFFIX_" + affix, 1., "tagDependent");

			if (AFFIX_NEIGHBORS && suffixNeighbor.containsKey(affix)) {
				int i = 0;
				for (String neighbor : suffixNeighbor.get(affix).keySet()) {
					if (word2Cnt.containsKey(parent + neighbor)) {
						Tools.addFeature(features, "COR_S_" + affix, 1., "tagDependent");
						break;
					}
					i++;
					if (i == TOP_AFFIX_NEIGHBORS)
						break;
				}
			}

			// context features
			if (AFFIX_CONTEXT) {
				Tools.addFeature(features,
						inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 1), 1.,
						"tagDependent");
				if (parent.length() >= 2)
					Tools.addFeature(features,
							inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 2), 1.,
							"tagDependent");
			}

			// MODIFY specific features
			int parentLen = parent.length();
			Tools.addFeature(features,
					inVocab + "MODIFY_" + parent.charAt(parentLen - 1) + "_" + word.charAt(parentLen - 1), 1., "other");
		}

		else if (type == DELETE) { // add last letter of parent
			// assuming affix is only on the right side
			affix = word.substring(parent.length());
			if (!suffixes.contains(affix))
				affix = "UNK";
			if (!affix.equals("UNK"))
				Tools.addFeature(features, inVocab + "SUFFIX_" + affix, 1., "tagDependent");

			if (AFFIX_NEIGHBORS && suffixNeighbor.containsKey(affix)) {
				int i = 0;
				for (String neighbor : suffixNeighbor.get(affix).keySet()) {
					if (word2Cnt.containsKey(parent + neighbor)) {
						Tools.addFeature(features, "COR_S_" + affix, 1., "tagDependent");
						break;
					}
					i++;
					if (i == TOP_AFFIX_NEIGHBORS)
						break;
				}
			}

			// context features
			if (AFFIX_CONTEXT) {
				Tools.addFeature(features,
						inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 1), 1.,
						"tagDependent");
				if (parent.length() >= 2)
					Tools.addFeature(features,
							inVocab + "SUFFIX_" + affix + "_BOUNDARY_" + parent.substring(parent.length() - 2), 1.,
							"tagDependent");
			}

			// DELETE specific features
			int parentLen = parent.length();
			Tools.addFeature(features, inVocab + "DELETE_" + parent.charAt(parentLen - 1), 1., "other");
		}

		else if (type == PREFIX) {
			assert word.length() != parent.length();
			affix = word.substring(0, word.length() - parent.length());
			if (affix.length() > MAX_AFFIX_LENGTH || !prefixes.contains(affix))
				affix = "UNK";
			if (!affix.equals("UNK")) {
				Tools.addFeature(features, inVocab + "PREFIX_" + affix, 1, "tagDependent");
			}

			if (AFFIX_NEIGHBORS && prefixNeighbor.containsKey(affix)) {
				int i = 0;
				for (String neighbor : prefixNeighbor.get(affix).keySet()) {
					if (word2Cnt.containsKey(neighbor + parent)) {
						Tools.addFeature(features, "COR_P_" + affix, 1., "tagDependent");
						break;
					}
					i++;
					if (i == TOP_AFFIX_NEIGHBORS)
						break;
				}
			}

			// context features
			if (AFFIX_CONTEXT) {
				Tools.addFeature(features, inVocab + "PREFIX_" + affix + "_BOUNDARY_" + parent.substring(0, 1), 1.,
						"tagDependent");
				if (parent.length() >= 2)
					Tools.addFeature(features, inVocab + "PREFIX_" + affix + "_BOUNDARY_" + parent.substring(0, 2), 1.,
							"tagDependent");
			}
		}

		// BIAS feature
		Tools.addFeature(features, "BIAS", 1., "other");

		// cache features
		cacheFeatures(word, parent, type, features);

		return features;
	}

	// get most frequent affixes
	static void selectMostFrequentAffixes() throws IOException {
		HashMap<String, Integer> suffixCnt = new HashMap<String, Integer>();
		HashMap<String, Integer> prefixCnt = new HashMap<String, Integer>();
		for (String word : word2Cnt.keySet()) {
			if (word2Cnt.get(word) < FREQ_THRESHOLD)
				continue;
			for (int i = 1; i < word.length(); i++) {
				String left = word.substring(0, i);
				String right = word.substring(i);

				// suffix case
				Integer cnt = word2Cnt.get(left);
				if (cnt != null && cnt > AFFIX_FREQ_THRESHOLD && right.length() <= MAX_AFFIX_LENGTH)
					if (suffixCnt.containsKey(right))
						suffixCnt.put(right, suffixCnt.get(right) + 1);
					else
						suffixCnt.put(right, 1);

				// prefix case
				cnt = word2Cnt.get(right);
				if (cnt != null && cnt > AFFIX_FREQ_THRESHOLD && left.length() <= MAX_AFFIX_LENGTH)
					if (prefixCnt.containsKey(left))
						prefixCnt.put(left, prefixCnt.get(left) + 1);
					else
						prefixCnt.put(left, 1);
			}
		}

		// sort and take top
		Map<String, Integer> sortedSuffixes = Tools.sortByValue(suffixCnt);
		Map<String, Integer> sortedPrefixes = Tools.sortByValue(prefixCnt);
		int i = 0;
		for (String suffix : sortedSuffixes.keySet()) {
			suffixes.add(suffix);
			i++;
			if (i == TOP_AFFIX_SELECT)
				break;
		}
		i = 0;
		for (String prefix : sortedPrefixes.keySet()) {
			prefixes.add(prefix);
			i++;
			if (i == TOP_AFFIX_SELECT)
				break;
		}

		// NEW: add suffixNeighbors
		if (AFFIX_NEIGHBORS) {
			suffixNeighbor = Tools.computeAffixCorrelation(suffixes, 's');
			prefixNeighbor = Tools.computeAffixCorrelation(prefixes, 'p');
		}
	}

	// Functions for Caching features
	static void cacheFeatures(String word, String parent, int type, HashMap<Integer, Double> features) {
		if (!w2P2TypeFeatures.containsKey(word))
			w2P2TypeFeatures.put(word, new HashMap<String, HashMap<Integer, HashMap<Integer, Double>>>());
		if (!w2P2TypeFeatures.get(word).containsKey(parent))
			w2P2TypeFeatures.get(word).put(parent, new HashMap<Integer, HashMap<Integer, Double>>());
		if (!w2P2TypeFeatures.get(word).get(parent).containsKey(type))
			w2P2TypeFeatures.get(word).get(parent).put(type, features);

		// else just ignore
	}

	static HashMap<Integer, Double> getStopFeatures(String word) {
		HashMap<Integer, Double> features = new HashMap<Integer, Double>();

		if (word.length() >= 2) { // check is only to avoid null exception
			Tools.addFeature(features, "STP_E_" + word.substring(word.length() - 2), 1., "tagDependent");
			Tools.addFeature(features, "STP_B_" + word.substring(0, 2), 1., "tagDependent");
		}

		// max dot feature
		if (DOT && word.length() >= 2) {
			double maxDot = word2MaxDot.get(word);
			Tools.addFeature(features, "STP_COS_" + (int) (10 * maxDot), 1., "other");
		}

		// length feature
		Tools.addFeature(features, "STP_LEN_" + word.length(), 1., "other");

		// BIAS feature
		Tools.addFeature(features, "BIAS", 1., "other");

		return features;
	}

	/* Check Methods */

	// Heuristic function to prune out
	static boolean checkHeuristic(String word, String parent, int type) {
		double dot_threshold = 0.0;
		if (type == MODIFY)
			dot_threshold = 0.5;
		Integer wordCnt = word2Cnt.get(parent);
		if (wordCnt != null && wordCnt > HEURISTIC_FREQ_THRESHOLD)
			if (2 * parent.length() >= word.length())
				if (Tools.dot(word, parent) > dot_threshold) // the dot takes care of the non-exist case by returning 0.
					return true;
		return false;
	}

	static boolean checkCacheExists(String word, String parent, int type) {
		if (w2P2TypeFeatures.containsKey(word))
			if (w2P2TypeFeatures.get(word).containsKey(parent))
				if (w2P2TypeFeatures.get(word).get(parent).containsKey(type))
					return true;

		return false;
	}

	/* Return Methods */

	public static HashMap<Integer, Double> returnTagFeatures(String tag, double[] weights,
			HashMap<String, Integer> gradFeature2Index, String decide) {

		HashMap<Integer, Double> returnFeatures = new HashMap<>();

		for (String s : gradFeature2Index.keySet()) {
			int index = gradFeature2Index.get(s);

			if (decide.equals("transition")) {
				if (transitionProbabilities.containsKey(s) && s.startsWith(tag + "|")) {
					returnFeatures.put(index, weights[index]);
				}
			} else if (decide.equals("emission")) {
				if (coarseProbabilities.containsKey(s) && s.startsWith(tag + "|")) {
					returnFeatures.put(index, weights[index]);
				}
			} else {
				if (initialProbabilities.containsKey(s) && s.startsWith(tag + "|")) {
					returnFeatures.put(index, weights[index]);
				}
			}
		}

		return returnFeatures;
	}

	/* Update Methods */

	public static void updateProbabilities(double[] weights, HashMap<String, Integer> gradFeature2Index) {
		for (String s : gradFeature2Index.keySet()) {
			int index = gradFeature2Index.get(s);

			if (transitionProbabilities.containsKey(s)) {
				transitionProbabilities.put(s, weights[index]);
			} else if (initialProbabilities.containsKey(s)) {
				initialProbabilities.put(s, weights[index]);
			} else if (coarseProbabilities.containsKey(s)) {
				coarseProbabilities.put(s, weights[index]);
			}
		}
	}

	public static void normalizeFeatures(double[] point, String decide) {
		for (int i = 0; i < tagSize; i++) {
			HashMap<Integer, Double> features = returnTagFeatures(tagList.get(i), point, gradFeature2Index, decide);

			double sum = Tools.sumDoubleListValues(features.values());

			for (int featureIndex : features.keySet()) {
				// point[featureIndex] = divide(Math.exp(point[featureIndex]), sum);
				point[featureIndex] = Tools.divide((point[featureIndex]), sum);
			}
		}
	}

}
