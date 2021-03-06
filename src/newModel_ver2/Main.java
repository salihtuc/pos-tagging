import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Properties;
import java.util.Random;

import org.apache.commons.lang3.tuple.ImmutablePair;

public class Main {

	static HashMap<String, Integer> transitionBigramFrequencies = new HashMap<>();
	static HashMap<String, Integer> transitionUnigramFrequencies = new HashMap<>();
	// static HashMap<String, Integer> transitionTrigramFrequencies = new
	// HashMap<>();
	static HashMap<String, Integer> emissionFrequencies = new HashMap<>();

	static ArrayList<String> wordList = new ArrayList<>();

	static ArrayList<ImmutablePair<ArrayList<String>, ArrayList<String>>> sentences = new ArrayList<>();

	static ArrayList<String> tagList = new ArrayList<>();
	static int tagSize = 12;
	static int iteration = 500;
	static String inputFile = "input.txt";
	static String outputFile = "output.txt";
	static String randomOutputFile = "randomOut.txt";

	static double prob = 0.0;

	public static void main(String[] args) {

		Properties prop = new Properties();
		InputStream input = null;
		String paramsFile;
		if (args.length > 0) {
			paramsFile = args[0];
		} else {
			paramsFile = "params.properties";
		}

		// get params
		try {

			input = new FileInputStream(paramsFile);

			// load a properties file
			prop.load(input);

			// get the property values
			inputFile = prop.getProperty("inputFile");
			outputFile = prop.getProperty("outputFile");
			randomOutputFile = prop.getProperty("randomOutputFile");
			tagSize = Integer.parseInt(prop.getProperty("tagSize"));
			iteration = Integer.parseInt(prop.getProperty("iteration"));

		} catch (IOException ex) {
			ex.printStackTrace();
		} finally {
			if (input != null) {
				try {
					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		fillTagList();
		readFile(inputFile);
		writeFile(randomOutputFile);

		int counter = 0;
		prob = calculateProbability1();
		while (counter < iteration) {
			iterate();

			System.out.println("At iteration: " + (counter + 1) + " Probability " + prob);
			counter++;
		}

		writeFile(outputFile);

	}

	/**
	 * This method is using for reading the input file and filling frequencies
	 * maps and sentences list.
	 * 
	 * @input: String file: The input file
	 * 
	 */
	private static void readFile(String file) {
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			// String start = "<start> ";
			String end = " <end>";

			String line;
			String prevTag = "";

			while ((line = br.readLine()) != null) {

				ArrayList<String> sentenceWords = new ArrayList<>();
				ArrayList<String> sentenceTags = new ArrayList<>();

				line = line.trim().toLowerCase();
				line = line + end;

				prevTag = "<s>";

				for (String word : line.split(" ")) {
					if (!wordList.contains(word))
						wordList.add(word);

					String tag = "";
					if (word.equals(end.trim()))
						tag = "</s>";
					else
						tag = returnRandomTag();

					sentenceWords.add(word);
					sentenceTags.add(tag);

					update(prevTag, tag, word, "", 1);

					prevTag = tag;
				}

				// update(prevTag, "</s>", "<end>", "", 1);

				sentences.add(new ImmutablePair<ArrayList<String>, ArrayList<String>>(sentenceWords, sentenceTags));

			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		// transitionFrequencies.put("<s>", sentences.size());
		// wordList.add("<end>");
	}

	private static void writeFile(String file) {
		PrintWriter pw = null;

		try {
			pw = new PrintWriter(new File(file));

			for (int i = 0; i < sentences.size(); i++) {
				ImmutablePair<ArrayList<String>, ArrayList<String>> sentence = sentences.get(i);

				for (int j = 0; j < sentence.left.size()-1; j++) {
					String word = sentence.left.get(j);
					String tag = sentence.right.get(j);

					pw.print(word + "/" + tag + " ");
				}
				pw.println();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} finally {
			pw.close();
		}
	}

	/**
	 * This method is using for creating tags depend on tagSize above and
	 * filling the tag list.
	 * 
	 */
	private static void fillTagList() {

		for (int i = 0; i < tagSize; i++) {
			tagList.add("t" + i);
		}

		tagSize = tagList.size();
	}

	/**
	 * This method is using for getting a random tag from the tag list.
	 * 
	 */
	private static String returnRandomTag() {
		Random r = new Random();

		return tagList.get(r.nextInt(tagSize));
	}

	/**
	 * This method is using for updating frequencies maps depend on the current
	 * word, current tag, previous tag and next tag.
	 * 
	 * @input: String prevTag: Previous word's tag.
	 * @input: String tag: Current tag
	 * @input: String word: Current word
	 * @input: String nextTag: Next word's tag
	 * 
	 */
	private static void update(String prevTag, String tag, String word, String nextTag, int value) {

		// Unigram Probabilities
		if (transitionUnigramFrequencies.containsKey(tag)) {
			transitionUnigramFrequencies.put(tag, transitionUnigramFrequencies.get(tag) + value);
		} else {
			transitionUnigramFrequencies.put(tag, value);
		}

		// Bigram Probabilities
		String bigramFeature = tag + "|" + prevTag;
		if (transitionBigramFrequencies.containsKey(bigramFeature)) {
			transitionBigramFrequencies.put(bigramFeature, transitionBigramFrequencies.get(bigramFeature) + value);
		} else {
			transitionBigramFrequencies.put(bigramFeature, value);
		}

		if (!nextTag.equals("")) {
			bigramFeature = nextTag + "|" + tag;
			if (transitionBigramFrequencies.containsKey(bigramFeature)) {
				transitionBigramFrequencies.put(bigramFeature, transitionBigramFrequencies.get(bigramFeature) + value);
			} else {
				transitionBigramFrequencies.put(bigramFeature, value);
			}
		}

		// Emission Probabilities
		String emissionFeature = word + "|" + tag;
		if (emissionFrequencies.containsKey(emissionFeature)) {
			emissionFrequencies.put(emissionFeature, emissionFrequencies.get(emissionFeature) + value);
		} else {
			emissionFrequencies.put(emissionFeature, value);
		}
	}

	private static boolean sampling(double prob) {
		// XXX: >= or >
		if (prob >= 1) {
			return true;
		} else {
			Random r = new Random();
			double randomProb = r.nextDouble();

			if (prob < randomProb) {
				return false;
			} else {
				return true;
			}

		}

	}

	private static void iterate() {
		for (int i = 0; i < sentences.size(); i++) {
			ImmutablePair<ArrayList<String>, ArrayList<String>> sentence = sentences.get(i);
			String prevTag = "<s>";
			int j = 0;

			for (j = 0; j < sentence.left.size() - 1; j++) {
				String word = sentence.left.get(j);
				String tag = sentence.right.get(j);
				String nextTag = sentence.right.get(j + 1);

				double pOld = prob;

				update(prevTag, tag, word, nextTag, -1);
				calculateChanges(tag, "negative");
				String newTag = returnRandomTag();

				update(prevTag, newTag, word, nextTag, 1);
				calculateChanges(newTag, "positive");

				double pNew = prob;

				boolean isAccepted = sampling(divide(pNew, pOld));

				if (!isAccepted) {

					update(prevTag, newTag, word, nextTag, -1);
					calculateChanges(newTag, "negative");

					update(prevTag, tag, word, nextTag, 1);
					calculateChanges(tag, "positive");

					prevTag = tag;
				} else {
					sentence.right.set(j, newTag);
					prevTag = newTag;
				}

			}
			// double pOld = prob;
			//
			// update(prevTag, sentence.right.get(j), sentence.left.get(j),
			// "</s>", -1);
			// calculateChanges(sentence.right.get(j), "negative");
			// String newTag = returnRandomTag();
			//
			// update(prevTag, newTag, sentence.left.get(j), "</s>", 1);
			// calculateChanges(newTag, "positive");
			//
			// double pNew = prob;
			//
			// boolean isAccepted = sampling(divide(pNew, pOld));
			//
			// if (!isAccepted) {
			// update(prevTag, newTag, sentence.left.get(j), "</s>", -1);
			// calculateChanges(newTag, "negative");
			// update(prevTag, sentence.right.get(j), sentence.left.get(j),
			// "</s>", 1);
			// calculateChanges(sentence.right.get(j), "positive");
			// } else {
			// sentence.right.set(j, newTag);
			// }

		}
	}

	private static double calculateProbabilityTransition(HashMap<String, Integer> transitionbigramFrequencies, HashMap<String, Integer> transitionunigramFrequencies, String tag,
			String prevTag, int value) {
		String key = tag + "|" + prevTag;

		if (transitionbigramFrequencies.containsKey(key) && transitionunigramFrequencies.containsKey(prevTag))
			return divide(transitionbigramFrequencies.get(key), transitionunigramFrequencies.get(prevTag) + value);
		else
			return 0;
	}

	private static double calculateProbabilityEmission(HashMap<String, Integer> transitionunigramFrequencies,
			HashMap<String, Integer> emissionFrequencies, String tag, String word, int value) {
		String key = word + "|" + tag;

		if (emissionFrequencies.containsKey(key) && transitionunigramFrequencies.containsKey(tag))
			return divide(emissionFrequencies.get(key), transitionunigramFrequencies.get(tag) + value);
		else
			return 0;
	}

	/*private static double calculateProbability() {
		double prob = 0;
		for (int i = 0; i < sentences.size(); i++) {
			ImmutablePair<ArrayList<String>, ArrayList<String>> sentence = sentences.get(i);
			String prevTag = "<s>";

			for (int j = 0; j < sentence.left.size() - 1; j++) {
				String word = sentence.left.get(j);
				String tag = sentence.right.get(j);

				prob += calculateProbabilityTransition(transitionFrequencies, tag, prevTag);
				prob += calculateProbabilityEmission(transitionFrequencies, emissionFrequencies, tag, word);

				prevTag = tag;
			}
			prob += calculateProbabilityTransition(transitionFrequencies, "</s>", prevTag);
		}

		return (prob);
	}
*/
	private static double calculateProbability1() {
		double prob = 0;
		
		for (String tag : tagList){
			for (String prevTag : tagList) {
				prob += getFrequency(transitionBigramFrequencies, tag + "|" + prevTag) * calculateProbabilityTransition(transitionBigramFrequencies,transitionUnigramFrequencies, tag, prevTag, 0);
			}
			for (String word : wordList) {
				prob += getFrequency(emissionFrequencies, word + "|" + tag) * calculateProbabilityEmission(transitionUnigramFrequencies, emissionFrequencies, tag, word, 0);
			}
		}

		return (prob);
	}
	
	private static int getFrequency(HashMap<String, Integer> map, String key) {
		if(map.containsKey(key))
			return map.get(key);
		else
			return 0;
	}

	private static void calculateChanges(String tag, String decide) {
		for (int i = 0; i < tagList.size(); i++) {
			String tagPrev = tagList.get(i);

			if (decide.equals("negative")) {
				prob -= (transitionBigramFrequencies.get(tag + "|" + tagPrev))
						* calculateProbabilityTransition(transitionBigramFrequencies,transitionUnigramFrequencies, tagPrev, tag, 1);
				prob -= (transitionBigramFrequencies.get(tag + "|" + tagPrev)) * calculateProbabilityTransition(transitionBigramFrequencies,transitionUnigramFrequencies, tag, tagPrev, 1);
				// prob -= (transitionFrequencies.get(tag)+1) *
				// calculateProbabilityEmission(transitionFrequencies,
				// emissionFrequencies, tag, word);
			} else {
				prob += (transitionBigramFrequencies.get(tag + "|" + tagPrev))
						* calculateProbabilityTransition(transitionBigramFrequencies,transitionUnigramFrequencies, tagPrev, tag, -1);
				
				prob += (transitionBigramFrequencies.get(tag + "|" + tagPrev)) * calculateProbabilityTransition(transitionBigramFrequencies,transitionUnigramFrequencies, tag, tagPrev, -1);
				// prob += calculateProbabilityEmission(transitionFrequencies,
				// emissionFrequencies, tag, word);
			}

		}

		for (int i = 0; i < wordList.size(); i++) {
			String word = wordList.get(i);

			if (decide.equals("negative")) {
				prob -= getFrequency(emissionFrequencies, (word + "|" + tag))
						* calculateProbabilityEmission(transitionUnigramFrequencies, emissionFrequencies, tag, word, 1);
			} else {
				prob += getFrequency(emissionFrequencies, (word + "|" + tag))
						* calculateProbabilityEmission(transitionUnigramFrequencies, emissionFrequencies, tag, word, -1);
			}
		}

	}

	public static double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		else
			return n / d;
	}
}
