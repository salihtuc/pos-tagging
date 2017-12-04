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

	static HashMap<String, Integer> transitionFrequencies = new HashMap<>();
	// static HashMap<String, Integer> transitionTrigramFrequencies = new
	// HashMap<>();
	static HashMap<String, Integer> emissionFrequencies = new HashMap<>();

	static double prob = 0;
	static int sentenceSize = 0;

	static ArrayList<ImmutablePair<ArrayList<String>, ArrayList<String>>> sentences = new ArrayList<>();

	static ArrayList<String> tagList = new ArrayList<>();
	static int tagSize = 12;
	static int iteration = 500;
	static String inputFile = "input.txt";
	static String outputFile = "output.txt";
	static String randomOutputFile = "randomOut.txt";

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
		transitionFrequencies.put("<s>", sentenceSize);
		transitionFrequencies.put("</s>", sentenceSize);
		
		writeFile(randomOutputFile);

		int counter = 0;
		prob = calculateProbability();

		while (counter < iteration) {
			iterate();

			System.out.println("At iteration: " + (counter + 1));
			counter++;
		}

		writeFile(outputFile);

	}

	/**
	 * This method is using for reading the input file and filling frequencies maps
	 * and sentences list.
	 * 
	 * @input: String file: The input file
	 * 
	 */
	private static void readFile(String file) {
		try (BufferedReader br = new BufferedReader(new FileReader(file))) {
			// String start = "<start> ";
			// String end = " <end>";

			String line;
			String prevTag = "";

			while ((line = br.readLine()) != null) {
				sentenceSize++;

				ArrayList<String> sentenceWords = new ArrayList<>();
				ArrayList<String> sentenceTags = new ArrayList<>();

				line = line.trim().toLowerCase();
				// line = start + line + end;

				prevTag = "<s>";

				for (String word : line.split(" ")) {
					String tag = returnRandomTag();

					sentenceWords.add(word);
					sentenceTags.add(tag);

					update(prevTag, tag, word, "", 1, "read");

					prevTag = tag;
				}

				update(prevTag, "</s>", "<end>", "", 1, "read");

				sentences.add(new ImmutablePair<ArrayList<String>, ArrayList<String>>(sentenceWords, sentenceTags));

			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void writeFile(String file) {
		PrintWriter pw = null;

		try {
			pw = new PrintWriter(new File(file));

			for (int i = 0; i < sentences.size(); i++) {
				ImmutablePair<ArrayList<String>, ArrayList<String>> sentence = sentences.get(i);

				for (int j = 0; j < sentence.left.size(); j++) {
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
	 * This method is using for creating tags depend on tagSize above and filling
	 * the tag list.
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
	private static void update(String prevTag, String tag, String word, String nextTag, int value, String decide) {

		double newProb = 0;
		// Unigram Probabilities
		if (transitionFrequencies.containsKey(tag)) {
			transitionFrequencies.put(tag, transitionFrequencies.get(tag) + value);
		} else {
			transitionFrequencies.put(tag, value);
		}

		// Bigram Probabilities
		String bigramFeature = tag + "|" + prevTag;
		if (transitionFrequencies.containsKey(bigramFeature)) {
			transitionFrequencies.put(bigramFeature, transitionFrequencies.get(bigramFeature) + value);
		} else {
			transitionFrequencies.put(bigramFeature, value);
		}

		if (!decide.equals("read"))
			if (value > 0)
				newProb += (transitionFrequencies.get(bigramFeature)-1) * divide(transitionFrequencies.get(bigramFeature)-1, transitionFrequencies.get(prevTag));
			else
				newProb -= (transitionFrequencies.get(bigramFeature)+1) * divide(transitionFrequencies.get(bigramFeature)+1, transitionFrequencies.get(prevTag));

		if (!nextTag.equals("")) {
			bigramFeature = nextTag + "|" + tag;
			if (transitionFrequencies.containsKey(bigramFeature)) {
				transitionFrequencies.put(bigramFeature, transitionFrequencies.get(bigramFeature) + value);
			} else {
				transitionFrequencies.put(bigramFeature, value);
			}

			if (!decide.equals("read"))
				if (value > 0)
					newProb += (transitionFrequencies.get(bigramFeature)-1) * divide(transitionFrequencies.get(bigramFeature)-1, transitionFrequencies.get(tag)-1);
				else
					newProb -= (transitionFrequencies.get(bigramFeature)+1) * divide(transitionFrequencies.get(bigramFeature)+1, transitionFrequencies.get(tag)+1);
		}

		// Emission Probabilities
		String emissionFeature = word + "|" + tag;
		if (emissionFrequencies.containsKey(emissionFeature)) {
			emissionFrequencies.put(emissionFeature, emissionFrequencies.get(emissionFeature) + value);
		} else {
			emissionFrequencies.put(emissionFeature, value);
		}

		if (!decide.equals("read"))
			if (value > 0)
				newProb += (emissionFrequencies.get(emissionFeature)-1) * divide(emissionFrequencies.get(emissionFeature)-1, transitionFrequencies.get(tag)-1);
			else
				newProb -= (emissionFrequencies.get(emissionFeature)+1) * divide(emissionFrequencies.get(emissionFeature)+1, transitionFrequencies.get(tag)+1);

		if (!decide.equals("read"))
			prob = Math.exp(Math.log(prob) + newProb);
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

				update(prevTag, tag, word, nextTag, -1, "update");
				String newTag = returnRandomTag();
				update(prevTag, newTag, word, nextTag, 1, "update");

				// double pNew = calculateProbability();

				boolean isAccepted = sampling(divide(prob, pOld));

				if (!isAccepted) {
					update(prevTag, newTag, word, nextTag, -1, "update");
					update(prevTag, tag, word, nextTag, 1, "update");
					prevTag = tag;
					prob = pOld;
				} else {
					sentence.right.set(j, newTag);
					prevTag = newTag;
				}

			}

			update(prevTag, sentence.right.get(j), sentence.left.get(j), "</s>", -1, "update");
			String newTag = returnRandomTag();
			update(prevTag, newTag, sentence.left.get(j), "</s>", 1, "update");

			double pOld = prob;
			// double pNew = calculateProbability();

			boolean isAccepted = sampling(divide(prob, pOld));

			if (!isAccepted) {
				update(prevTag, newTag, sentence.left.get(j), "</s>", -1, "update");
				update(prevTag, sentence.right.get(j), sentence.left.get(j), "</s>", 1, "update");
				prob = pOld;
			} else {
				sentence.right.set(j, newTag);
			}

		}
	}

	private static double calculateProbabilityTransition(HashMap<String, Integer> transitionFrequencies, String tag,
			String prevTag) {
		String key = tag + "|" + prevTag;

		if (transitionFrequencies.containsKey(key) && transitionFrequencies.containsKey(prevTag))
			return divide(transitionFrequencies.get(key), transitionFrequencies.get(prevTag));
		else
			return 0;
	}

	private static double calculateProbabilityEmission(HashMap<String, Integer> transitionFrequencies,
			HashMap<String, Integer> emissionFrequencies, String tag, String word) {
		String key = word + "|" + tag;

		if (emissionFrequencies.containsKey(key) && transitionFrequencies.containsKey(tag))
			return divide(emissionFrequencies.get(key), transitionFrequencies.get(tag));
		else
			return 0;
	}

	private static double calculateProbability() {
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

		return Math.exp(prob);
	}

	public static double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		else
			return n / d;
	}
}
