import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

public class Viterbi {

	private static double[][] viterbi;
	private static int[][] vitPrev;
	static int smoothed = 0;
	static double sumScore = 0;

	public static HashMap<String, String> tagWord2Morph = new HashMap<>();
	public static HashMap<String, Double> tagWord2Score = new HashMap<>();

	@SuppressWarnings("resource")
	protected static void tagging(String sysoutput) throws IOException {
		int countWord = 0;
		int cs = 0;

		BufferedWriter output = null;
		BufferedWriter outputMorph = null;

		File file = new File(sysoutput + ".txt");
		output = new BufferedWriter(new FileWriter(file));

		File fileMorph = new File(sysoutput + "Morph.txt");
		outputMorph = new BufferedWriter(new FileWriter(fileMorph));

		JointModel.tagList.add("<s>");
		JointModel.tagList.add("</s>");
		
		
		for(int i = 0; i < JointModel.tagList.size(); i++) {
			String emitStart = JointModel.tagList.get(i) + "|<start>";
			tagWord2Morph.put(emitStart, "<start>");
			
			String emit = JointModel.tagList.get(i) + "|<end>";
			tagWord2Morph.put(emit, "<end>");
			
			if(emit.equals("</s>|<end>")) {
				tagWord2Score.put("</s>|<end>", 1.0);
			}
			else if(emitStart.equals("<s>|<start>")) {
				tagWord2Score.put("<s>|<start>", 1.0);
			}
			else {
				tagWord2Score.put("</s>|<end>", 0.0);
				tagWord2Score.put("<s>|<start>", 0.0);
			}
		}

		for (String sentence : JointModel.sentences) {
			sentence = sentence.replace("<start> ", "") + " <end>";
			String[] words = sentence.toLowerCase().split("\\s+");
			viterbi = new double[JointModel.tagList.size()][words.length];
			vitPrev = new int[JointModel.tagList.size()][words.length];

			for (int i = 0; i < words.length; i++) {

				String word = words[i];
				// if (!JointModel.allWords.contains(word)) { //XXX
				// //// System.out.println(word+" smoothed");
				// leplaceSmoothing(word, String.valueOf(i));
				// }

				for (int j = 0; j < JointModel.tagList.size(); j++) {
					//// System.out.println(statess.get(j)+" "+word);
					double score = vitScore(JointModel.tagList.get(j), word, i, j);
					//// System.out.println(statess.get(j)+" "+word+"
					//// "+score);
					viterbi[j][i] = score;
				}

			}

			String tag = "";
			String tagMorph = "";
			int idx = 0;
			double max = 0;
			int n = words.length - 1;
			for (int i = 0; i < JointModel.tagList.size(); i++) {
				if (viterbi[i][n] > max) {
					max = viterbi[i][n];
					idx = i;
				}
			}
			tag = "<end>" + "/" + JointModel.tagList.get(idx) + " " + tag;

			idx = vitPrev[idx][n];
			for (int i = n - 1; i >= 0; i--) {

				String cek = JointModel.tagList.get(idx);
				countWord++;

				String key = cek + "|" + words[i];

				tag = words[i] + "/" + cek + " " + tag;
				tagMorph = words[i] + "/" + tagWord2Morph.get(key) + " " + tagMorph;

				idx = vitPrev[idx][i];

			}
			output.write(tag + "\n");
			outputMorph.write(tagMorph + "\n");
			cs++;

		}

		System.out.println("Total evaluate sentence: " + cs);
		System.out.println("Total smoothed words: " + smoothed);
		System.out.println("Total evaluate words: " + countWord);

		output.close();
		outputMorph.close();

	}

	// predict the parent of a word
	static Pair<Pair<String, Integer>, Double> predict(String word, String tag) {
		ArrayList<Pair<String, Integer>> candidateParents = JointModel.getCandidates(word); // no need to restrict to
																							// heuristics
		double bestScore = -Double.MAX_VALUE, score;
		Pair<String, Integer> bestParentAndType = null;

		// for multinomial
		ArrayList<Sample.MultinomialObject> multinomial = new ArrayList<Sample.MultinomialObject>();
		double Z = 0.;

		for (Pair<String, Integer> parentAndType : candidateParents) {
			String parent = parentAndType.getKey();
			int type = parentAndType.getValue();
			score = scoreParent(word, parent, type, tag);
			if (type == JointModel.STOP)
				score *= JointModel.STOP_FACTOR;
			if (score > bestScore) {
				bestScore = score;
				bestParentAndType = parentAndType;
			}
			multinomial.add(new Sample.MultinomialObject(parent, type, Math.exp(score)));
			Z += Math.exp(score);
		}

		// normalize the multinomial
		for (Sample.MultinomialObject obj : multinomial)
			obj.score /= Z;

		// Map<String, Double> f2W = new HashMap<String, Double>();
		// for (int featureIndex : JointModel.getFeatures(word,
		// bestParentAndType.getKey(), bestParentAndType.getValue()).keySet())
		// f2W.put(JointModel.index2Feature.get(featureIndex),
		// JointModel.feature2Weight.get(JointModel.index2Feature.get(featureIndex)));

		// System.out.println(JointModel.feature2Index);

		return new ImmutablePair<Pair<String, Integer>, Double>(bestParentAndType, divide(Math.exp(bestScore), Z));
	}

	// function to choose predictions from the sampled points
	static ArrayList<Integer> predictSampledPoints(ArrayList<Integer> sampledPoints) {
		ArrayList<Integer> points = new ArrayList<Integer>();
		int lastPoint = -1;
		for (int point : sampledPoints) {
			if (point == -1)
				break;
			if (point < lastPoint)
				break;
			points.add(point);
			lastPoint = point;
		}
		return points;
	}

	static String segment(String word, String tag) {

		// produces a segmentation
		if (word.length() < JointModel.MIN_SEG_LENGTH) {
			tagWord2Morph.put(tag + "|" + word, word);
			tagWord2Score.put(tag + "|" + word, 1.0); // FIXME

			return word;
		}

		if (word.contains("-")) {
			String[] parts = word.split("-");
			String seg = segment(parts[0], tag);
			for (int i = 1; i < parts.length; i++)
				seg += ("-" + segment(parts[i], tag));
			return seg;
		}

		// TODO
		if (word.contains("'")) {
			String[] parts = word.split("'");
			// assert parts.length<=2; //just a check
			String suffix = "'";
			// if (parts.length==2)
			suffix += parts[parts.length - 1];
			if (JointModel.suffixes.contains(suffix))
				Tools.incrementMap(JointModel.suffixDist, suffix);
			String seg = segment(parts[0], tag);
			for (int i = 1; i < parts.length; i++)
				seg += ("-'" + segment(parts[i], tag));
			return seg;
		}

		ImmutablePair<Pair<String, Integer>, Double> parentTypeAndScore = (ImmutablePair<Pair<String, Integer>, Double>) predict(
				word, tag);
		Pair<String, Integer> parentAndType = parentTypeAndScore.getLeft();
		sumScore += parentTypeAndScore.getRight();

		if (parentAndType == null) {
			tagWord2Morph.put(tag + "|" + word, word);
			tagWord2Score.put(tag + "|" + word, 1.0); // FIXME

			return word;
		}

		String parent = parentAndType.getKey();
		int type = parentAndType.getValue();

		if (type == JointModel.PUNCTUATION) {
			tagWord2Morph.put(tag + "|" + word, word);
			tagWord2Score.put(tag + "|" + word, sumScore); // FIXME

			return word;
		}

		if (type == JointModel.STOP) {
			tagWord2Morph.put(tag + "|" + word, word);
			tagWord2Score.put(tag + "|" + word, sumScore); // FIXME

			return word;
		}

		int parentLen = parent.length();

		// IMP: cases for REPEAT, MODIFY, etc here
		// else segment
		// suffix case
		if (type == JointModel.SUFFIX) {
			String suffix = word.substring(parentLen);
			if (JointModel.suffixes.contains(suffix))
				Tools.incrementMap(JointModel.suffixDist, suffix);
			return segment(parent, tag) + "-" + suffix;
		} else if (type == JointModel.REPEAT)
			return segment(parent, tag) + word.charAt(parentLen) + "-" + word.substring(parentLen + 1);
		else if (type == JointModel.MODIFY) {
			// TODO : check
			String parentSeg = segment(parent, tag);
			return parentSeg.substring(0, parentSeg.length() - 1) + word.charAt(parentLen - 1) + "-"
					+ word.substring(parentLen);
		} else if (type == JointModel.DELETE) {
			String parentSeg = segment(parent, tag);
			int parentSegLen = parentSeg.length();
			if (parentSeg.charAt(parentSegLen - 2) == '-')
				return parentSeg.substring(0, parentSegLen - 1) + word.substring(parentLen - 1);
			else
				return parentSeg.substring(0, parentSegLen - 1) + "-" + word.substring(parentLen - 1);
		}
		// prefix case
		else if (type == JointModel.PREFIX) {
			String prefix = word.substring(0, word.length() - parentLen);
			if (JointModel.prefixes.contains(prefix))
				Tools.incrementMap(JointModel.prefixDist, prefix);
			return prefix + "-" + segment(parent, tag);
		} else if (type == JointModel.COMPOUND)
			return segment(word.substring(0, parentLen), tag) + "-" + segment(word.substring(parentLen), tag);

		// null should not be returned at all. Having this to debug, instead of an
		// assert
		return null;
	}

	// returns the logScore
	static double scoreParent(String word, String parent, int type, String tag) {
		return Tools.featureWeightProduct(JointModel.getFeatures(word, parent, type));
	}

	private static void leplaceSmoothing(String word, String i) {
		Iterator it = JointModel.generalEmissionProbabilities.entrySet().iterator();
		double sum = 0;
		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			String wordd = (String) pair.getKey();
			String tag = wordd.split("\\|", 2)[0];
			double value = (double) pair.getValue();
			if (tag.equals(i)) {
				sum += value + Double.MIN_VALUE;
				// sum += value + -1*Double.MAX_VALUE;
			}
		}

		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			String wordd = (String) pair.getKey();
			String tag = wordd.split("\\|", 2)[0];
			if (tag.equals(i)) {
				JointModel.generalEmissionProbabilities.replace(wordd,
						JointModel.generalEmissionProbabilities.get(wordd) / sum);
				//// sum += value + 0.00001;
			}
		}

		JointModel.generalEmissionProbabilities.put(i + "|" + word, Double.MIN_VALUE / sum);
		// emissions.put(i + "|" + word, -1*Double.MAX_VALUE / sum);
		// JointModel.allWords.add(word);

		System.out.println("\"" + word + "\" has been smoothed");
		smoothed++;
	}

	private static double vitScore(String tag, String word, int idxWord, int idxTag) {
		double max = 0.0;
		double emProb = Double.MIN_VALUE;

		String emit = tag + "|" + word;

		if (tagWord2Score.keySet().contains(emit))
			emProb = tagWord2Score.get(emit);
		else {
			String segmentedWord = segment(word, tag);
			// System.out.println(word + "::" + tag);

			tagWord2Morph.put(emit, segmentedWord);
			tagWord2Score.put(emit, sumScore);

			sumScore = 0;
			emProb = tagWord2Score.get(emit);
		}

		if (idxWord == 0) {
			double transProb = Double.MIN_VALUE;
			String trans = "<s>|" + tag;
			if (JointModel.transitionProbabilities.keySet().contains(trans))
				transProb = JointModel.transitionProbabilities.get(trans);
			max = (1 * transProb * emProb);
			vitPrev[idxTag][idxWord] = 0;
		} else {

			for (int i = 0; i < JointModel.tagList.size(); i++) {
				String prevTag = JointModel.tagList.get(i);
				double prevVitS = viterbi[i][idxWord - 1];
				double transProb = Double.MIN_VALUE;
				// double transProb = -1*Double.MAX_VALUE;
				String trans = prevTag + "|" + tag;
				if (JointModel.transitionProbabilities.keySet().contains(trans))
					transProb = JointModel.transitionProbabilities.get(trans);
				else if(JointModel.initialProbabilities.keySet().contains(trans))
					transProb = JointModel.initialProbabilities.get(trans);

				double score = (prevVitS * transProb * emProb);
				if (score > max) {
					max = score;
					vitPrev[idxTag][idxWord] = i;
				}
			}
		}

		return max;
	}

	/** divides two doubles. (0 / 0 = 0!) && (1 / 0 = 0!) */
	public static double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		else
			return n / d;
	}

}
