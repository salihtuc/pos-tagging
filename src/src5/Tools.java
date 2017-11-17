import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class Tools {

	static double dot(String a, String b) {
		double sum = 0.;
		double norm1 = 0, norm2 = 0;

		ArrayList<Double> vec1 = JointModel.wordVec.get(a);
		ArrayList<Double> vec2 = JointModel.wordVec.get(b);

		if (vec1 == null || vec2 == null)
			return -0.5; // FIXME - make sure not using just plain DOT

		for (int i = 0; i < vec1.size(); i++) {
			sum += vec1.get(i) * vec2.get(i);

			norm1 += Math.pow(vec1.get(i), 2);
			norm2 += Math.pow(vec2.get(i), 2);
		}

		return divide(sum, (Math.sqrt(norm1) * Math.sqrt(norm2)));
	}

	static int getFeatureIndex(String feature, String decideAffix) {
		if (!JointModel.feature2Index.containsKey(feature)) {
			if (JointModel.TEST)
				return -1; // if in testing phase, and feature does not exist already, do not create new

			int index = JointModel.feature2Index.size();
			JointModel.feature2Index.put(feature, index);
			JointModel.index2Feature.add(feature);

			if(JointModel.generalEmissionProbabilities.containsKey(feature))
				JointModel.coarseProbabilities.put(feature, JointModel.generalEmissionProbabilities.get(feature));
			else if(JointModel.suffixProbabilities.containsKey(feature) && decideAffix.equals("suffix"))
				JointModel.coarseProbabilities.put(feature, JointModel.suffixProbabilities.get(feature));
			else if(JointModel.prefixProbabilities.containsKey(feature) && decideAffix.equals("prefix"))
				JointModel.coarseProbabilities.put(feature, JointModel.prefixProbabilities.get(feature));
			else if(decideAffix.equals("stopSuffix")) {
				String[] parts = feature.split("\\|STP_E_");
				String tag = parts[0];
				String affix = parts[1];
				
				String suffixFeature = tag + "|SUFFIX_" + affix;
				
				if(JointModel.suffixProbabilities.containsKey(suffixFeature)) {
					JointModel.coarseProbabilities.put(feature, JointModel.suffixProbabilities.get(suffixFeature));
				}
				else {
					JointModel.coarseProbabilities.put(feature, 0.00000000001);
				}
			}
			else if(decideAffix.equals("stopPrefix")) {
				String[] parts = feature.split("\\|STP_B_");
				String tag = parts[0];
				String affix = parts[1];
				
				String prefixFeature = tag + "|PREFIX_" + affix;
				
				if(JointModel.prefixProbabilities.containsKey(prefixFeature)) {
					JointModel.coarseProbabilities.put(feature, JointModel.prefixProbabilities.get(prefixFeature));
				}
				else {
					JointModel.coarseProbabilities.put(feature, 0.00000000001);
				}
			}
			else {
				JointModel.coarseProbabilities.put(feature, 0.00000000001);
			}

			return index;
		}
		return JointModel.feature2Index.get(feature);

	}

	static void addFeature(HashMap<Integer, Double> features, String newFeature, double value, String decide, String decideAffix) {
		int featureIndex;
		if (decide.equals("tagDependent")) {
			for (String tagDependentFeature : returnTagDependentFeatures(newFeature)) {
				featureIndex = getFeatureIndex(tagDependentFeature, decideAffix);
				if (featureIndex != -1)
					features.put(featureIndex, value);
			}
		} else {
			featureIndex = getFeatureIndex(newFeature, decideAffix);
			if (featureIndex != -1)
				features.put(featureIndex, value);
		}

	}
	
	static ArrayList<String> returnTagDependentFeatures(String feature) {
		ArrayList<String> tagDependentFeatures = new ArrayList<>();

		for (String tag : JointModel.tagList) {
			String newFeature = tag + "|" + feature;
			tagDependentFeatures.add(newFeature);
		}

		return tagDependentFeatures;
	}

	protected static void printDoubleArray(double[] array) {
		for (double d : array) {
			System.out.print(d + " ");
			// Main.pw.print(d + " ");
		}
		System.out.println();
		// Main.pw.println();
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

	// dot product of feature and weights(global)
	static double featureWeightProduct(HashMap<Integer, Double> features) {
		double sum = 0.;
		if (features == null || features.size() == 0)
			return 0.;
		for (int i : features.keySet())
			if (i < JointModel.coarseProbabilities.size()) {// check if weight exists for the feature
				sum += features.get(i) * JointModel.coarseProbabilities.get(JointModel.index2Feature.get(i));
			}
		return sum;
	}

	// dot product of feature and weights(global)
	static double featureWeightProductWithTag(HashMap<Integer, Double> features, String tag) {
		double sum = 0.;
		if (features == null || features.size() == 0)
			return 0.;
		for (int i : features.keySet()) {
			
			String feature = JointModel.index2Feature.get(i);

			if (feature.startsWith(tag + "|")) {	// Tag dependent
				if (i < JointModel.coarseProbabilities.size()) { // check if weight exists for the feature
					sum += features.get(i) * JointModel.coarseProbabilities.get(feature);
				}
			}
			else if(!feature.contains("|")) {	// Not dependent
				if (i < JointModel.coarseProbabilities.size()) { // check if weight exists for the feature
					sum += features.get(i) * JointModel.coarseProbabilities.get(feature);
				}
			}
		}
		return sum;
	}

	/** divides two doubles. (0 / 0 = 0!) && (1 / 0 = 0!) */
	public static double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		else
			return n / d;
	}

	public static double sumDoubleListValues(Collection<Double> list) {
		double sum = 0.0;

		for (double d : list) {
			sum += d;
		}

		return sum;
	}

	static HashMap<String, Map<String, Double>> computeAffixCorrelation(LinkedHashSet<String> affixes, char type)
			throws IOException {
		System.out.print("Computing affix correlation - " + type + " ...");
		String[] affixArray = new String[affixes.size()];
		affixArray = affixes.toArray(affixArray);
		double[][] correlationMatrix = new double[affixArray.length][affixArray.length];
		HashMap<String, HashSet<String>> affix2Word = new HashMap<String, HashSet<String>>();
		for (String affix : affixArray) {
			affix2Word.put(affix, new HashSet<String>());
			for (String word : JointModel.word2Cnt.keySet())
				if (type == 's' && word.endsWith(affix))
					affix2Word.get(affix).add(word.substring(0, word.length() - affix.length()));
				else if (type == 'p' && word.startsWith(affix))
					affix2Word.get(affix).add(word.substring(affix.length()));
		}

		HashMap<String, Map<String, Double>> affixNeighbor = new HashMap<String, Map<String, Double>>();
		for (int i = 0; i < affixArray.length; i++) {
			int bestJ = 0;
			double bestCorrelation = 0.;
			HashMap<String, Double> neighbor2Score = new HashMap<String, Double>();
			for (int j = 0; j < affixArray.length; j++)
				if (i != j) {
					HashSet<String> tmp = Tools.clone(affix2Word.get(affixArray[i]));
					tmp.retainAll(affix2Word.get(affixArray[j]));
					correlationMatrix[i][j] = ((double) tmp.size()) / affix2Word.get(affixArray[i]).size();
					neighbor2Score.put(affixArray[j], correlationMatrix[i][j]);
					if (correlationMatrix[i][j] > bestCorrelation) {
						bestCorrelation = correlationMatrix[i][j];
						bestJ = j;
					}
				}

			affixNeighbor.put(affixArray[i], Tools.sortByValue(neighbor2Score));
		}

		System.out.println("done.");
		return affixNeighbor;
	}

	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
		List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>(map.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<K, V>>() {
			public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2) {
				return -(o1.getValue()).compareTo(o2.getValue()); // change sign to make ascending
			}
		});

		Map<K, V> result = new LinkedHashMap<K, V>();
		for (Map.Entry<K, V> entry : list) {
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}

	static HashSet<String> clone(HashSet<String> map) {
		HashSet<String> newMap = new HashSet<String>();
		for (String key : map)
			newMap.add(key);
		return newMap;
	}

	static void incrementMap(Map<String, Double> map, String key) {
		Double value = map.get(key);
		if (value == null)
			map.put(key, 1.);
		else
			map.put(key, value + 1);
	}

}
