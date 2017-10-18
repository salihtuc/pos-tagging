
/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import lbfgsb.LBFGSBException;
import lbfgsb.Minimizer;
//import lbfgsb.Minimizer;
import lbfgsb.Result;

public class Main {

	static PrintWriter pw = null;
	
	public static void main(String[] args) {

		try {
			pw = new PrintWriter(new File("outp.txt"));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}

//		for (String sentence : sentences)
//			fillEmissionMap(sentence); // Creating emissionProbabilities map

//		fillInitialFromFile("initialProb.txt");
		
		JointModel.initialize();

		
		// TODO
		JointModel.tagFeatureWeights = createWeightsArray(JointModel.transitionProbabilities, JointModel.initialProbabilities, JointModel.generalEmissionProbabilities,
				JointModel.gradFeature2Index);

		
		/* LBFGS-B part */
		Minimizer alg = new Minimizer();
		alg.getStopConditions().setMaxIterations(10000);
		alg.setDebugLevel(1);

		Result ret;
		try {
			ret = alg.run(new Function(), JointModel.tagFeatureWeights);

			double finalValue = ret.functionValue;
			double[] finalGradient = ret.gradient;

			System.out.println("Final Value: " + finalValue);
			Main.pw.println("Final Value: " + finalValue);
			System.out.println("Gradients:");
			Tools.printDoubleArray(finalGradient);

			JointModel.updateProbabilities(finalGradient, JointModel.gradFeature2Index);
			JointModel.tagFeatureWeights = finalGradient;

			for (String s : JointModel.transitionProbabilities.keySet()) {
				Main.pw.println(s + " " + JointModel.transitionProbabilities.get(s));
			}
			for (String s : JointModel.initialProbabilities.keySet()) {
				Main.pw.println(s + " " + JointModel.initialProbabilities.get(s));
			}
			for (String s : JointModel.generalEmissionProbabilities.keySet()) {
				Main.pw.println(s + " " + JointModel.generalEmissionProbabilities.get(s));
			}

		} catch (LBFGSBException e) {
			e.printStackTrace();
		}
		pw.close();
	}


//	private static void fillEmissionMap(String sentence) {
//		words = new ArrayList<String>(Arrays.asList(sentence.split(" ")));
//		// double value = 1.0 / (allWords.size()); // For uniform values
//		double value = 0.000000001; // For zero values
//		
//		value = generalInitialValue;
//		Random r = new Random();
//
//		for (int i = 0; i < tagSize; i++) {
//			for (String word : words) {
////				value = r.nextGaussian();
////				value = (r.nextInt(100) / 10000.0);
//
//				String key = tagList.get(i) + "|" + word;
//				if (!emissionProbabilities.containsKey(key)) {
//					emissionProbabilities.put(key, value);
//				}
//			}
//		}
//	}

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
	
//	private static void fillInitialFromFile(String fileName) {
//		try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
//			String line;
//			while ((line = br.readLine()) != null) {
//				String[] values = line.split(" ");
//				String key = values[0];
//				double prob = Double.parseDouble(values[1]);
//
//				if (key.endsWith("_t")) { // Transition or Initial
//
//					key = key.substring(0, key.length() - 2);
//
//					if (key.contains("<s>") || key.contains("</s>")) {
//						initialProbabilities.put(key, prob);
//					} else {
//						transitionProbabilities.put(key, prob);
//					}
//				} else {
//					if (emissionProbabilities.containsKey(key)) {
//						emissionProbabilities.put(key, prob);
//					}
//				}
//			}
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}

	// This function is using for creating lattices.
	static int sentenceCount = 0;

	
	public static Node returnEndState(HashMap<Integer, Node> lattice) {
		Node returnNode = null;
		for (Node n : lattice.values()) {
			if (n.next.isEmpty()) {
				returnNode = n;
			}
		}
		return returnNode;
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

}
