
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

		JointModel.initialize();

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
			
			JointModel.normalizeFeatures(finalGradient, "transition");
			JointModel.normalizeFeatures(finalGradient, "initial");
			JointModel.normalizeFeatures(finalGradient, "generalEmission");

			JointModel.updateProbabilities(finalGradient, JointModel.gradFeature2Index);
//			JointModel.updateGeneralEmissions();
			JointModel.tagFeatureWeights = finalGradient;
			
			Tools.printDoubleArray(finalGradient);

			for (String s : JointModel.transitionProbabilities.keySet()) {
				Main.pw.println(s + " " + JointModel.transitionProbabilities.get(s));
			}
			for (String s : JointModel.initialProbabilities.keySet()) {
				Main.pw.println(s + " " + JointModel.initialProbabilities.get(s));
			}
			for (String s : JointModel.coarseProbabilities.keySet()) {
				Main.pw.println(s + " " + JointModel.coarseProbabilities.get(s));
			}
			
			for (String s : JointModel.generalEmissionProbabilities.keySet()) {
				System.out.println(s + " " + JointModel.generalEmissionProbabilities.get(s));
			}

			Viterbi.tagging("outViterbi.txt");
			
		} catch (LBFGSBException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		pw.close();
	}


	protected static double[] createGradArray(HashMap<String, Double> transitionMap, HashMap<String, Double> initialMap, HashMap<String, Double> generalEmissionMap,
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
		
		for (String s : generalEmissionMap.keySet()) {
			index = gradFeature2Index.get(s);
			weights[index] = generalEmissionMap.get(s);
		}

//		for (String s : emissionMap.keySet()) {
//			index = gradFeature2Index.get(s);
//			weights[index] = emissionMap.get(s);
//		}

		return weights;
	}
	
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
