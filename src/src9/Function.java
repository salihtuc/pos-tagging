
/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

public class Function implements DifferentiableFunction {

	// -------------------------------------- LBFGS-B
	// ---------------------------------------------------------

	public double functionValue = 0.0;
	public static final double LAMBDA_EM = 0.0;

	public int iterCount = 0;
	ArrayList<Double> originalList = new ArrayList<>();
	ArrayList<Double> negativeList = new ArrayList<>();

	@Override
	public FunctionValues getValues(double[] point) {
		if (iterCount != 0) {
			Main.updateProbabilities(point, Main.gradFeature2Index);

			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
				for (int j = 0; j < globalMap.size(); j++) {
					HashMap<Integer, Node> targetMap = globalMap.get(j);

					// Iterate over targetMap
					Main.iterateFromStartToEnd(targetMap);
					Main.iterateFromEndToStart(targetMap);
				}
			}
		}
		iterCount++;

		FunctionValues fv = new FunctionValues(functionValue(point), gradient(point));

		return fv;
	}

	public double functionValue(double[] iterWeights) {
		// originalList.clear();
		// negativeList.clear();
		functionValue = 0.0;
		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) { // Iterate
																							// for
																							// all
																							// sentences
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> lattice : globalMap.values()) { // Each
																		// sentence's
																		// lattices
				for (Node endState : Main.returnEndStates(lattice)) {
					if (counter == 0) {
						sumNum += sumListValues(endState.tagScores); // Original
																		// lattice's
																		// score
						counter++;
					} else {
						sumDenom += sumListValues(endState.tagScores); // Negative
																		// samples'
																		// scores
					}
				}
			}

			double div = divide(sumNum, sumDenom);
			if (Double.isFinite(div) && div != 0 && div > 0)
				functionValue += Math.log(div); // Each sentence's scores added.

		}

		// Regularization
		for (double weight : iterWeights) {
			functionValue -= LAMBDA_EM * Math.pow(weight, 2);
		}
		functionValue *= -1; // -1 because of minimizing

		return functionValue;
	}

	/**
	 * @param iterWeights
	 * @return
	 */
	public double[] gradient(double[] iterWeights) {

		double[] grad = new double[iterWeights.length];

		HashMap<String, Double> gradTransitionProbabilities = new HashMap<>();
		HashMap<String, Double> gradInitialProbabilities = new HashMap<>();
		HashMap<String, Double> gradEmissionProbabilities = new HashMap<>();

		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.tagSize; j++) {

				/* Variables for transition */
				double transOriginalNum = 0.0;
				double transNeighborNum = 0.0;

				double transOriginalDenom = 0.0;
				double transNeighborDenom = 0.0;

				double orgScore = 0.0;
				double negScore = 0.0;

				ArrayList<Double> values = new ArrayList<>();
				ArrayList<Double> valuesDenom = new ArrayList<>();

				/* Variables for initials */
				ArrayList<Double> valuesInit = new ArrayList<>();
				ArrayList<Double> valuesInitDenom = new ArrayList<>();

				ArrayList<Double> valuesInitEnd = new ArrayList<>();
				ArrayList<Double> valuesInitEndDenom = new ArrayList<>();

				double initOriginalNum = 0.0;
				double initNeighborNum = 0.0;

				double initOriginalDenom = 0.0;
				double initNeighborDenom = 0.0;

				double initOriginalEndNum = 0.0;
				double initNeighborEndNum = 0.0;

				double initOriginalEndDenom = 0.0;
				double initNeighborEndDenom = 0.0;

				double initOrgScore = 0.0;
				double initNegScore = 0.0;

				double initOrgScoreEnd = 0.0;
				double initNegScoreEnd = 0.0;

				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {

					for (int a = 0; a < globalMap.size(); a++) {

						HashMap<Integer, Node> lattice = globalMap.get(a);
						ArrayList<Integer> startStates = Main.returnStartStates(lattice);

						for (int k = 1; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
								values.add(p(i, j, n, n2)); // calculate transition i-j

								for (int k2 = 0; k2 < Main.tagSize; k2++) {
									valuesDenom.add(p(i, k2, n, n2)); // calculate transition  i-k
								}

								if (n2.isEndState) { // if node end state calculate initial i-</s>
									valuesInitEnd.add(Main.initialProbabilities.get(Main.tagList.get(i) + "|</s>")); // calculate  initial  i-</s>
									for (int k2 = 0; k2 < Main.tagSize; k2++) {
										valuesInitEndDenom.add(Main.initialProbabilities.get(Main.tagList.get(k2) + "|</s>")); // calculate initial  all  tags-</s>
									}

								}
							}
							if (startStates.contains(n.stateNum)) { // if node  start  state  calculate initial  <s>-i

								double emission = 0;

								if (Main.emissionProbabilities.containsKey(Main.tagList.get(i) + "|" + n.word)) {
									emission = Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + n.word);
								}

								valuesInit.add(emission + Main.initialProbabilities.get("<s>|" + Main.tagList.get(i))); // calculate  initial  <s>-i
								for (int k2 = 0; k2 < Main.tagSize; k2++) {

									double emission2 = 0;

									if (Main.emissionProbabilities.containsKey(Main.tagList.get(k2) + "|" + n.word)) {
										emission2 = Main.emissionProbabilities.get(Main.tagList.get(k2) + "|" + n.word);
									}

									valuesInitDenom.add(emission2 + Main.initialProbabilities.get("<s>|" + Main.tagList.get(k2))); // calculate  initial  <s>-all tags
								}
							}

						}

						if (a == 0) {

							transOriginalNum += Math.exp(Main.logSumOfExponentials(values));
							transOriginalDenom += Math.exp(Main.logSumOfExponentials(valuesDenom));

							initOriginalNum += Math.exp(Main.logSumOfExponentials(valuesInit));
							initOriginalDenom += Math.exp(Main.logSumOfExponentials(valuesInitDenom));

							initOriginalEndNum += Math.exp(Main.logSumOfExponentials(valuesInitEnd));
							initOriginalEndDenom += Math.exp(Main.logSumOfExponentials(valuesInitEndDenom));
							
							values.clear();
							valuesDenom.clear();

							valuesInit.clear();
							valuesInitDenom.clear();

							valuesInitEnd.clear();
							valuesInitEndDenom.clear();

						} else {

							transNeighborNum += Math.exp(Main.logSumOfExponentials(values));
							transNeighborDenom += Math.exp(Main.logSumOfExponentials(valuesDenom));

							initNeighborNum += Math.exp(Main.logSumOfExponentials(valuesInit));
							initNeighborDenom += Math.exp(Main.logSumOfExponentials(valuesInitDenom));

							initNeighborEndNum += Math.exp(Main.logSumOfExponentials(valuesInitEnd));
							initNeighborEndDenom += Math.exp(Main.logSumOfExponentials(valuesInitEndDenom));

							values.clear();
							valuesDenom.clear();
							
							valuesInit.clear();
							valuesInitDenom.clear();

							valuesInitEnd.clear();
							valuesInitEndDenom.clear();

						}
					}

				}

				orgScore += divide(transOriginalNum, transOriginalDenom); // divide for all sentences transOriginalNum / transOriginalDenom
				negScore += divide(transNeighborNum, transNeighborDenom); // divide for all negative sentences transNeighborNum / transNeighborDenom

				/* Calculation for initial */

				initOrgScore += divide(initOriginalNum, initOriginalDenom);
				initNegScore += divide(initNeighborNum, initNeighborDenom);

				initOrgScoreEnd += divide(initOriginalEndNum, initOriginalEndDenom);
				initNegScoreEnd += divide(initNeighborEndNum, initNeighborEndDenom);

				String key = (Main.tagList.get(i) + "|" + Main.tagList.get(j));
				double score = orgScore - negScore;
				gradTransitionProbabilities.put(key, score);

				String initKey = "<s>|" + Main.tagList.get(i);
				double initScore = initOrgScore - initNegScore;
				gradInitialProbabilities.put(initKey, initScore);

				String initKeyEnd = Main.tagList.get(i) + "|</s>";
				double initScoreEnd = initOrgScoreEnd - initNegScoreEnd;
				gradInitialProbabilities.put(initKeyEnd, initScoreEnd);

			}

			/* re-estimation of emission probabilities NEW VERSION */

			HashMap<String, Double> originalwordsNum = new HashMap<String, Double>();
			HashMap<String, Double> neighborsWordNum = new HashMap<String, Double>();

			ArrayList<Double> listOrgDenom = new ArrayList<>();
			ArrayList<Double> listNegDenom = new ArrayList<>();

			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {

				for (int a = 0; a < globalMap.size(); a++) {
					HashMap<Integer, Node> lattice = globalMap.get(a);
					for (int k = 1; k < lattice.size(); k++) {
						Node node = lattice.get(k);

						if (Main.taggingDictionary.containsKey(node.word)
								&& Main.taggingDictionary.get(node.word).contains(i)) {

							double g = gamma(i, node);

							if (a == 0) {
								if (originalwordsNum.containsKey(node.word)) {
									ArrayList<Double> list = new ArrayList<>();

									list.add(originalwordsNum.get(node.word));
									list.add(g);

									originalwordsNum.put(node.word, Main.logSumOfExponentials(list));

								} else {
									originalwordsNum.put(node.word, g);
								}
								listOrgDenom.add(g);

							} else {
								if (neighborsWordNum.containsKey(node.word)) {
									ArrayList<Double> list = new ArrayList<>();

									list.add(neighborsWordNum.get(node.word));
									list.add(g);

									neighborsWordNum.put(node.word, Main.logSumOfExponentials(list));

								} else {
									neighborsWordNum.put(node.word, g);
								}

								listNegDenom.add(g);
							}

						}
					}
				}
			}

			for (String s : originalwordsNum.keySet()) {
				
				
				String key = Main.tagList.get(i) + "|" + s;
				double valOrg = originalwordsNum.get(s);
				double valNeg = neighborsWordNum.get(s);

				
				gradEmissionProbabilities.put(key, Math.exp(valOrg -  Main.logSumOfExponentials(listOrgDenom)) - Math.exp(valNeg - Main.logSumOfExponentials(listNegDenom)));

			}
			listOrgDenom.clear();
			listNegDenom.clear();

		}

		grad = Main.createGradArray(gradTransitionProbabilities, gradInitialProbabilities, gradEmissionProbabilities,
				Main.gradFeature2Index);

		for (int i = 0; i < grad.length; i++) {
			grad[i] *= -1;
		}

		double dSum = 0.0;
		double dSumPre = 0.0;
		for (int i = 0; i < grad.length; i++) {

			if (!Double.isFinite(grad[i]))
				grad[i] = 0.0;

			dSumPre += grad[i];
			grad[i] += 2 * LAMBDA_EM * iterWeights[i];
			dSum += grad[i];
		}

		Main.pw.println(
				"***********************************************NEW ITERATION*********************************************");
		Main.pw.println("Pre Total: " + dSumPre);
		System.out.println("Pre Total: " + dSumPre);
		Main.pw.println("Post Total: " + dSum);
		System.out.println("Post Total: " + dSum);

		return grad;

	}

	/**
	 * @param i
	 *            the number of state s_i
	 * @param j
	 *            the number of state s_j
	 * @return P
	 */
	public double p(int i, int j, Node node, Node next) {
		double num = 0.0;

		double emission = 0;

		if (Main.emissionProbabilities.containsKey(Main.tagList.get(i) + "|" + next.word)) {
			emission = Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + next.word);
		}

		num = Main.transitionProbabilities.get(Main.tagList.get(i) + "|" + Main.tagList.get(j)) + emission
				+ Math.log(next.beta.get(j));

		return num;
	}

	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num = 0.0;
		num = Math.log(node.alpha.get(i)) + Math.log(node.beta.get(i));

		return num;
	}

	/** divides two doubles. (0 / 0 = 0!) && (1 / 0 = 0!) */
	public double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		else
			return n / d;
	}

	public double[] updateWeights(double[] weights, double[] gradients) {

		double lambdaValue = 0.0;

		for (int i = 0; i < weights.length; i++) {
			weights[i] += (gradients[i] - lambdaValue);
		}

		return weights;
	}

	public double sumListValues(List<Double> list) {

		return Math.exp(Main.logSumOfExponentials((ArrayList<Double>) list));
	}

	public double sumListValuesLog(List<Double> list) {

		return Main.logSumOfExponentials((ArrayList<Double>) list);
	}

}
