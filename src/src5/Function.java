
/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.LongSummaryStatistics;
import java.util.Map.Entry;

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

public class Function implements DifferentiableFunction {

	// -------------------------------------- LBFGS-B
	// ---------------------------------------------------------

	public double functionValue = 0.0;
	public static final double LAMBDA_EM = 1;
	public int iterCount = 0;
	ArrayList<Double> originalList = new ArrayList<>();
	ArrayList<Double> negativeList = new ArrayList<>();

	@Override
	public FunctionValues getValues(double[] point) {

		// System.out.println("Points:");
		// Main.printDoubleArray(point);

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
			// originalList.add(sumNum);
			// negativeList.add(sumDenom);
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

			/* re-estimation of transition probabilities NEW VERSION */
			for (int j = 0; j < Main.tagSize; j++) {
				double transOriginalNum = 0.0;
				double transNeighborNum = 0.0;

				double transOriginalDenom = 0.0;
				double transNeighborDenom = 0.0;

				double orgScore = 0.0;
				double negScore = 0.0;

				ArrayList<Double> values = new ArrayList<>();
				ArrayList<Double> valuesDenom = new ArrayList<>();

				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					transOriginalNum = 0.0;
					transNeighborNum = 0.0;

					transOriginalDenom = 0.0;
					transNeighborDenom = 0.0;
					for (int a = 0; a < globalMap.size(); a++) {

						HashMap<Integer, Node> lattice = globalMap.get(a);

						for (int k = 1; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
								values.add(p(i, j, n, n2));  //calculate transition i-j
								for (int k2 = 0; k2 < Main.tagSize; k2++) {
									valuesDenom.add(p(i, k2, n, n2));  //calculate transition i-k
								}
							}

						}

						if (a == 0) {

							transOriginalNum += Math.exp(Main.logSumOfExponentials(values));
							transOriginalDenom += Math.exp(Main.logSumOfExponentials(valuesDenom));

							values.clear();
							valuesDenom.clear();
						} else {

							transNeighborNum += Math.exp(Main.logSumOfExponentials(values));
							transNeighborDenom += Math.exp(Main.logSumOfExponentials(valuesDenom));

							values.clear();
							valuesDenom.clear();

						}
					}

					orgScore += divide(transOriginalNum, transOriginalDenom);  //divide for all sentences transOriginalNum / transOriginalDenom
					negScore += divide(transNeighborNum, transNeighborDenom);  //divide for all negative sentences transNeighborNum / transNeighborDenom

				}

				String key = (Main.tagList.get(i) + "|" + Main.tagList.get(j));
				double score = orgScore - negScore;
				gradTransitionProbabilities.put(key, score);

			}

			/* re-estimation of initial probabilities NEW VERSION */

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

			double initorgScore = 0.0;
			double initnegScore = 0.0;

			double initorgScoreEnd = 0.0;
			double initnegScoreEnd = 0.0;

			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
				initOriginalNum = 0.0;
				initNeighborNum = 0.0;

				initOriginalDenom = 0.0;
				initNeighborDenom = 0.0;

				initOriginalEndNum = 0.0;
				initNeighborEndNum = 0.0;

				initOriginalEndDenom = 0.0;
				initNeighborEndDenom = 0.0;

				for (int a = 0; a < globalMap.size(); a++) {

					HashMap<Integer, Node> lattice = globalMap.get(a);
					ArrayList<Integer> startStates = Main.returnStartStates(lattice);

					for (int k = 1; k < lattice.size() - 1; k++) {
						Node n = lattice.get(k);
						for (int k1 : n.next) {
							Node n2 = lattice.get(k1);

							if (n2.isEndState) {  // if node end state calculate initial i-</s>

								valuesInitEnd.add(Main.initialProbabilities.get(Main.tagList.get(i) + "|</s>"));  // calculate initial i-</s>
								for (int k2 = 0; k2 < Main.tagSize; k2++) {
									valuesInitEndDenom
											.add(Main.initialProbabilities.get(Main.tagList.get(k2) + "|</s>"));  // calculate initial all tags-</s>
								}

							}

						}

						if (startStates.contains(n.stateNum)) {   // if node start state calculate initial <s>-i

							valuesInit.add(Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + n.word)
									+ Main.initialProbabilities.get("<s>|" + Main.tagList.get(i)));  // calculate initial <s>-i
							for (int k2 = 0; k2 < Main.tagSize; k2++) {
								valuesInitDenom.add(Main.emissionProbabilities.get(Main.tagList.get(k2) + "|" + n.word)
										+ Main.initialProbabilities.get("<s>|" + Main.tagList.get(k2)));  // calculate initial <s>- all tags
							}
						}

					}

					if (a == 0) {

						initOriginalNum += Math.exp(Main.logSumOfExponentials(valuesInit));
						initOriginalDenom += Math.exp(Main.logSumOfExponentials(valuesInitDenom));

						initOriginalEndNum += Math.exp(Main.logSumOfExponentials(valuesInitEnd));
						initOriginalEndDenom += Math.exp(Main.logSumOfExponentials(valuesInitEndDenom));

						valuesInit.clear();
						valuesInitDenom.clear();

						valuesInitEnd.clear();
						valuesInitEndDenom.clear();

					} else {

						initNeighborNum += Math.exp(Main.logSumOfExponentials(valuesInit));
						initNeighborDenom += Math.exp(Main.logSumOfExponentials(valuesInitDenom));

						initNeighborEndNum += Math.exp(Main.logSumOfExponentials(valuesInitEnd));
						initNeighborEndDenom += Math.exp(Main.logSumOfExponentials(valuesInitEndDenom));

						valuesInit.clear();
						valuesInitDenom.clear();

						valuesInitEnd.clear();
						valuesInitEndDenom.clear();

					}
				}

				initorgScore += divide(initOriginalNum, initOriginalDenom);
				initnegScore += divide(initNeighborNum, initNeighborDenom);

				initorgScoreEnd += divide(initOriginalEndNum, initOriginalEndDenom);
				initnegScoreEnd += divide(initNeighborEndNum, initNeighborEndDenom);
			}

			String initKey = "<s>|" + Main.tagList.get(i);
			double initScore = initorgScore - initnegScore;
			gradInitialProbabilities.put(initKey, initScore);

			String initKeyEnd = Main.tagList.get(i) + "|</s>";
			double initScoreEnd = initorgScoreEnd - initnegScoreEnd;
			gradInitialProbabilities.put(initKeyEnd, initScoreEnd);

			/* re-estimation of emission probabilities NEW VERSION */
			HashMap<String, Double> originalMap = new HashMap<>();
			HashMap<String, Double> negativeMap = new HashMap<>();

			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
				HashMap<String, Double> originalwordsNum = new HashMap<String, Double>();
				HashMap<String, Double> neighborsWord = new HashMap<String, Double>();

				ArrayList<Double> listNum = new ArrayList<>();
				ArrayList<Double> listDenom = new ArrayList<>();

				double emissionOriginalDenom = 0.0;
				double emissionNeighborDenom = 0.0;

				for (int a = 0; a < globalMap.size(); a++) {
					HashMap<Integer, Node> lattice = globalMap.get(a);
					for (int k = 1; k < lattice.size(); k++) {
						Node node = lattice.get(k);

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
							listNum.add(g);

						} else {
							if (neighborsWord.containsKey(node.word)) {
								ArrayList<Double> list = new ArrayList<>();

								list.add(neighborsWord.get(node.word));
								list.add(g);

								neighborsWord.put(node.word, Main.logSumOfExponentials(list));

							} else {
								neighborsWord.put(node.word, g);
							}

							listDenom.add(g);
						}
					}
				}
				emissionOriginalDenom += Main.logSumOfExponentials(listNum);
				emissionNeighborDenom += Main.logSumOfExponentials(listDenom);
				listNum.clear();
				listDenom.clear();
				for (String s : originalwordsNum.keySet()) {
					double valOrg = originalwordsNum.get(s);

					if (originalMap.containsKey(s)) {

						ArrayList<Double> list = new ArrayList<>();

						list.add(originalMap.get(s));
						list.add(divide(valOrg, emissionOriginalDenom));

						originalMap.put(s, Main.logSumOfExponentials(list));

					} else {
						originalMap.put(s, divide(valOrg, emissionOriginalDenom));
					}

					double valNeg = neighborsWord.get(s);

					if (negativeMap.containsKey(s)) {

						ArrayList<Double> list = new ArrayList<>();

						list.add(negativeMap.get(s));
						list.add(divide(valNeg, emissionNeighborDenom));

						negativeMap.put(s, Main.logSumOfExponentials(list));

					} else {
						negativeMap.put(s, divide(valNeg, emissionNeighborDenom));
					}
				}
				originalwordsNum.clear();
				neighborsWord.clear();
			}

			Iterator<Entry<String, Double>> originals = originalMap.entrySet().iterator();
			Iterator<Entry<String, Double>> neighbors = negativeMap.entrySet().iterator();
			while (originals.hasNext()) {
				Entry<String, Double> pairs = originals.next();
				Double firstVal = (Double) pairs.getValue();
				Entry<String, Double> pairs2 = neighbors.next();
				Double secondVal = (Double) pairs2.getValue();

				String key = Main.tagList.get(i) + "|" + pairs.getKey();

				gradEmissionProbabilities.put(key, Math.exp(firstVal) - Math.exp(secondVal));

			}
		}

		grad = Main.createGradArray(gradTransitionProbabilities, gradInitialProbabilities, gradEmissionProbabilities,
				Main.gradFeature2Index);

		for (int i = 0; i < grad.length; i++) {
			grad[i] *= -1;
			// grad[i] -= iterWeights[i];
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

		// num = Main.transitionProbabilities.get(Main.tagList.get(i) + "|" +
		// Main.tagList.get(j)) *
		// Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + next.word)
		// *
		// Math.exp(node.alpha.get(i) + next.beta.get(j)) ;

		num = (Main.transitionProbabilities.get(Main.tagList.get(i) + "|" + Main.tagList.get(j))
				+ Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + next.word) + node.alpha.get(i)
				+ next.beta.get(j));

		double denom = sumListValues(node.tagScores);

		return divide(num, denom);
	}

	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num = 0.0;
		num = (node.alpha.get(i) + node.beta.get(i));
		double denom = sumListValues(node.tagScores);

		return divide(num, denom);
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

		// double sum = 0.0;
		//
		// for (double d : list) {
		// sum += (d);
		// }

		return Math.exp(Main.logSumOfExponentials((ArrayList<Double>) list));
	}

}
