
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

	double allSumNum = 0;
	double allSumDenom = 0;

	public double functionValue(double[] iterWeights) {
		originalList.clear();
		negativeList.clear();
		functionValue = 0.0;
		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) { // Iterate for all sentences
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> lattice : globalMap.values()) { // Each sentence's lattices
				for (Node endState : Main.returnEndStates(lattice)) {
					if (counter == 0) {
						sumNum += sumListValues(endState.tagScores); // Original lattice's score
						allSumNum += sumNum;
						counter++;
					} else {
						sumDenom += sumListValues(endState.tagScores); // Negative samples' scores
						allSumDenom += sumDenom;
					}
				}
			}
			originalList.add(sumNum);
			negativeList.add(sumDenom);
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

		double totalTransition = 0;
		double totalEmission = 0;
		double[] grad = new double[iterWeights.length];

		HashMap<String, Double> gradTransitionProbabilities = new HashMap<>();
		HashMap<String, Double> gradInitialProbabilities = new HashMap<>();
		HashMap<String, Double> gradEmissionProbabilities = new HashMap<>();

		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.tagSize; j++) {
				double transOriginalNum = 0.0;
				double transNeighborNum = 0.0;
				double orgScore = 0.0;
				double negScore = 0.0;
				
				double initialOriginalNum = 0.0;
				double initialNeighborNum = 0.0;
				double initialOrgScore = 0.0;
				double initialNegScore = 0.0;
				
				int sentenceCounter = 0;
				
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					transOriginalNum = 0.0;
					transNeighborNum = 0.0;
					for (int a = 0; a < globalMap.size(); a++) {
						double num = 0;
						double initNum = 0.0;
						HashMap<Integer, Node> lattice = globalMap.get(a);
						ArrayList<Integer> startStates = Main.returnStartStates(lattice);

						for (int k = 1; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
								num += p(i, j, n, n2);
							}
							
							if(startStates.contains(n.stateNum)) {
								
								ArrayList<Double> list = new ArrayList<>();
								
								list.add(Main.emissionProbabilities.get("t" + (i+1) + "|" + n.word));
								list.add(Main.initialProbabilities.get("<s>|t" + (i+1)));
								
								initNum += Math.exp(Main.logSumOfExponentials(list));
								
							}
						}

						if (a == 0) {
							transOriginalNum += num;
							initialOriginalNum += initNum;
						} else {
							transNeighborNum += num;
							initialNeighborNum += initNum;
						}
					}

					orgScore += divide(transOriginalNum, originalList.get(sentenceCounter));
					negScore += divide(transNeighborNum, negativeList.get(sentenceCounter));
					
					initialOrgScore += divide(initialOriginalNum, originalList.get(sentenceCounter));
					initialNegScore += divide(initialNeighborNum, negativeList.get(sentenceCounter));
					
					sentenceCounter++;
				}

				String key = (Main.tagList.get(i) + "|" + Main.tagList.get(j));
				double score = orgScore - negScore;
				gradTransitionProbabilities.put(key, score);
				
				String initKey = "<s>|" + Main.tagList.get(i);
				double initScore = initialOrgScore - initialNegScore;
				gradInitialProbabilities.put(initKey, initScore);
				// }

				/* re-estimation of emission probabilities NEW VERSION */
				HashMap<String, Double> originalMap = new HashMap<>();
				HashMap<String, Double> negativeMap = new HashMap<>();
				double emissionOriginalDenom = 0.0;
				double emissionNeighborDenom = 0.0;
				sentenceCounter = 0;
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					HashMap<String, Double> originalwordsNum = new HashMap<String, Double>();
					HashMap<String, Double> neighborsWord = new HashMap<String, Double>();
					for (int a = 0; a < globalMap.size(); a++) {
						HashMap<Integer, Node> lattice = globalMap.get(a);
						for (int k = 1; k < lattice.size() - 1; k++) {
							Node node = lattice.get(k);
							double num = 0;
							for (int k1 : node.next) {
								Node n2 = lattice.get(k1);
								num += p(i, j, node, n2);
							}

							if (a == 0) {
								if (originalwordsNum.containsKey(node.word)) {
									originalwordsNum.put(node.word, originalwordsNum.get(node.word) + num);

								} else {
									originalwordsNum.put(node.word, num);
								}
								// emissionOriginalDenom += g;

							} else {
								if (neighborsWord.containsKey(node.word)) {
									neighborsWord.put(node.word, neighborsWord.get(node.word) + num);

								} else {
									neighborsWord.put(node.word, num);
								}
								// emissionNeighborDenom += g;

							}

						}
					}
					for (String s : originalwordsNum.keySet()) {
						// double valOrg = divide(divide(originalwordsNum.get(s),
						// emissionOriginalDenom), originalList.get(sentenceCounter));
						double valOrg = divide(originalwordsNum.get(s), originalList.get(sentenceCounter));

						if (originalMap.containsKey(s)) {
							originalMap.put(s, originalwordsNum.get(s) + valOrg);

						} else {
							originalMap.put(s, valOrg);
						}

						// double valNeg = divide(divide(neighborsWord.get(s), emissionNeighborDenom),
						// negativeList.get(sentenceCounter));
						double valNeg = divide(neighborsWord.get(s), negativeList.get(sentenceCounter));

						if (negativeMap.containsKey(s)) {
							negativeMap.put(s, negativeMap.get(s) + valNeg);

						} else {
							negativeMap.put(s, valNeg);
						}
					}
					sentenceCounter++;
				}

				Iterator<Entry<String, Double>> originals = originalMap.entrySet().iterator();
				Iterator<Entry<String, Double>> neighbors = negativeMap.entrySet().iterator();
				while (originals.hasNext()) {
					Entry<String, Double> pairs = originals.next();
					Double firstVal = (Double) pairs.getValue();
					Entry<String, Double> pairs2 = neighbors.next();
					Double secondVal = (Double) pairs2.getValue();

					key = Main.tagList.get(i) + "|" + pairs.getKey();

					gradEmissionProbabilities.put(key, firstVal - secondVal);
					totalEmission += (firstVal - secondVal);

				}
			}
		}

		grad = Main.createGradArray(gradTransitionProbabilities, gradInitialProbabilities, gradEmissionProbabilities, Main.gradFeature2Index);

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

		ArrayList<Double> list = new ArrayList<>();

		list.add(Main.transitionProbabilities.get(Main.tagList.get(i) + "|" + Main.tagList.get(j)));
		list.add(Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + node.word));
		list.add(Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + next.word));

		num = Math.exp(Main.logSumOfExponentials(list));

		return num;
	}

	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num, denom = 0.0;
		num = (node.alpha.get(i)) * (node.beta.get(i));

		// for (int k = 0; k < Main.tagSize; k++) {
		// denom += (node.alpha.get(k)) * (node.beta.get(k));
		// }

		// return divide(num, denom);

		return num;
	}

	/** divides two doubles. (0 / 0 = 0!) && (1 / 0 = 0!) */
	public double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		// else if(!Double.isFinite(n) || !Double.isFinite(d))
		// return 0;
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

		double sum = 0.0;

		for (double d : list) {
			sum += (d);
		}

		return (sum);
	}

}
