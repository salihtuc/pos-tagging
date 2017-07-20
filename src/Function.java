/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

public class Function implements DifferentiableFunction {

	// -------------------------------------- LBFGS-B
	// ---------------------------------------------------------

	public double functionValue = 0.0;
	public static final double LAMBDA_EM = 1;

	@Override
	public FunctionValues getValues(double[] point) {

//		System.out.println("Points:");
		//Main.printDoubleArray(point);

		return new FunctionValues(functionValue(point), gradient(point));
	}

	public double functionValue(double[] iterWeights) {
		
		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {	// Iterate for all sentences
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> map : globalMap.values()) {	// Each sentence's lattices
				for (Node endState : Main.returnEndStates(map)) {	// Nodes of each lattice
					if (counter == 0) {
						sumNum += addListValues(endState.alpha);	// Original lattice's score
					} else {
						sumDenom += addListValues(endState.alpha);	// Negative samples' scores
					}
				}
			}
			functionValue += (sumNum - sumDenom);	// Each sentence's scores added.
		}

		// Regularization
		for (double weight : iterWeights) {
			functionValue -= LAMBDA_EM * Math.pow(weight, 2);
		}
		functionValue *= -1;	// -1 because of minimizing
		
		return functionValue;
	}

	public double[] gradient(double[] iterWeights) {

				double[] grad = new double[iterWeights.length];
		HashMap<String, Double> gradtransitionProbabilities = new HashMap<>();
		HashMap<String, Double> grademissionProbabilities = new HashMap<>();

			for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.tagSize; j++) {
				double transOriginal = 0.0;
				double transNeighbor = 0.0;
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {

					for (int a = 0; a < globalMap.size(); a++) {
						double num = 0;
						double denom = 0;

						HashMap<Integer, Node> lattice = globalMap.get(a);

						for (int k = 0; k < lattice.size() - 1; k++) {
							for (int k1 : lattice.get(k).next) {
								Node n2 = lattice.get(k1);
								for (int k2 = 0; k2 < Main.tagSize; k2++) {
									denom += p(i, k2, lattice.get(k), n2);
								}
								num += p(i, j, lattice.get(k), n2);
							}

						}
						if (a == 0) {
							transOriginal += divide(num, denom);
						} else {
							transNeighbor += divide(num, denom);
						}
					}
				}
				gradtransitionProbabilities.put((Main.tagList.get(i) + "-" + Main.tagList.get(j)),transOriginal-transNeighbor);
			}

		}

		// /* re-estimation of emission probabilities */
		//
		// for (int i = 0; i < Main.tagSize; i++) {
		// double num = 0;
		// double denom = 0;
		//
		// for (HashMap<Integer, HashMap<Integer, Node>> globalMap :
		// Main.sentenceGlobals) {
		// for (int a = 0; a < globalMap.size(); a++) {
		// ArrayList<String> uniqueWords =
		// Main.returnUniqueWords(globalMap.get(a));
		// for (int k = 0; k < globalMap.get(a).size(); k++) {
		// for (int j = 0; j < uniqueWords.size(); j++) {
		//
		// double g = gamma(i, globalMap.get(a).get(k));
		// num += g * (uniqueWords.get(j).equals(globalMap.get(a).get(k).word) ?
		// 1 : 0);
		// denom += g;
		// grademissionProbabilities.put(Main.tagList.get(i) + "-" +
		// uniqueWords.get(j),
		// divide(num, denom));
		// }
		// }
		// }
		//
		// }
		// }

		/* re-estimation of emission probabilities */

		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.allWords.size(); j++) {
				double emissionOriginal = 0.0;
				double emissionNeighbor = 0.0;

				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					for (int a = 0; a < globalMap.size(); a++) {
						double num = 0;
						double denom = 0;
						// TODO ArrayList<String> uniqueWords =
						// Main.returnUniqueWords(globalMap.get(a));
						for (int k = 0; k < globalMap.get(a).size(); k++) {

							double g = gamma(i, globalMap.get(a).get(k));
							num += g * (Main.allWords.get(j).equals(globalMap.get(a).get(k).word) ? 1 : 0);
							denom += g;

						}
						if (a == 0) {
							emissionOriginal += divide(num, denom);
						} else {
							emissionNeighbor += divide(num, denom);
						}
					}
				}
				grademissionProbabilities.put(Main.tagList.get(i) + "-" + Main.allWords.get(j), emissionOriginal-emissionNeighbor);

			}
		}

		grad = Main.createWeightsArray(gradtransitionProbabilities, grademissionProbabilities);

		for (int i = 0; i < grad.length; i++)
			grad[i] *= -1;
			
		for (int i = 0; i < grad.length; i++)
			grad[i] += 2 * LAMBDA_EM * iterWeights[i];

		return grad;

	}

	/**
	 * @param i
	 *            the number of state s_i
	 * @param j
	 *            the number of state s_j
	 * @return P
	 */
	public double p(int i, int j, Node node,Node next) {
		double num, denom = 0.0;
		if(node.isEndState)
			num=node.alpha.get(i)*Main.transitionProbabilities.get(Main.tagList.get(i) + "-" + Main.tagList.get(j));
		else
			num = node.alpha.get(i) * Main.transitionProbabilities.get(Main.tagList.get(i) + "-" + Main.tagList.get(j))
				* Main.emissionProbabilities.get(Main.tagList.get(j) + "-" + next.word) * next.beta.get(j);

		for (int k = 0; k < Main.tagSize; k++) {
			denom += node.alpha.get(k) * node.beta.get(k);
		}

		return divide(num, denom);
	}
	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num, denom = 0.0;
		num = node.alpha.get(i) * node.beta.get(i);

		for (int k = 0; k < Main.tagSize; k++) {
			denom += node.alpha.get(k) * node.beta.get(k);
		}

		return divide(num, denom);
	}

	/** divides two doubles. 0 / 0 = 0! */
	public double divide(double n, double d) {
		if (n == 0)
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

	public double addListValues(List<Double> list) {

		double sum = 0.0;

		for (double d : list) {
			sum += d;
		}

		return sum;
	}

}
