/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

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
//		Main.printDoubleArray(point);

		return new FunctionValues(functionValue(point), gradient(point));
	}

	public double functionValue(double[] iterWeights) {

		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> map : globalMap.values()) {
				for (Node endState : Main.returnEndStates(map)) {
					if (counter == 0) {
						sumNum += addListValues(endState.tagScores);
					} else {
						sumDenom += addListValues(endState.tagScores);
					}
				}
			}

			functionValue += (sumNum - sumDenom);
		}

		for (double weight : iterWeights) {
			functionValue -= LAMBDA_EM * Math.pow(weight, 2);
		}
		functionValue *= -1;

		return functionValue;
	}

	public double[] gradient(double[] iterWeights) {

		double[] grad = new double[iterWeights.length];
		HashMap<String, Double> gradtransitionProbabilities = new HashMap<>();
		HashMap<String, Double> grademissionProbabilities = new HashMap<>();

		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.tagSize; j++) {
				double num = 0;
				double denom = 0;
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					for (int a = 0; a < globalMap.size(); a++) {
						for (int k = 0; k < globalMap.get(a).size() - 1; k++) {
							for (int k1 : globalMap.get(a).get(k).next) {
								Node n2 = globalMap.get(a).get(k1);
								for (int k2 = 0; k2 < Main.tagSize; k2++) {
									denom += p(i, k2, globalMap.get(a).get(k), n2);
								}
								num += p(i, j, globalMap.get(a).get(k), n2);
							}

						}
					}
				}
				gradtransitionProbabilities.put((Main.tagList.get(i) + "-" + Main.tagList.get(j)), divide(num, denom));
			}

		}

		/* re-estimation of emission probabilities */
		
		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.allWords.size(); j++) {
				double num = 0;
				double denom = 0;
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					for (int a = 0; a < globalMap.size(); a++) {
						//TODO ArrayList<String> uniqueWords = Main.returnUniqueWords(globalMap.get(a));
						for (int k = 0; k < globalMap.get(a).size(); k++) {

							double g = gamma(i, globalMap.get(a).get(k));
							num += g * (Main.allWords.get(j).equals(globalMap.get(a).get(k).word) ? 1 : 0);
							denom += g;

						}
					}
				}
				grademissionProbabilities.put(Main.tagList.get(i) + "-" + Main.allWords.get(j), divide(num, denom));

			}
		}

		grad = Main.createWeightsArray(gradtransitionProbabilities, grademissionProbabilities);

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
