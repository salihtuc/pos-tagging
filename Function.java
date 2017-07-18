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

	@Override
	public FunctionValues getValues(double[] point) {

		System.out.println("Points:");
		Main.printDoubleArray(point);
		
		return new FunctionValues(functionValue(), gradient(point));
	}
	
	public double functionValue() {

		double sumNum = 0.0;
		double sumDenom = 0.0;
		
		int counter = 0;
        for(HashMap<Integer, Node> map : Main.globalMap.values()) {
        		for(Node endState : Main.returnEndStates(map)) {
        			if(counter == 0) {
        				sumNum += addListValues(endState.tagScores);
        			}
        			else {
        				sumDenom += addListValues(endState.tagScores);
        			}
        		}
        }
        
        return functionValue + ((sumNum - sumDenom));
    }

	public double[] gradient(double[] iterWeights) {

		double[] grad = new double[iterWeights.length];
		HashMap<String, Double> gradtransitionProbabilities = new HashMap<>();
		HashMap<String, Double> grademissionProbabilities = new HashMap<>();

		/* re-estimation of transition probabilities */
		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.tagSize; j++) {
				double num = 0;
				double denom = 0;
				for (int a = 0; a < Main.globalMap.size(); a++) {
					for (int k = 0; k < Main.globalMap.get(a).size(); k++) {
						num += p(i, j, Main.globalMap.get(a).get(k));
						denom += gamma(i, Main.globalMap.get(a).get(k));

					}
				}
				gradtransitionProbabilities.put((Main.tagList.get(i) + "-" + Main.tagList.get(j)), divide(num, denom));
			}
		}

		/* re-estimation of emission probabilities */
		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.words.size(); j++) {
				double num = 0;
				double denom = 0;
				for (int a = 0; a < Main.globalMap.size(); a++) {
					for (int k = 0; k < Main.globalMap.get(a).size(); k++) {
						double g = gamma(i, Main.globalMap.get(a).get(k));
						num += g;
						denom += g;
					}
				}
				grademissionProbabilities.put(Main.tagList.get(i) + "-" + Main.words.get(j), divide(num, denom)); // NUM/DENOM
																													// is
																													// always
																													// 1.0
																													// !!
			}
		}
		grad = Main.createWeightsArray(gradtransitionProbabilities, grademissionProbabilities);
		return grad;

	}

	/**
	 * @param i
	 *            the number of state s_i
	 * @param j
	 *            the number of state s_j
	 * @return P
	 */
	public double p(int i, int j, Node node) {
		double num, denom = 0.0;
		num = node.alpha.get(i) * Main.transitionProbabilities.get(Main.tagList.get(i) + "-" + Main.tagList.get(j))
				* Main.emissionProbabilities.get(Main.tagList.get(i) + "-" + node.word) * node.beta.get(j);

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
		
		for(int i = 0; i < weights.length; i++) {
			weights[i] += (gradients[i] - lambdaValue);
		}
		
		return weights;
	}
	
	public double addListValues(List<Double> list) {
		
		double sum = 0.0;
		
		for(double d : list) {
			sum += d;
		}
		
		return sum;
	}

}
