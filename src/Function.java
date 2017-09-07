
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

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

public class Function implements DifferentiableFunction {

	// -------------------------------------- LBFGS-B
	// ---------------------------------------------------------

	public double functionValue = 0.0; //TODO
	public static final double LAMBDA_EM = 1;
	public int cc = 0;
	ArrayList<Double> originalList = new ArrayList<>();
	ArrayList<Double> negativeList = new ArrayList<>();

	@Override
	public FunctionValues getValues(double[] point) {

		// System.out.println("Points:");
		// Main.printDoubleArray(point);

		FunctionValues fv = new FunctionValues(functionValue(point), gradient(point));
		
		// Update the feature probabilities!!
		//Main.updateProbabilities(point);
		
		for(HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
			for(int j = 0; j < globalMap.size(); j++) {
				HashMap<Integer, Node> targetMap = globalMap.get(j);
//				System.out.println("Lattice: " + j + "\n" + targetMap);
//				System.out.println("----------------------");
				
				// Iterate over targetMap
				Main.iterateFromStartToEnd(targetMap);
				Main.iterateFromEndToStart(targetMap);
			}
		}
		
		return fv;
	}

	double allSumNum = 0;
	double allSumDenom = 0;
	
	public double functionValue(double[] iterWeights) {
		functionValue = 0.0;
		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) { // Iterate for all sentences
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> lattice : globalMap.values()) { // Each sentence's lattices
//				Main.pw.println(lattice);
				for (Node endState : Main.returnEndStates(lattice)) { // Nodes of each lattice
					if (counter == 0) {
						sumNum += sumListValues(endState.alpha); // Original lattice's score	//exp or not?!!!
						allSumNum += sumNum;
						counter++;
					} else {
						sumDenom += sumListValues(endState.alpha); // Negative samples' scores
						allSumDenom += sumDenom;
					}
				}
			}
			originalList.add(sumNum);
			negativeList.add(sumDenom);
			double div =  divide(sumNum, sumDenom);
			if(Double.isFinite(div) && div != 0)
				functionValue += Math.log(div); // Each sentence's scores added.
		}
		//functionValue = Math.log(functionValue);
		//functionValue = Math.exp(functionValue);
		// Regularization
		for (double weight : iterWeights) {
			functionValue -= LAMBDA_EM * Math.pow(weight, 2);
		}
		functionValue *= -1; // -1 because of minimizing
		
		

		return functionValue;
	}
	public static ArrayList<Integer> passedList = new ArrayList<>();

	public double[] gradient(double[] iterWeights) {

		double totalTransition = 0;
		double totalEmission = 0;
		double[] grad = new double[iterWeights.length];
		
		HashMap<String, Double> gradtransitionProbabilities = new HashMap<>();
		HashMap<String, Double> grademissionProbabilities = new HashMap<>();

		for (int i = 0; i < Main.tagSize; i++) {
			for (int j = 0; j < Main.tagSize; j++) {
				double transOriginalNum = 0.0;
				double transOriginalDenom = 0.0;
				double transNeighborNum = 0.0;
				double transNeighborDenom = 0.0;
				double orgScore = 0.0;
				double negScore = 0.0;
				int sentenceCounter = 0;
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					transOriginalNum = 0.0;
					transOriginalDenom = 0.0;
					transNeighborNum = 0.0;
					transNeighborDenom = 0.0;
					for (int a = 0; a < globalMap.size(); a++) {
						double num = 0;
						double denom = 0;

						HashMap<Integer, Node> lattice = globalMap.get(a);

						for (int k = 0; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
//								for (int k2 = 0; k2 < Main.tagSize; k2++) {
//									denom += p(i, k2, n, n2);
//								}
								num += p(i, j, n, n2);
							}

						}
						
						if (a == 0) {
							transOriginalNum += num;
//							transOriginalDenom += denom;
						} else {
							transNeighborNum += num;
//							transNeighborDenom += denom;
						}
					}
//					orgScore += divide(divide(transOriginalNum, transOriginalDenom), originalList.get(sentenceCounter));
//					negScore += divide(divide(transNeighborNum, transNeighborDenom), negativeList.get(sentenceCounter));
					
					orgScore += divide(transOriginalNum, originalList.get(sentenceCounter));
					negScore += divide(transNeighborNum, negativeList.get(sentenceCounter));
					sentenceCounter++;
				}
				cc++;
//				gradtransitionProbabilities.put((Main.tagList.get(i) + "-" + Main.tagList.get(j)),
//						orgScore - negScore);
				
				String key = (Main.tagList.get(i) + "-" + Main.tagList.get(j));
//				gradtransitionProbabilities.put(key,
//						(orgScore - negScore));
				
				if(key.endsWith("<s>")) {
					gradtransitionProbabilities.put(key, 0.0);
				}
				else {
					gradtransitionProbabilities.put(key,
							(orgScore - negScore));
					totalTransition += orgScore - negScore;
				}
				
				
				
						//divide(transOriginalNum, transOriginalDenom));
			}

			/* re-estimation of emission probabilities NEW VERSION */
			// for (int i = 0; i < Main.tagSize; i++) {
			HashMap<String, Double> originalMap = new HashMap<>();
			HashMap<String, Double> negativeMap = new HashMap<>();
			double emissionOriginalDenom = 0.0;
			double emissionNeighborDenom = 0.0;
			int sentenceCounter = 0;
			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
				HashMap<String, Double> originalwordsNum = new HashMap<String, Double>();
				HashMap<String, Double> neighborsWord = new HashMap<String, Double>();
				for (int a = 0; a < globalMap.size(); a++) {
					for (int k = 0; k < globalMap.get(a).size(); k++) {
						Node node = globalMap.get(a).get(k);

						double g = gamma(i, node);

						if (a == 0) {
							if (originalwordsNum.containsKey(node.word)) {
								originalwordsNum.put(node.word, originalwordsNum.get(node.word) + g);

							} else {
								originalwordsNum.put(node.word, g);
							}
//							emissionOriginalDenom += g;

						} else {
							if (neighborsWord.containsKey(node.word)) {
								neighborsWord.put(node.word, neighborsWord.get(node.word) + g);

							} else {
								neighborsWord.put(node.word, g);
							}
//							emissionNeighborDenom += g;

						}

					}
				}
				for(String s : originalwordsNum.keySet()) {
//					double valOrg = divide(divide(originalwordsNum.get(s), emissionOriginalDenom), originalList.get(sentenceCounter));
					double valOrg = divide(originalwordsNum.get(s), originalList.get(sentenceCounter));

					if (originalMap.containsKey(s)) {
						originalMap.put(s, originalwordsNum.get(s) + valOrg);

					} else {
						originalMap.put(s, valOrg);
					}
					
//					double valNeg = divide(divide(neighborsWord.get(s), emissionNeighborDenom), negativeList.get(sentenceCounter));
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
				
				String key = Main.tagList.get(i) + "-" + pairs.getKey();
				
//				grademissionProbabilities.put(key,
//						firstVal - secondVal);
				
				if(key.startsWith("<s>") && key.endsWith("<s>")) {
					grademissionProbabilities.put(key,
							1.0);
					
					totalEmission += 1.0;
				}
				else if(!key.startsWith("<s>") && !key.endsWith("<s>")) {
					grademissionProbabilities.put(key,
							firstVal - secondVal);
					
					totalEmission += firstVal - secondVal;
				}
				else {
					grademissionProbabilities.put(key, 0.0);
				}
				
						//divide(firstVal, emissionOriginalDenom));

			}
		}
		
		passedList.clear();
		for(String s : gradtransitionProbabilities.keySet()) {
			double d = gradtransitionProbabilities.get(s);
			gradtransitionProbabilities.put(s, divide(d, totalTransition));
		}
		for(String s : grademissionProbabilities.keySet()) {
			double d = grademissionProbabilities.get(s);
			grademissionProbabilities.put(s, divide(d, totalEmission));
		}
		grad = Main.createWeightsArray(gradtransitionProbabilities, grademissionProbabilities);
		
//		System.out.println("The list: " + Main.featureList);
//		System.out.println("The Map: " + gradtransitionProbabilities.keySet());

		for (int i = 0; i < grad.length; i++) {
			grad[i] *= -1;
//			grad[i] -= iterWeights[i];
		}

		
//		Main.pw.println("Transition-Emission Old::::");
//		for(String s : Main.transitionProbabilities.keySet()) {
//			System.out.println(s + " " + Main.transitionProbabilities.get(s)); 
//			Main.pw.println(s + " " + Main.transitionProbabilities.get(s)); 
//		}
//		
//		for(String s : Main.emissionProbabilities.keySet()) {
//			System.out.println(s + " " + Main.emissionProbabilities.get(s)); 
//			Main.pw.println(s + " " + Main.emissionProbabilities.get(s)); 
//		}
//		
//		Main.pw.println("Grads:::.");
		
		double dSum = 0.0;
		double dSumPre = 0.0;
		for (int i = 0; i < grad.length; i++) {
			
			if(!Double.isFinite(grad[i]))
				grad[i] = 0.0;
			
			dSumPre += grad[i];
			grad[i] += 2 * LAMBDA_EM * iterWeights[i];
			dSum += grad[i];
			//Main.pw.println(grad[i]);
		}
//		
//		System.out.println("Length: " + grad.length);
		Main.updateProbabilities(grad);
		Main.pw.println("***********************************************NEW ITERATION*********************************************");
		Main.pw.println("Pre Total: " + dSumPre);
		System.out.println("Pre Total: " + dSumPre);
		Main.pw.println("Post Total: " + dSum);
		System.out.println("Post Total: " + dSum);
//		for(String s : Main.transitionProbabilities.keySet()) {
//			Main.pw.println(s + " " + Main.transitionProbabilities.get(s)); 
//		}
//		
//		for(String s : Main.emissionProbabilities.keySet()) {
//			Main.pw.println(s + " " + Main.emissionProbabilities.get(s)); 
//		}
		
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
		double num, denom = 0.0;
		if (node.next.isEmpty())
			num = (node.alpha.get(i)) * Main.transitionProbabilities.get(Main.tagList.get(i) + "-" + Main.tagList.get(j));
		else
			num = (node.alpha.get(i)) * Main.transitionProbabilities.get(Main.tagList.get(i) + "-" + Main.tagList.get(j))
					* Main.emissionProbabilities.get(Main.tagList.get(j) + "-" + next.word) * (next.beta.get(j));

//		for (int k = 0; k < Main.tagSize; k++) {
//			denom += (node.alpha.get(k)) * (node.beta.get(k));
//		}

//		return divide(num, denom);
		
		return num;
	}

	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num, denom = 0.0;
		num = (node.alpha.get(i)) * (node.beta.get(i));

//		for (int k = 0; k < Main.tagSize; k++) {
//			denom += (node.alpha.get(k)) * (node.beta.get(k));
//		}

//		return divide(num, denom);
		
		return num;
	}

	/** divides two doubles. (0 / 0 = 0!) && (1 / 0 = 0!) */
	public double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
//		else if(!Double.isFinite(n) || !Double.isFinite(d))
//			return 0;
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
	
	private void updateFeatures(HashMap<String, Double> transitionMap, HashMap<String, Double> emissionMap, double[] updatedValues) {
		for(int i = 0; i < updatedValues.length; i++) {
			
		}
	}

}
