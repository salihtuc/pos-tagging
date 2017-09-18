
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
		//originalList.clear();
		//negativeList.clear();
		functionValue = 0.0;
		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) { // Iterate for all sentences
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> lattice : globalMap.values()) { // Each sentence's lattices
				for (Node endState : Main.returnEndStates(lattice)) {
					if (counter == 0) {
						sumNum += sumListValues(endState.tagScores); // Original lattice's score
						counter++;
					} else {
						sumDenom += sumListValues(endState.tagScores); // Negative samples' scores
					}
				}
			}
			//originalList.add(sumNum);
			//negativeList.add(sumDenom);
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
				double transOriginalNum = 0.0;
				double transNeighborNum = 0.0;
				double orgScore = 0.0;
				double negScore = 0.0;
				
				double initialOriginalNum = 0.0;
				double initialNeighborNum = 0.0;
				double initialOrgScore = 0.0;
				double initialNegScore = 0.0;
				
				double initialOriginalNumEnd = 0.0;
				double initialNeighborNumEnd = 0.0;
				double initialOrgScoreEnd = 0.0;
				double initialNegScoreEnd = 0.0;
				
				//int sentenceCounter = 0;
				
				ArrayList<Double> values = new ArrayList<>();
//				ArrayList<Double> valuesNeg = new ArrayList<>();
				
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					transOriginalNum = 0.0;
					transNeighborNum = 0.0;
					for (int a = 0; a < globalMap.size(); a++) {
						double num = 0;
						double initNum = 0.0;
						double initNumEnd = 0.0;
						
						HashMap<Integer, Node> lattice = globalMap.get(a);
						ArrayList<Integer> startStates = Main.returnStartStates(lattice);

						for (int k = 1; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
//								num += p(i, j, n, n2);
								values.add(p(i, j, n, n2));
								
								
								if(n2.isEndState) {
//									ArrayList<Double> list = new ArrayList<>();
//									
//									list.add(Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + n2.word));
//									list.add(Main.initialProbabilities.get(Main.tagList.get(j) + "|</s>"));
									
									initNumEnd += Math.exp(Main.initialProbabilities.get(Main.tagList.get(j) + "|</s>"));
								}

							}
							
							if(startStates.contains(n.stateNum)) {
								
//								ArrayList<Double> list = new ArrayList<>();
//								
//								list.add(Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + n.word));
//								list.add(Main.initialProbabilities.get("<s>|" + Main.tagList.get(i)));
								
								initNum += Math.exp(Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + n.word) + Main.initialProbabilities.get("<s>|" + Main.tagList.get(i)));
								
							}
						}

						if (a == 0) {
//							transOriginalNum += Math.log(num);
//							initialOriginalNum += Math.log(initNum);
//							initialOriginalNumEnd += Math.log(initNumEnd);
							
							transOriginalNum += Math.exp(Main.logSumOfExponentials(values));
							initialOriginalNum += (initNum);
							initialOriginalNumEnd += (initNumEnd);
							
							values.clear();
						} else {
							
							transNeighborNum += Math.exp(Main.logSumOfExponentials(values));
							initialNeighborNum += (initNum);
							initialNeighborNumEnd += (initNumEnd);
							
							values.clear();
							
//							transNeighborNum += Math.log(num);
//							initialNeighborNum += Math.log(initNum);
//							initialNeighborNumEnd += Math.log(initNumEnd);
						}
					}

					orgScore += transOriginalNum;
					negScore += transNeighborNum;
					
					initialOrgScore += initialOriginalNum;
					initialNegScore += initialNeighborNum;
					
					initialOrgScoreEnd += initialOriginalNumEnd;
					initialNegScoreEnd += initialNeighborNumEnd;
					
					//sentenceCounter++;
				}

				String key = (Main.tagList.get(i) + "|" + Main.tagList.get(j));
				double score = orgScore - negScore;
				gradTransitionProbabilities.put(key, score);
				
				String initKey = "<s>|" + Main.tagList.get(i);
				double initScore = Math.log(initialOrgScore / initialNegScore);
				gradInitialProbabilities.put(initKey, initScore);
				
				String initKeyEnd = Main.tagList.get(j) + "|</s>";
				double initScoreEnd = Math.log(initialOrgScoreEnd / initialNegScoreEnd);
				gradInitialProbabilities.put(initKeyEnd, initScoreEnd);
				
//				String key = (Main.tagList.get(i) + "|" + Main.tagList.get(j));
//				double score = orgScore - negScore;
//				gradTransitionProbabilities.put(key, score);
//				
//				String initKey = "<s>|" + Main.tagList.get(i);
//				double initScore = initialOrgScore - initialNegScore;
//				gradInitialProbabilities.put(initKey, initScore);
//				
//				String initKeyEnd = Main.tagList.get(j) + "|</s>";
//				double initScoreEnd = initialOrgScoreEnd - initialNegScoreEnd;
//				gradInitialProbabilities.put(initKeyEnd, initScoreEnd);
				
			}

				/* re-estimation of emission probabilities NEW VERSION */
				HashMap<String, Double> originalMap = new HashMap<>();
				HashMap<String, Double> negativeMap = new HashMap<>();
				
				//sentenceCounter = 0;
				ArrayList<Double> values = new ArrayList<>();
				
				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : Main.sentenceGlobals) {
					HashMap<String, Double> originalwordsNum = new HashMap<String, Double>();
					HashMap<String, Double> neighborsWord = new HashMap<String, Double>();
					for (int a = 0; a < globalMap.size(); a++) {
						HashMap<Integer, Node> lattice = globalMap.get(a);
						for (int k = 1; k < lattice.size(); k++) {
							Node node = lattice.get(k);
							
								double g = gamma(i,  node);
							

							if (a == 0) {
								if (originalwordsNum.containsKey(node.word)) {
									ArrayList<Double> list = new ArrayList<>();
									
									list.add(originalwordsNum.get(node.word));
									list.add(g);
									
									originalwordsNum.put(node.word, Main.logSumOfExponentials(list));

								} else {
									originalwordsNum.put(node.word, g);
								}
								// emissionOriginalDenom += g;

							} else {
								if (neighborsWord.containsKey(node.word)) {
									ArrayList<Double> list = new ArrayList<>();
									
									list.add(neighborsWord.get(node.word));
									list.add(g);
									
									neighborsWord.put(node.word, Main.logSumOfExponentials(list));


								} else {
									neighborsWord.put(node.word, g);
								}
								// emissionNeighborDenom += g;

							}

						}
					}
					for (String s : originalwordsNum.keySet()) {
						// double valOrg = divide(divide(originalwordsNum.get(s),
						// emissionOriginalDenom), originalList.get(sentenceCounter));
						double valOrg = originalwordsNum.get(s);

						if (originalMap.containsKey(s)) {
//							originalMap.put(s, originalMap.get(s) + valOrg);
							
							ArrayList<Double> list = new ArrayList<>();
							
							list.add(originalMap.get(s));
							list.add(valOrg);
							
							originalMap.put(s, Main.logSumOfExponentials(list));

						} else {
							originalMap.put(s, valOrg);
						}

						// double valNeg = divide(divide(neighborsWord.get(s), emissionNeighborDenom),
						// negativeList.get(sentenceCounter));
						double valNeg = neighborsWord.get(s);

						if (negativeMap.containsKey(s)) {
//							negativeMap.put(s, negativeMap.get(s) + valNeg);
							
							ArrayList<Double> list = new ArrayList<>();
							
							list.add(negativeMap.get(s));
							list.add(valNeg);
							
							negativeMap.put(s, Main.logSumOfExponentials(list));

						} else {
							negativeMap.put(s, valNeg);
						}
					}
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
		
//		num = Main.transitionProbabilities.get(Main.tagList.get(i) + "|" + Main.tagList.get(j)) * Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + next.word) * 
//				Math.exp(node.alpha.get(i) + next.beta.get(j)) ;
		
		num = (Main.transitionProbabilities.get(Main.tagList.get(i) + "|" + Main.tagList.get(j)) + Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + next.word) + 
		node.alpha.get(i)+next.beta.get(j));	

		double denom = sumListValues(node.tagScores);
		
//		list.add(Main.transitionProbabilities.get(Main.tagList.get(i) + "|" + Main.tagList.get(j)));
//		list.add(Main.emissionProbabilities.get(Main.tagList.get(i) + "|" + node.word));
//		list.add(Main.emissionProbabilities.get(Main.tagList.get(j) + "|" + next.word));

//		num = Math.exp(num);

		return divide(num, denom);
	}

	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num = 0.0;
//		num = (node.alpha.get(i)) + (node.beta.get(i));
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

//		double sum = 0.0;
//
//		for (double d : list) {
//			sum += (d);
//		}

		return Math.exp(Main.logSumOfExponentials((ArrayList<Double>)list));
	}

}
