
/*
 * @author: Necva Bolucu (@necvabolucu)
 * @author: Salih Tuc (@salihtuc)
 * 
 * */

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

public class Function implements DifferentiableFunction {

	// -------------------------------------- LBFGS-B
	// ---------------------------------------------------------

	public double functionValue = 0.0;
	public static final double LAMBDA_EM = 1.0;
	
			
	public int iterCount = 0;
	ArrayList<Double> originalList = new ArrayList<>();
	ArrayList<Double> negativeList = new ArrayList<>();

	@Override
	public FunctionValues getValues(double[] point) {
		
		System.out.println("After lbfgs");
		findMinAndMax(point);

//		JointModel.normalizeFeatures(point, "transition");
//		JointModel.normalizeFeatures(point, "initial");
//		JointModel.normalizeFeatures(point, "emission");
//		
//		System.out.println("After 1st normalization");
//		findMinAndMax(point);
		
		if(iterCount == 0) {
			JointModel.updateGeneralEmissions();
		}
		else if (iterCount != 0) {

			JointModel.updateProbabilities(point, JointModel.gradFeature2Index);
			JointModel.updateGeneralEmissions();

			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : JointModel.sentenceGlobals) {
				for (int j = 0; j < globalMap.size(); j++) {
					HashMap<Integer, Node> targetMap = globalMap.get(j);

					String decideOriginal;
					
					if(j == 0)
						decideOriginal = "original";
					else
						decideOriginal = "negative";
					// Iterate over targetMap
					JointModel.iterateFromStartToEnd(targetMap, decideOriginal);
					JointModel.iterateFromEndToStart(targetMap, decideOriginal);

				}
			}
		}
		iterCount++;

		FunctionValues fv = new FunctionValues(functionValue(point), gradient(point));
		
		System.out.println("After gradient");
		findMinAndMax(fv.gradient);
		
//		JointModel.normalizeFeatures(fv.gradient, "transition");
//		JointModel.normalizeFeatures(fv.gradient, "initial");
//		JointModel.normalizeFeatures(fv.gradient, "emission");
//		
//		System.out.println("After 2nd normalization"); 
//		findMinAndMax(fv.gradient);
		
		System.out.println("Inside-Iteration: " + iterCount);

		return fv;
	}

	public double functionValue(double[] iterWeights) {
		// originalList.clear();
		// negativeList.clear();
		functionValue = 0.0;
		for (HashMap<Integer, HashMap<Integer, Node>> globalMap : JointModel.sentenceGlobals) { // Iterate
																							// for
																							// all
																							// sentences
			int counter = 0;
			double sumNum = 0.0;
			double sumDenom = 0.0;
			for (HashMap<Integer, Node> lattice : globalMap.values()) { // Each sentence's lattices
				for (Node endState : Main.returnEndStates(lattice)) {
					if (counter == 0) {
						sumNum += sumListValuesExp(endState.tagScores); // Original lattice's score
						counter++;
					} else {
						sumDenom += sumListValuesExp(endState.tagScores); // Negative samples' scores
					}
				}
			}

			double div = divide(sumNum, sumDenom);
			if (Double.isFinite(div) && div != 0 && div > 0)
				functionValue += Math.log(div); // Each sentence's scores added.

		}

		// Regularization
//		for (double weight : iterWeights) {
//			functionValue -= LAMBDA_EM * Math.pow(weight, 2);
//		}
//		functionValue *= -1; // -1 because of minimizing

		System.out.println("f = " + functionValue);
		
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
		
		HashMap<String, Double> gradCoarseProbabilities = new HashMap<>();
		HashMap<String, Double> gradCoarseProbabilitiesNegative = new HashMap<>();
		HashMap<String, Double> coarseWordProbabilities = new HashMap<>();

		for (int i = 0; i < JointModel.tagSize; i++) {
			for (int j = 0; j < JointModel.tagSize; j++) {
				
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
				String decideOriginal;

				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : JointModel.sentenceGlobals) {
					
					for (int a = 0; a < globalMap.size(); a++) {
						
						if(a == 0)
							decideOriginal = "original";
						else
							decideOriginal = "negative";
					

						HashMap<Integer, Node> lattice = globalMap.get(a);
						ArrayList<Integer> startStates = JointModel.returnStartStates(lattice);

						for (int k = 1; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
								values.add(p(i, j, n, n2, decideOriginal));  //calculate transition i-j
								
								for (int k2 = 0; k2 < JointModel.tagSize; k2++) {
									valuesDenom.add(p(i, k2, n, n2, decideOriginal));  //calculate transition i-k
								}
								
								if (n2.isEndState) {  // if node end state calculate initial i-</s>
									valuesInitEnd.add(JointModel.initialProbabilities.get(JointModel.tagList.get(i) + "|</s>"));  // calculate initial i-</s>
									for (int k2 = 0; k2 < JointModel.tagSize; k2++) {
										valuesInitEndDenom
												.add(JointModel.initialProbabilities.get(JointModel.tagList.get(k2) + "|</s>"));  // calculate initial all tags-</s>
									}

								}
							}
							if (startStates.contains(n.stateNum)) {   // if node start state calculate initial <s>-i

								double emission = 0;
								
								String key = JointModel.tagList.get(i) + "|" + n.word;
								
								if(decideOriginal.equals("original") && JointModel.generalEmissionProbabilities.containsKey(key)) {
									emission = JointModel.generalEmissionProbabilities.get(key);
								}
								else if(JointModel.generalEmissionProbabilitiesNegative.containsKey(key)) {
									emission = JointModel.generalEmissionProbabilitiesNegative.get(key);
								}
								
								valuesInit.add(emission
										+ JointModel.initialProbabilities.get("<s>|" + JointModel.tagList.get(i)));  // calculate initial <s>-i
								
								for (int k2 = 0; k2 < JointModel.tagSize; k2++) {
									
									double emission2 = 0;
									
									String key2 = JointModel.tagList.get(k2) + "|" + n.word;
									
									if(decideOriginal.equals("original") && JointModel.generalEmissionProbabilities.containsKey(key2)) {
										emission = JointModel.generalEmissionProbabilities.get(key2);
									}
									else if(JointModel.generalEmissionProbabilitiesNegative.containsKey(key2)) {
										emission = JointModel.generalEmissionProbabilitiesNegative.get(key2);
									}
									
									valuesInitDenom.add(emission2
											+ JointModel.initialProbabilities.get("<s>|" + JointModel.tagList.get(k2)));  // calculate initial <s>- all tags
								}
							}

						}
						
						if (a == 0) {

							transOriginalNum += Math.exp(Tools.logSumOfExponentials(values));
							transOriginalDenom += Math.exp(Tools.logSumOfExponentials(valuesDenom));

							values.clear();
							valuesDenom.clear();
							

							initOriginalNum += Math.exp(Tools.logSumOfExponentials(valuesInit));
							initOriginalDenom += Math.exp(Tools.logSumOfExponentials(valuesInitDenom));

							initOriginalEndNum += Math.exp(Tools.logSumOfExponentials(valuesInitEnd));
							initOriginalEndDenom += Math.exp(Tools.logSumOfExponentials(valuesInitEndDenom));

							valuesInit.clear();
							valuesInitDenom.clear();

							valuesInitEnd.clear();
							valuesInitEndDenom.clear();
							
						} else {

							transNeighborNum += Math.exp(Tools.logSumOfExponentials(values));
							transNeighborDenom += Math.exp(Tools.logSumOfExponentials(valuesDenom));

							values.clear();
							valuesDenom.clear();
							
							initNeighborNum += Math.exp(Tools.logSumOfExponentials(valuesInit));
							initNeighborDenom += Math.exp(Tools.logSumOfExponentials(valuesInitDenom));

							initNeighborEndNum += Math.exp(Tools.logSumOfExponentials(valuesInitEnd));
							initNeighborEndDenom += Math.exp(Tools.logSumOfExponentials(valuesInitEndDenom));

							valuesInit.clear();
							valuesInitDenom.clear();

							valuesInitEnd.clear();
							valuesInitEndDenom.clear();

						}
					}

				}
				
				orgScore += divide(transOriginalNum, transOriginalDenom);  //divide for all sentences transOriginalNum / transOriginalDenom
				negScore += divide(transNeighborNum, transNeighborDenom);  //divide for all negative sentences transNeighborNum / transNeighborDenom

				/* Calculation for initial */
				
				initOrgScore += divide(initOriginalNum, initOriginalDenom);
				initNegScore += divide(initNeighborNum, initNeighborDenom);
				
				initOrgScoreEnd += divide(initOriginalEndNum, initOriginalEndDenom);
				initNegScoreEnd += divide(initNeighborEndNum, initNeighborEndDenom);

				String key = (JointModel.tagList.get(i) + "|" + JointModel.tagList.get(j));
				double score = orgScore - negScore;
				gradTransitionProbabilities.put(key, score);
				
				String initKey = "<s>|" + JointModel.tagList.get(i);
				double initScore = initOrgScore - initNegScore;
				gradInitialProbabilities.put(initKey, initScore);

				String initKeyEnd = JointModel.tagList.get(i) + "|</s>";
				double initScoreEnd = initOrgScoreEnd - initNegScoreEnd;
				gradInitialProbabilities.put(initKeyEnd, initScoreEnd);

			}

			/* re-estimation of emission probabilities */
			
			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : JointModel.sentenceGlobals) {

				for (int a = 0; a < globalMap.size(); a++) {
					HashMap<Integer, Node> lattice = globalMap.get(a);
					for (int k = 1; k < lattice.size(); k++) {
						Node node = lattice.get(k);
						
						if(a == 0) {
							
							double logScore;
							for (Pair<String, Integer> candidate : JointModel.getCandidates(node.word, false)) {
								HashMap<Integer, Double> features = JointModel.getFeatures(node.word,
										candidate.getKey(), candidate.getValue());
								logScore = Tools.featureWeightProduct(features);

								for (int featureIndex : features.keySet()) {
									String feature = JointModel.index2Feature.get(featureIndex);
									
									if (feature.startsWith(JointModel.tagList.get(i) + "|")) {	// Tag-dependent
										if (coarseWordProbabilities.containsKey(feature)) {
											double[] logArray = new double[2];
											logArray[0] = coarseWordProbabilities.get(feature);
											logArray[1] = logScore;

											coarseWordProbabilities.put(feature, Tools.logSumOfExponentials(logArray));
										} else {
											coarseWordProbabilities.put(feature, logScore);
										}
									}
									
									if(!feature.contains("|") && i == 0) {
										if (coarseWordProbabilities.containsKey(feature)) {
											double[] logArray = new double[2];
											logArray[0] = coarseWordProbabilities.get(feature);
											logArray[1] = logScore;

											coarseWordProbabilities.put(feature, Tools.logSumOfExponentials(logArray));
										} else {
											coarseWordProbabilities.put(feature, logScore);
										}
									}
								}
							}

							for (String feature : coarseWordProbabilities.keySet()) {
								String key = JointModel.tagList.get(i) + "|" + node.word;
//								double prob = divide(Math.exp(coarseWordProbabilities.get(feature)),
//										JointModel.generalEmissionProbabilities.get(key));
								
								// XXX
								double prob = Math.exp(coarseWordProbabilities.get(feature) - JointModel.emissionDenominator.get(node.word)); //* JointModel.generalEmissionProbabilitiesBackup.get(key);

//								double featureWeightProd = coarseWordProbabilities.get(feature) * JointModel.weights.get(feature);
//								double prob = Math.exp(featureWeightProd - JointModel.generalEmissionProbabilities.get(key));
								
								if (gradCoarseProbabilities.containsKey(feature)) {
									double val = gradCoarseProbabilities.get(feature);

									gradCoarseProbabilities.put(feature, prob * val);
								} else {
									gradCoarseProbabilities.put(feature, prob);
								}
							}
							
							coarseWordProbabilities.clear();
						}
						else {	// Negatives
							double logScore;
							ArrayList<String> neighbors = JointModel.getNeighbors(node.word);
							
					        for(String neighbor : neighbors) {
					            for(Pair<String, Integer> candidate : JointModel.getCandidates(neighbor, false)) {
					                HashMap<Integer, Double> features = JointModel.getFeatures(neighbor, candidate.getKey(), candidate.getValue());
					                logScore = Tools.featureWeightProduct(features);
					                
									for (int featureIndex : features.keySet()) {
										String feature = JointModel.index2Feature.get(featureIndex);
										if (feature.startsWith(JointModel.tagList.get(i) + "|")) {
											if (coarseWordProbabilities.containsKey(feature)) {
												double[] logArray = new double[2];
												logArray[0] = coarseWordProbabilities.get(feature);
												logArray[1] = logScore;

												coarseWordProbabilities.put(feature,
														Tools.logSumOfExponentials(logArray));
											} else {
												coarseWordProbabilities.put(feature, logScore);
											}
										}
										if(!feature.contains("|") && i == 0) {
											if (coarseWordProbabilities.containsKey(feature)) {
												double[] logArray = new double[2];
												logArray[0] = coarseWordProbabilities.get(feature);
												logArray[1] = logScore;

												coarseWordProbabilities.put(feature,
														Tools.logSumOfExponentials(logArray));
											} else {
												coarseWordProbabilities.put(feature, logScore);
											}
										}
									}
					            }
					        }
					        
//					        for(String neighborWord : JointModel.getNeighbors(node.word)) {
						        for (String feature : coarseWordProbabilities.keySet()) {
									String key = JointModel.tagList.get(i) + "|" + node.word;
//									double prob = divide(Math.exp(coarseWordProbabilities.get(feature)),
//											JointModel.generalEmissionProbabilitiesNegative.get(key));
									
									// XXX
									double prob = Math.exp(coarseWordProbabilities.get(feature) - JointModel.emissionDenominatorNegative.get(node.word)); //* JointModel.generalEmissionProbabilitiesBackup.get(key);
									
//									double featureWeightProd = coarseWordProbabilities.get(feature) * JointModel.weights.get(feature);
//									double prob = Math.exp(featureWeightProd - JointModel.generalEmissionProbabilitiesNegative.get(key));
	
									if (gradCoarseProbabilitiesNegative.containsKey(feature)) {
										double val = gradCoarseProbabilitiesNegative.get(feature);
	
										gradCoarseProbabilitiesNegative.put(feature, prob * val);
									} else {
										gradCoarseProbabilitiesNegative.put(feature, prob);
									}
								}
//					        }
							
							coarseWordProbabilities.clear();
						}
					}
				}
			}
		}

//		System.out.println("For gradCoarse:");
//		System.out.println("Min: " + Collections.min(gradCoarseProbabilities.values()) + " Max: " + Collections.max(gradCoarseProbabilities.values()));
//		
//		System.out.println("For gradCoarseNegative:");
//		System.out.println("Min: " + Collections.min(gradCoarseProbabilitiesNegative.values()) + " Max: " + Collections.max(gradCoarseProbabilitiesNegative.values()));
		
		for(String feature : gradCoarseProbabilities.keySet()) {
			double prob;
			
			if(gradCoarseProbabilitiesNegative.containsKey(feature)) {
				prob = gradCoarseProbabilities.get(feature) - gradCoarseProbabilitiesNegative.get(feature);
			}
			else {
				prob = gradCoarseProbabilities.get(feature);
			}
			gradEmissionProbabilities.put(feature, prob);
		}
		
		grad = Main.createGradArray(gradTransitionProbabilities, gradInitialProbabilities, gradEmissionProbabilities,
				JointModel.gradFeature2Index);

		for (int i = 0; i < grad.length; i++) {
			grad[i] *= -1;
		}

		for (int i = 0; i < grad.length; i++) {

			if (!Double.isFinite(grad[i]))
				grad[i] = 0.0;

			grad[i] += 2 * LAMBDA_EM * iterWeights[i];
		}

		Main.pw.println(
				"***********************************************NEW ITERATION*********************************************");

		return grad;

	}

	/**
	 * @param i
	 *            the number of state s_i
	 * @param j
	 *            the number of state s_j
	 * @return P
	 */
	public double p(int i, int j, Node node, Node next, String decideOriginal) {
		double num = 0.0;

		double emission = 0;
		
		String key = JointModel.tagList.get(i) + "|" + next.word;
		
		if(decideOriginal.equals("original") && JointModel.generalEmissionProbabilities.containsKey(key)) {
			emission = JointModel.generalEmissionProbabilities.get(key);
		}
		else if(JointModel.generalEmissionProbabilitiesNegative.containsKey(key)) {
			emission = JointModel.generalEmissionProbabilitiesNegative.get(key);
		}
		
		num = JointModel.transitionProbabilities.get(JointModel.tagList.get(i) + "|" + JointModel.tagList.get(j))
				+ emission
				+ Math.log(next.beta.get(j));

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

	public double sumListValuesExp(List<Double> list) {

		return Math.exp(Tools.logSumOfExponentials((ArrayList<Double>) list));
	}

	public double sumListValuesLog(List<Double> list) {

		return Tools.logSumOfExponentials((ArrayList<Double>) list);
	}
	
	public void findMinAndMax(double[] array) {
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (int i = 0; i < array.length; i++) {
			if (array[i] < min)
				min = array[i];
			if (array[i] > max)
				max = array[i];
		}
		
		System.out.println("Min: " + min + " Max: " + max);
	}
	
}
