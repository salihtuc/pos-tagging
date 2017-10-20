
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

	public double functionValue = 0.0;
	public static final double LAMBDA_EM = 1.0;
	
			
	public int iterCount = 0;
	ArrayList<Double> originalList = new ArrayList<>();
	ArrayList<Double> negativeList = new ArrayList<>();

	@Override
	public FunctionValues getValues(double[] point) {

		JointModel.normalizeFeatures(point, "transition");
		JointModel.normalizeFeatures(point, "emission");
		JointModel.normalizeFeatures(point, "initial");
		
		if (iterCount != 0) {

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
			for (HashMap<Integer, Node> lattice : globalMap.values()) { // Each
																		// sentence's
																		// lattices
				for (Node endState : Main.returnEndStates(lattice)) {
					if (counter == 0) {
						sumNum += sumListValuesExp(endState.tagScores); // Original
																		// lattice's
																		// score
						counter++;
					} else {
						sumDenom += sumListValuesExp(endState.tagScores); // Negative
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

				for (HashMap<Integer, HashMap<Integer, Node>> globalMap : JointModel.sentenceGlobals) {
					
					for (int a = 0; a < globalMap.size(); a++) {

						HashMap<Integer, Node> lattice = globalMap.get(a);
						ArrayList<Integer> startStates = JointModel.returnStartStates(lattice);

						for (int k = 1; k < lattice.size() - 1; k++) {
							Node n = lattice.get(k);
							for (int k1 : n.next) {
								Node n2 = lattice.get(k1);
								values.add(p(i, j, n, n2));  //calculate transition i-j
								
								for (int k2 = 0; k2 < JointModel.tagSize; k2++) {
									valuesDenom.add(p(i, k2, n, n2));  //calculate transition i-k
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
								
								if(JointModel.generalEmissionProbabilities.containsKey(JointModel.tagList.get(i) + "|" + n.word)) {
									emission = JointModel.generalEmissionProbabilities.get(JointModel.tagList.get(i) + "|" + n.word);
								}
								
								valuesInit.add(emission
										+ JointModel.initialProbabilities.get("<s>|" + JointModel.tagList.get(i)));  // calculate initial <s>-i
								for (int k2 = 0; k2 < JointModel.tagSize; k2++) {
									
									double emission2 = 0;
									
									if(JointModel.generalEmissionProbabilities.containsKey(JointModel.tagList.get(k2) + "|" + n.word)) {
										emission2 = JointModel.generalEmissionProbabilities.get(JointModel.tagList.get(k2) + "|" + n.word);
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

			/* re-estimation of emission probabilities NEW VERSION */
			
			HashMap<String, Double> originalMap = new HashMap<>();
			HashMap<String, Double> negativeMap = new HashMap<>();
			
			HashMap<String, Double> originalwordsNum = new HashMap<String, Double>();
			HashMap<String, Double> neighborsWord = new HashMap<String, Double>();

			ArrayList<Double> listNum = new ArrayList<>();
			ArrayList<Double> listDenom = new ArrayList<>();

			double emissionOriginalDenom = 0.0;
			double emissionNeighborDenom = 0.0;

			for (HashMap<Integer, HashMap<Integer, Node>> globalMap : JointModel.sentenceGlobals) {

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

									originalwordsNum.put(node.word, Tools.logSumOfExponentials(list));

								} else {
									originalwordsNum.put(node.word, g);
								}
								listNum.add(g);

							} else {
								if (neighborsWord.containsKey(node.word)) {
									ArrayList<Double> list = new ArrayList<>();

									list.add(neighborsWord.get(node.word));
									list.add(g);

									neighborsWord.put(node.word, Tools.logSumOfExponentials(list));

								} else {
									neighborsWord.put(node.word, g);
								}

								listDenom.add(g);
							}

					}
				}
			}
			
			emissionOriginalDenom += Tools.logSumOfExponentials(listNum);
			emissionNeighborDenom += Tools.logSumOfExponentials(listDenom);
			listNum.clear();
			listDenom.clear();
			for (String s : originalwordsNum.keySet()) {
				double valOrg = originalwordsNum.get(s);

				if (originalMap.containsKey(s)) {

					ArrayList<Double> list = new ArrayList<>();

					list.add(originalMap.get(s));
					list.add((valOrg - emissionOriginalDenom));

					originalMap.put(s, Tools.logSumOfExponentials(list));

				} else {
					originalMap.put(s, (valOrg - emissionOriginalDenom));
				}

				double valNeg = neighborsWord.get(s);

				if (negativeMap.containsKey(s)) {

					ArrayList<Double> list = new ArrayList<>();

					list.add(negativeMap.get(s));
					list.add((valNeg - emissionNeighborDenom));

					negativeMap.put(s, Tools.logSumOfExponentials(list));

				} else {
					negativeMap.put(s, (valNeg - emissionNeighborDenom));
				}
			}
			originalwordsNum.clear();
			neighborsWord.clear();

			Iterator<Entry<String, Double>> originals = originalMap.entrySet().iterator();
			Iterator<Entry<String, Double>> neighbors = negativeMap.entrySet().iterator();
			while (originals.hasNext()) {
				Entry<String, Double> pairs = originals.next();
				Double firstVal = (Double) pairs.getValue();
				Entry<String, Double> pairs2 = neighbors.next();
				Double secondVal = (Double) pairs2.getValue();

				String key = JointModel.tagList.get(i) + "|" + pairs.getKey();

				gradEmissionProbabilities.put(key, Math.exp(firstVal) - Math.exp(secondVal));

			}
			
			originalMap.clear();
			negativeMap.clear();
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
	public double p(int i, int j, Node node, Node next) {
		double num = 0.0;

		double emission = 0;
		
		if(JointModel.generalEmissionProbabilities.containsKey(JointModel.tagList.get(i) + "|" + next.word)) {
			emission = JointModel.generalEmissionProbabilities.get(JointModel.tagList.get(i) + "|" + next.word);
		}
		
		num = JointModel.transitionProbabilities.get(JointModel.tagList.get(i) + "|" + JointModel.tagList.get(j))
				+ emission
				+ Math.log(next.beta.get(j));

		return num;
	}

	/** computes gamma(i, node) */

	public double gamma(int i, Node node) {
		double num = 0.0;
		num = Math.log(node.alpha.get(i)) + Math.log(node.beta.get(i));

		return num;
	}
	
	public double gammaEmission() {
		
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
	
}
