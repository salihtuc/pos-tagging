import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

public class Tools {
	
	static double dot(String a, String b) {
        double sum = 0.;
        ArrayList<Double> vec1 = JointModel.wordVec.get(a);
        ArrayList<Double> vec2 = JointModel.wordVec.get(b);
        if(vec1==null || vec2==null) return -0.5; //FIXME - make sure not using just plain DOT
        for(int i=0; i < vec1.size(); i++)
            sum += vec1.get(i) * vec2.get(i);
        return sum;
    }

	static int getFeatureIndex(String feature) {
		if (!JointModel.feature2Index.containsKey(feature)) {
			if (JointModel.TEST)
				return -1; // if in testing phase, and feature does not exist already, do not create new

			int index = JointModel.feature2Index.size();
			JointModel.feature2Index.put(feature, index);
			JointModel.index2Feature.add(feature);
			JointModel.weights.add(0.);
			return index;
		}
		return JointModel.feature2Index.get(feature);

	}
	
	static void addFeature(HashMap<Integer, Double> features, String newFeature, double value) {
        int featureIndex = getFeatureIndex(newFeature);
        if(featureIndex!=-1)
            features.put(featureIndex,value);
    }
	
	protected static void printDoubleArray(double[] array) {
		for (double d : array) {
			System.out.print(d + " ");
			// Main.pw.print(d + " ");
		}
		System.out.println();
		// Main.pw.println();
	}
	
	public static double max(double[] values) {
		double max = -Double.MAX_VALUE;
		for (double value : values) {
			if (value > max)
				max = value;
		}
		return max;
	}
	
	public static double logSumOfExponentials(double[] xs) {
		if (xs.length == 1)
			return xs[0];
		double max = max(xs);
		double sum = 0.0;
		for (double x : xs)
			if (x != Double.NEGATIVE_INFINITY)
				sum += Math.exp(x - max);
		return max + java.lang.Math.log(sum);
	}

	public static double logSumOfExponentials(ArrayList<Double> x) {
		double[] xs = new double[x.size()];
		for (int i = 0; i < x.size(); i++)
			xs[i] = x.get(i);
		return logSumOfExponentials(xs);
	}
	
	//	dot product of feature and weights(global)
    static double featureWeightProduct(HashMap<Integer, Double> features) {
        double sum = 0.;
        if(features==null || features.size()==0) return 0.;
        for(int i : features.keySet())
            if( i < JointModel.weights.size())  //check if weight exists for the feature
                sum += features.get(i) * JointModel.weights.get(i);
        return sum;
    }
    
    /** divides two doubles. (0 / 0 = 0!) && (1 / 0 = 0!) */
	public static double divide(double n, double d) {
		if (n == 0 || d == 0)
			return 0;
		else
			return n / d;
	}
	
	public static double sumDoubleListValues(Collection<Double> list) {
		double sum = 0.0;
		
		for(double d : list) {
			sum += d;
		}
		
		return sum;
	}
}
