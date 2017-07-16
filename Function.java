
import java.util.HashMap;

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

public class Function implements DifferentiableFunction {

    // -------------------------------------- LBFGS-B ---------------------------------------------------------


    @Override
    public FunctionValues getValues(double[] point) {

    	//TODO	
    	
        return new FunctionValues(0.0, gradient(point));
    }


    double[] gradient(double[] iterWeights) {
        return new double[1];
    }


}

