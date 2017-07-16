
import java.util.HashMap;

import lbfgsb.DifferentiableFunction;
import lbfgsb.FunctionValues;

/**
 * Created by ghostof2007 on 5/6/14.
 *
 * Function to be used for LBFGS
 */
public class Function implements DifferentiableFunction {

    // -------------------------------------- LBFGS-B ---------------------------------------------------------


    @Override
    public FunctionValues getValues(double[] point) {

    		
    	
        return new FunctionValues(0.0, gradient(point));
    }


    double[] gradient(double[] iterWeights) {
        return new double[1];
    }


}

