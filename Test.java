import java.util.ArrayList;

public class Test {
	
	public static void main(String[] args) {

		ArrayList<String> verboseFlags = getVerboseFlags(args);

		double[][] matrix = 
			{{3, 2, -2},
			 {-2, -1, 3},
			 {6, 0, -12}};

		double[] b = 
			{-5,
			  6,
			 -18};

		GaussianArithmetic ga = new GaussianArithmetic(verboseFlags);

		ga.printMatrixWithRHS(matrix, b);

		double[] solutionVector = ga.solveWithPartialPivoting(matrix, b);
		ga.printSolutionVector(solutionVector);
	}

	public static ArrayList<String> getVerboseFlags(String[] args) {
		ArrayList<String> verboseFlags = new ArrayList<String>();

		for(int arg=0; arg<args.length; arg++) {
			String currArg = args[arg];

			if(currArg.equals("-verbose")) {
				verboseFlags.add("verbose");
			}
			else if(currArg.equals("-partial_pivoting")) {
				verboseFlags.add("partial_pivoting");
			}
			else {
				System.out.println("Unrecognized flag: " + currArg);
			}
		}

		return verboseFlags;
	}
}