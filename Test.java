import java.util.ArrayList;

public class Test {
	
	public static void main(String[] args) {

		ArrayList<String> verboseFlags = getVerboseFlags(args);

		double[][] matrix = 
			{{3, 2, -2},
			 {-2, -1, 3},
			 {6, 0, -12}};

		GaussianArithmetic ga = new GaussianArithmetic(verboseFlags);

		ga.printMatrix(matrix);

		double[] solutionVector = ga.solveWithPartialPivoting(matrix);
		ga.printSolutionVector(solutionVector);
	}

	public static ArrayList<String> getVerboseFlags(String[] args) {
		ArrayList<String> verboseFlags = new ArrayList<String>();

		for(int arg=0; arg<args.length; arg++) {
			String currArg = args[arg];

			if(currArg.equals("-verbose")) {
				verboseFlags.add("verbose");
			}
			else if(currArg.equals("-forward_elimination")) {
				verboseFlags.add("forward_elimination");
			}
			else {
				System.out.println("Unrecognized flag: " + currArg);
			}
		}

		return verboseFlags;
	}
}