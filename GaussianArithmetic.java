import java.util.ArrayList;

public class GaussianArithmetic {

	private boolean VERBOSE=false, FORWARD_ELIMINATION=false;

	public GaussianArithmetic() {
		// If you don't care about verbose flags
	}

	public GaussianArithmetic(ArrayList<String> verboseFlags) {
		if(verboseFlags.contains("verbose")) {
			VERBOSE = true;
		}
		if(verboseFlags.contains("forward_elimination")) {
			FORWARD_ELIMINATION = true;
		}
	}

	public double[] solveWithPartialPivoting(double[][] matrix) {

		// FORWARD ELIMINATION
		int[] ip = new int[matrix[0].length];
		for(int row=0; row<ip.length; row++) {
			ip[row] = -1;
		}

		if(VERBOSE || FORWARD_ELIMINATION) { printMatrixWithIpVector(matrix, ip); }

		// The outer loop - this eliminates the variable k
		for(int k=0; k<matrix.length; k++) {

			if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("---------------FORWARD ELIMINATION-------------\n"); }
			if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("Current column to pivot on is " + (k+1) + "\n"); }

			// Find the largest of the candidate pivots
			int maxPivotRow = maxPivotRow(matrix, k);
			if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("The row with the maximum pivot is row " + (maxPivotRow+1) + "." + "\n"); }

			// If the largest is zero, no possible pivots. Get out of here
			if(matrix[maxPivotRow][k]==0) {
				if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("Warning: Pivot in Gaussain Elimination is 0!" + "\n"); }
				return null;
			}

			// Swap the rows to get the pivot into position
			if(k!=maxPivotRow) {
				if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("Swapping row " + (k+1) + " with row " + (maxPivotRow+1) + "\n"); }

				for(int col=k; col<matrix[0].length; col++) {
					double topVal = matrix[k][col];
					double bottomVal = matrix[maxPivotRow][col];
					matrix[k][col] = bottomVal;
					matrix[maxPivotRow][col] = topVal;
				}

				ip[k] = maxPivotRow;

				if(VERBOSE || FORWARD_ELIMINATION) { printMatrixWithIpVector(matrix, ip); }
			}
			else {
				ip[k] = k;

				if(VERBOSE || FORWARD_ELIMINATION) { printMatrixWithIpVector(matrix, ip); }
			}

			// Loop on the rows
			for(int i=k+1; i<matrix.length; i++) {
				// Get the multiplier for the row i
				String numerator = String.format("%-4.6s", matrix[i][k]);
				String denominator = String.format("%-4.6s", matrix[k][k]);
				String multiplier = String.format("%-4.6s", (matrix[i][k]/matrix[k][k]));
				if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("Performing SAXPY on each item in row " + (i+1) + "\n"); }
				if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("The multiplier for this SAXPY is A_" + (i+1) + "," + (k+1) + " / A_" + (k+1) + "," + (k+1) + " = " + numerator + " / " + denominator + " = " + multiplier); }
				
				if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("Setting A_" + (i+1) + "," + (k+1) + " = multiplier = " + multiplier + "\n"); }
				matrix[i][k] = matrix[i][k] / matrix[k][k];
				
				if(VERBOSE || FORWARD_ELIMINATION) { System.out.println("SAXPY: r_" + (i+1) + " = r_" + (i+1) + " - ( multiplier * r_" + (k+1) + " )" + " = r_" + (i+1) + " - ( " + multiplier + " * r_" + (k+1) + " )" + "\n"); }

				// Loop on the columns - innermost loop
				for(int j=k+1; j<matrix.length; j++) {
					matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
				}

				if(VERBOSE || FORWARD_ELIMINATION) { printMatrixWithIpVector(matrix, ip); }
			}

		}

		// SOLVING

		return new double[]{1,-2,-4};
	}

	private int maxPivotRow(double[][] matrix, int k) {
		int maxPivotRow = k;
		double record = matrix[maxPivotRow][k];

		for(int row = k; row<matrix.length; row++) {
			if(matrix[row][k]>record) {
				record = matrix[row][k];
				maxPivotRow = row;
			}
		}

		return maxPivotRow;
	}

	public void printMatrix(double[][] matrix) {
		for(int row=0; row<matrix.length; row++) {
			System.out.print("[ ");
			for(int col=0; col<matrix[0].length; col++) {
				System.out.printf(" %-6.6s  ", matrix[row][col]);
			}
			System.out.println("]");
		}
		System.out.println();
	}

	private void printMatrixWithIpVector(double[][] matrix, int[] ip) {
		for(int row=0; row<matrix.length; row++) {
			System.out.print("[ ");
			for(int col=0; col<matrix[0].length; col++) {
				System.out.printf(" %-6.6s  ", matrix[row][col]);
			}
			System.out.print("]");

			if(ip[row]!=-1) {
				System.out.print(" [ " + (ip[row]+1) + " ]");
			}
			else {
				System.out.print(" [ ? ]");	
			}

			System.out.println("\n");
		}
	}

	public void printSolutionVector(double[] solutionVector) {
		if(solutionVector == null) {
			System.out.println("Solution vector is null");
			return;
		}

		System.out.println("Solution vector");
		System.out.println("---------------");

		for(int row = 0; row < solutionVector.length; row++) {
			System.out.println("x_" + (row+1) + " = " + solutionVector[row]);
		}
	}
}