import java.util.ArrayList;

public class GaussianArithmetic {

	private boolean VERBOSE=false, PARTIAL_PIVOTING=false;

	public GaussianArithmetic() {
		// If you don't care about verbose flags
	}

	public GaussianArithmetic(ArrayList<String> verboseFlags) {
		if(verboseFlags.contains("verbose")) {
			VERBOSE = true;
		}
		if(verboseFlags.contains("partial_pivoting")) {
			PARTIAL_PIVOTING = true;
		}
	}

	public double[] solveWithPartialPivoting(double[][] matrix, double[] b) {

		// FORWARD ELIMINATION
		int[] ip = new int[matrix[0].length];
		for(int row=0; row<ip.length; row++) {
			ip[row] = -1;
		}

		if(VERBOSE || PARTIAL_PIVOTING) { printMatrixWithIpVector(matrix, ip); }

		if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("---------------FORWARD ELIMINATION-------------\n"); }

		// The outer loop - this eliminates the variable k
		for(int k=0; k<matrix[0].length; k++) {

			if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("Current column to pivot on is " + (k+1) + "\n"); }

			// Find the largest of the candidate pivots
			int maxPivotRow = maxPivotRow(matrix, k);
			if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("The row with the maximum pivot is row " + (maxPivotRow+1) + "." + "\n"); }

			// If the largest is zero, no possible pivots. Get out of here
			if(matrix[maxPivotRow][k]==0) {
				if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("Warning: Pivot in Gaussain Elimination is 0!" + "\n"); }
				return null;
			}

			// Swap the rows to get the pivot into position
			if(k!=maxPivotRow) {
				if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("Swapping row " + (k+1) + " with row " + (maxPivotRow+1) + "\n"); }

				for(int col=k; col<matrix[0].length; col++) {
					double topVal = matrix[k][col];
					double bottomVal = matrix[maxPivotRow][col];
					matrix[k][col] = bottomVal;
					matrix[maxPivotRow][col] = topVal;
				}

				ip[k] = maxPivotRow;

				if(VERBOSE || PARTIAL_PIVOTING) { printMatrixWithIpVector(matrix, ip); }
			}
			else {
				ip[k] = k;

				if(VERBOSE || PARTIAL_PIVOTING) { printMatrixWithIpVector(matrix, ip); }
			}

			// Loop on the rows
			for(int i=k+1; i<matrix.length; i++) {
				// Get the multiplier for the row i
				String numerator = String.format("%-4.6s", matrix[i][k]);
				String denominator = String.format("%-4.6s", matrix[k][k]);
				String multiplier = String.format("%-4.6s", (matrix[i][k]/matrix[k][k]));
				if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("Performing SAXPY on each item in row " + (i+1) + "\n"); }
				if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("The multiplier for this SAXPY is A_" + (i+1) + "," + (k+1) + " / A_" + (k+1) + "," + (k+1) + " = " + numerator + " / " + denominator + " = " + multiplier); }
				
				if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("Setting A_" + (i+1) + "," + (k+1) + " = multiplier = " + multiplier + "\n"); }
				matrix[i][k] = matrix[i][k] / matrix[k][k];
				
				if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("SAXPY: r_" + (i+1) + " = r_" + (i+1) + " - ( multiplier * r_" + (k+1) + " )" + " = r_" + (i+1) + " - ( " + multiplier + " * r_" + (k+1) + " )" + "\n"); }

				// Loop on the columns - innermost loop
				for(int j=k+1; j<matrix[0].length; j++) {
					matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
				}

				if(VERBOSE || PARTIAL_PIVOTING) { printMatrixWithIpVector(matrix, ip); }
			}

		}

		// SOLVING PART 1

		if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("---------------SOLVING PART ONE-------------\n"); } 

		if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("State of b vector from input.\n"); }
		if(VERBOSE || PARTIAL_PIVOTING) { printVector(b); }

		if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("Beginning solving of b vector from stored ip vector and modified A matrix\n"); }

		for(int k=0; k<matrix[0].length; k++) {

			// Perform swapping if necessary
			if(ip[k]!=k) {
				double topVal = b[k];
				double bottomVal = b[ip[k]];
				b[k] = bottomVal;
				b[ip[k]] = topVal;

				if(VERBOSE || PARTIAL_PIVOTING) { 
					System.out.println("Swap b_" + (k+1) + " with b_" + (ip[k]+1) + "\n");
					printVector(b); 
				}

				for(int i=k+1; i<matrix[0].length; i++) {
					String bIStr = String.format("%-4.6s", b[i]);
					String bkStr = String.format("%-4.6s", b[k]);
					String matrixPosStr = String.format("%-4.6s", matrix[i][k]);
					if(VERBOSE || PARTIAL_PIVOTING) { 
						System.out.print("b_" + (i+1) + " = b_" + (i+1) + " - ( A_" + (i+1) + "," + (k+1) + " * b_" + (k+1) + " ) = "); 
					}
					
					b[i] = b[i] - (matrix[i][k]*b[k]);

					String bIModStr = String.format("%-4.6s", b[i]);
					if(VERBOSE || PARTIAL_PIVOTING) { 
						System.out.println(bIStr + " - ( " + matrixPosStr + " * " + bkStr + " ) = " + bIModStr + "\n");
						printVector(b); 
					}
				}
			}


		}

		// SOLVING PART 2

		if(VERBOSE || PARTIAL_PIVOTING) { System.out.println("---------------SOLVING PART TWO-------------\n"); }

		double[] solutionVector = new double[matrix.length];
		for(int row=0; row<solutionVector.length; row++) {
			solutionVector[row] = 0;
		}

		// Solve with the upper triangular system
		for(int i=matrix.length-1; i>=0; i--) {

			// Solve for b
			for(int j=i+1; j<matrix.length; j++) {
				if(VERBOSE || PARTIAL_PIVOTING) { 
					String bIStr = String.format("%-4.6s", b[i]);
					String solVectorStr = String.format("%-4.6s", solutionVector[j]);
					String matrixPosStr = String.format("%-4.6s", matrix[i][j]);
					System.out.print("b_" + (i+1) + " = b_" + (i+1) + " - ( A_" + (i+1) + "," + (j+1) + " * x_" + (j+1) + " = ");
					System.out.print(bIStr + " - ( " + matrixPosStr + " * " + solVectorStr + " = "); 
				}

				b[i] = b[i] - matrix[i][j] * solutionVector[j];
			
				if(VERBOSE || PARTIAL_PIVOTING) { 
					String bIModStr = String.format("%-4.6s", b[i]);
					System.out.println(bIModStr + "\n"); 
				}
			}

			if(VERBOSE || PARTIAL_PIVOTING) { 
				System.out.print("x_" + (i+1) + " = b_" + (i+1) + " / A_" + (i+1) + "," + (i+1) + " = "); 
			}
			
			// Set x = A / b
			solutionVector[i] = b[i] / matrix[i][i];

			if(VERBOSE || PARTIAL_PIVOTING) { 
				String vectorSolStr = String.format("%-4.6s", solutionVector[i]);
				System.out.println(vectorSolStr + "\n"); 
			}
		}

		return solutionVector;
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
				System.out.printf(" %-4.6s  ", matrix[row][col]);
			}
			System.out.println("]");
		}
		System.out.println();
	}

	public void printMatrixWithRHS(double[][] matrix, double[] rhs) {
		for(int row=0; row<matrix.length; row++) {
			System.out.print("[ ");
			for(int col=0; col<matrix[0].length; col++) {
				System.out.printf(" %-4.6s  ", matrix[row][col]);
			}
			System.out.print("]");

			
			System.out.printf(" [ %-4.6s ]\n", rhs[row]);
		}

		System.out.println();
	}

	public void printVector(double[] vector) {
		for(int row=0; row<vector.length; row++) {
				System.out.printf("[ %-4.6s ]\n", vector[row]);
		}
		System.out.println();
	}

	private void printMatrixWithIpVector(double[][] matrix, int[] ip) {
		for(int row=0; row<matrix.length; row++) {
			System.out.print("[ ");
			for(int col=0; col<matrix[0].length; col++) {
				System.out.printf(" %-4.6s  ", matrix[row][col]);
			}
			System.out.print("]");

			if(ip[row]!=-1) {
				System.out.println(" [ " + (ip[row]+1) + " ]");
			}
			else {
				System.out.println(" [ ? ]");	
			}
		}

		System.out.println();
	}

	public void printSolutionVector(double[] solutionVector) {
		if(solutionVector == null) {
			System.out.println("Solution vector is null");
			return;
		}

		System.out.println("Solution vector");
		System.out.println("---------------");

		for(int row = 0; row < solutionVector.length; row++) {
			String currXStr = String.format("%-4.6s", solutionVector[row]);
			System.out.println("[ x_" + (row+1) + " = " + currXStr + " ]");
		}
	}
}