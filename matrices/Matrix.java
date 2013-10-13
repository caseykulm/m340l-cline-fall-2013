package matrices;

import java.text.DecimalFormat;

public class Matrix {

	/**
	 * Contents of the matrix.
	 */
	private double[][] matrix;
	
	/**
	 * Number of rows in the matrix.
	 */
	private int m;
	
	
	/**
	 * Number of columns in the matrix.
	 */
	private int n;
	
	/**
	 * Creates an m x n matrix with zeros every 
	 * spot.
	 * 
	 * @param m Number of rows
	 * @param n Number of columns
	 */
	public Matrix(int m, int n) {
		this.m = m;
		this.n = n;
		matrix = new double[m][n];
		
		for(int i=0; i<m; i++) {
			for(int j=0; j<n; j++) {
				matrix[i][j] = 0;
			}
		}
	}
	
	/**
	 * Creates an m x n matrix object from the 
	 * provided matrix.
	 * 
	 * @param matrix Input matrix
	 */
	public Matrix(double[][] matrix) {
		this.m = matrix.length;
		this.n = matrix[0].length;
		this.matrix = matrix;
	}
	
	/**
	 * @return The number of rows in the matrix
	 */
	public int getRows() {
		return m;
	}
	
	/**
	 * @return The number of columns in the matrix
	 */
	public int getCols() {
		return n;
	}
	
	/**
	 * @param i Row index
	 * @param j Column index
	 * @return The A_ij element in the matrix
	 */
	public double getVal(int i, int j) {
		return matrix[i][j];
	}
	
	// ELEMENTARY ROW OPERATIONS
	
    /**
     * Swap the i'th and j'th rows in the matrix.
     * 
     * @param rowI i'th row in the matrix
     * @param rowJ j'th row in the matrix
     */
    private void swap(int rowI, int rowJ) {
        double[] temp = matrix[rowI];
        matrix[rowI] = matrix[rowJ];
        matrix[rowJ] = temp;
    }
    
    /**
     * Replaces the rowNum'th row by scaling it by scaleFactor.
     * 
     * @param rowNum Index of the row to be scaled
     * @param scaleFactor Factor to scale the rows by
     */
    private void scale(int rowNum, int scaleFactor) {
    	for(int i=0; i<n; i++) {
    		matrix[rowNum][i] *= scaleFactor;
    	}
    }
    
    /**
     * Performs the SAXPY operation on a matrix.
     * EX: r_2 = r_2 + 3*r_1
     * 
     * @param changeRow Row to be changed
     * @param scaleRow Row to be scaled and added
     * @param scaleFactor Factor to scale the scaleRow by
     */
    private void saxpy(int changeRow, int scaleRow, int scaleFactor) {
    	for(int i=0; i<n; i++) {
    		matrix[changeRow][i] += scaleFactor*matrix[scaleRow][i];
    	}
    }
    
    /* (non-Javadoc)
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(Object obj) {
    	if(!(obj instanceof Matrix))
    		return false;
    	
    	Matrix otherMatrix = (Matrix) obj;
    			
    	if(this.m != otherMatrix.getRows() || this.n != otherMatrix.getCols())
    		return false;
    	
    	for(int i=0; i<m; i++) {
    		for(int j=0; j<n; j++) {
    			if(matrix[i][j] != otherMatrix.getVal(i, j))
    				return false;
    		}
    	}
    	
    	return true;
    }
    
    /* (non-Javadoc)
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
    	String matrixString = "";
    	
    	for(int i=0; i<m; i++) {
    		matrixString += "[ ";
    		for(int j=0; j<n; j++) {
    			matrixString +=  String.format("%9.3f ", matrix[i][j]);
    		}
    		matrixString += "]\n";
    	}
    	
    	return matrixString;
    }
    
    /**
     * Displays a formatted matrix to System.out of the format 
     * String.format("%9.3f", matrix[i][j]) for each element.
     */
    public void display() {
    	System.out.println(this.toString());
    }
}
