using System;

public static class Cholesky {

    // choldc from numerical recipes p97, except that a new matrix is created and returned instead of placing back into a.
    public static void Choldc(double[,] a, double[] p) {
        int n = a.GetLength(0);
        if (a.GetLength(1) != n) {
            throw new InvalidOperationException("matrix is not square");
        }
        
        for (int i=0; i<n; i++) {
            for (int j=i; j<n; j++) {
                double sum = a[i,j];
                for (int k=i-1; k>=0; k--) sum -= a[i,k]*a[j,k];
                if (i == j) {
                    if (sum <= 0) {
                        throw new InvalidOperationException("matrix is not positive definite");
                    }
                    p[i] = Math.Sqrt(sum); 
                } else {
                    a[j,i] = sum/p[i];
                }
            }
        }
    }

    // cholsl from numerical recipes p97. Solves A.x = b, given A is a LU decomposition.
    // p is the diagonal elements
    public static void BackSubstitution(double[,] a, double[] p, double[] b, double[] x) {
        int n = b.Length;
        if (a.GetLength(0) != n || a.GetLength(1) != n || x.Length != n) {
            throw new InvalidOperationException("matrix and vector have different lengths");
        }

        // var x = new double[n];
        for (int i=0; i<n; i++) {
            double sum = b[i];
            for (int k=i-1; k>=0; k--) {
                sum -= a[i,k]*x[k];
            }
            x[i] = sum/p[i];
        }
        for (int i=n-1; i>=0; i--) {
            double sum = x[i];
            for (int k=i+1; k<n; k++) {
                sum -= a[k,i]*x[k];
            }
            x[i] = sum/p[i];
        }
    }
}
