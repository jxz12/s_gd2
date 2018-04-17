using System;
using System.Collections.Generic;
using GraphStuff;

public static class Majorization {
    public static IEnumerable<double> Chol(int[,] d, Vector2[] positions, double eps=0.00001, int maxIter=1000) {
        int n = positions.Length;

        // first find the laplacian for the left hand side
        var laplacian_w = new double[n,n];
        WeightLaplacian(d, laplacian_w, n);

        // cut out the first row and column
        var cholesky_Lw = new double[n-1,n-1];
        for (int i=1; i<n; i++) {
            for (int j=1; j<n; j++) {
                cholesky_Lw[i-1,j-1] = laplacian_w[i,j];
            }
        }
        // p used to store diagonal values, as in Numerical Recipes
        var cholesky_p = new double[n-1];
        Cholesky.Choldc(cholesky_Lw, cholesky_p);

        // delta = w_ij * d_ij as in Majorization, Gansner et al.
        var deltas = new double[n,n];
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                double dist = d[i,j];
                deltas[i, j] = 1.0 / dist;
            }
        }

        var LXt = new double[n,n]; // the laplacian for the right hand side
        var LXt_Xt = new double[n-1]; // skip the first position as it's fixed to (0,0)
        var Xt1 = new double[n-1]; // X(t+1)

        double prevStress = GraphIO.CalculateStress(d, positions, n);
        // majorize
        for (int k=0; k<maxIter; k++) {
            PositionLaplacian(deltas,positions,LXt,n);

            // solve for x axis
            Multiply_x(LXt,positions,LXt_Xt);
            Cholesky.BackSubstitution(cholesky_Lw,cholesky_p,LXt_Xt,Xt1);
            for (int i=1; i<n; i++) positions[i].x = Xt1[i-1];

            // solve for y axis
            Multiply_y(LXt,positions,LXt_Xt);
            Cholesky.BackSubstitution(cholesky_Lw,cholesky_p,LXt_Xt,Xt1);
            for (int i=1; i<n; i++) positions[i].y = Xt1[i-1];

            double stress = GraphIO.CalculateStress(d, positions, n);
            yield return stress;
            if ((prevStress - stress) / prevStress < eps)
                yield break;
            prevStress = stress;
        }
    }

    
    public static IEnumerable<double> Conj(int[,] d, Vector2[] positions, double eps=0.0001, int maxIter=1000) {
        int n = positions.Length;

        // first find the laplacian for the left hand side
        var laplacian_w = new double[n,n];
        Majorization.WeightLaplacian(d, laplacian_w, n);

        // cut out the first row and column
        var Lw = new double[n-1,n-1];
        for (int i=1; i<n; i++) {
            for (int j=1; j<n; j++) {
                Lw[i-1,j-1] = laplacian_w[i,j];
            }
        }

        // delta = w_ij * d_ij as in Majorization, Gansner et al.
        var deltas = new double[n,n];
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                double dist = d[i,j];
                deltas[i, j] = 1.0 / dist;
            }
        }

        var LXt = new double[n,n]; // the laplacian for the right hand side
        var LXt_Xt = new double[n-1]; // skip the first position as it's fixed to (0,0)
        var Xt1 = new double[n-1]; // X(t+1)

        // temporary variables to speed up conjugate gradient
        var r = new double[n-1];
        var p = new double[n-1];
        var Ap = new double[n-1];

        double prevStress = GraphIO.CalculateStress(d, positions, n);
        // majorize
        for (int k=0; k<maxIter; k++) {
            PositionLaplacian(deltas, positions, LXt, n);

            // solve for x axis
            Multiply_x(LXt, positions, LXt_Xt);
            ConjugateGradient.Cg(Lw, Xt1, LXt_Xt, r, p, Ap, .1, 10);
            for (int i=1; i<n; i++) positions[i].x = Xt1[i-1];

            // solve for y axis
            Multiply_y(LXt, positions, LXt_Xt);
            ConjugateGradient.Cg(Lw, Xt1, LXt_Xt, r, p, Ap, .1, 10);
            for (int i=1; i<n; i++) positions[i].y = Xt1[i-1];

            double stress = GraphIO.CalculateStress(d, positions, n);
            yield return stress;
            if ((prevStress - stress) / prevStress < eps)
                yield break;
            prevStress = stress;
        }
    }

    public static IEnumerable<double> Local(int[,] d, Vector2[] positions, double eps=0.00001, int maxIter=100) {
        int n = positions.Length;

        double prevStress = GraphIO.CalculateStress(d, positions, n);
        // majorize
        for (int k=0; k<maxIter; k++) {
            for (int i=0; i<n; i++) {

                double topSumX=0, topSumY=0, botSum=0;
                for (int j=0; j<n; j++) {
                    if (i!=j) {
                        double d_ij = d[i,j];
                        double w_ij = 1/(d_ij*d_ij);
                        double magnitude = (positions[i] - positions[j]).Magnitude();

                        topSumX += w_ij * (positions[j].x + d_ij*(positions[i].x - positions[j].x)/(magnitude));
                        topSumY += w_ij * (positions[j].y + d_ij*(positions[i].y - positions[j].y)/(magnitude));
                        botSum += w_ij;
                    }
                }

                double newX = topSumX/botSum;
                double newY = topSumY/botSum;
                positions[i] = new Vector2(newX, newY);
            }

            double stress = GraphIO.CalculateStress(d, positions, n);
            yield return stress;
            if ((prevStress - stress) / prevStress < eps)
                yield break;
            prevStress = stress;
        }
    }


    // weight = w_ij
    public static void WeightLaplacian(int[,] d, double[,] result, int n) {
        CreateLaplacian((i, j) => -1.0 / (d[i,j] * d[i,j]), result, n);
    }

    // delta = w_ij*d_ij^2
    public static void PositionLaplacian(double[,] deltas, Vector2[] positions, double[,] result, int n) {
        CreateLaplacian((i, j) => -deltas[i, j] * (1.0 / (positions[i] - positions[j]).Magnitude()), result, n);
    }

    public static void CreateLaplacian(Func<int,int,double> Entry, double[,] result, int n) {
        if (result.GetLength(0) != n || result.GetLength(1) != n) {
            throw new System.InvalidOperationException("wrong size of result matrix");
        }
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (i != j) {
                    result[i,j] = Entry(i,j);
                }
            }
            double sum=0;
            for (int k=0; k<n; k++) {
                if (k != i) {
                    sum -= result[i,k];
                }
            }
            result[i,i] = sum;
        }
    }

    // multiplies a and b and places result in c, skipping the first row for majorization
    public static void Multiply_x(double[,] a, Vector2[] b, double[] c) {
        int n = b.Length;
        if (a.GetLength(0) != n || a.GetLength(1) != n || c.Length != n-1) {
            throw new System.InvalidOperationException("matrix and vector have different lengths");
        }
        for (int i=1; i<n; i++) {
            double sum = 0;
            for (int j=0; j<n; j++) {
                sum += b[j].x * a[i,j];
            }
            c[i-1] = sum;
        }
    }
    // the same, but for the y axis
    public static void Multiply_y(double[,] a, Vector2[] b, double[] c) {
        int n = b.Length;
        if (a.GetLength(0) != n || a.GetLength(1) != n || c.Length != n-1) {
            throw new System.InvalidOperationException("matrix and vector have different lengths");
        }
        for (int i=1; i<n; i++) {
            double sum = 0;
            for (int j=0; j<n; j++) {
                sum += b[j].y * a[i,j];
            }
            c[i-1] = sum;
        }
    }
    // ^the above could use delegates, but they're slooow
}
