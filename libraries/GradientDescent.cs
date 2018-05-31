using System;
using System.Collections.Generic;
using GraphStuff;

public static class GradientDescent
{
    static Vector2 Gradient(Vector2 x_i, Vector2 x_j, double d_ij)
    {
        Vector2 ij = x_i - x_j;
        double mag = ij.Magnitude();
        double w_ij = 1.0 / (d_ij * d_ij);
        
        Vector2 m = 2 * w_ij * (d_ij - mag) * (ij / mag);
        return m;
    }
    public static IEnumerable<double> Full(int[,] d, Vector2[] X, double eta = .1, int numIterations = 15)
    {
        int n = X.Length;

        // relax
        for (int k=0; k<numIterations; k++) {
            var gradients = new Vector2[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    Vector2 gradient = Gradient(X[i], X[j], d[i, j]);

                    gradients[i] += gradient;
                    gradients[j] -= gradient;
                }
            }
            for (int i = 0; i < n; i++)
            {
                X[i] += eta * gradients[i];
            }
                
            double stress = GraphIO.CalculateStress(d, X, n);
            yield return stress;
        }
    }

    public static IEnumerable<double> Momentum(int[,] d, Vector2[] X, double eta = .1, double beta = .9, int numIterations = 15)
    {
        int n = X.Length;

        // relax
        for (int k=0; k<numIterations; k++) {
            var gradients = new Vector2[n];
            var momentums = new Vector2[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    Vector2 gradient = Gradient(X[i], X[j], d[i, j]);

                    gradients[i] += gradient;
                    gradients[j] -= gradient;
                }
            }
            for (int i = 0; i < n; i++)
            {
                momentums[i] = momentums[i] * beta + gradients[i] * eta;
                X[i] += momentums[i];
            }

            double stress = GraphIO.CalculateStress(d, X, n);
            yield return stress;
        }
    }
}
