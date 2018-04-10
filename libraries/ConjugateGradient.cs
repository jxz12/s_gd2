using System;

public static class ConjugateGradient {
    // solves Ax = b, from initial guess x, and overwrites it.
    // r, p, Ap are temporary variables, but should be supplied so that they can be reused over iterations
    public static void Cg(
        double[,] a, double[] x, double[] b,
        double[] r, double[] p, double[] Ap, double tolerance=.1f, int maxIterations=10
     ) {
        int n = x.GetLength(0);
        if (a.GetLength(0) != n || a.GetLength(1) != n || b.GetLength(0) != n || r.GetLength(0) != n || p.GetLength(0) != n) {
            throw new System.ArgumentException("wrong dimensions for multiplication");
        }
        double rsold, rsnew, alpha;

        Multiply(a,x,r);
        Minus(b,r,r);
        Assign(p,r);
        rsold = Multiply(r,r);

        for (int i=0; i<maxIterations; i++) {
            Multiply(a,p,Ap);
            alpha = rsold / Multiply(p,Ap);
            Multiply(p,alpha,p);
            Plus(x,p,x);
            Multiply(Ap,alpha,Ap);
            Minus(r,Ap,r);
            rsnew = Multiply(r,r);
            // if (Math.Sqrt(rsnew) < tolerance) {
            if (rsnew < tolerance) {
                // Debug.Log(i);
                break;
            }
            Multiply(p,rsnew/rsold/alpha,p);
            Plus(r,p,p);
            rsold = rsnew;
        }
    }

// function [x] = conjgrad(A, b, x)
//     r = b - A * x;
//     p = r;
//     rsold = r' * r;

//     for i = 1:length(b)
//         Ap = A * p;
//         alpha = rsold / (p' * Ap);
//         x = x + alpha * p;
//         r = r - alpha * Ap;
//         rsnew = r' * r;
//         if sqrt(rsnew) < 1e-10
//               break;
//         end
//         p = r + (rsnew / rsold) * p;
//         rsold = rsnew;
//     end
// end

    static void Multiply(double[,] a, double[] b, double[] ans) {
        int n = ans.GetLength(0);
        // if (a.GetLength(0) != n || a.GetLength(1) != n || b.GetLength(0) != n) {
        //     throw new System.ArgumentException("wrong dimensions for multiplication");
        // }
        for (int i=0; i<n; i++) {
            double temp=0;
            for (int j=0; j<n; j++) {
                temp += a[i,j] * b[j];
            }
            ans[i] = temp;
        }
    }
    static void Multiply(double[] a, double b, double[] ans) {
        int n = a.GetLength(0);
        for (int i=0; i<n; i++) {
            ans[i] = a[i] * b;
        }
    }
    static double Multiply(double[] a, double[] b) {
        int n = a.GetLength(0);
        double temp=0;
        for (int i=0; i<n; i++) {
            temp += a[i] * b[i];
        }
        return temp;
    }
    static void Plus(double[] a, double[] b, double[] ans) {
        int n = a.GetLength(0);
        for (int i=0; i<n; i++) {
            ans[i] = a[i] + b[i];
        }
    }
    static void Minus(double[] a, double[] b, double[] ans) {
        int n = a.GetLength(0);
        for (int i=0; i<n; i++) {
            ans[i] = a[i] - b[i];
        }
    }
    static void Assign(double[] a, double[] b) {
        int n = a.GetLength(0);
        for (int i=0; i<n; i++) {
            a[i] = b[i];
        }
    }
}
