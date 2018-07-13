using System;
using System.Collections.Generic;
using System.Linq;
using GraphStuff;

public static class Sgd2
{
    const double limit = 1;

    //[System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    static double Satisfy(ref Vector2 x_i, ref Vector2 x_j, double d, double c)
    {
        Vector2 pq = x_i - x_j;
        // mag is |p-q|
        double mag = pq.Magnitude();
        // r is minimum distance each vertex has to move to satisfy the constraint
        double r = (mag - d) / 2;

        double w_ij = 1.0 / (d * d);
        // weight by a maximum of 2
        double wc = w_ij * c;
        if (wc > limit)
            wc = limit;
        r = wc * r;

        Vector2 m = pq * r / mag;
        x_i -= m;
        x_j += m;

        return r;
    }
    public static IEnumerable<double> Schedule(int[,] d, double cMin, int numIter)
    {
        int maxIter = numIter;
        double cMax = EtaMax(d);
        double[] eta = new double[maxIter];
        var f = ExpDecay(cMax, cMin, numIter);
        //var f = Reciprocal(cMax, cMin, numIter);
        //var f = Sqrt(cMax, cMin, numIter);
        for (int t = 0; t < maxIter; t++)
            eta[t] = f(t);
        return eta;
    }

    static int[,] CreatePairs(int n)
    {
        int nn = (n*(n-1))/2;
        var indexPairs = new int[nn, 2];
        int idx = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < i; j++)
            {
                indexPairs[idx, 0] = i;
                indexPairs[idx, 1] = j;
                idx++;
            }
        }
        return indexPairs;
    }
	
    public static IEnumerable<double> Full(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        var pairs = CreatePairs(n);
        var rnd = new Random();

        // relax
        foreach (double c in eta)
        {
            GraphIO.FYShuffle2(pairs, rnd);
            for (int ij = 0; ij < nn; ij++)
            {
                int i = pairs[ij, 0], j = pairs[ij, 1];
                Satisfy(ref positions[i], ref positions[j], d[i, j], c);
            }
            double stress = GraphIO.CalculateStress(d, positions, n);
            yield return stress;
        }
    }

    public static IEnumerable<double> Adaptive(int[,] d, Vector2[] positions, double eps=0.001, int maxIter=1000) {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        var pairs = CreatePairs(n);
        var rnd = new Random();

        double prevStress = GraphIO.CalculateStress(d, positions, n);
        int t = 0, numSteps = 15;
        double cMin = 0.1;
        Func<int, double> Eta = ExpDecay(EtaMax(d), cMin, numSteps);
        double c = Eta(t);

        // relax
        for (int it=0; it<maxIter; it++) {
            GraphIO.FYShuffle2(pairs, rnd);

            var prevPositions = (Vector2[])positions.Clone();
            double maxMovement = 0;
            for (int ij=0; ij<nn; ij++)
            {
                int i = pairs[ij, 0], j = pairs[ij, 1];
                double r = Satisfy(ref positions[i], ref positions[j], d[i, j], c);

                maxMovement = Math.Max(maxMovement, r*r);
            }

            double stress = GraphIO.CalculateStress(d, positions, n);
            yield return stress;

            //Console.Error.WriteLine("                      " + c + " " + maxMovement);
            if (maxMovement < eps)
                yield break;

            if (prevStress - stress < 0)
                c = Eta(t++);
			
            prevStress = stress;
        }

    }


    public static IEnumerable<double> Ordered(int[,] d, Vector2[] positions, IEnumerable<double> eta)
    {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        var pairs = CreatePairs(n);

        var indices = new int[nn];
        for (int i = 0; i < nn; i++)
            indices[i] = i;
        GraphIO.FYShuffle(indices, new Random());

        var result = indices.OrderByDescending(ij => d[pairs[ij, 0], pairs[ij, 1]]);

        // relax
        foreach (double c in eta) {
            foreach (int ij in result)
            {
                int i = pairs[ij, 0], j = pairs[ij, 1];
                Satisfy(ref positions[i], ref positions[j], d[i,j], c);
            }
            yield return GraphIO.CalculateStress(d, positions, n);
        }

    }
    public static IEnumerable<double> Once(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        var pairs = CreatePairs(n);
        var rnd = new Random();
        GraphIO.FYShuffle2(pairs, rnd);

        // relax
        foreach (double c in eta) { 
            for (int ij=0; ij<nn; ij++)
            {
                int i = pairs[ij, 0], j = pairs[ij, 1];
                Satisfy(ref positions[i], ref positions[j], d[i,j], c);
            }
            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }
    public static IEnumerable<double> NoRand(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;

        // relax
        foreach (double c in eta) {
            for (int i=0; i<n; i++)
                for (int j=0; j<i; j++)
                    Satisfy(ref positions[i], ref positions[j], d[i,j], c);

            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }
    public static IEnumerable<double> Indices(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;

        var indices = new int[n];
        for (int i=0; i<n; i++)
            indices[i] = i;
        var rnd = new Random();
        //GraphIO.FYShuffle(indices, rnd);

        // relax
        foreach (double c in eta) { 
            GraphIO.FYShuffle(indices, rnd);
            for (int i=0; i<n; i++)
                for (int j=0; j<i; j++)
                    Satisfy(ref positions[indices[i]], ref positions[indices[j]], d[indices[i],indices[j]], c);

            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }
    public static IEnumerable<double> EveryPair(int[,] d, Vector2[] positions, IEnumerable<double> eta)
    {
        int n = positions.Length;
        int nn = (n*(n-1))/2;
        var rnd = new Random();

        //var indices = CreatePairs(n);

        // relax
        foreach (double c in eta) { 
            for (int ij=0; ij<nn; ij++) {
                int idx = rnd.Next(nn);
                int i = (int)((1 + Math.Sqrt(8 * idx + 1)) / 2);
                int j = idx - (i * (i - 1)) / 2;
                //int i = indices[idx, 0], j = indices[idx, 1];

                Satisfy(ref positions[i], ref positions[j], d[i,j], c);
            }
            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }

    public static IEnumerable<double> Alternating(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;
        int nn = (n*(n-1))/2;


        int[,] pairs1 = CreatePairs(n), pairs2 = CreatePairs(n);
        var rnd = new Random();
        GraphIO.FYShuffle2(pairs1, rnd);
        GraphIO.FYShuffle2(pairs2, rnd);

        // relax
        int k = 0;
        foreach (double c in eta) {
            int[,] pairs = k++ % 2 == 0 ? pairs1 : pairs2;
            for (int ij=0; ij<nn; ij++)
            {
                int i = pairs[ij, 0], j = pairs[ij, 1];
                Satisfy(ref positions[i], ref positions[j], d[i,j], c);
            }

            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }

    public static IEnumerable<double> Modulo(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        // relax

        int prime = 646957;
        int[] primitives = new int[] {5, 6, 7, 17, 18, 20, 21, 24, 26, 28, 45, 46, 50, 53, 55, 58, 66, 68, 72, 73};

        int k = 0;
        foreach (double c in eta) { 
            int primitive = primitives[k++];
            int modulo = 1;
            for (int ij=0; ij<prime; ij++)
            {
                modulo = (modulo*primitive) % prime;
                if (modulo > nn)
                    continue;

                int idx = modulo - 1;
                int i = (int)( (1+Math.Sqrt(8*idx+1))/2 );
                int j = idx - (i*(i-1))/2;

                Satisfy(ref positions[i], ref positions[j], d[i,j], c);
            }
            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }
    public static IEnumerable<double> Sliced(int[,] d, Vector2[] positions, IEnumerable<double> eta)
    {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        int numSlices = 2000;
        int blockSize = (nn + numSlices - 1) / numSlices;
        var slices = new int[numSlices][];

        int ij = 0;
        for (int i=0; i<numSlices-1; i++)
        {
            slices[i] = new int[blockSize];
            for (int j = 0; j < blockSize; j++)
                slices[i][j] = ij++;
        }
        int lastBlock = nn - ij;
        slices[numSlices - 1] = new int[lastBlock];
        for (int j = 0; j < lastBlock; j++)
            slices[numSlices - 1][j] = ij++;

        //for (int s = 0; s < numSlices; s++)
        //{
        //    Console.ReadLine();
        //    for (int idx = 0; idx < slices[s].Length; idx++)
        //    {
        //        Console.Error.WriteLine(slices[s][idx]);
        //    }
        //}


        var rnd = new Random();
        foreach (double c in eta)
        {
            GraphIO.FYShuffle(slices, rnd);
            for (int s=0; s<numSlices; s++)
            {
                GraphIO.FYShuffle(slices[s], rnd);
                for (int b=0; b<slices[s].Length; b++)
                {
                    int idx = slices[s][b];
                    int i = (int)( (1+Math.Sqrt(8*idx+1))/2 );
                    int j = idx - (i*(i-1))/2;
                    
                    Satisfy(ref positions[i], ref positions[j], d[i,j], c);
                }
            }
            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }
    private class Mom
    {
        public int i;
        public int j;
        public Vector2 m;
        public Mom(int i, int j, Vector2 m)
        {
            this.i = i;
            this.j = j;
            this.m = m;
        }
    }
    public static IEnumerable<double> Momentum(int[,] d, Vector2[] positions, IEnumerable<double> eta) {
        int n = positions.Length;
        int nn = (n*(n-1))/2;

        var pairs = CreatePairs(n);
        var rnd = new Random();
        //Func<int, double> Cooling = ExpDecay(CMax(d), 1, numIterations);
        //Func<int, double> Cooling = Reciprocal(d, cMin, numIterations);
        //Func<int, double> Cooling = Cosine(d, cMin, numIterations);

        double beta = .5;
        //int mo = 3;
        //var momentum = new Queue<Mom>();
        //for (int i = 0; i < mo; i++)
        //    momentum.Enqueue(new Mom(0, 0, new Vector2(0, 0)));

        var momentum = new Vector2[n];

        // relax
        foreach (double c in eta) {
            GraphIO.FYShuffle2(pairs, rnd);

            for (int ij=0; ij<nn; ij++)
            {
                int i = pairs[ij, 0], j = pairs[ij, 1];
                Vector2 pq = positions[i] - positions[j];
                // mag is |p-q|
                double mag = pq.Magnitude();
                // r is minimum distance each vertex has to move to satisfy the constraint
                double r = (mag-d[i,j]) / 2;

                double w_ij = 1.0 / (d[i,j] * d[i,j]);
                // weight by a maximum of 2
                double wc = w_ij * c;
                if (wc > limit)
                    wc = limit;
                r = wc * r;

                Vector2 m = pq * r / mag;
                positions[i] -= m;
                positions[j] += m;

                momentum[i] -= m;
                momentum[j] += m;

                //momentum.Enqueue(new Mom(i, j, m));
                //momentum.Dequeue();

                //foreach (var a in momentum)
                //{
                //    positions[a.i] -= a.m;
                //    positions[a.j] += a.m;
                //    a.m *= beta;
                //}
            }
            for (int i=0; i<n; i++)
            {
                momentum[i] *= beta;
                positions[i] += momentum[i];
            }
            yield return GraphIO.CalculateStress(d, positions, n);
        }
    }

    static void OverSatisfy(ref Vector2 x_i, ref Vector2 x_j, double d_ij, double c)
    {
        Vector2 pq = x_i - x_j;
        // mag is |p-q|
        double mag = pq.Magnitude();
        // r is minimum distance each vertex has to move to satisfy the constraint
        double om = Math.Min(1.0 / (d_ij * d_ij) * c, 1);
        double r = om*(mag-d_ij)/2;

        Vector2 m = pq * r / mag;
        x_i -= 2*m;
        x_j += 2*m;
    }
    static double EtaMax(int[,] d)
    {
        int n = d.GetLength(0), m = d.GetLength(1);
        int max = 0;
        for (int i=0; i<n; i++)
            for (int j=0; j<m; j++)
                max = d[i,j] > max ? d[i,j] : max;

        double cMax = limit * max * max;
        return cMax;
    }
    public static Func<int, double> ExpDecay(double cMax, double cMin, int numIter)
    {
        double lambda = Math.Log(cMin / cMax) / (numIter - 1);

        Console.Error.WriteLine("cMax=" + cMax + " rate=" + lambda);
        return k=> cMax * Math.Exp(lambda * k);
    }
    public static Func<int, double> Reciprocal(double cMax, double cMin, int numIter)
    {
        double a = cMax;
        double b = (a - cMin) / (cMin * (numIter-1));

        Console.Error.WriteLine("a=" + a + " b=" + b);
        return k =>a / (1 + k * b);
    }
    public static Func<int, double> Cosine(double cMax, double cMin, int numIter)
    {
        return k => cMin + .5*(cMax - cMin) * (1 + Math.Cos(((double)k / (numIter-1)) * Math.PI));
    }
    public static Func<int, double> Quadratic(double cMax, double cMin, int numIter)
    {
        double b = (-Math.Sqrt(cMin / cMax) + 1) / (numIter-1);
        return k => cMax * Math.Pow(b * k - 1, 2);
    }
    public static Func<int, double> Linear(double cMax, double cMin, int numIter)
    {
        double b = (cMin - cMax) / (numIter - 1);
        return k => cMax + b * k;
    }
    public static Func<int, double> Sqrt(double cMax, double cMin, int numIter)
    {
        double b = (Math.Pow(cMax/cMin, 2) - 1) / (numIter - 1);
        return k => cMax / Math.Sqrt(1+b*k);
    }
    
}
