using System;
using System.Diagnostics;
using System.Collections.Generic;

using SparseCollections;
using GraphStuff;

public class Ortmann {

    public static void Main(string[] args) {
        string name = args[0];
        SparseMatrix<bool> adjMatrix = GraphIO.Read("data/" + name + ".txt");
        int n = adjMatrix.MaxIdx()+1;
        int[,] d_full = ShortestPaths.Bacon(adjMatrix, n);

        int[] nums = new int[] { 10, 50, 100, 200 };
        foreach (int numPivots in nums)
        {
            double bestStress = double.MaxValue;
            for (int i=0; i<25; i++)
            {
                var positions = SparseStress(adjMatrix, numPivots);
                Console.Error.WriteLine("finished run " + i);
                var stress = GraphIO.CalculateStress(d_full, positions, n);

                Console.Error.WriteLine("stress=" + stress);
                Console.WriteLine(numPivots + " " + stress);

                if (stress < bestStress)
                    GraphIO.Write("svg/ortmann/" + name + "_" + numPivots + ".svg", positions, adjMatrix, 5, 50);

            }
            Console.WriteLine();
        }
        Console.ReadLine();
    }


    public static Vector2[] SparseStress(SparseMatrix<bool> adjMatrix, int numPivots)
    {
        var sw = new Stopwatch();
        sw.Start();

        int n = adjMatrix.MaxIdx()+1;
        Console.Error.WriteLine("vertices="+n);
        // create sparse distance matrix
        var d = new SparseMatrix<int>();
        Dictionary<int, List<int>> regions = ShortestPaths.SparseBacon(adjMatrix, d, 1, numPivots);

        var w = new SparseMatrix<double>();

        //////////////
        // 1-stress //
        //////////////

        var constraints = new List<Tuple<int,int>>();
        foreach (var ij in adjMatrix.IndexPairs) {
            int i=ij.Item1, j=ij.Item2;
            if (i>j) {
                constraints.Add(ij);
                d[i, j] = 1;
                w[i, j] = w[j, i] = 1;
            }
        }
        Console.Error.WriteLine("edges="+constraints.Count);

        ///////////////////
        // sparse stress //
        ///////////////////

        // for every pivot
        foreach (int pivot in regions.Keys) {
            // calculate the number of vertices in the region at each distance
            var regionDistances = new Dictionary<int, int>();
            int maxDistance = 0;
            foreach (int regionMember in regions[pivot]) {
                if (regionMember != pivot) {
                    int distance = d[pivot, regionMember];
                    maxDistance = Math.Max(distance, maxDistance);
                    if (!regionDistances.ContainsKey(distance)) regionDistances[distance] = 1;
                    else regionDistances[distance] += 1;
                }
            }
            var cumulDistances = new List<int>();
            cumulDistances.Add(1);
            for (int i=1; i<=maxDistance; i++) {
                cumulDistances.Add(cumulDistances[i-1] + regionDistances[i]);
            }
            // for every other vertex
            for (int thisVertex=0; thisVertex<n; thisVertex++) {
                int distance = d[pivot, thisVertex];

                // cut out paths to itself and its direct neighbours
                if (thisVertex != pivot && distance > 1) {
                    int s = regions[pivot].Count;
                    if (distance/2 < cumulDistances.Count-1)
                        s = cumulDistances[distance/2];

                    double weight = 1/Math.Pow(distance,2);
                    w[pivot, thisVertex] = weight;
                    w[thisVertex, pivot] = weight * s;

                    d[pivot, thisVertex] = d[thisVertex, pivot] = distance;
                    constraints.Add(Tuple.Create(pivot, thisVertex));
                }
            }
        }
        Console.Error.WriteLine("constraints="+constraints.Count);

        // initialise positions
        var positions = new Vector2[n];
        var rnd = new Random();
        for (int i=0; i<n; i++) {
            positions[i] = new Vector2(rnd.NextDouble(), rnd.NextDouble());
        }

        // relax
        //double c = 1000;
        //double rate = 0.75;
        //int iterations = 30;

        double wMax = 0, wMin = double.MaxValue;
        foreach (var ij in w.IndexPairs)
        {
            wMax = w[ij] > wMax ? w[ij] : wMax;
            wMin = w[ij] < wMin ? w[ij] : wMin;
        }
        double cMax = 1.0 / wMin;
        double cMin = 0.1 / wMax;
        int iterations = 15;
        double lambda = Math.Log(cMin / cMax) / (iterations - 1);
        Console.Error.WriteLine("cMax=" + cMax + " rate=" + lambda);

        for (int k=0; k<iterations+5; k++) {
            FYShuffle(constraints, rnd);
            double c = cMax * Math.Exp(lambda * k);

            foreach (var ij in constraints) {
                int i=ij.Item1, j=ij.Item2;

                Vector2 pq = positions[i] - positions[j];
                // mag is |p-q|
                double mag = pq.Magnitude();
                // r is minimum distance each vertex has to move to satisfy the constraint
                double r = (d[i,j] - mag) / 2;

                // the weighting for i and j are different
                double wc = w[i,j] * c;
                wc = Math.Min(wc, 1);
                double r_i = wc * r;

                wc = w[j,i] * c;
                wc = Math.Min(wc, 1);
                double r_j = wc * r;

                positions[i] += pq * (r_i/mag);
                positions[j] -= pq * (r_j/mag);
            }
            Console.Error.Write(k);
        }
        Console.Error.WriteLine("done!");

        Console.Error.WriteLine("time="+sw.ElapsedMilliseconds+"ms");
        return positions;
    }

    static void FYShuffle<T>(List<T> l, Random rnd)
    {
        int n = l.Count;
        for (int i = 0; i < n; i++)
        {
            int j = rnd.Next(i, n);
            T temp = l[i];
            l[i] = l[j];
            l[j] = temp;
        }
    }
}
