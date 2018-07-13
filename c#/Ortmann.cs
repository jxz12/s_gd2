using System;
using System.Diagnostics;
using System.Collections.Generic;

using SparseCollections;
using GraphStuff;

public class Ortmann {

    public static void Main(string[] args) {
        string name = args[0];
        SparseMatrix<bool> adjMatrix = GraphIO.Read("data/" + name + ".txt");

        int numPivots = int.Parse(args[1]);
        var positions = SparseSGD(adjMatrix, numPivots, 15);
        GraphIO.Write("svg/ortmann/" + name + "_" + numPivots + ".svg", positions, adjMatrix, 10, 50);

        //int n = adjMatrix.MaxIdx()+1;
        //int[,] d_full = ShortestPaths.Bacon(adjMatrix, n);
        //int[] nums = new int[] { 10, 50, 100, 200 };
        //foreach (int numPivots in nums)
        //{
        //    double bestStress = double.MaxValue;
        //    for (int i=0; i<25; i++)
        //    {
        //        var positions = SparseMajorization(adjMatrix, numPivots, 100);
        //        Console.Error.WriteLine("finished run " + i);
        //        var stress = GraphIO.CalculateStress(d_full, positions, n);

        //        Console.Error.WriteLine("stress=" + stress);
        //        Console.WriteLine(numPivots + " " + stress);

        //        if (stress < bestStress)
        //            GraphIO.Write("svg/ortmann/" + name + "_" + numPivots + ".svg", positions, adjMatrix, 10, 50);

        //    }
        //    Console.WriteLine();
        //}
        //Console.ReadLine();
    }

    public static Vector2[] SparseSGD(SparseMatrix<bool> adjMatrix, int numPivots, int numIter)
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
                d[i, j] = d[j, i] = 1;
                w[i, j] = w[j, i] = 1.0;
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
            foreach (int regionMember in regions[pivot])
            {
                if (regionMember != pivot)
                {
                    int distance = d[pivot, regionMember];
                    maxDistance = Math.Max(distance, maxDistance);
                    if (!regionDistances.ContainsKey(distance)) regionDistances[distance] = 1;
                    else regionDistances[distance] += 1;
                }
            }
            var cumulDistances = new List<int>();
            cumulDistances.Add(1);
            for (int i = 1; i <= maxDistance; i++)
            {
                cumulDistances.Add(cumulDistances[i - 1] + regionDistances[i]);
            }
            var terms = new List<int>();
            // for every other vertex
            for (int thisVertex = 0; thisVertex < n; thisVertex++)
            {
                int distance = d[pivot, thisVertex];

                // cut out paths to itself and its direct neighbours
                if (thisVertex != pivot && distance > 1)
                {
                    int s = regions[pivot].Count;
                    if (distance / 2 < cumulDistances.Count - 1)
                        s = cumulDistances[distance / 2];

                    double weight = 1f / (distance * distance);
                    //double w_pi = weight * numPivots / terms.Count;
                    double w_ip = weight * s;

                    d[pivot, thisVertex] = d[thisVertex, pivot] = distance;
                    w[thisVertex, pivot] = w_ip;
                }
            }
        }

        foreach (var ij in d.IndexPairs)
        {
            int i = ij.Item1, j = ij.Item2;
            if (i < j)
                constraints.Add(ij);
        }
        Console.Error.WriteLine("constraints="+constraints.Count);

        // initialise positions
        var positions = new Vector2[n];
        var rnd = new Random();
        for (int i=0; i<n; i++) {
            positions[i] = new Vector2(rnd.NextDouble(), rnd.NextDouble());
        }

		
        // calculate annealing
        double wMax = 0, wMin = double.MaxValue;
        foreach (var w_ij in w)
        {
            wMax = w_ij > wMax ? w_ij : wMax;
            wMin = w_ij < wMin ? w_ij : wMin;
        }
        double cMax = 1 / wMin;
        double cMin = 0.1 / wMax;

        double lambda = Math.Log(cMin / cMax) / (numIter - 1);
        Console.Error.WriteLine("cMax=" + cMax + " cMin=" + cMin + " rate=" + lambda);

        // SGD
        for (int k=0; k<numIter+5; k++) {
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

    public static Vector2[] SparseMajorization(SparseMatrix<bool> adjMatrix, int numPivots, int numIter)
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

        foreach (var ij in adjMatrix.IndexPairs) {
            int i=ij.Item1, j=ij.Item2;
            if (i>j) {
                d[i, j] = d[j, i] = 1;
                w[i, j] = w[j, i] = 1.0;
            }
        }

        ///////////////////
        // sparse stress //
        ///////////////////

        // for every pivot
        foreach (int pivot in regions.Keys) {

            // calculate the number of vertices in the region at each distance
            var regionDistances = new Dictionary<int, int>();
            int maxDistance = 0;
            foreach (int regionMember in regions[pivot])
            {
                if (regionMember != pivot)
                {
                    int distance = d[pivot, regionMember];
                    maxDistance = Math.Max(distance, maxDistance);
                    if (!regionDistances.ContainsKey(distance)) regionDistances[distance] = 1;
                    else regionDistances[distance] += 1;
                }
            }
            var cumulDistances = new List<int>();
            cumulDistances.Add(1);
            for (int i = 1; i <= maxDistance; i++)
            {
                cumulDistances.Add(cumulDistances[i - 1] + regionDistances[i]);
            }
            var terms = new List<int>();
            // for every other vertex
            for (int thisVertex = 0; thisVertex < n; thisVertex++)
            {
                int distance = d[pivot, thisVertex];

                // cut out paths to itself and its direct neighbours
                if (thisVertex != pivot && distance > 1)
                {
                    int s = regions[pivot].Count;
                    if (distance / 2 < cumulDistances.Count - 1)
                        s = cumulDistances[distance / 2];

                    double weight = 1f / (distance * distance);
                    //double w_pi = weight * numPivots / terms.Count;
                    double w_ip = weight * s;

                    d[pivot, thisVertex] = d[thisVertex, pivot] = distance;
                    w[thisVertex, pivot] = w_ip;
                }
            }
        }

        // initialise positions
        var X = new Vector2[n];
        var rnd = new Random();
        for (int i=0; i<n; i++) {
            X[i] = new Vector2(rnd.NextDouble(), rnd.NextDouble());
        }

        var wBot = new Dictionary<int, double>();
        foreach (int i in adjMatrix.GetRowIndices())
        {
            double bot = 0;
            foreach (int j in adjMatrix.GetColumnIndicesInRow(i))
                if (i != j)
                    bot += w[i, j];
            foreach (int p in regions.Keys)
                if (p!= i && adjMatrix[i,p] == false)
                    bot += w[i, p];
            wBot[i] = bot;
        }


        for (int k = 0; k < numIter; k++)
        {
            foreach (int i in adjMatrix.GetRowIndices())
            {
                double xTop = 0, yTop = 0;
                // do neighbours
                foreach (int j in adjMatrix.GetColumnIndicesInRow(i))
                {
                    if (i != j)
                    {
                        double dist = (X[i] - X[j]).Magnitude();
                        xTop += w[i, j] * (X[j].x + (d[i, j] * (X[i].x - X[j].x)) / dist);
                        yTop += w[i, j] * (X[j].y + (d[i, j] * (X[i].y - X[j].y)) / dist);
                    }
                }
                // do pivots
                foreach (int p in regions.Keys)
                {
                    if (p != i && adjMatrix[i, p] == false)
                    {
                        double w_ip = w[i, p];
                        double dist = (X[i] - X[p]).Magnitude();
                        xTop += w[i, p] * (X[p].x + (d[i, p] * (X[i].x - X[p].x)) / dist);
                        yTop += w[i, p] * (X[p].y + (d[i, p] * (X[i].y - X[p].y)) / dist);
                    }
                }
                X[i] = new Vector2(xTop / wBot[i], yTop / wBot[i]);
            }
            Console.Error.Write(k);
        }

        Console.Error.WriteLine("done!");

        Console.Error.WriteLine("time="+sw.ElapsedMilliseconds+"ms");
        return X;
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
