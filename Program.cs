using System;
using System.IO;
using System.Linq;
using System.Diagnostics;
using System.Collections.Generic;

using SparseCollections;
using GraphStuff;

public class Program {
    public static void Main(string[] args) {

		string graphDir = args[0];
        SparseMatrix<bool> adjMatrix = GraphIO.Read(graphDir);
			
        using (StreamWriter file = new StreamWriter(args[1]))
        {
            Vector2[] X = Stress((d, x) => Sgd2.Full(d, x, Sgd2.Schedule(d, 0.1, 15)), adjMatrix, file);
            //Vector2[] X = Stress((d, x) => Majorization.Chol(d, x), adjMatrix, file);
            //Vector2[] X = Stress((d, x) => Sgd2.Adaptive(d, x, 0.001), adjMatrix, file);
            //Vector2[] X = Stress((d, x) => GradientDescent.Full(d, x, .02, 100))

            GraphIO.Write(args[2], X, adjMatrix);
        }
    }

    public static Vector2[] Stress(Func<int[,], Vector2[], IEnumerable<double>> Layout, SparseMatrix<bool> adjMatrix, StreamWriter file)
    {
        int n = adjMatrix.MaxIdx() + 1;
        Console.Error.WriteLine("vertices="+n);
        Console.Error.WriteLine("edges="+adjMatrix.Count());

        var positions = new Vector2[n];
        var rnd = new Random();

        int[,] shortestPaths = ShortestPaths.Bacon(adjMatrix, n);

        var sw = new Stopwatch();

        //for (int k = 0; k < 25; k++)
        //{
        //    GC.Collect();

            // initialise positions
            positions[0] = new Vector2(0,0); // for majorization
            for (int i=1; i<n; i++)
                positions[i] = new Vector2(rnd.NextDouble(), rnd.NextDouble());
            
			sw.Restart();
			
            int iteration = 0;
            foreach (double stress in Layout(shortestPaths, positions))
            {
                double time = sw.ElapsedMilliseconds;
                string info = iteration + " " + stress + " " + time;
                Console.Error.WriteLine(info);
                //Console.WriteLine(info);
                file.WriteLine(info);
                iteration += 1;
            }

        //    Console.Error.WriteLine("finished run " + (k + 1) + "\n");
        //    Console.WriteLine();
        //    file.WriteLine();
        //}

        return positions;
    }

}
