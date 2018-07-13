using SparseCollections;
using System;
using System.IO;
using System.Collections.Generic;

namespace GraphStuff {
    public struct Vector2 {
        public double x;
        public double y;
        public Vector2(double x, double y) {
            this.x = x;
            this.y = y;
        }
        public double Magnitude() {
            return Math.Sqrt(x*x + y*y);
        }
        public static Vector2 operator +(Vector2 v1, Vector2 v2) {
            return new Vector2(v1.x + v2.x, v1.y + v2.y);
        }
        public static Vector2 operator -(Vector2 v1, Vector2 v2) {
            return new Vector2(v1.x - v2.x, v1.y - v2.y);
        }
        public static Vector2 operator *(Vector2 v, double a) {
            return new Vector2(v.x * a, v.y * a);
        }
        public static Vector2 operator *(double a, Vector2 v) {
            return new Vector2(v.x * a, v.y * a);
        }
        public static Vector2 operator /(Vector2 v, double a) {
            return new Vector2(v.x / a, v.y / a);
        }
    }

    public static class GraphIO {
        public static SparseMatrix<bool> Read(string path) {

            var adjMatrix = new SparseMatrix<bool>();

            // run through once to squash indices down to contiguous integers numbered from 0
            var indices = new HashSet<int>();
            var mapping = new Dictionary<int,int>();
            int squashedIndex = 0;
            foreach (var ij in EnumerateIndices(path)) {
                int i = ij.Item1, j = ij.Item2;
                if (!indices.Contains(i)) {
                    indices.Add(i);
                    mapping[i] = squashedIndex;
                    squashedIndex++;
                }
                if (!indices.Contains(j)) {
                    indices.Add(j);
                    mapping[j] = squashedIndex;
                    squashedIndex++;
                }
            }
            // fill in the adjacency matrix
            foreach (var ij in EnumerateIndices(path)) {
                int i = ij.Item1, j = ij.Item2;
                int i2 = mapping[i], j2 = mapping[j];
                adjMatrix[i2,j2] = adjMatrix[j2,i2] = true; // make the matrix unweighted
            }
            return adjMatrix;
        }

        static IEnumerable<Tuple<int,int>> EnumerateIndices(string s) {
            using (var sreader = new StreamReader(s)) {
                string line = sreader.ReadLine();
                while (line != null) {
                    var numbers = line.Split(' ');
                    int i = int.Parse(numbers[0]);
                    int j = int.Parse(numbers[1]);
                    yield return Tuple.Create(i,j);
                    line = sreader.ReadLine();
                }
            }
        }

        // writes an svg file to stdout
        public static void Write(string path, Vector2[] positions, SparseMatrix<bool> adjMatrix, double scale=5, int border=50) {
            int n = positions.Length;
            for (int i=0; i<n; i++) {
                positions[i] *= scale;
            }
            // find min and max values
            double xMin, yMin, xMax, yMax;
            xMin = xMax = positions[0].x;
            yMin = yMax = positions[0].y;
            for (int i=1; i<n; i++) {
                if (positions[i].x < xMin) xMin = positions[i].x;
                if (positions[i].x > xMax) xMax = positions[i].x;
                if (positions[i].y < yMin) yMin = positions[i].y;
                if (positions[i].y > yMax) yMax = positions[i].y;
            }
            double xRange = xMax - xMin;
            double yRange = yMax - yMin;
            // scale and offset so that the whole thing fits within a box of 1000x1000 with a border of 50
            Vector2 offset = new Vector2(xMin-border, yMin-border);

            for (int i=0; i<n; i++) {
                positions[i] -= offset;
            }

            List<string> lines = new List<string>();
            lines.Add(String.Format("<svg width=\"{0:.0}\" height=\"{1:.0}\" xmlns=\"http://www.w3.org/2000/svg\">", xRange+2*border, yRange+2*border));
            lines.Add("<style type=\"text/css\">line{stroke:black;stroke-width:1;stroke-opacity:0.5;stroke-linecap:round;}</style>");
            foreach (var ij in adjMatrix.IndexPairs) { 
                int i=ij.Item1, j=ij.Item2;
                if (i > j) {
                    lines.Add(String.Format("<line x1=\"{0:.###}\" x2=\"{1:.###}\" y1=\"{2:.###}\" y2=\"{3:.###}\"/>", positions[i].x, positions[j].x, positions[i].y, positions[j].y));
                }
            }
            lines.Add("</svg>");
            File.WriteAllLines(path, lines);
        }

        public static void RandomiseIndices<T>(ref SparseMatrix<T> adjMatrix) {
            int n = adjMatrix.MaxIdx()+1;
            //Console.Error.WriteLine(n);
            int[] indices = new int[n];
            for (int i=0; i<n; i++)
                indices[i] = i;
            var rnd = new Random();
            FYShuffle(indices, rnd);
            var randomised = new SparseMatrix<T>();
            foreach (var ij in adjMatrix.IndexPairs) {
                randomised[indices[ij.Item1], indices[ij.Item2]] = adjMatrix[ij];
            }
            adjMatrix = randomised;
        }
        public static void RandomiseIndices(ref int[,] adjMatrix) {
            int n = adjMatrix.GetLength(0), m = adjMatrix.GetLength(1);
            //Console.Error.WriteLine(n);
            int[] indices = new int[n];
            for (int i=0; i<n; i++)
                indices[i] = i;
            var rnd = new Random();
            FYShuffle(indices, rnd);
            var randomised = new int[n,n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    randomised[indices[i], indices[j]] = adjMatrix[i, j];
            adjMatrix = randomised;
        }

        public static void FYShuffle<T>(T[] l, Random rnd) {
            int n = l.GetLength(0);
            T temp;
            int j;
            for (int i=0; i<n; i++) {
                j = rnd.Next(i,n);
                temp = l[i];
                l[i] = l[j];
                l[j] = temp;
            }
        }
        public static void FYShuffle2(int[,] l, Random rnd) {
            int n = l.GetLength(0);
            int temp0, temp1;
            int j;
            for (int i=0; i<n; i++) {
                j = rnd.Next(i,n);
                temp0 = l[i,0];
                temp1 = l[i,1];
                l[i,0] = l[j,0];
                l[i,1] = l[j,1];
                l[j,0] = temp0;
                l[j,1] = temp1;
            }
        }

        public static string MakeBtree(string outputPath, int depth) {
            var q = new Queue<int>();
            var sb = new System.Text.StringBuilder();

            q.Enqueue(0);

            int newIdx = 1;
            for (int i=0; i<depth; i++) {
                int current = q.Peek();
                int size = q.Count;
                for (int j=0; j<size; j++) {
                    current = q.Dequeue();
                    q.Enqueue(newIdx);
                    q.Enqueue(newIdx+1);
                    // sb.Append(current.ToString() + ":" + newIdx + "," + (newIdx+1) + "\n");
                    sb.Append(current+" "+newIdx+"\n"+current+" "+(newIdx+1)+"\n");
                    newIdx += 2;
                }
            }
            return sb.ToString();
        }

        public static double CalculateStress(int[,] d, Vector2[] positions, int n) {
            double stress = 0;
            for (int i=0; i<n; i++) {
                for (int j=0; j<i; j++) {
                    double d_ij = d[i,j];
                    stress += 1.0/(d_ij*d_ij) * Math.Pow((positions[i]-positions[j]).Magnitude() - d_ij, 2);
                }
            }
            return stress;
        }
    }
}
