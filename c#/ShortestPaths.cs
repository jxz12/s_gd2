using SparseCollections;
using System;
using System.Collections.Generic;
using System.Linq;

public static class ShortestPaths {
    public static int[,] Bacon(SparseMatrix<bool> adjMatrix, int n) {
        // calculate the distance using Kevin Bacon number (shortest path)
        var distMatrix = new int[n,n];
        foreach (int sourceIdx in adjMatrix.GetRowIndices()) {
            // keep a dictionary of previously visited nodes to backtrack the shortest paths to get there
            var prevVisits = new Dictionary<int, int>();

            // create a queue for the BFS
            var currentSearchNodes = new Queue<int>();
            currentSearchNodes.Enqueue(sourceIdx);

            // then take breadth steps from the node in order to fill up previousVisitedFlags
            while (currentSearchNodes.Count > 0) {
                int currentNode = currentSearchNodes.Dequeue();
                // try hopping to each connected node
                foreach (int nextNode in adjMatrix.GetColumnIndicesInRow(currentNode)) {
                    // if we havent seen it before
                    if (prevVisits.ContainsKey(nextNode) == false) {
                        // add it to the queue to explore from there
                        currentSearchNodes.Enqueue(nextNode);
                        // and set its flag to the node we hopped from
                        prevVisits[nextNode] = currentNode;
                    }
                }
            }

            // backtrack for each node to get the distance
            foreach (var kvp in prevVisits) {
                int destinationIdx = kvp.Key;
                int current = kvp.Value;
                int distance = 1;
                while (current != sourceIdx) {
                    current = prevVisits[current];
                    distance++;
                }
                distMatrix[sourceIdx, destinationIdx] = distance;
            }
        }
        return distMatrix;
    }


    // choose landmark points by drawing from a probability density of the MaxMin strategy from
    // "Sparse multidimensional scaling using landmark points, Vin de Silva & Joshua B. Tenenbaum, 2004"
    // and returns the list of closest vertices for each pivot
    public static Dictionary<int, List<int>> SparseBacon(SparseMatrix<bool> adjMatrix, SparseMatrix<int> distMatrix, int numSeeds, int numPivots) {

        // P(1:N) <- randperm(N)
        var allIdxs = adjMatrix.GetRowIndices().ToArray();
        FYShuffle(allIdxs);
        var pivots = new int[numPivots];

        // l(1:s) <- P(1:s)
        int n = allIdxs.Length;
        for (int i=0; i<numSeeds; i++) {
            int pivot = allIdxs[i];
            pivots[i] = pivot;
            SingleBacon(adjMatrix, distMatrix, pivot); 
        }

        // for (j = 1:N)
        int[,] mins = new int[n,2];
        for (int j=0; j<n; j++) {
            int min = int.MaxValue, argmin = -1;
            for (int i=0; i<numSeeds; i++) {
                int temp = distMatrix[pivots[i],j];
                if (temp < min) {
                    min = temp;
                    argmin = pivots[i];
                }
            }
            mins[j,0] = min;
            mins[j,1] = argmin;
        }

        int[] cumulativeProb = new int[n];
        Random rnd = new Random();
        // for (i = s+1:n)
        for (int i=numSeeds; i<numPivots; i++) {
            // int argmax = 0;
            // for (int k=1; k<n; k++) {
            //     if (mins[k,0] > mins[argmax,0]) argmax = k;
            // }
            // pivots[i] = argmax;
            // SingleBacon(adjMatrix, distMatrix, argmax);


            // use cumulative probability to select the next pivot
            int totalProb = 0;
            for (int k=0; k<n; k++) {
                totalProb += mins[k,0];
                // totalProb += mins[k,0] * adjMatrix.GetRowDataCount(k);
                // totalProb += adjMatrix.GetRowDataCount(k);
                cumulativeProb[k] = totalProb;
            }
            int sample = rnd.Next(0,totalProb);
            for (int k=0; k<n; k++) {
                if (sample < cumulativeProb[k]) {
                    pivots[i] = k;
                    SingleBacon(adjMatrix, distMatrix,k);
                    break;
                }
            }

            for (int j=0; j<n; j++) {
                int temp = distMatrix[pivots[i],j];
                if (temp < mins[j,0]) {
                    mins[j,0] = temp;
                    mins[j,1] = pivots[i];
                }
            }
        }

        // work out the closest pivot for each other vertex
        var closestVertices = new Dictionary<int, List<int>>();
        foreach (int i in pivots) {
            closestVertices[i] = new List<int>();
            mins[i,0] = 0;
            mins[i,1] = i;
        }
        for (int i=0; i<mins.GetLength(0); i++) {
            int closestPivot = mins[i,1];

            closestVertices[closestPivot].Add(i);
        }
        // foreach (var kvp in closestVertices) {
        //     for (int i=0; i<kvp.Value.Count; i++) print(kvp.Key + " "+kvp.Value[i] +" "+ distMatrix[kvp.Key,kvp.Value[i]]);
        // }
        return closestVertices;
    }

    public static void SingleBacon(SparseMatrix<bool> adjMatrix, SparseMatrix<int> distMatrix, int sourceIdx) {

        // keep a dictionary of previously visited nodes to backtrack the shortest paths to get there
        var prevVisits = new Dictionary<int, int>();

        // create a queue for the BFS
        var currentSearchNodes = new Queue<int>();
        currentSearchNodes.Enqueue(sourceIdx);

        // then take breadth steps from the node in order to fill up previousVisitedFlags
        while (currentSearchNodes.Count > 0) {
            int currentNode = currentSearchNodes.Dequeue();
            // try hopping to each connected node
            foreach (int nextNode in adjMatrix.GetColumnIndicesInRow(currentNode)) {
                // if we havent seen it before
                if (prevVisits.ContainsKey(nextNode) == false) {
                    // add it to the queue to explore from there
                    currentSearchNodes.Enqueue(nextNode);
                    // and set its flag to the node we hopped from
                    prevVisits[nextNode] = currentNode;
                }
            }
        }

        // backtrack for each node to get the distance
        foreach (var kvp in prevVisits) {
            int destinationIdx = kvp.Key;
            int current = kvp.Value;
            int distance = 1;
            while (current != sourceIdx) {
                current = prevVisits[current];
                distance++;
                // if (distance > vertices.Count) throw new InvalidOperationException("shortest path is bigger than longest possible path somehow lol");
            }
            // ignore distances of 1 because we have them already, and make the matrix symmetric.
            distMatrix[sourceIdx,destinationIdx] = distMatrix[destinationIdx,sourceIdx] = distance;
        }
    }

    static void FYShuffle(int[] l) {
        Random rnd = new Random();
        int n = l.Length;
        for (int i=0; i<n; i++) {
            int j = rnd.Next(i,n);
            int temp = l[i];
            l[i] = l[j];
            l[j] = temp;
        }
    }

}
