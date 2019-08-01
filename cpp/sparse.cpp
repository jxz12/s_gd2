#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <limits>

#include <chrono>
#include <iostream>

// I like C#
template<typename T>
using Set = std::unordered_set<T>;
template<typename K, typename V>
using Dict = std::unordered_map<K, V>;
template<typename T>
using List = std::vector<T>;

using Graph = List<List<int>>;

struct term
{
    int i, j;
    double d, w_ij, w_ji;
    term(int i, int j, double d) : i(i), j(j), d(d), w_ij(0), w_ji(0) {}
};


Graph build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    std::vector<std::unordered_set<int>> undirected(n);
    Graph graph(n);

    for (int ij = 0; ij < m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw "i or j bigger than n";

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if edge not seen
        {
            undirected[i].insert(j);
            undirected[j].insert(i);
            graph[i].push_back(j);
            graph[j].push_back(i);
        }
    }
    return graph;
}

void maxmin_bfs(const Graph& graph, const int p, List<int>& mins, List<int>& argmins)
{
    int n = graph.size();
    std::queue<int> q;
    List<int> d(n, -1);

    q.push(p);
    d[p] = 0;
    while (!q.empty())
    {
        int current = q.front();
        q.pop();
        for (int next : graph[current])
        {
            if (d[next] == -1)
            {
                q.push(next);
                d[next] = d[current] + 1;
                if (d[next] < mins[next])
                {
                    mins[next] = d[next];
                    argmins[next] = p;
                }
            }
        }
    }
}

// returns closest pivot for each vertex
List<int> maxmin_random_sp(const Graph& graph, int n_pivots, int p0=0)
{
    int n = graph.size();
    List<int> pivots;
    pivots.push_back(p0);

    List<int> mins(n, std::numeric_limits<int>::max());
    List<int> argmins(n, -1);

    // first pivot
    mins[p0] = 0;
    argmins[p0] = p0;
    maxmin_bfs(graph, p0, mins, argmins);

    // remaining pivots
    while (pivots.size() < n_pivots)
    {
        int max = mins[0], argmax = 0;
        for (int i = 1; i < n; i++)
        {
            if (mins[i] > max)
            {
                max = mins[i];
                argmax = i;
            }
        }
        pivots.push_back(argmax);
        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs(graph, argmax, mins, argmins);
    }
    return argmins;
}

List<term> MSSP(const Graph& graph, const List<int> pivots)
{
    int n = pivots.size();

    // get list of pivots, set of region
    Dict<int, Set<int>> regions;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(pivots[i]) == regions.end())
        {
            regions[pivots[i]] = Set<int>();
        }
        regions[pivots[i]].insert(i);
    }

    Dict<int, Dict<int, term>> termsDict;
    for (auto const& region : regions)
    {
        // q contains next to visit
        std::queue<int> q;
        List<int> d(n, -1);

        int p = region.first;
        q.push(p);
        d[p] = 0;

        // q2 contains visited vertices, for s calculation
        std::queue<int> q2;
        int s = 0;
        q2.push(0);

        while (!q.empty())
        {
            int current = q.front();
            q.pop();

            for (int next : graph[current])
            {
                if (d[next] == -1)
                {
                    q.push(next);
                    d[next] = d[current] + 1;

                    // empty the second queue enough to calculate s
                    while (!q2.empty() && q2.front() <= d[next]/2)
                    {
                        q2.pop();
                        s += 1;
                    }
                    if (region.second.find(next) != region.second.end())
                    {
                        q2.push(d[next]);
                    }

                    //std::cerr << next << " " << d[next] << " " << s << std::endl;
                    int i = next;
                    if (p < i)
                    {
                        if (termsDict.find(p) == termsDict.end())
                            termsDict[p] = Dict<int, term>();
                        if (termsDict[p].find(i) == termsDict[p].end())
                            termsDict[p].insert({ i, term(p, i, d[next]) });

                        termsDict[p].at(i).w_ji = (double)s / (d[next] * d[next]);
                    }
                    else
                    {
                        if (termsDict.find(i) == termsDict.end())
                            termsDict[i] = Dict<int, term>();
                        if (termsDict[i].find(p) == termsDict[i].end())
                            termsDict[i].insert({ p, term(i, p, d[next]) });

                        termsDict[i].at(p).w_ij = (double)s / (d[next] * d[next]);
                    }
                }
            }
        }
    }
    // 1-stress
    for (int i=0; i<n; i++)
    {
        for (int j : graph[i])
        {
            if (i < j)
            {
                if (termsDict.find(i) == termsDict.end())
                    termsDict[i] = Dict<int, term>();
                if (termsDict[i].find(j) == termsDict[i].end())
                    termsDict[i].insert({ j, term(i, j, 1) });
                else
                    termsDict[i].at(j).d = 1;

                termsDict[i].at(j).w_ij = termsDict[i].at(j).w_ji = 1;
            }
        }
    }
    List<term> terms;
    for (auto const& i : termsDict)
    {
        for (auto const& term : i.second)
        {
            terms.push_back(term.second);
        }
    }
    return terms;
}

void sgd(double* X, List<term> &terms, const List<double> &etas)
{
    // iterate through step sizes
    int iteration = 0;
    for (double eta : etas)
    {
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());

        for (term const& t : terms)
        {
            // cap step size
            double mu_i = eta * t.w_ij;
            if (mu_i > 1)
                mu_i = 1;

            // cap step size
            double mu_j = eta * t.w_ji;
            if (mu_j > 1)
                mu_j = 1;

            double d_ij = t.d;
            int i = t.i, j = t.j;

            double dx = X[i*2]-X[j*2], dy = X[i*2+1]-X[j*2+1];
            double mag = sqrt(dx*dx + dy*dy);

            double r = (mag-d_ij) / (2*mag);
            double r_x = r * dx;
            double r_y = r * dy;

            X[i*2] -= mu_i * r_x;
            X[i*2+1] -= mu_i * r_y;
            X[j*2] += mu_j * r_x;
            X[j*2+1] += mu_j * r_y;
        }
    }
}


int main()
{
    int I[9] = { 0,0,1,2,3,4,4,5,6 };
    int J[9] = { 1,2,3,3,4,5,6,7,7 };

    Graph g = build_graph_unweighted(8, 9, I, J);

    auto pivots = maxmin_random_sp(g, 3, 0);
    auto terms = MSSP(g, pivots);
    //for (auto const& term : terms)
    //{
    //    std::cout << term.i << " " << term.j << " " << term.d << " " << term.w_ij << " " << term.w_ji << std::endl;
    //}
}
