//#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <cmath>

// for testing
#include <chrono>
#include <iostream>

#include "layout.hpp"

using std::vector;

void sgd(double* X, vector<term> &terms, const vector<double> &etas, double delta)
{
    // iterate through step sizes
    for (double eta : etas)
    {
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());

        double Delta_max = 0;
        for (const term &t : terms)
        {
            // cap step size
            double w_ij = t.w;
            double mu = eta * w_ij;
            if (mu > 1)
                mu = 1;

            double d_ij = t.d;
            int i = t.i, j = t.j;

            double dx = X[i*2]-X[j*2], dy = X[i*2+1]-X[j*2+1];
            double mag = sqrt(dx*dx + dy*dy);

            // check distances for early stopping
            double Delta = mu * (mag-d_ij) / 2;
            if (Delta > Delta_max)
                Delta_max = Delta;

            double r = Delta / mag;
            double r_x = r * dx;
            double r_y = r * dy;
            
            X[i*2] -= r_x;
            X[i*2+1] -= r_y;
            X[j*2] += r_x;
            X[j*2+1] += r_y;
        }
        //std::cerr << ++iteration << ", eta: " << eta << ", Delta: " << Delta_max << std::endl;
        if (Delta_max < delta)
            return;
    }
}

vector<vector<int>> build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    vector<std::unordered_set<int>> undirected(n);
    vector<vector<int>> graph(n);

    for (int ij=0; ij<m; ij++)
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

// calculates the unweighted shortest paths between indices I and J
// using a breadth-first search, returning a vector of terms
vector<term> bfs(int n, int m, int* I, int* J)
{
    auto graph = build_graph_unweighted(n, m, I, J);

    int nC2 = (n*(n-1))/2;
    vector<term> terms;
    terms.reserve(nC2);

    int terms_size_goal = 0; // to keep track of when to stop searching i<j

    for (int source=0; source<n-1; source++) // no need to do final vertex because i<j
    {
        vector<int> d(n, -1); // distances from source
        std::queue<int> q;
        
        d[source] = 0;
        q.push(source);

        terms_size_goal += n-source-1; // this is how many terms exist for i<j

        while (!q.empty() && terms.size() <= terms_size_goal)
        {
            int current = q.front();
            q.pop();
            for (int next : graph[current])
            {
                if (d[next] == -1)
                {
                    q.push(next);
                    int d_ij = d[current] + 1;
                    d[next] = d_ij;

                    if (source < next) // only add terms for i<j
                    {
                        double w_ij = 1.0 / ((double)d_ij*d_ij);
                        terms.push_back(term(source, next, d_ij, w_ij));
                    }
                }
            }
        }
        if (terms.size() != terms_size_goal)
        {
            throw "graph is not strongly connected, or is not indexed from zero";
        }
    }
    return terms;
}


vector<vector<edge>> build_graph_weighted(int n, int m, int* I, int* J, double* V)
{
    // used to make graph undirected, in case graph is not already
    vector<std::unordered_map<int, double>> undirected(n);
    vector<vector<edge>> graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw "i or j bigger than n";

        double v = V[ij];
        if (v <= 0)
            throw "v less or equal 0";

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if key not there
        {
            undirected[i].insert({j, v});
            undirected[j].insert({i, v});
            graph[i].push_back(edge(j, v));
            graph[j].push_back(edge(i, v));
        }
        else
        {
            if (undirected[j][i] != v)
                throw "graph weights not symmetric";
        }
    }
    return graph;
}

// calculates the unweighted shortest paths between indices I and J
// using Dijkstra's algorithm, returning a vector of terms
vector<term> dijkstra(int n, int m, int* I, int* J, double* V)
{
    auto graph = build_graph_weighted(n, m, I, J, V);

    int nC2 = (n*(n-1))/2;
    vector<term> terms;
    terms.reserve(nC2);

    int terms_size_goal = 0; // to keep track of when to stop searching i<j

    for (int source=0; source<n-1; source++) // no need to do final vertex because i<j
    {
        vector<bool> visited(n, false);
        vector<double> d(n, std::numeric_limits<double>::max()); // init 'tentative' distances to infinity

        // I am not using a fibonacci heap. I AM NOT USING A FIBONACCI HEAP
        // edges are used 'invisibly' here
        std::priority_queue<edge, vector<edge>, edge_comp> pq;

        d[source] = 0;
        pq.push(edge(source,0));

        terms_size_goal += n-source-1; // this is how many terms exist for i<j

        while (!pq.empty() && terms.size() <= terms_size_goal)
        {
            int current = pq.top().target;
            double d_ij = pq.top().weight;
            pq.pop();
            
            if (!visited[current]) // ignore redundant elements in queue
            {
                visited[current] = true;

                if (source < current) // only add terms for i<j
                {
                    double w_ij = 1.0 / (d_ij*d_ij);
                    terms.push_back(term(source, current, d_ij, w_ij));
                }
                for (edge e : graph[current])
                {
                    // here the edge is not 'invisible'
                    int next = e.target;
                    double weight = e.weight;

                    if (d_ij + weight < d[next]) // update tentative value of d 
                    {
                        d[next] = d_ij + weight;
                        pq.push(edge(next, d[next]));
                    }
                }
            }
        }
        if (terms.size() != terms_size_goal)
        {
            throw "graph is not strongly connected, or is not indexed from zero";
        }
    }
    return terms;
}


vector<double> schedule(const vector<term> &terms, int t_max, double eps)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (int i=1; i<terms.size(); i++)
    {
        double w = terms[i].w;
        if (w < w_min) w_min = w;
        if (w > w_max) w_max = w;
    }
    double eta_max = 1.0 / w_min;
    double eta_min = eps / w_max;

    double lambda = log(eta_max/eta_min) / ((double)t_max-1);

    // initialize step sizes
    vector<double> etas;
    etas.reserve(t_max);
    for (int t=0; t<t_max; t++)
        etas.push_back(eta_max * exp(-lambda * t));

    return etas;
}

vector<double> schedule_convergent(const vector<term> &terms, int t_max, double eps, int t_maxmax)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (int i=1; i<terms.size(); i++)
    {
        double w = terms[i].w;
        if (w < w_min) w_min = w;
        if (w > w_max) w_max = w;
    }
    double eta_max = 1.0 / w_min;
    double eta_min = eps / w_max;

    double lambda = log(eta_max/eta_min) / ((double)t_max-1);

    // initialize step sizes
    vector<double> etas;
    etas.reserve(t_maxmax);
    double eta_switch = 1.0 / w_max;
    for (int t=0; t<t_maxmax; t++)
    {
        double eta = eta_max * exp(-lambda * t);
        if (eta < eta_switch)
            break;

        etas.push_back(eta);
    }
    int tau = etas.size();
    for (int t=tau; t<t_maxmax; t++)
    {
        double eta = eta_switch / (1 + lambda*((double)t-tau));
        etas.push_back(eta);
    }

    return etas;
}

void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps)
{
    try
    {
        //auto start = std::chrono::steady_clock::now();

        vector<term> terms = bfs(n, m, I, J);
        //auto end = std::chrono::steady_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;

        vector<double> etas = schedule(terms, t_max, eps);
        sgd(X, terms, etas);
        //end = std::chrono::steady_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}

void layout_weighted(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps)
{
    try
    {
        vector<term> terms = dijkstra(n, m, I, J, V);
        vector<double> etas = schedule(terms, t_max, eps);
        sgd(X, terms, etas);
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}
void layout_unweighted_convergent(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax)
{
    try
    {
        vector<term> terms = bfs(n, m, I, J);
        vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
        sgd(X, terms, etas, delta);
    } 
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}
void layout_weighted_convergent(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax)
{
    try
    {
        vector<term> terms = dijkstra(n, m, I, J, V);
        vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
        sgd(X, terms, etas, delta);
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}

// d and w should be condensed distance matrices
void mds_direct(int n, double* X, double* d, double* w, int t_max, double* eta)
{
    // initialize SGD
    int nC2 = (n*(n-1))/2;
    vector<term> terms;
    terms.reserve(nC2);
    int ij=0;
    for (int i=0; i<n; i++) // unpack the condensed distance matrices
    {
        for (int j=i+1; j<n; j++)
        {
            terms.push_back(term(i, j, d[ij], w[ij]));
            ij += 1;
        }
    }

    // initialize step sizes
    vector<double> etas;
    etas.reserve(t_max);
    for (int t=0; t<t_max; t++)
    {
        etas.push_back(eta[t]);
    }
    
    sgd(X, terms, etas, 0);
}


