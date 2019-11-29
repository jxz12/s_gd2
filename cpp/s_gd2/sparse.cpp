//#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <random>
#include <exception>

// for debug and test
#include <chrono>
#include <iostream>

#include "layout.hpp"

using std::vector;
using std::unordered_map;
using std::unordered_set;

void sgd(double* X, vector<term_sparse> &terms, const vector<double> &etas, const int seed)
{
    // seed random number generator
    std::minstd_rand rng(seed);
    // iterate through step sizes
    for (double eta : etas)
    {
        // shuffle terms
        std::shuffle(terms.begin(), terms.end(), rng);

        for (const term_sparse& t : terms)
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

// returns closest pivot for each vertex, not the pivots themselves
vector<int> maxmin_random_sp_unweighted(const vector<vector<int>>& graph, int n_pivots, int p0, int seed)
{
    int n = graph.size();

    vector<int> mins(n, std::numeric_limits<int>::max());
    vector<int> argmins(n, -1);

    // first pivot
    mins[p0] = 0;
    argmins[p0] = p0;
    maxmin_bfs_unweighted(graph, p0, mins, argmins);
    for (int i = 0; i < n; i++)
    {
        if (argmins[i] == -1)
            throw std::invalid_argument("graph has multiple connected components");
    }

    // remaining pivots
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> uniform(0, 1);
    for (int i = 1; i < n_pivots; i++)
    {
        // int max = mins[0], argmax = 0;
        // for (int i = 1; i < n; i++)
        // {
        //     if (mins[i] > max)
        //     {
        //         max = mins[i];
        //         argmax = i;
        //     }
        // }
        // maxmin non-random above

        // choose pivots with probability min
        int min_total = 0;
        for (int i = 0; i < n; i++)
        {
            min_total += mins[i];
        }
        double rn = uniform(rng) * min_total;
        int cumul = 0;
        int argmax = 0;
        for (int i = 1; i < n; i++)
        {
            cumul += mins[i];
            if (cumul >= rn)
            {
                argmax = i;
                break;
            }
        }

        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs_unweighted(graph, argmax, mins, argmins);
    }
    return argmins;
}
void maxmin_bfs_unweighted(const vector<vector<int>>& graph, const int p, vector<int>& mins, vector<int>& argmins)
{
    int n = graph.size();
    std::queue<int> q;
    vector<int> d(n, -1);

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

// is not actually a multi-source shortest path, because regions come for free with maxmin_random_sp
vector<term_sparse> MSSP_unweighted(const vector<vector<int>>& graph, const vector<int>& closest_pivots)
{
    int n = graph.size();

    // get pivots and their regions, but in sets
    unordered_map<int, unordered_set<int>> regions;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(closest_pivots[i]) == regions.end())
        {
            regions[closest_pivots[i]] = unordered_set<int>();
        }
        regions[closest_pivots[i]].insert(i);
    }

    unordered_map<int, unordered_map<int, term_sparse>> termsDict;
    for (const auto& region : regions)
    {
        // q contains next to visit
        std::queue<int> q;
        vector<int> d(n, -1);

        int p = region.first;
        q.push(p);
        d[p] = 0;

        // q2 contains visited vertices' distances for s calculation
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

                    int i = next;
                    if (i < p)
                    {
                        if (termsDict.find(i) == termsDict.end())
                            termsDict[i] = unordered_map<int, term_sparse>();
                        if (termsDict[i].find(p) == termsDict[i].end())
                            termsDict[i].insert({ p, term_sparse(i, p, d[next]) });

                        termsDict[i].at(p).w_ij = s / ((double)d[next] * d[next]);
                    }
                    else
                    {
                        if (termsDict.find(p) == termsDict.end())
                            termsDict[p] = unordered_map<int, term_sparse>();
                        if (termsDict[p].find(i) == termsDict[p].end())
                            termsDict[p].insert({ i, term_sparse(p, i, d[next]) });

                        termsDict[p].at(i).w_ji = s / ((double)d[next] * d[next]);
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
                    termsDict[i] = unordered_map<int, term_sparse>();
                if (termsDict[i].find(j) == termsDict[i].end())
                    termsDict[i].insert({ j, term_sparse(i, j, 1) });
                else
                    termsDict[i].at(j).d = 1;

                termsDict[i].at(j).w_ij = termsDict[i].at(j).w_ji = 1;
            }
        }
    }
    vector<term_sparse> terms;
    for (const auto& i : termsDict)
    {
        for (const auto& j : i.second)
        {
            terms.push_back(j.second);
        }
    }
    return terms;
}

// returns closest pivot for each vertex, not pivots themselves
vector<int> maxmin_random_sp_weighted(const vector<vector<edge>>& graph, int n_pivots, int p0, int seed)
{
    int n = graph.size();

    vector<double> mins(n, std::numeric_limits<double>::max());
    vector<int> argmins(n, -1);

    // first pivot
    mins[p0] = 0;
    argmins[p0] = p0;
    maxmin_bfs_weighted(graph, p0, mins, argmins);
    for (int i = 0; i < n; i++)
    {
        if (argmins[i] == -1)
            throw std::invalid_argument("graph has multiple connected components");
    }

    // remaining pivots
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> uniform(0, 1);
    for (int i = 1; i < n_pivots; i++)
    {
        // choose pivots with probability min
        double min_total = 0;
        for (int i = 0; i < n; i++)
        {
            min_total += mins[i];
        }
        double rn = uniform(rng) * min_total;
        double cumul = 0;
        int argmax = n-1;
        for (int i = 1; i < n; i++)
        {
            cumul += mins[i];
            if (cumul >= rn)
            {
                argmax = i;
                break;
            }
        }

        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs_weighted(graph, argmax, mins, argmins);
    }
    return argmins;
}
void maxmin_bfs_weighted(const vector<vector<edge>>& graph, const int p, vector<double>& mins, vector<int>& argmins)
{
    int n = graph.size();
    vector<bool> visited(n, false);
    vector<double> d(n, std::numeric_limits<double>::max()); // init 'tentative' distances to infinity

    // I am not using a fibonacci heap. I AM NOT USING A FIBONACCI HEAP
    // edges are used 'invisibly' here
    std::priority_queue<edge, vector<edge>, edge_comp> pq;

    d[p] = 0;
    pq.push(edge(p,0));

    while (!pq.empty())
    {
        int current = pq.top().target;
        double d_pi = pq.top().weight;
        pq.pop();
        
        if (!visited[current]) // ignore redundant elements in queue
        {
            visited[current] = true;

            if (d_pi < mins[current])
            {
                mins[current] = d_pi;
                argmins[current] = p;
            }
            for (edge e : graph[current])
            {
                // here the edge is not 'invisible'
                int next = e.target;
                double weight = e.weight;

                if (d_pi + weight < d[next]) // update tentative value of d 
                {
                    d[next] = d_pi + weight;
                    pq.push(edge(next, d[next]));
                }
            }
        }
    }
}

// again, not a proper MSSP because we get regions for free with maxmin_random_sp
vector<term_sparse> MSSP_weighted(const vector<vector<edge>>& graph, const vector<int>& closest_pivots)
{
    int n = graph.size();

    // get pivots and their regions, but in sets
    unordered_map<int, unordered_set<int>> regions;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(closest_pivots[i]) == regions.end())
        {
            regions[closest_pivots[i]] = unordered_set<int>();
        }
        regions[closest_pivots[i]].insert(i);
    }

    unordered_map<int, unordered_map<int, term_sparse>> termsDict;
    for (const auto& region : regions)
    {
        int p = region.first;

        vector<bool> visited(n, false);
        vector<double> d(n, std::numeric_limits<double>::max()); // init 'tentative' distances to infinity

        // edges are used 'invisibly' in this queue
        std::priority_queue<edge, vector<edge>, edge_comp> pq;

        // init initial edges so that pivot-pivot term is avoided
        for (edge e : graph[p])
        {
            // here the edge is not 'invisible'
            int next = e.target;
            double weight = e.weight;

            d[next] = weight; // init tentative value of d
            pq.push(edge(next, d[next]));
        }
        d[p] = 0;
        visited[p] = true;

        // q2 contains visited vertices' distances for s calculation
        std::queue<double> q2;
        int s = 1;

        while (!pq.empty())
        {
            int current = pq.top().target;
            double d_pi = pq.top().weight;
            pq.pop();
            
            if (!visited[current]) // ignore redundant elements in queue
            {
                visited[current] = true;

                // empty the second queue enough to calculate s
                while (!q2.empty() && q2.front() <= d_pi/2)
                {
                    q2.pop();
                    s += 1;
                }
                if (region.second.find(current) != region.second.end())
                {
                    q2.push(d_pi);
                }

                int i = current;
                if (i < p)
                {
                    if (termsDict.find(i) == termsDict.end())
                        termsDict[i] = unordered_map<int, term_sparse>();
                    if (termsDict[i].find(p) == termsDict[i].end())
                        termsDict[i].insert({ p, term_sparse(i, p, d_pi) });

                    termsDict[i].at(p).w_ij = s / ((double)d_pi * d_pi);
                }
                else
                {
                    if (termsDict.find(p) == termsDict.end())
                        termsDict[p] = unordered_map<int, term_sparse>();
                    if (termsDict[p].find(i) == termsDict[p].end())
                        termsDict[p].insert({ i, term_sparse(p, i, d_pi) });

                    termsDict[p].at(i).w_ji = s / ((double)d_pi * d_pi);
                }

                // update tentative distances
                for (edge e : graph[current])
                {
                    // here the edge is not 'invisible'
                    int next = e.target;
                    double weight = e.weight;

                    if (d_pi + weight < d[next]) 
                    {
                        d[next] = d_pi + weight; // update tentative value of d
                        pq.push(edge(next, d[next]));
                    }
                }
            }
        }
    }
    // 1-stress
    for (int i=0; i<n; i++)
    {
        for (edge e : graph[i])
        {
            int j = e.target;
            double d_ij = e.weight;
            if (i < j)
            {
                if (termsDict.find(i) == termsDict.end())
                    termsDict[i] = unordered_map<int, term_sparse>();
                if (termsDict[i].find(j) == termsDict[i].end())
                    termsDict[i].insert({ j, term_sparse(i, j, d_ij) });
                else
                    termsDict[i].at(j).d = d_ij;

                termsDict[i].at(j).w_ij = termsDict[i].at(j).w_ji = 1/(d_ij*d_ij);
            }
        }
    }
    vector<term_sparse> terms;
    for (const auto& i : termsDict)
    {
        for (const auto& j : i.second)
        {
            terms.push_back(j.second);
        }
    }
    return terms;
}

vector<double> schedule(const vector<term_sparse> &terms, int t_max, double eps)
{
    double w_min = std::numeric_limits<double>::max();
    double w_max = std::numeric_limits<double>::min();
    for (const auto& term : terms)
    {
        if (term.w_ij < w_min && term.w_ij != 0) w_min = term.w_ij;
        if (term.w_ji < w_min && term.w_ji != 0) w_min = term.w_ji;

        if (term.w_ij > w_max) w_max = term.w_ij;
        if (term.w_ji > w_max) w_max = term.w_ji;
    }
    double eta_max = 1.0 / w_min;
    double eta_min = eps / w_max;

    double lambda = log(eta_max/eta_min) / ((double)t_max-1);

    // initialize step sizes
    std::vector<double> etas;
    etas.reserve(t_max);
    for (int t=0; t<t_max; t++)
        etas.push_back(eta_max * exp(-lambda * t));

    return etas;
}


void layout_sparse_unweighted(int n, double* X, int m, int* I, int* J, int p, int t_max, double eps, int seed)
{
    vector<vector<int>> g = build_graph_unweighted(n, m, I, J);

    auto closest_pivots = maxmin_random_sp_unweighted(g, p, 0);
    auto terms = MSSP_unweighted(g, closest_pivots);
    auto etas = schedule(terms, t_max, eps);
    sgd(X, terms, etas, seed);
}
void layout_sparse_weighted(int n, double* X, int m, int* I, int* J, double* V, int p, int t_max, double eps, int seed)
{
    vector<vector<edge>> g = build_graph_weighted(n, m, I, J, V);

    auto closest_pivots = maxmin_random_sp_weighted(g, p, 0);
    auto terms = MSSP_weighted(g, closest_pivots);
    auto etas = schedule(terms, t_max, eps);
    sgd(X, terms, etas, seed);
}