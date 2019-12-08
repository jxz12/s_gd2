#include <vector>
#include <queue>
#include <set>
#include <map>
#include <limits>
#include <cmath>
#include <stdexcept>

// for debug and test
#include <iostream>

#include "layout.hpp"
#include "randomkit.h"

using std::vector;

void sgd(double* X, vector<term_sparse> &terms, const vector<double> &etas, const int seed)
{
    // seed random number generator
    // std::minstd_rand rng(seed);
    rk_state rstate;
    rk_seed(seed, &rstate);

    // iterate through step sizes
    for (unsigned i_eta=0; i_eta<etas.size(); i_eta++)
    {
        const double eta = etas[i_eta];
        // shuffle terms
        // std::shuffle(terms.begin(), terms.end(), rng);
        // fisheryates_shuffle(terms);
        fisheryates_shuffle(terms, rstate);

        unsigned n_terms = terms.size();
        for (unsigned i_term=0; i_term<n_terms; i_term++)
        {
            const term_sparse &t = terms[i_term];
            const int &i = t.i, &j = t.j;
            const double &d_ij = t.d;

            // cap step size
            double mu_i = eta * t.w_ij;
            if (mu_i > 1)
                mu_i = 1;

            // cap step size
            double mu_j = eta * t.w_ji;
            if (mu_j > 1)
                mu_j = 1;

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
void fisheryates_shuffle(vector<term_sparse> &terms, rk_state &rstate)
{
    int n = terms.size();
    for (int i=n-1; i>=1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        term_sparse temp = terms[i];
        terms[i] = terms[j];
        terms[j] = temp;
    }
}

// returns closest pivot for each vertex, not the pivots themselves
vector<int> maxmin_random_sp_unweighted(const vector<vector<int> >& graph, int n_pivots, int p0, int seed)
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
    // std::mt19937 rng(seed);
    // std::uniform_real_distribution<> uniform(0, 1);
    rk_state rstate;
    rk_seed(seed, &rstate);
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
        int sample = rk_interval(min_total, &rstate);
        int cumul = 0;
        int argmax = -1;
        for (int i = 0; i < n; i++)
        {
            cumul += mins[i];
            if (cumul >= sample)
            {
                argmax = i;
                break;
            }
        }
        if (argmax == -1)
            throw std::invalid_argument("unweighted pivot sampling failed");

        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs_unweighted(graph, argmax, mins, argmins);
    }
    return argmins;
}
void maxmin_bfs_unweighted(const vector<vector<int> >& graph, const int p, vector<int>& mins, vector<int>& argmins)
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

        for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
        {
            const int &next = graph[current][i_edge];
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
vector<term_sparse> MSSP_unweighted(const vector<vector<int> >& graph, const vector<int>& closest_pivots)
{
    int n = graph.size();

    // get pivots and their regions, but in sets
    std::map<int, std::set<int> > regions;
    std::map<int, std::map<int, term_sparse> > termsDict;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(closest_pivots[i]) == regions.end())
        {
            regions[closest_pivots[i]] = std::set<int>();
        }
        regions[closest_pivots[i]].insert(i);
    }

    std::map<int, std::set<int> >::iterator region;
    for (region=regions.begin(); region!=regions.end(); ++region)
    {
        // q contains next to visit
        std::queue<int> q;
        vector<int> d(n, -1);

        int p = region->first;
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

            for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
            {
                const int &next = graph[current][i_edge];
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
                    if (region->second.find(next) != region->second.end())
                    {
                        q2.push(d[next]);
                    }

                    int i = next;
                    if (i < p)
                    {
                        if (termsDict.find(i) == termsDict.end())
                            termsDict[i] = std::map<int, term_sparse>();
                        if (termsDict[i].find(p) == termsDict[i].end())
                            termsDict[i].insert(std::pair<int, term_sparse>(p, term_sparse(i, p, d[next])));

                        // termsDict[i].at(p).w_ij = s / ((double)d[next] * d[next]);
                        termsDict[i].find(p)->second.w_ij = s / ((double)d[next] * d[next]);
                    }
                    else
                    {
                        if (termsDict.find(p) == termsDict.end())
                            termsDict[p] = std::map<int, term_sparse>();
                        if (termsDict[p].find(i) == termsDict[p].end())
                            termsDict[p].insert(std::pair<int, term_sparse>(i, term_sparse(p, i, d[next])));

                        // termsDict[p].at(i).w_ji = s / ((double)d[next] * d[next]);
                        termsDict[p].find(i)->second.w_ji = s / ((double)d[next] * d[next]);
                    }
                }
            }
        }
    }
    // 1-stress
    for (int i=0; i<n; i++)
    {
        for (unsigned i_edge=0; i_edge<graph[i].size(); i_edge++)
        {
            const int j = graph[i][i_edge];
            if (i < j)
            {
                if (termsDict.find(i) == termsDict.end())
                    termsDict[i] = std::map<int, term_sparse>();
                if (termsDict[i].find(j) == termsDict[i].end())
                    termsDict[i].insert(std::pair<int, term_sparse>(j, term_sparse(i, j, 1)));
                else
                    // termsDict[i].at(j).d = 1;
                    termsDict[i].find(j)->second.d = 1;

                // termsDict[i].at(j).w_ij = termsDict[i].at(j).w_ji = 1;
                termsDict[i].find(j)->second.w_ij = termsDict[i].find(j)->second.w_ji = 1;
            }
        }
    }
    vector<term_sparse> terms;
    std::map<int, std::map<int, term_sparse> >::iterator it;
    for (it=termsDict.begin(); it!=termsDict.end(); ++it)
    {
        std::map<int, term_sparse>::iterator jt;
        for (jt=it->second.begin(); jt!=it->second.end(); ++jt)
        {
            terms.push_back(jt->second);
        }
    }
    return terms;
}

// returns closest pivot for each vertex, not pivots themselves
vector<int> maxmin_random_sp_weighted(const vector<vector<edge> >& graph, int n_pivots, int p0, int seed)
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
    // std::mt19937 rng(seed);
    // std::uniform_real_distribution<> uniform(0, 1);
    rk_state rstate;
    rk_seed(seed, &rstate);
    for (int i = 1; i < n_pivots; i++)
    {
        // choose pivots with probability min
        double min_total = 0;
        for (int i = 0; i < n; i++)
        {
            min_total += mins[i];
        }
        double sample = rk_double(&rstate) * min_total;
        double cumul = 0;
        int argmax = -1;
        for (int i = 0; i < n; i++)
        {
            cumul += mins[i];
            if (cumul >= sample)
            {
                argmax = i;
                break;
            }
        }
        if (argmax == -1)
            throw std::invalid_argument("weighted pivot sampling failed");

        mins[argmax] = 0;
        argmins[argmax] = argmax;
        maxmin_bfs_weighted(graph, argmax, mins, argmins);
    }
    return argmins;
}
void maxmin_bfs_weighted(const vector<vector<edge> >& graph, const int p, vector<double>& mins, vector<int>& argmins)
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
            for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
            {
                const edge &e = graph[current][i_edge];
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
vector<term_sparse> MSSP_weighted(const vector<vector<edge> >& graph, const vector<int>& closest_pivots)
{
    int n = graph.size();

    // get pivots and their regions, but in sets
    std::map<int, std::set<int> > regions;
    std::map<int, std::map<int, term_sparse> > termsDict;
    for (int i = 0; i < n; i++)
    {
        if (regions.find(closest_pivots[i]) == regions.end())
        {
            regions[closest_pivots[i]] = std::set<int>();
        }
        regions[closest_pivots[i]].insert(i);
    }

    std::map<int, std::set<int> >::iterator region;
    for (region=regions.begin(); region!=regions.end(); ++region)
    {
        int p = region->first;

        vector<bool> visited(n, false);
        vector<double> d(n, std::numeric_limits<double>::max()); // init 'tentative' distances to infinity

        // edges are used 'invisibly' in this queue
        std::priority_queue<edge, vector<edge>, edge_comp> pq;

        // init initial edges so that pivot-pivot term is avoided
        for (unsigned i_edge=0; i_edge<graph[p].size(); i_edge++)
        {
            const edge &e = graph[p][i_edge];

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
                if (region->second.find(current) != region->second.end())
                {
                    q2.push(d_pi);
                }

                int i = current;
                if (i < p)
                {
                    if (termsDict.find(i) == termsDict.end())
                        termsDict[i] = std::map<int, term_sparse>();
                    if (termsDict[i].find(p) == termsDict[i].end())
                        termsDict[i].insert(std::pair<int, term_sparse>(p, term_sparse(i, p, d_pi)));

                    // termsDict[i].at(p).w_ij = s / ((double)d_pi * d_pi);
                    termsDict[i].find(p)->second.w_ij = s / ((double)d_pi * d_pi);
                }
                else
                {
                    if (termsDict.find(p) == termsDict.end())
                        termsDict[p] = std::map<int, term_sparse>();
                    if (termsDict[p].find(i) == termsDict[p].end())
                        termsDict[p].insert(std::pair<int, term_sparse>(i, term_sparse(p, i, d_pi)));

                    // termsDict[p].at(i).w_ji = s / ((double)d_pi * d_pi);
                    termsDict[p].find(i)->second.w_ji = s / ((double)d_pi * d_pi);
                }

                // update tentative distances
                for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
                {
                    const edge &e = graph[current][i_edge];

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
        for (unsigned i_edge=0; i_edge<graph[i].size(); i_edge++)
        {
            const edge &e = graph[i][i_edge];

            int j = e.target;
            double d_ij = e.weight;
            if (i < j)
            {
                if (termsDict.find(i) == termsDict.end())
                    termsDict[i] = std::map<int, term_sparse>();
                if (termsDict[i].find(j) == termsDict[i].end())
                    termsDict[i].insert(std::pair<int, term_sparse>(j, term_sparse(i, j, d_ij)));
                else
                    // termsDict[i].at(j).d = d_ij;
                    termsDict[i].find(j)->second.d = d_ij;

                // termsDict[i].at(j).w_ij = termsDict[i].at(j).w_ji = 1/(d_ij*d_ij);
                termsDict[i].find(j)->second.w_ij = termsDict[i].find(j)->second.w_ji = 1/(d_ij*d_ij);
            }
        }
    }
    vector<term_sparse> terms;
    std::map<int, std::map<int, term_sparse> >::iterator it;
    for (it=termsDict.begin(); it!=termsDict.end(); ++it)
    {
        std::map<int, term_sparse>::iterator jt;
        for (jt=it->second.begin(); jt!=it->second.end(); ++jt)
        {
            terms.push_back(jt->second);
        }
    }
    return terms;
}

vector<double> schedule(const vector<term_sparse> &terms, int t_max, double eps)
{
    double w_min = std::numeric_limits<double>::max();
    double w_max = std::numeric_limits<double>::min();
    for (unsigned i=0; i<terms.size(); i++)
    {
        const term_sparse &term = terms[i];
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
    vector<vector<int> > g = build_graph_unweighted(n, m, I, J);

    vector<int> closest_pivots = maxmin_random_sp_unweighted(g, p, 0, seed);
    vector<term_sparse> terms = MSSP_unweighted(g, closest_pivots);
    vector<double> etas = schedule(terms, t_max, eps);
    sgd(X, terms, etas, seed);
}
void layout_sparse_weighted(int n, double* X, int m, int* I, int* J, double* V, int p, int t_max, double eps, int seed)
{
    vector<vector<edge> > g = build_graph_weighted(n, m, I, J, V);

    vector<int> closest_pivots = maxmin_random_sp_weighted(g, p, 0, seed);
    vector<term_sparse> terms = MSSP_weighted(g, closest_pivots);
    vector<double> etas = schedule(terms, t_max, eps);
    sgd(X, terms, etas, seed);
}