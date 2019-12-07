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

void sgd(double* X, vector<term> &terms, const vector<double> &etas, const int seed)
{
    // seed random number generator
    // std::minstd_rand rng(seed);
    // std::mt19937 rng(seed);
    // srand(seed);
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
            const term &t = terms[i_term];
            const int &i = t.i, &j = t.j;
            const double &w_ij = t.w;
            const double &d_ij = t.d;

            // cap step size
            double mu = eta * w_ij;
            if (mu > 1)
                mu = 1;

            double dx = X[i*2]-X[j*2], dy = X[i*2+1]-X[j*2+1];
            double mag = sqrt(dx*dx + dy*dy);

            // check distances for early stopping
            double r = (mu * (mag-d_ij)) / (2*mag);
            double r_x = r * dx;
            double r_y = r * dy;
            
            X[i*2] -= r_x;
            X[i*2+1] -= r_y;
            X[j*2] += r_x;
            X[j*2+1] += r_y;
        }
    }
}
void fisheryates_shuffle(vector<term> &terms, rk_state &rstate)
{
    int n = terms.size();
    for (unsigned i=n-1; i>=1; i--)
    {
        unsigned j = rk_interval(i, &rstate);
        term temp = terms[i];
        terms[i] = terms[j];
        terms[j] = temp;
    }
}
double calculate_stress(double* X, const vector<term> &terms)
{
    double stress = 0;
    unsigned n_terms = terms.size();
    for  (unsigned i_term=0; i_term<n_terms; i_term++)
    {
        const term &t = terms[i_term];
        const double &w_ij = t.w;
        const double &d_ij = t.d;
        const int &i = t.i, &j = t.j;

        double dx = X[i*2]-X[j*2], dy = X[i*2+1]-X[j*2+1];
        double stretch = d_ij - sqrt(dx*dx + dy*dy);
        stress += w_ij*stretch*stretch;
    }
    return stress;
}

vector<vector<int> > build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    vector<std::set<int> > undirected(n);
    vector<vector<int> > graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

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
    vector<vector<int> > graph = build_graph_unweighted(n, m, I, J);

    int nC2 = (n*(n-1))/2;
    vector<term> terms;
    terms.reserve(nC2);

    unsigned terms_size_goal = 0; // to keep track of when to stop searching i<j

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

            for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
            {
                const int &next = graph[current][i_edge];
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
            throw std::invalid_argument("graph is not strongly connected, or is not indexed from zero");
        }
    }
    return terms;
}


vector<vector<edge> > build_graph_weighted(int n, int m, int* I, int* J, double* V)
{
    // used to make graph undirected, in case graph is not already
    vector<std::map<int, double> > undirected(n);
    vector<vector<edge> > graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw std::invalid_argument("i or j bigger than n");

        double v = V[ij];
        if (v <= 0)
            throw std::invalid_argument("edge length less than or equal to 0");

        if (i != j && undirected[j].find(i) == undirected[j].end()) // if key not there
        {
            // undirected[i].insert({j, v});
            // undirected[j].insert({i, v});
            undirected[i].insert(std::pair<int,double>(j, v));
            undirected[j].insert(std::pair<int,double>(i, v));
            graph[i].push_back(edge(j, v));
            graph[j].push_back(edge(i, v));
        }
        else
        {
            if (undirected[j][i] != v)
                throw std::invalid_argument("graph edge lengths not symmetric");
        }
    }
    return graph;
}

// calculates the unweighted shortest paths between indices I and J
// using Dijkstra's algorithm, returning a vector of terms
vector<term> dijkstra(int n, int m, int* I, int* J, double* V)
{
    vector<vector<edge> > graph = build_graph_weighted(n, m, I, J, V);

    int nC2 = (n*(n-1))/2;
    vector<term> terms;
    terms.reserve(nC2);

    unsigned terms_size_goal = 0; // to keep track of when to stop searching i<j

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

                for (unsigned i_edge=0; i_edge<graph[current].size(); i_edge++)
                {
                    const edge &e = graph[current][i_edge];
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
            throw std::invalid_argument("graph is not strongly connected, or is not indexed from zero");
        }
    }
    return terms;
}


vector<double> schedule(const vector<term> &terms, int t_max, double eps)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (unsigned i=1; i<terms.size(); i++)
    {
        const double &w = terms[i].w;
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
    for (unsigned i=1; i<terms.size(); i++)
    {
        const double &w = terms[i].w;
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

void sgd_threshold(double* X, vector<term> &terms, const vector<double> &etas, double delta, const int seed)
{
    // seed random number generator
    // std::minstd_rand rng(seed);
    // srand(seed);
    rk_state rstate;
    rk_seed(seed, &rstate);

    // iterate through step sizes
    for (unsigned i_eta=0; i_eta<etas.size(); i_eta++)
    {
        // shuffle terms
        // std::shuffle(terms.begin(), terms.end(), rng);
        // fisheryates_shuffle(terms);
        fisheryates_shuffle(terms, rstate);

        const double &eta = etas[i_eta];
        unsigned n_terms = terms.size();
        double Delta_max = 0;
        for (unsigned i_term=0; i_term<n_terms; i_term++)
        {
            const term &t = terms[i_term];
            const int &i = t.i, &j = t.j;
            const double &w_ij = t.w;
            const double &d_ij = t.d;

            // cap step size
            double mu = eta * w_ij;
            if (mu > 1)
                mu = 1;

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
        if (Delta_max < delta)
            return;
    }
}
void sgd3D(double* X, vector<term> &terms, const vector<double> &etas, const int seed)
{
    // seed random number generator
    // std::minstd_rand rng(seed);
    // srand(seed);
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
            const term &t = terms[i_term];
            const int &i = t.i, &j = t.j;
            const double &w_ij = t.w;
            const double &d_ij = t.d;

            // cap step size
            double mu = eta * w_ij;
            if (mu > 1)
                mu = 1;

            double dx = X[i*3]-X[j*3], dy = X[i*3+1]-X[j*3+1], dz = X[i*3+2]-X[j*3+2];
            double mag = sqrt(dx*dx + dy*dy + dz*dz);

            double r = mu * (mag-d_ij) / (2*mag);
            double r_x = r * dx;
            double r_y = r * dy;
            double r_z = r * dz;
            
            X[i*3] -= r_x;
            X[i*3+1] -= r_y;
            X[i*3+2] -= r_z;
            X[j*3] += r_x;
            X[j*3+1] += r_y;
            X[j*3+2] += r_z;
        }
    }
}

void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps, int seed)
{
    vector<term> terms = bfs(n, m, I, J);
    vector<double> etas = schedule(terms, t_max, eps);
    sgd(X, terms, etas, seed);
}

void layout_weighted(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, int seed)
{
    vector<term> terms = dijkstra(n, m, I, J, V);
    vector<double> etas = schedule(terms, t_max, eps);
    sgd(X, terms, etas, seed);
}
void layout_unweighted_convergent(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax, int seed)
{
    vector<term> terms = bfs(n, m, I, J);
    vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
    sgd_threshold(X, terms, etas, delta, seed);
}
void layout_weighted_convergent(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax, int seed)
{
    vector<term> terms = dijkstra(n, m, I, J, V);
    vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
    sgd_threshold(X, terms, etas, delta, seed);
}

// d and w should be condensed distance matrices
void mds_direct(int n, int kd, double* X, double* d, double* w, int t_max, double* eta, int seed)
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
    
    if (kd == 2)
        sgd(X, terms, etas, seed);
    else if (kd == 3)
        sgd3D(X, terms, etas, seed);
    else
        throw std::invalid_argument("only 2 or 3 dimensional layouts are supported");
}