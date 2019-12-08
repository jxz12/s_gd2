#ifndef LAYOUT_HPP
#define LAYOUT_HPP

#include "randomkit.h"
#include <vector>
using std::vector;

///////////////////////
// visible to python //
///////////////////////
void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps, int seed);
void layout_weighted(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, int seed);
void layout_unweighted_convergent(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax, int seed);
void layout_weighted_convergent(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax, int seed);

void layout_sparse_unweighted(int n, double* X, int m, int* I, int* J, int p, int t_max, double eps, int seed);
void layout_sparse_weighted(int n, double* X, int m, int* I, int* J, double* V, int p, int t_max, double eps, int seed);

void mds_direct(int n, int kd, double* X, double* d, double* w, int t_max, double* etas, int seed);


//////////////
// cpp only //
//////////////
struct term
{
    int i, j;
    double d, w;
    term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}
};
void sgd(double* X, vector<term> &terms, const vector<double> &etas, const int seed);
void sgd_threshold(double* X, vector<term> &terms, const vector<double> &etas, const double delta, const int seed);
void sgd3D(double* X, vector<term> &terms, const vector<double> &etas, const int seed);
 
void fisheryates_shuffle(vector<term> &terms, rk_state &rstate);
double calculate_stress(double* X, const vector<term> &terms);

// for Dijkstra
struct edge
{
    // NOTE: this will be used for 'invisible' edges in the Dijkstra priority queue
    int target;
    double weight;
    edge(int target, double weight) : target(target), weight(weight) {}
};
struct edge_comp
{
    bool operator() (const edge &lhs, const edge &rhs) const
    {
        return lhs.weight > rhs.weight;
    }
};

vector<vector<int> > build_graph_unweighted(int n, int m, int* I, int* J);
vector<vector<edge> > build_graph_weighted(int n, int m, int* I, int* J, double* V);
vector<term> bfs(int n, int m, int* I, int* J);
vector<term> dijkstra(int n, int m, int* I, int* J, double* V);

vector<double> schedule(const vector<term> &terms, int t_max, double eps);
vector<double> schedule_convergent(const vector<term> &terms, int t_max, double eps, int t_maxmax);


// Ortmann et al. stuff
struct term_sparse
{
    int i, j;
    double d, w_ij, w_ji;
    term_sparse(int i, int j, double d) : i(i), j(j), d(d), w_ij(0), w_ji(0) {}
};
void sgd(double* X, vector<term_sparse>& terms, const vector<double>& etas, const int seed);

void fisheryates_shuffle(vector<term_sparse> &terms, rk_state &rstate);

vector<int> maxmin_random_sp_unweighted(const vector<vector<int> >& graph, int n_pivots, int p0, int seed);
vector<int> maxmin_random_sp_weighted(const vector<vector<edge> >& graph, int n_pivots, int p0, int seed);
void maxmin_bfs_unweighted(const vector<vector<int> >& graph, const int p, vector<int>& mins, vector<int>& argmins);
void maxmin_bfs_weighted(const vector<vector<edge> >& graph, const int p, vector<double>& mins, vector<int>& argmins);
vector<term_sparse> MSSP_unweighted(const vector<vector<int> >& graph, const vector<int>& pivots);
vector<term_sparse> MSSP_weighted(const vector<vector<edge> >& graph, const vector<int>& pivots);

#endif
