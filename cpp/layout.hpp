#ifndef LAYOUT_HPP
#define LAYOUT_HPP

#include <vector>

struct term
{
	int i, j;
	double d, w;
	term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}
};
struct edge
{
    // NOTE: this will be used for 'invisible' edges too!
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

void sgd(double* X, std::vector<term*> &terms, std::vector<double> &etas, double delta=0);
void sgd_horizontal(double* X, std::vector<term*> &terms, std::vector<double> &etas, double delta=0); // only changes x coords

std::vector<double> schedule(std::vector<term*> &terms, double t_max=15, double eps=.1);
// std::vector<double> schedule_convergent(std::vector<term*> terms, double t_max=30, double eps=.1, double delta=.03);

std::vector<term*> bfs(int n, int m, int* I, int* J);
void sgd_unweighted(int n, double* X, int m, double* I, double* J, int t_max, double eps, bool horizontal=false);
// void sgd_unweighted_convergent(int n, double* X, int m, double* I, double* J, int t_max, double eps, double delta, bool horizontal=false);

std::vector<term*> dijkstra(int n, int m, int* I, int* J, double* V);
void sgd_weighted(int n, double* X, int m, double* I, double* J, double* V, int t_max, double eps, bool horizontal=false);

void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* etas, bool horizontal=false);

#endif
