#include "layout.hpp"

#include <cmath>
#include <algorithm>
#include <queue>
#include <vector>

void sgd(double* X, std::vector<term*> &terms, std::vector<double> &etas, double delta)
{
	// iterate through step sizes
	for (double eta : etas)
	{
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());

        double mag_max = 0;
		for (term* ijdw : terms)
		{
			// cap step size
            double w_ij = ijdw->w;
			double mu = eta * w_ij;
			if (mu > 1)
				mu = 1;

			double d_ij = ijdw->d;

			int i = ijdw->i, j = ijdw->j;
			double del_x = X[i*2] - X[j*2], del_y = X[i*2+1] - X[j*2+1];
			double mag = sqrt(del_x*del_x + del_y*del_y);
            if (mag > mag_max)
                mag_max = mag;
			
			double r = mu * (mag - d_ij) / (2 * mag);
			double r_x = r * del_x;
			double r_y = r * del_y;
			
			X[i*2] -= r_x;
			X[i*2+1] -= r_y;
			X[j*2] += r_x;
			X[j*2+1] += r_y;
		}
        if (mag_max < delta)
            return;
	}
}

void sgd_horizontal(double* X, std::vector<term*> &terms, std::vector<double> &etas, double delta)
{
	// iterate through step sizes
	for (double eta : etas)
	{
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());

        double mag_max = 0;
		for (term* ijdw : terms)
		{
			// cap step size
            double w_ij = ijdw->w;
			double mu = eta * w_ij;
			if (mu > 1)
				mu = 1;

			double d_ij = ijdw->d;

			int i = ijdw->i, j = ijdw->j;
			double del_x = X[i*2] - X[j*2], del_y = X[i*2+1] - X[j*2+1];
			double mag = sqrt(del_x*del_x + del_y*del_y);
            if (mag > mag_max)
                mag_max = mag;
			
			double r = mu * (mag - d_ij) / (2 * mag);
			double r_x = r * del_x;
			//double r_y = r * del_y;
			
			X[i*2] -= r_x;
			//X[i*2+1] -= r_y;
			X[j*2] += r_x;
			//X[j*2+1] += r_y;
            // do not change y coord positions, but still take them into account
		}
        if (mag_max < delta)
            return;
	}
}

// calculates the unweighted shortest paths between indices I and J
// using a breadth-first search, returning a vector of term*s
std::vector<term*> bfs(int n, int m, int* I, int* J)
{
    std::vector<std::vector<int>> graph(n);
    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw "i or j bigger than n";

        graph[i].push_back(j);
    }

	int nC2 = (n*(n-1))/2;
	std::vector<term*> terms;
	terms.reserve(nC2);

    int terms_size_goal = 0; // to keep track of when to stop searching i<j

    for (int source=0; source<n-1; source++) // no need to do final vertex because i<j
    {
        std::vector<int> d(n, -1); // distances from source
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
                        double w_ij = 1.0 / (d_ij*d_ij);
                        terms.push_back(new term(source, next, d_ij, w_ij));
                    }
                }
            }
        }
        if (terms.size() != terms_size_goal)
        {
            throw "graph is not strongly connected";
        }
    }
    return terms;
}


std::vector<double> schedule(std::vector<term*> &terms, double t_max, double eps)
{
    double w_min = terms[0]->w, w_max = terms[0]->w;
    for (int i=1; i<terms.size(); i++)
    {
        double w = terms[i]->w;
        if (w < w_min)
            w_min = w;

        if (w > w_max)
            w_max = w;
    }
    double eta_max = 1.0 / w_min;
    double eta_min = eps / w_max;

    double lambda = log(eta_max/eta_min) / (t_max-1);

    // initialize step sizes
    std::vector<double> etas;
    etas.reserve(t_max);
    for (int t=0; t<t_max; t++)
        etas.push_back(eta_max * exp(-lambda * t));

    return etas;
}

void sgd_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps, bool horizontal)
{
    std::vector<term*> terms = bfs(n, m, I, J);
    std::vector<double> etas = schedule(terms, t_max, eps);

    if (!horizontal)
        sgd(X, terms, etas, 0);
    else
        sgd_horizontal(X, terms, etas, 0);
}

// calculates the unweighted shortest paths between indices I and J
// using dijkstra's algorithm, returning a vector of term*s

std::vector<term*> dijkstra(int n, int m, int* I, int* J, double* V)
{
    std::vector<std::vector<edge>> graph(n);
    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw "i or j bigger than n";

        double v = V[ij];
        graph[i].push_back(edge(j, v));
    }

	int nC2 = (n*(n-1))/2;
	std::vector<term*> terms;
	terms.reserve(nC2);

    int terms_size_goal = 0; // to keep track of when to stop searching i<j

    for (int source=0; source<n-1; source++) // no need to do final vertex because i<j
    {
        std::vector<bool> visited(n, false);
        std::vector<double> d(n, n); // distances from source, n is effectively infinity

        // I am not using a fibonacci heap. I AM NOT USING A FIBONACCI HEAP
        // edges are used 'invisibly' here
        std::priority_queue<edge, std::vector<edge>, edge_comp> pq;

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
                    terms.push_back(new term(source, current, d_ij, w_ij));
                }
                for (edge e : graph[current])
                {
                    // here the edge is not 'invisible'
                    int next = e.target;
                    double weight = e.weight;

                    if (d[next] > d_ij + weight) 
                    {
                        d[next] = d_ij + weight; // update tentative value of d
                        pq.push(edge(next, d[next]));
                    }
                }
            }
        }
        if (terms.size() != terms_size_goal)
        {
            throw "graph is not strongly connected";
        }
    }
    return terms;
}

// d and w should be condensed distance matrices
void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta, bool horizontal)
{
	// initialize SGD
	int nC2 = (n*(n-1))/2;
	std::vector<term*> terms;
	terms.reserve(nC2);
    int ij=0;
	for (int i=0; i<n; i++) // unpack the condensed distance matrices
	{
		for (int j=i+1; j<n; j++)
		{
			terms.push_back(new term(i, j, d[ij], w[ij]));
            ij += 1;
		}
	}

    // initialize step sizes
    std::vector<double> etas;
    etas.reserve(t_max);
    for (int t=0; t<t_max; t++)
    {
        etas.push_back(etas[t]);
    }

    // perform optimisation
    if (!horizontal)
        sgd(X, terms, etas, 0);
    else
        sgd_horizontal(X, terms, etas, 0);
}


