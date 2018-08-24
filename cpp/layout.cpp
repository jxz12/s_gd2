#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>

struct term {
	term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}
	int i, j;
	double d, w;
};

void sgd(double* X, std::vector<term> terms, std::vector<double> eta);
void sgd_horizontal(double* X, std::vector<term> terms, std::vector<double> eta); // only changes x coords
void bacon(int n, int m, int* I, int* J, int* d_out); // BFS shortest paths algorithm

void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta);
void sgd_direct_horizontal(int n, double* X, double* d, double* w, int t_max, double* eta);

void sgd(double* X, std::vector<term> terms, std::vector<double> eta)
{
	// iterate through step sizes
	for (auto it = eta.begin(); it != eta.end(); ++it)
	{
        double step = *it;
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());
		for (auto ij = terms.begin(); ij != terms.end(); ++ij)
		{
			// cap step size
			double mu = step * ij->w;
			if (mu > 1)
				mu = 1;

			double d_ij = ij->d;

			int i = ij->i, j = ij->j;
			double del_x = X[i*2] - X[j*2], del_y = X[i*2+1] - X[j*2+1];
			double mag = sqrt(del_x*del_x + del_y*del_y);
			
			double r = mu * (mag - d_ij) / (2 * mag);
			double r_x = r * del_x;
			double r_y = r * del_y;
			
			X[i*2] -= r_x;
			X[i*2+1] -= r_y;
			X[j*2] += r_x;
			X[j*2+1] += r_y;
		}
	}
}

void sgd_horizontal(double* X, std::vector<term> terms, std::vector<double> eta)
{
	// iterate through step sizes
	for (auto it = eta.begin(); it != eta.end(); ++it)
	{
        double step = *it;
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());
		for (auto ij = terms.begin(); ij != terms.end(); ++ij)
		{
			// cap step size
			double mu = step * ij->w;
			if (mu > 1)
				mu = 1;

			double d_ij = ij->d;

			int i = ij->i, j = ij->j;
			double del_x = X[i*2] - X[j*2], del_y = X[i*2+1] - X[j*2+1];
			double mag = sqrt(del_x*del_x + del_y*del_y);
			
			double r = mu * (mag - d_ij) / (2 * mag);
			double r_x = r * del_x;
			//double r_y = r * del_y;
			
			X[i*2] -= r_x;
			//X[i*2+1] -= r_y;
			X[j*2] += r_x;
			//X[j*2+1] += r_y;
            // do not change y coord positions, but still take them into account
		}
	}
}

// d and w should be condensed distance matrices
void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta)
{
	// initialize SGD
	int nC2 = (n*(n-1))/2;
	std::vector<term> terms;
	terms.reserve(nC2);
    int ij = 0;
	for (int i=0; i<n; i++)
	{
		for (int j=i+1; j<n; j++)
		{
			terms.push_back(term(i, j, d[ij], w[ij]));
            ij += 1;
		}
	}

    // initialize step sizes
    std::vector<double> eta_vec;
    eta_vec.reserve(t_max);
    for (int t=0; t<t_max; t++)
    {
        eta_vec.push_back(eta[t]);
    }

    // perform optimisation
    sgd(X, terms, eta_vec);
}

// d and w should be condensed distance matrices
void sgd_direct_horizontal(int n, double* X, double* d, double* w, int t_max, double* eta)
{
	// initialize SGD
	int nC2 = (n*(n-1))/2;
	std::vector<term> terms;
	terms.reserve(nC2);
    int ij = 0;
	for (int i=0; i<n; i++)
	{
		for (int j=i+1; j<n; j++)
		{
			terms.push_back(term(i, j, d[ij], w[ij]));
            ij += 1;
		}
	}

    // initialize step sizes
    std::vector<double> eta_vec;
    eta_vec.reserve(t_max);
    for (int t=0; t<t_max; t++)
    {
        eta_vec.push_back(eta[t]);
    }

    // perform optimisation
    sgd_horizontal(X, terms, eta_vec);
}

// calculates the unweighted shortest paths between indices I and J
// using a breadth-first search
void bacon(int n, int m, int* I, int* J, int* d_out)
{
    std::vector<std::vector<int>> graph(n);
    for (int ij=0; ij<m; ij++)
        graph[I[ij]].push_back(J[ij]);

    for (int source=0; source<n; source++)
    {
        int depth=0;
        std::vector<bool> visited(n, false);
        std::queue<int> to_visit;
        
        visited[source] = true;
        to_visit.push(source);
        to_visit.push(-1);

        while (1)
        {
            int next = to_visit.front();
            to_visit.pop();
            if (next == -1)
            {
                if (to_visit.empty())
                    break;

                to_visit.push(-1);
                depth++;
            }
            else
            {
                d_out[source*n + next] = depth;
                for (auto it = graph[next].begin(); it != graph[next].end(); ++it)
                {
                    int nextnext = *it;
                    if (visited[nextnext] == false)
                    {
                        to_visit.push(nextnext);
                        visited[nextnext] = true;
                    }
                }
            }
        }
    }
}

