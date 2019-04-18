#include <cmath>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <limits>

#include <chrono>
#include <iostream>

struct term
{
	int i, j;
	double d, w;
	term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}
};

void sgd(double* X, std::vector<term> &terms, const std::vector<double> &etas);
void sgd(double* X, std::vector<term> &terms, const std::vector<double> &etas, double delta);
void sgd_horizontal(double* X, std::vector<term> &terms, const std::vector<double> &etas, double delta); // only changes x coords

std::vector<term> bfs(int n, int m, int* I, int* J);
std::vector<term> dijkstra(int n, int m, int* I, int* J, double* V);

std::vector<double> schedule(const std::vector<term> &terms, int t_max, double eps);
std::vector<double> schedule_convergent(const std::vector<term> &terms, int t_max, double eps, int t_maxmax);

void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps);
void layout_weighted(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps);
void layout_unweighted_convergent(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax);
void layout_weighted_convergent(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax);

// focus is never convergent
void layout_unweighted_focus(int n, double* X, int m, int* I, int* J, int f, int t_max, double eps, int t_maxmax);
void layout_weighted_focus(int n, double* X, int m, int* I, int* J, int f, double* V, int t_max, double eps, int t_maxmax);
// horizontal is always convergent
void layout_unweighted_horizontal(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax);
void layout_weighted_horizontal(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax);

void mds_direct(int n, double* X, double* d, double* w, int t_max, double* etas);


void sgd(double* X, std::vector<term> &terms, const std::vector<double> &etas)
{
	// iterate through step sizes
    int iteration = 0;
	for (double eta : etas)
	{
        // shuffle terms
        std::random_shuffle(terms.begin(), terms.end());

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

			double r = mu * (mag-d_ij) / (2*mag);
			double r_x = r * dx;
			double r_y = r * dy;
			
			X[i*2] -= r_x;
			X[i*2+1] -= r_y;
			X[j*2] += r_x;
			X[j*2+1] += r_y;
		}
        std::cerr << ++iteration << ", eta: " << eta << std::endl;
	}
}


void sgd(double* X, std::vector<term> &terms, const std::vector<double> &etas, double delta)
{
	// iterate through step sizes
    int iteration = 0;
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
        std::cerr << ++iteration << ", eta: " << eta << ", Delta: " << Delta_max << std::endl;
        if (Delta_max < delta)
            return;
	}
}

void sgd_horizontal(double* X, std::vector<term> &terms, const std::vector<double> &etas, double delta)
{
	// iterate through step sizes
    int iteration = 0;
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

            double Delta = mu * (mag-d_ij) / 2.0;
            if (Delta > Delta_max)
                Delta_max = Delta;

			double r = Delta / mag;
			double r_x = r * dx;
			//double r_y = r * dy;
			
			X[i*2] -= r_x;
			//X[i*2+1] -= r_y;
			X[j*2] += r_x;
			//X[j*2+1] += r_y;
		}
        std::cerr << ++iteration << ", eta: " << eta << ", Delta: " << Delta_max << std::endl;
        if (Delta_max < delta)
            return;
	}
}

std::vector<std::vector<int>> build_graph_unweighted(int n, int m, int* I, int* J)
{
    // used to make graph undirected, in case it is not already
    std::vector<std::unordered_set<int>> undirected(n);
    std::vector<std::vector<int>> graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw "i or j bigger than n";

        if (undirected[j].find(i) == undirected[j].end()) // if edge not seen
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
std::vector<term> bfs(int n, int m, int* I, int* J)
{
    auto graph = build_graph_unweighted(n, m, I, J);

	int nC2 = (n*(n-1))/2;
	std::vector<term> terms;
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
                        terms.push_back(term(source, next, d_ij, w_ij));
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

// calculates the unweighted shortest paths between indices I and J
// using dijkstra's algorithm, returning a vector of terms
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

std::vector<std::vector<edge>> build_graph_weighted(int n, int m, int* I, int* J, double* V)
{
    // used to make graph undirected, in case graph is not already
    std::vector<std::unordered_map<int, double>> undirected(n);
    std::vector<std::vector<edge>> graph(n);

    for (int ij=0; ij<m; ij++)
    {
        int i = I[ij], j = J[ij];
        if (i >= n || j >= n)
            throw "i or j bigger than n";

        double v = V[ij];
        if (v <= 0)
            throw "v less or equal 0";

        if (undirected[j].find(i) == undirected[j].end()) // if key not there
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

std::vector<term> dijkstra(int n, int m, int* I, int* J, double* V)
{
    auto graph = build_graph_weighted(n, m, I, J, V);

	int nC2 = (n*(n-1))/2;
	std::vector<term> terms;
	terms.reserve(nC2);

    int terms_size_goal = 0; // to keep track of when to stop searching i<j

    for (int source=0; source<n-1; source++) // no need to do final vertex because i<j
    {
        std::vector<bool> visited(n, false);
        std::vector<double> d(n, std::numeric_limits<double>::max()); // init 'tentative' distances to infinity

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
                    terms.push_back(term(source, current, d_ij, w_ij));
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


std::vector<double> schedule(const std::vector<term> &terms, int t_max, double eps)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (int i=1; i<terms.size(); i++)
    {
        double w = terms[i].w;
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

std::vector<double> schedule_convergent(const std::vector<term> &terms, int t_max, double eps, int t_maxmax)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (int i=1; i<terms.size(); i++)
    {
        double w = terms[i].w;
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
        double eta = eta_switch / (1 + lambda*(t-tau));
        etas.push_back(eta);
    }

    return etas;
}
std::vector<double> schedule_focus(const std::vector<term> &terms, int t_max, double eps, int t_maxmax)
{
    double w_min = terms[0].w, w_max = terms[0].w;
    for (int i=1; i<terms.size(); i++)
    {
        double w = terms[i].w;
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
    for (int t=0; t<t_maxmax; t++)
        etas.push_back(eta_max * exp(-lambda * t));

    return etas;
}

void focus_terms(std::vector<term> &terms, const std::vector<double> &etas, int focus_idx)
{
    double min_eta = etas[0];
    for (int i=1; i<etas.size(); i++)
    {
        if (etas[i] < min_eta)
            min_eta = etas[i];
    }
    double min_w = terms[0].w;
    for (int i=1; i<terms.size(); i++)
    {
        if (terms[i].w < min_w)
            min_w = terms[i].w;
    }

    // this is necessary because we do not know the behaviour of c++ infinity*infinity
    double my_infinity = 1.0 / min_eta / min_w;

    for (term &t: terms)
    {
        if (t.i == focus_idx || t.j == focus_idx)
        {
            t.w = my_infinity;
        }
    }
}

void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps)
{
    try
    {
        //auto start = std::chrono::steady_clock::now();

        std::vector<term> terms = bfs(n, m, I, J);
        //auto end = std::chrono::steady_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << "ms" << std::endl;

        std::vector<double> etas = schedule(terms, t_max, eps);
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
        std::vector<term> terms = dijkstra(n, m, I, J, V);
        std::vector<double> etas = schedule(terms, t_max, eps);
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
        std::vector<term> terms = bfs(n, m, I, J);
        std::vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
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
        std::vector<term> terms = dijkstra(n, m, I, J, V);
        std::vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
        sgd(X, terms, etas, delta);
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}
// focus is never convergent
void layout_unweighted_focus(int n, double* X, int m, int* I, int* J, int f, int t_max, double eps, int t_maxmax)
{
    try
    {
        std::vector<term> terms = bfs(n, m, I, J);
        std::vector<double> etas = schedule_focus(terms, t_max, eps, t_maxmax);
        focus_terms(terms, etas, f);
        sgd(X, terms, etas);
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}
void layout_weighted_focus(int n, double* X, int m, int* I, int* J, int f, double* V, int t_max, double eps, int t_maxmax)
{
    try
    {
        std::vector<term> terms = dijkstra(n, m, I, J, V);
        std::vector<double> etas = schedule_focus(terms, t_max, eps, t_maxmax);
        focus_terms(terms, etas, f);
        sgd(X, terms, etas);
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}
// horizontal is always convergent
void layout_unweighted_horizontal(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax)
{
    try
    {
        std::vector<term> terms = bfs(n, m, I, J);
        std::vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
        sgd_horizontal(X, terms, etas, delta);
    }
    catch (const char* msg)
    {
        std::cerr << "Error: " << msg << std::endl;
    }
}
void layout_weighted_horizontal(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax)
{
    try
    {
        std::vector<term> terms = dijkstra(n, m, I, J, V);
        std::vector<double> etas = schedule_convergent(terms, t_max, eps, t_maxmax);
        sgd_horizontal(X, terms, etas, delta);
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
	std::vector<term> terms;
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
    std::vector<double> etas;
    etas.reserve(t_max);
    for (int t=0; t<t_max; t++)
    {
        etas.push_back(etas[t]);
    }
    
    sgd(X, terms, etas, 0);
}


