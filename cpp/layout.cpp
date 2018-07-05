#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "layout.hpp"

struct term {
	term(int i, int j, int d, int w) : i(i), j(j), d(d), w(w) {}
	int i, j;
	double d, w;
};

std::vector<std::vector<double> > sgd(const std::vector<double> &d, const std::vector<double> &w, const std::vector<double> &eta)
{
	// check whether the input condensed matrix is the correct length
	int nC2 = d.size();
	if (w.size() != nC2)
		throw "d and w do not have the same length";

	double n = (1 + sqrt((double)(1+8*nC2))) / 2;
	if (n != floor(n))
		throw "d is not the correct length to be a distance matrix";

	// initialize SGD
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

	// initialize positions
	srand(time(NULL));
	std::vector<std::vector<double> > X(n);
	for (int i=0; i<n; i++)
	{
		X[i].resize(2);
		// initialize randomly within 1x1 square
		X[i][0] = (double)rand() / RAND_MAX;
		X[i][1] = (double)rand() / RAND_MAX;
	}
 
	// perform SGD
	for (int k=0; k<eta.size(); k++)
	{
		fisher_yates<term>(terms);
		for (int t=0; t<nC2; t++)
		{
			// cap step size
			double mu = eta[k] * terms[t].w;
			if (mu > 1)
				mu = 1;

			double d_ij = terms[t].d;

			int i = terms[t].i, j = terms[t].j;
			double del_x = X[i][0] - X[j][0], del_y = X[i][1] - X[j][1];
			double mag = sqrt(del_x*del_x + del_y*del_y);
			
			double r = mu * (mag - d_ij) / (2 * mag);
			double r_x = r * del_x;
			double r_y = r * del_y;
			
			X[i][0] -= r_x;
			X[i][1] -= r_y;
			X[j][0] += r_x;
			X[j][1] += r_y;
		}
	}
	return X;
}

