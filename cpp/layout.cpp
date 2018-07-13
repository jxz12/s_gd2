#include <vector>
#include <ctime>
#include <cstdlib>
#include <cmath>

void sgd(int n, double* X, double* d, double* w, int t_max, double* eta);

template <class T>
void fisher_yates(std::vector<T> &list)
{
	// srand(time(NULL)); // assume it has been seeded
	int n = list.size();
	for (int i=n-1; i>0; i--)
	{
		int j = rand() % (i+1);
		T temp = list[i];
		list[i] = list[j];
		list[j] = temp;
	}
}

struct term {
	term(int i, int j, double d, double w) : i(i), j(j), d(d), w(w) {}
	int i, j;
	double d, w;
};

void sgd(int n, double* X, double* d, double* w, int t_max, double* eta)
{
	int nC2 = (n*(n-1))/2;

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
 
	// perform SGD
	for (int t=0; t<t_max; t++)
	{
		fisher_yates<term>(terms);
		for (int ij=0; ij<nC2; ij++)
		{
			// cap step size
			double mu = eta[t] * terms[ij].w;
			if (mu > 1)
				mu = 1;

			double d_ij = terms[ij].d;

			int i = terms[ij].i, j = terms[ij].j;
			double del_x = X[i*2] - X[j*2], del_y = X[i*2+1] - X[j*2+1];
			double mag = sqrt(del_x*del_x + del_y*del_y);
			
			double r = mu * (mag - d_ij) / (2 * mag);
			double r_x = r * del_x;
			double r_y = r * del_y;
			
			// X[i][0] -= r_x;
			// X[i][1] -= r_y;
			// X[j][0] += r_x;
			// X[j][1] += r_y;
			X[i*2] -= r_x;
			X[i*2+1] -= r_y;
			X[j*2] += r_x;
			X[j*2+1] += r_y;
		}
	}
}

