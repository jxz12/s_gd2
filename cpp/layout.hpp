#ifndef LAYOUT_HPP
#define LAYOUT_HPP

#include <vector>
#include <ctime>
#include <cstdlib>

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

#endif
