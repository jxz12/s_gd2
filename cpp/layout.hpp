#ifndef LAYOUT_HPP
#define LAYOUT_HPP

#include <vector>
#include <ctime>
#include <cstdlib>

std::vector<std::vector<double> > sgd(const std::vector<double> &d, const std::vector<double> &w, const std::vector<double> &eta);

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
