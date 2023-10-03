#include <iostream>
#include "../manseglib.hpp"

int main()
{
	constexpr int size = 1000;
	ManSeg::ManSegArray a(size);
	ManSeg::ManSegArray b(size);
	// double a[size];
	// double b[size];

	std::cout << "filling\n";
	for(int i = 0; i < size; ++i)
	{
		b.pairs[i] = 0;
		a.pairs[i] = i;
	}

	std::cout << "now copying\n";
	for(int i = 0; i < size; ++i)
	{
		for(int j = 0; j < size; ++j) b.pairs[i] += a.pairs[j];
	}
	std::cout << "finished copying\n";

	for(int i = 0; i < size; ++i)
		std::cout << b.pairs[i] << " ";
	std::cout << std::endl;

	return 0;
}
