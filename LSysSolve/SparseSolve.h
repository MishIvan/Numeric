#pragma once
#include <vector>
#include <algorithm>
#include <iterator>

struct SparseElement
{
	int row;
	int column;
	double value;
};

double FullnessDegree(const std::vector<SparseElement> matrix);