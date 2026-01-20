#pragma once
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

struct SparseElement
{
	int row;
	int column;
	double value;
};

void PrintMatrix(std::vector<SparseElement> matrix, int n);
int SparseRotationSolve(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	std::vector<SparseElement>& x, int n);
double ErrorMeasure(std::vector<SparseElement>& A,
	std::vector<SparseElement>& b,
	std::vector<SparseElement>& x, int n);
