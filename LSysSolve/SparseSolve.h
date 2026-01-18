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

void PrintMatrix(std::vector<SparseElement> matrix);
double FullnessDegree(const std::vector<SparseElement> matrix);
int SparseRotationSolve(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	std::vector<SparseElement>& x);
double ErrorMeasure(std::vector<SparseElement>& A,
	std::vector<SparseElement>& b,
	std::vector<SparseElement>& x);
