#pragma once
#include <map>
#include <algorithm>
#define MAX_ITERATION_NUMBER 60000

double SparseNormv(const std::map<int, double>& vect);
double ErrorMeasure(const std::map<int, std::map<int, double>>& A,
	const std::map<int, double>& b,
	const std::map<int, double>& x);

int SparseRotationSolve(const std::map<int, std::map<int, double>>& A,
	const std::map<int, double>& b,
	std::map<int, double>& x,
	int n);
int SparseRelaxation(const std::map<int, std::map<int, double>>& A,
	const std::map<int, double>& b,
	std::map<int, double>& x,
	double omega);
int SparseGradientDescent(const std::map<int, std::map<int, double>>& A,
	const std::map<int, double>& b,
	std::map<int, double>& x);
