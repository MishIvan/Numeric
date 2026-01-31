#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

#define MAX_ITERATION_NUMBER 60000
#define EPS 1.0e-12

// элемент матрицы
struct SparseElement
{
	int row; // строка
	int column; // столбец
	double value; // значение
};

void PrintMatrix(std::vector<SparseElement> matrix, int n);
std::vector<SparseElement> SparseTranspose(const std::vector<SparseElement>& matrix);
std::vector <SparseElement> SparseMultiply(const std::vector<SparseElement>& first,
	const std::vector<SparseElement>& second);
double ErrorMeasure(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	const std::vector<SparseElement>& x);


int SparseRotationSolve(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	std::vector<SparseElement>& x, int n);
int SparseRelaxation(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	std::vector<SparseElement>& x,
	double omega);
int SparseGradientDescent(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	std::vector<SparseElement>& x);

