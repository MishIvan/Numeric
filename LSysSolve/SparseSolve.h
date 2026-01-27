#pragma once
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

// элемент матрицы
struct SparseElement
{
	int row; // строка
	int column; // столбец
	double value; // значение
};

double FindElement(const std::vector<SparseElement>& matrix, int _row, int _column = 1);
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
int SparseRelaxation(std::vector<SparseElement>& A,
	std::vector<SparseElement>& b,
	std::vector<SparseElement>& x,
	int n, double omega);
int SparseGradientDescent(const std::vector<SparseElement>& A,
	const std::vector<SparseElement>& b,
	std::vector<SparseElement>& x);

