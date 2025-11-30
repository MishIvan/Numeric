#pragma once
#include "VECTOR.h"
/// <summary>
/// Класс "матрица" матрица размера M x N
/// </summary>
class MATRIX
{
protected:
	const int MIN_SIZE_FOR_THREAD = 100; // минимальный размер матрицы, начиная с которого запускаются потоки для вычисления обратной матрицы
	int m_rows; //  число строк
	int m_columns; // число столбцов
	double* m_data; // данные матрицы
	double FormMatrixCompactScheme(MATRIX& alpha);
public:
	MATRIX(int M, int N, double val = 0.0);
	MATRIX(const MATRIX& src);
	inline int rows() { return m_rows; }
	inline int columns() { return m_columns; }

	MATRIX& operator=(const MATRIX& src);
	inline double& operator()(int i, int j) { return *(m_data + i * m_columns + j); };
	friend MATRIX operator*(const MATRIX& matr1, const MATRIX& matr2);
	friend ostream& operator<<(ostream& s, MATRIX& matr);
	friend istream& operator>>(istream& s, MATRIX& matr);
	friend VECTOR operator*(const MATRIX& matr, const VECTOR& v);
	MATRIX Transpose();
	MATRIX Reverse();
	MATRIX Invert();

	double Determinant();
	double Minor(int i, int j);
	bool IsSymmetric();

	static bool readFromFile(const char* fileName, MATRIX& matr);
	static bool writeToFile(const char* fileName, MATRIX& matr);
	
	friend bool Gauss(const MATRIX& a, const VECTOR& b, VECTOR& x);
	friend void CompactSchemeSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void QRDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void LLTDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void TriangleSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void CopyColumn(MATRIX& A, VECTOR& v, int j);
	

	bool QRDecomposition(MATRIX& Q, MATRIX& R);
	bool CholeskyDecomposition(MATRIX& L);

	~MATRIX();
};

class MATRIXEXT : public MATRIX
{
public:
	MATRIXEXT(int M, int N) : MATRIX(M, N) {};
	int getCountNotNumsUnderMD(double val);
	
};

