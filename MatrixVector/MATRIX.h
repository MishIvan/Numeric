#pragma once
#include "VECTOR.h"
#include <complex>
#include <vector>
#include <cmath>
#include <thread>
/// <summary>
/// Класс "матрица" матрица размера M x N
/// </summary>
template <typename T>
class MATRIX
{
private:
	const int MIN_SIZE_FOR_THREAD = 100; // минимальный размер матрицы, начиная с которого запускаются потоки для вычисления обратной матрицы
	int m_rows; //  число строк
	int m_columns; // число столбцов
	T* m_data; // данные матрицы
	double FormMatrixCompactScheme(MATRIX& alpha);
public:
	MATRIX(int M, int N);
	MATRIX(const MATRIX& src);
	inline int rows() { return m_rows; }
	inline int columns() { return m_columns; }

	MATRIX& operator=(const MATRIX& src);
	inline T& operator()(int i, int j) { return *(m_data + i * m_columns + j); };
	friend MATRIX<T> operator*(const MATRIX<T>& matr1, const MATRIX<T>& matr2);
	friend ostream& operator<<(ostream& s, MATRIX<T>& matr);
	friend istream& operator>>(istream& s, MATRIX<T>& matr);
	friend VECTOR<T> operator*(const MATRIX<T>& matr, const VECTOR<T>& v);
	friend MATRIX<T> operator*(T alf, const MATRIX<T>& matr);
	friend MATRIX<T> operator+(const MATRIX<T>& matr1, const MATRIX<T>& matr2);
	friend MATRIX<T> operator-(const MATRIX<T>& matr1, const MATRIX<T>& matr2);
	MATRIX& operator/=(T val);
	MATRIX& operator*=(const MATRIX& matr);

	MATRIX Transpose();
	MATRIX Reverse();
	MATRIX Invert();
	MATRIX InvertFaddev();
	T Determinant();
	T Sp();
	T Minor(int i, int j);
	bool IsSymmetric();

	static bool readFromFile(const char* fileName, MATRIX& matr);
	static bool writeToFile(const char* fileName, MATRIX& matr);
	
	friend bool Gauss(const MATRIX<T>& a, const VECTOR<T>& b, VECTOR<T>& x);
	friend void CompactSchemeSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	friend void QRDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	friend void LLTDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	friend void TriangleSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

	void CopyColumn(VECTOR<T>& v, int j);
	VECTOR<T> CopyColumn2Vector(int j);
	
	bool QRDecomposition(MATRIX& Q, MATRIX& R);
	bool CholeskyDecomposition(MATRIX& L);

	void EigenvaluesAndVectorsKrylov(complex<T>* lambda, complex<T>** vect);
	void EigenvaluesAndVectorsLeVerrierFaddeev(complex<T>* lambda, complex<T>** vect);
	void Rotate(T* lambda, T** vect);

	~MATRIX();
};


