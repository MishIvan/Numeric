#pragma once
#include <iostream> 
#include <fstream>
#include <math.h>

using namespace std;
class MATRIX;

/// <summary>
/// Класс вектор - n-мерный вектор
/// </summary>
class VECTOR
{
	double* m_data;
	int m_size;
public:
	VECTOR(int n, double val = 0.0);
	VECTOR(const VECTOR& src);
	inline int size() { return m_size; } // возврат размерности вектора
	double norm();
	void Reverse();

	VECTOR& operator=(const VECTOR& src);
	inline const double* data() const { return m_data; };
	friend double operator*(const VECTOR& v1, const VECTOR& v2);
	friend VECTOR operator*(double a, const VECTOR& v);
	friend VECTOR operator+(const VECTOR& v1, const VECTOR& v2);
	friend VECTOR operator-(const VECTOR& v1, const VECTOR& v2);
	friend ostream& operator<<(ostream& s, VECTOR& v);
	friend istream& operator>>(istream& s, VECTOR& v);
	inline double& operator[](int i) { return *(m_data + i); }
	static bool readFromFile(const char* fileName, VECTOR& vect);
	static bool writeToFile(const char* fileName, VECTOR& vect);
	friend VECTOR operator*(const MATRIX& matr, const VECTOR& v);

	friend bool Gauss(const MATRIX& a, const VECTOR& b, VECTOR& x);
	friend void CompactSchemeSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void QRDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void LLTDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void TriangleSolve(MATRIX& A, VECTOR& b, VECTOR& x);
	friend void CopyColumn(MATRIX& A, VECTOR& v, int j);

	~VECTOR();
};

