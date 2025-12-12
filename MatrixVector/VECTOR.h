#pragma once
#include <iostream> 
#include <fstream>
#include <complex>
#include <math.h>

using namespace std;

template <typename T>
class MATRIX;

/// <summary>
/// Класс вектор - n-мерный вектор
/// </summary>
template <typename T>
class VECTOR
{
	T* m_data;
	int m_size;
public:
	VECTOR(int n);
	VECTOR(const VECTOR<T>& src);
	inline int size() { return m_size; } // возврат размерности вектора
	double norm();
	void Reverse();

	VECTOR& operator=(const VECTOR& src);
	friend double operator*(const VECTOR& v1, const VECTOR& v2);
	friend VECTOR operator*(T a, const VECTOR& v);
	friend VECTOR operator+(const VECTOR& v1, const VECTOR& v2);
	friend VECTOR operator-(const VECTOR& v1, const VECTOR& v2);
	friend ostream& operator<<(ostream& s, VECTOR& v);
	friend istream& operator>>(istream& s, VECTOR& v);
	inline T& operator[](int i) { return *(m_data + i); }
	static bool readFromFile(const char* fileName, VECTOR<T>& vect);
	static bool writeToFile(const char* fileName, VECTOR<T>& vect);

	friend bool Gauss(const MATRIX<T>& a, const VECTOR<T>& b, VECTOR<T>& x);
	friend void CompactSchemeSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	friend void QRDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	friend void LLTDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	friend void TriangleSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

	~VECTOR();
};

