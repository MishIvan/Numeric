#pragma once
#include <iostream> 
#include <fstream>
#include <complex>
#include <math.h>
#include <vector>
#include <thread>

#define EPS 1.0e-12
#define MAX_ITER_NUMBER 30000


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
	T norm();
	void Reverse();

	VECTOR& operator=(const VECTOR& src);
	template <typename T> friend T operator*(const VECTOR<T>& v1, const VECTOR<T>& v2);
	template <typename T> friend VECTOR<T> operator*(T a, const VECTOR<T>& v);

	template <typename T> friend VECTOR<T> operator+(const VECTOR<T>& v1, const VECTOR<T>& v2);
	template <typename T> friend VECTOR<T> operator-(const VECTOR<T>& v1, const VECTOR<T>& v2);
	template <typename T> friend ostream& operator<<(ostream& s, VECTOR<T>& v);
	template <typename T> friend istream& operator>>(istream& s, VECTOR<T>& v);
	inline T& operator[](int i) { return *(m_data + i); }
	static bool readFromFile(const char* fileName, VECTOR<T>& vect);
	static bool writeToFile(const char* fileName, VECTOR<T>& vect);

	template <typename T> friend bool Gauss(const MATRIX<T>& a, const VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void CompactSchemeSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void QRDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void LLTDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void TriangleSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void LUDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

	~VECTOR();
};

/// <summary>
/// Конструктор вектора: выделение памяти под указатель и инициализация компонент вектора значениями val
/// </summary>
/// <param name="n">размерность вектора</param>
template <typename T>
VECTOR<T>::VECTOR(int n)
{
	m_data = 0; m_size = 0;
	if (n < 1)
		return;
	m_size = n;
	m_data = new T[m_size];
	memset(m_data, 0, n*sizeof(T));

}
/// <summary>
/// Конструктор копирования вектора src
/// </summary>
/// <param name="src">исходный вектор для копирования</param>
template <typename T>
VECTOR<T>::VECTOR(const VECTOR& src)
{
	if (!this->m_data) delete[] this->m_data;
	this->m_data = new T[src.m_size];
	this->m_size = src.m_size;
	for (int i = 0; i < this->m_size; i++)
		*(this->m_data + i) = *(src.m_data + i);
}

/// <summary>
/// Возвращает евклидову норму вектора
/// </summary>
/// <returns>yорма вектора</returns>
template <typename T>
T VECTOR<T>::norm()
{
	T val = 0;
	for (int i = 0; i < m_size; i++)
		val += abs(*(m_data + i)) * abs(*(m_data + i));
	return sqrt(val);
}

/// <summary>
/// Оператор присвоения
/// </summary>
/// <param name="src">исходный вестор для присвоения</param>
/// <returns></returns>
template <typename T>
VECTOR<T>& VECTOR<T>::operator=(const VECTOR& src)
{
	if (!this->m_data) delete[] this->m_data;
	this->m_data = new T[src.m_size];
	this->m_size = src.m_size;
	for (int i = 0; i < this->m_size; i++)
		*(this->m_data + i) = *(src.m_data + i);
	return *this;
}
/// <summary>
/// Возвращает скалярное произведение векторов v1 и v2
/// Вектора должны быть одной размерности
/// </summary>
/// <param name="v1">первый вектор</param>
/// <param name="v2">второй вектор</param>
/// <returns>значение скалярного произведения</returns>
template <typename T>
T operator*(const VECTOR<T>& v1, const VECTOR<T>& v2)
{
	T prod = 0.0;
	if (v1.m_size != v2.m_size) return prod;
	for (int i = 0; i < v1.m_size; i++)
		prod += *(v1.m_data + i) * (v2.m_data + i);
	return prod;
}

/// <summary>
/// Перегрузка оператора умножения ветора на вещественное число
/// размерности векторов должны совпадать
/// </summary>
/// <param name="v">вектор</param>
/// <param name="a">вещественное число</param>
/// <returns>вектор - результат умножения </returns>
template <typename T>
VECTOR<T> operator*(T a, const VECTOR<T>& v)
{
	VECTOR<T> res(v.m_size);
	for (int i = 0; i < res.m_size; i++)
		*(res.m_data + i) = a * *(v.m_data + i);
	return res;
}

/// <summary>
/// Перегрузка операторов вычитания и сложения векторов v1 и v2
/// размерности векторов должны совпадать
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns>вектор - разность (сумму) векторов  </returns>
template <typename T>
VECTOR<T> operator-(const VECTOR<T>& v1, const VECTOR<T>& v2)
{
	VECTOR<T> res(v1.m_size);
	if (v1.m_size != v2.m_size)
	{
		throw "Векторы имеют разную размерность";
		return res;
	}
	for (int i = 0; i < res.m_size; i++)
		*(res.m_data + i) = *(v1.m_data + i) - *(v2.m_data + i);
	return res;
}
template <typename T>
VECTOR<T> operator+(const VECTOR<T>& v1, const VECTOR<T>& v2)
{
	VECTOR<T> res(v1.m_size);
	if (v1.m_size != v2.m_size)
	{
		throw "Векторы имеют разную размерность";
		return res;
	}
	for (int i = 0; i < res.m_size; i++)
		*(res.m_data + i) = *(v1.m_data + i) + *(v2.m_data + i);
	return res;
}

/// <summary>
/// Перегрузка оператора вывода на консоль или в файловый поток (файл нужно будет открыть для чтения)
/// </summary>
/// <param name="s"></param>
/// <param name="v">вектор</param>
/// <returns></returns>
template <typename T>
ostream& operator<<(ostream& s, VECTOR<T>& v)
{
	int n = v.m_size;
	for (int i = 0; i < n; i++)
		s << *(v.m_data + i) << ' ';
	return s;
}
/// <summary>
/// Перегрузка оператора ввода с консоли или в файловый поток (файл должен быть открыт для записи) 
/// </summary>
/// <param name="s"></param>
/// <param name="v">вектор</param>
/// <returns></returns>
template <typename T>
istream& operator>>(istream& s, VECTOR<T>& v)
{
	int n = v.m_size;
	for (int i = 0; i < n; i++)
		s >> *(v.m_data + i);
	return s;
}

/// <summary>
/// Расставляет элементы ветора в обратном порядке: последний становится первым, предпоследний вторым и т.д.
/// </summary>
template <typename T>
void VECTOR<T>::Reverse() {
	int size = m_size;
	for (int start = 0, end = size - 1; start < end; ++start, --end)
	{
		// Меняем местами первый и последний элемент
		swap(*(m_data + start), *(m_data + end));
	}
}

/// <summary>
/// Чтение вектора из текстового файла
/// </summary>
/// <param name="fileName">полное имя файла</param>
/// <param name="vect">вектор, создаваемый по данным файла</param>
/// <returns>true - в случае успешного считывания данных, false - в случае ошибки </returns>
template <typename T>
bool VECTOR<T>::readFromFile(const char* fileName, VECTOR<T>& vect)
{
	ifstream fs;
	fs.open(fileName);
	if (fs.is_open())
	{
		int n = 0;
		fs >> n;
		if (n < 1) return false;
		VECTOR v_out(n);
		fs >> v_out;
		vect = v_out;
		fs.close();
		return true;
	}
	return false;
}
/// <summary>
/// Запись ветора в текстовый файл
/// </summary>
/// <param name="fileName">полное имя текстового файла</param>
/// <returns>true - в случае успешного считывания данных, false - в случае ошибки</returns>
template <typename T>
bool VECTOR<T>::writeToFile(const char* fileName, VECTOR<T>& vect)
{
	ofstream fs;
	fs.open(fileName);
	if (fs.is_open())
	{
		int n = 0;
		fs << vect.m_size << endl;
		if (vect.m_size < 1) return false;
		fs << vect;
		fs.close();
		return true;
	}
	return false;

}
/// <summary>
/// Деструктор. Освобождение памяти
/// </summary>
template <typename T>
VECTOR<T>::~VECTOR()
{
	if (m_data) delete[] m_data;
}

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
	T FormMatrixCompactScheme(MATRIX& alpha);
public:
	MATRIX(int M, int N);
	MATRIX(const MATRIX& src);
	inline int rows() { return m_rows; }
	inline int columns() { return m_columns; }

	MATRIX& operator=(const MATRIX& src);
	inline T& operator()(int i, int j) { return *(m_data + i * m_columns + j); };
	template <typename T>friend MATRIX<T> operator*(const MATRIX<T>& matr1, const MATRIX<T>& matr2);
	template <typename T> friend ostream& operator<<(ostream& s, MATRIX<T>& matr);
	template <typename T> friend istream& operator>>(istream& s, MATRIX<T>& matr);
	template <typename T> friend VECTOR<T> operator*(const MATRIX<T>& matr, VECTOR<T>& v);
	template <typename T> friend MATRIX<T> operator*(T alf, const MATRIX<T>& matr);
	template <typename T> friend MATRIX<T> operator+(const MATRIX<T>& matr1, const MATRIX<T>& matr2);
	template <typename T> friend MATRIX<T> operator-(const MATRIX<T>& matr1, const MATRIX<T>& matr2);
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

	template <typename T> friend bool Gauss(const MATRIX<T>& a, const VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void CompactSchemeSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void QRDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void LLTDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void TriangleSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);
	template <typename T> friend void LUDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

	void CopyColumn(VECTOR<T>& v, int j);
	VECTOR<T> CopyColumn2Vector(int j);

	bool QRDecomposition(MATRIX& Q, MATRIX& R);
	bool CholeskyDecomposition(MATRIX& L);
	T LUDecomposition(MATRIX& alfa);

	void EigenvaluesAndVectorsKrylov(complex<T>* lambda, complex<T>** vect);
	void EigenvaluesAndVectorsLeVerrierFaddeev(complex<T>* lambda, complex<T>** vect);
	void Rotate(T* lambda, T** vect);

	~MATRIX();
};

/// <summary>
/// Конструктор матрицы M x N, всем элементам которой присваивается значение val
/// </summary>
/// <param name="M">число строк матрицы</param>
/// <param name="N">число столбцов матрицы</param>
/// <param name="val"></param>
template <typename T>
MATRIX<T>::MATRIX(int M, int N)
{
	m_columns = m_rows = 0;
	m_data = nullptr;
	if (M < 0 || N < 0)
	{
		throw invalid_argument("Задана отрицательная размерность матрицы");
		return;
	}
	m_rows = M; m_columns = N;
	this->m_data = new T[this->m_rows * this->m_columns];
	memset(m_data, 0, this->m_rows*this->m_columns*sizeof(T));
}
/// <summary>
/// Конструктор копирования матрицы из матрицы src 
/// </summary>
/// <param name="src">исходная матрица для копирования</param>
template <typename T>
MATRIX<T>::MATRIX(const MATRIX& src)
{
	if (!this->m_data) delete[] this->m_data;
	this->m_columns = src.m_columns;
	this->m_rows = src.m_rows;
	this->m_data = new T[this->m_rows * this->m_columns];
	for (int i = 0; i < this->m_rows; i++)
		for (int j = 0; j < this->m_columns; j++)
			*(this->m_data + i * this->m_columns + j) = *(src.m_data + i * src.m_columns + j);
}
/// <summary>
/// Перегрузка оператора присваивания
/// </summary>
/// <param name="src">матрица для присваивания</param>
/// <returns></returns>
template <typename T>
MATRIX<T>& MATRIX<T>::operator=(const MATRIX& src)
{
	if (this->m_rows != src.m_rows || this->m_columns != src.m_columns)
	{
		if (!this->m_data) delete[] this->m_data;
		this->m_columns = src.m_columns;
		this->m_rows = src.m_rows;
		this->m_data = new T[this->m_rows * this->m_columns];
	}
	for (int i = 0; i < this->m_rows; i++)
		for (int j = 0; j < this->m_columns; j++)
			*(this->m_data + i * this->m_columns + j) = *(src.m_data + i * src.m_columns + j);
	return *this;
}

/// <summary>
/// Перегрузка оператора умножения матриц matr1 и matr2
/// Число столбцов матрицы matr1 должно быть равно числу строк матрицы matr2
/// </summary>
/// <param name="matr1">первая матрица</param>
/// <param name="matr2">вторая матрица</param>
/// <returns>матрицу, произведение</returns>
template <typename T>
MATRIX<T> operator*(const MATRIX<T>& matr1, const MATRIX<T>& matr2)
{
	MATRIX<T> pmatr(matr1.m_rows, matr2.m_columns);
	if (matr1.m_columns != matr2.m_rows)
	{
		throw("Несоответствие размерности матриц");
		return pmatr;
	}
	int i, j;
	for (i = 0; i < matr1.m_rows; i++)
		for (j = 0; j < matr2.m_columns; j++)
		{
			T prod = 0;
			for (int k = 0; k < matr1.m_columns; k++)
			{
				if (*(matr1.m_data + i * matr1.m_columns + k) == 0 || *(matr2.m_data + k * matr2.m_columns + j) == 0)
					continue;
				prod += *(matr1.m_data + i * matr1.m_columns + k) * *(matr2.m_data + k * matr2.m_columns + j);
			}
			*(pmatr.m_data + i * pmatr.m_columns + j) = prod;
		}
	return pmatr;
}
/// <summary>
/// Перегрузка оператора вывода на консоль или в файловый поток (файл должен быть открыт для чтения) 
/// матрица выводится построчно
/// </summary>
/// <param name="s"></param>
/// <param name="matr">матрица</param>
/// <returns></returns>
template <typename T>
ostream& operator<<(ostream& s, MATRIX<T>& matr)
{
	int m = matr.m_rows;
	int n = matr.m_columns;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) // вывод строки в поток
			s << *(matr.m_data + i * n + j) << ' ';
		s << endl;
	}
	return s;
}
/// <summary>
/// Перегрузка оператора ввода с консоли или в файловый поток (файл должен быть открыт для записи) 
/// матрицы вводится построчно
/// </summary>
/// <param name="s"></param>
/// <param name="matr">матрица</param>
/// <returns>поток</returns>
template <typename T>
istream& operator>>(istream& s, MATRIX<T>& matr)
{
	int m = matr.m_rows;
	int n = matr.m_columns;
	T* buff = new T[n];
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++) // считывание строки матрицы
			s >> *(buff + j);
		for (int j = 0; j < n; j++) // запись строки матрицы
			*(matr.m_data + i * n + j) = *(buff + j);

	}
	delete[] buff;
	return s;
}
/// <summary>
/// Вычисление определителя квадратной матрицы
/// </summary>
/// <returns></returns>
template <typename T>
T MATRIX<T>::Determinant()
{
	if (m_rows != m_columns) return NAN;
	MATRIX alpha(m_rows, m_columns);
	return FormMatrixCompactScheme(alpha);
}
/// <summary>
/// Считывание данных их файла и загрузка их в матрицу matr
/// </summary>
/// <param name="fileName">полное имя файла</param>
/// <param name="matr">матрица</param>
/// <returns>true - в случае успешного считывания данных, false - в случае ошибки</returns>
template <typename T>
bool MATRIX<T>::readFromFile(const char* fileName, MATRIX& matr)
{
	ifstream fs;
	fs.open(fileName);
	if (fs.is_open())
	{
		fs >> matr.m_rows >> matr.m_columns;
		fs >> matr;
		fs.close();
		return true;
	}
	return false;
}
/// <summary>
///  Запись данных в файл
/// </summary>
/// <param name="fileName">полное имя файла для записи</param>
/// <returns>true - в случае успешного считывания данных, false - в случае ошибки</returns>
template <typename T>
bool MATRIX<T>::writeToFile(const char* fileName, MATRIX& matr)
{
	ofstream fs;
	fs.open(fileName);
	if (fs.is_open()) {
		fs << matr.m_rows << ' ' << matr.m_columns << endl;
		fs << matr;
		fs.close();
		return true;
	}
	return false;
}
/// <summary>
/// Является ли квадратная матрица симметричной
/// </summary>
/// <returns></returns>
template <typename T>
bool MATRIX<T>::IsSymmetric()
{
	if (m_rows != m_columns)
	{
		throw "Симметричной может быть только квадратная матрица";
		return false;
	}
	int n = m_rows;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (*(m_data + i * n + j) != *(m_data + j * n + i)) return false;
		}
	return true;
}
/// <summary>
/// Вычислить след матрицы, сумму её колонок по диагонали
/// </summary>
/// <returns>значение следа матрицы</returns>
template <typename T>
T MATRIX<T>::Sp()
{
	T val = 0.0;
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_columns; j++)
			if (i == j) val += *(m_data + i * m_columns + j);
	return val;
}

/// <summary>
/// Сложение мариц matr1 и matr2
/// </summary>
/// <param name="matr1">матрица</param>
/// <param name="matr2">матрица</param>
/// <returns>матрица - результат сложения</returns>
template <typename T>
MATRIX<T> operator+(const MATRIX<T>& matr1, const MATRIX<T>& matr2)
{
	MATRIX<T> prod(matr1.m_rows, matr1.m_columns);
	if (matr1.m_columns != matr2.m_columns && matr1.m_rows != matr2.m_rows)
	{
		throw("Размерности матриц не совпадают");
		return prod;
	}
	for (int i = 0; i < matr1.m_rows; i++)
		for (int j = 0; j < matr1.m_columns; j++)
			*(prod.m_data + i * prod.m_columns + j) = *(matr1.m_data + i * matr1.m_columns + j) + *(matr2.m_data + i * matr2.m_columns + j);
	return prod;
}

/// <summary>
/// Вычитание марицы matr2 из матрицы matr1
/// </summary>
/// <param name="matr1">матрица</param>
/// <param name="matr2">матрица</param>
/// <returns>матрица - результат вычитания</returns>
template <typename T>
MATRIX<T> operator-(const MATRIX<T>& matr1, const MATRIX<T>& matr2)
{
	MATRIX<T> prod(matr1.m_rows, matr1.m_columns);
	if (matr1.m_columns != matr2.m_columns && matr1.m_rows != matr2.m_rows)
	{
		throw("Размерности матриц не совпадают");
		return prod;
	}
	for (int i = 0; i < matr1.m_rows; i++)
		for (int j = 0; j < matr1.m_columns; j++)
			*(prod.m_data + i * prod.m_columns + j) = *(matr1.m_data + i * matr1.m_columns + j) - *(matr2.m_data + i * matr2.m_columns + j);
	return prod;
}


/// <summary>
/// Умножение марицы matr на вектор v
/// </summary>
/// <param name="matr">матрица</param>
/// <param name="v">вектор</param>
/// <returns>результат умножения, вектор</returns>
template <typename T>
VECTOR<T> operator*(const MATRIX<T>& matr, VECTOR<T>& v)
{
	VECTOR<T> prod(matr.m_rows);
	if (matr.m_columns != v.size()) return prod;
	for (int i = 0; i < matr.m_rows; i++)
	{
		T val = 0;
		for (int j = 0; j < matr.m_columns; j++)
		{
			if (*(matr.m_data + i * matr.m_columns + j) == 0.0 || v[j] == 0.0) continue;
			val += *(matr.m_data + i * matr.m_columns + j) * v[j];
		}
		prod[i] = val;
	}
	return prod;
}

/// <summary>
/// Умножение матрицы на матрицу matr справа. Результат записывается в матрицу
/// </summary>
/// <param name="matr">матрица, на которую умножают справа</param>
/// <returns>эта матрица как результат умножения</returns>
template <typename T>
MATRIX<T>& MATRIX<T>::operator*=(const MATRIX& matr)
{
	if (m_columns != matr.m_rows)
	{
		throw("Несоответсвие размерности матриц");
		return *this;
	}
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_columns; j++)
		{
			T val = 0;
			for (int k = 0; k < m_columns; k++)
			{
				if (*(m_data + i * m_columns + k) == 0.0 || *(matr.m_data + k * matr.m_columns + j) == 0.0)
					continue;
				val += *(m_data + i * m_columns + k) * *(matr.m_data + k * matr.m_columns + j);
			}
			*(m_data + i * m_columns + j) = val;
		}
	return *this;
}

/// <summary>
/// Деление элементов матрицы на число val
/// </summary>
/// <param name="val">число, на которое делятся элементы матрицы</param>
/// <returns>матрица, элементы которой поделены на число val</returns>
template <typename T>
MATRIX<T>& MATRIX<T>::operator/=(T val)
{
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_columns; j++)
			*(m_data + i * m_columns + j) /= val;
	return *this;
}
/// <summary>
/// Умножение марицы matr на скаляр alf
/// </summary>
/// <param name="matr">матрица</param>
/// <param name="alf"></param>
/// <returns>результат умножения, матрица</returns>
template <typename T>
MATRIX<T> operator*(T alf, const MATRIX<T>& matr)
{
	MATRIX<T> prod(matr.m_rows, matr.m_columns);
	for (int i = 0; i < matr.m_rows; i++)
		for (int j = 0; j < matr.m_columns; j++)
			*(prod.m_data + i * prod.m_columns + j) += *(matr.m_data + i * matr.m_columns + j) * alf;
	return prod;
}

/// <summary>
/// QR разложение квадратной матрицы
/// </summary>
/// <param name="Q">матрица Q</param>
/// <param name="R">матрица R</param>
template <typename T>
bool MATRIX<T>::QRDecomposition(MATRIX& Q, MATRIX& R)
{
	if (m_columns != m_rows) return false;
	int n = m_rows;
	T sum = 0;

	for (int j = 0; j < n; j++)
	{
		// q(j) = a(j)
		for (int k = 0; k < n; k++)
			*(Q.m_data + k * n + j) = *(m_data + k * n + j);

		for (int i = 0; i <= j - 1; i++)
		{
			// rij = q(i)^T*a(j)
			sum = 0.0;
			for (int k = 0; k < n; k++)
				sum += *(Q.m_data + k * n + i) * *(m_data + k * n + j);
			*(R.m_data + i * n + j) = sum;

			// r(i,j)*q(i)
			for (int k = 0; k < n; k++)
				*(Q.m_data + k * n + j) -= *(R.m_data + i * n + j) * *(Q.m_data + k * n + i);
		}

		// r(j,j) = || q(j)) ||2
		sum = 0.0;
		for (int k = 0; k < n; k++)
			sum += *(Q.m_data + k * n + j) * *(Q.m_data + k * n + j);
		*(R.m_data + j * n + j) = sqrt(sum);

		if (*(R.m_data + j * n + j) == 0.0)
		{
			return false;
		}
		for (int k = 0; k < n; k++)
			*(Q.m_data + k * n + j) /= *(R.m_data + j * n + j);

	}

	return true;
}
/// <summary>
///  Разложение Холецкого матрицы A = L*L^t
/// </summary>
/// <param name="L">матрица в разложении Холецкого</param>
/// <returns></returns>
template <typename T>
bool MATRIX<T>::CholeskyDecomposition(MATRIX& L)
{
	if (m_rows != m_columns) return false;
	int n = m_rows;
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++) {

			T sum = 0;
			for (int k = 0; k < j; k++)
				sum += *(L.m_data + i * n + k) * *(L.m_data + j * n + k);

			if (i == j)
			{
				*(L.m_data + i * n + j) = sqrt(*(m_data + i * n + i) - sum);
				if (*(L.m_data + i * n + j) <= 0.0 || isnan(*(L.m_data + i * n + j))) return false;
			}
			else
				*(L.m_data + i * n + j) = (*(m_data + i * n + j) - sum) / (*(L.m_data + j * n + j));
		}

	return true;
}
/// <summary>
/// Получение транспонированной матрицы 
/// </summary>
/// <returns>транспонированную матрицу</returns>
template <typename T>
MATRIX<T> MATRIX<T>::Transpose()
{
	MATRIX<T> tr(m_columns, m_rows);
	for (int i = 0; i < m_columns; i++)
		for (int j = 0; j < m_rows; j++)
			*(tr.m_data + i * tr.m_columns + j) = *(m_data + j * m_columns + i);
	return tr;
}

/// <summary>
/// Деструктор: освобождение памяти, занимаемой под матрицу
/// </summary>
template <typename T>
MATRIX<T>::~MATRIX()
{
	if (!m_data) delete[] m_data;
}


/// <summary>
/// Решение системы линейных алгебраических уравнений (СЛАУ) методом Гаусса с вычислением определителя матрицы системы 
/// </summary>
/// <param name="a">матрица коэффициентов СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">решение СЛАУ</param>
/// <param name="det">определитель матрицы a</param>
/// <returns></returns>
template <typename T>
bool Gauss(const MATRIX<T>& a, const VECTOR<T>& b, VECTOR<T>& x)
{
	int i, k, m;
	T amm, aim;

	// матрица должна быть квадратной и размерность вектора должна совпадать 
	// с размерностью матрицы
	if (a.m_columns != b.m_size || a.m_columns != a.m_rows) return false;
	int size = a.m_rows;

	// сведение исходной системы к системе с верхней треугольной матрицей
	MATRIX<T> alf(a.m_rows, a.m_columns);
	VECTOR<T> bet(b.m_size);
	alf = a;
	bet = b;
	for (m = 0; m <= size - 2; m++)
	{
		amm = *(alf.m_data + m * size + m);
		for (k = m; k <= size - 1; k++)
			*(alf.m_data + m * size + k) /= amm;
		*(bet.m_data + m) /= amm;
		for (i = m + 1; i <= size - 1; i++)
		{
			aim = *(alf.m_data + i * size + m);
			for (k = m; k <= size - 1; k++)
				*(alf.m_data + i * size + k) -= *(alf.m_data + m * size + k) * aim;
			*(bet.m_data + i) -= *(bet.m_data + m) * aim;
		}//end i 
	}//end m 

	// нахождение решения СЛАУ с верхней треугольной матрицей
	*(x.m_data + size - 1) = *(bet.m_data + size - 1) / *(alf.m_data + size * (size - 1) + size - 1);
	for (i = size - 2; i >= 0; i--)
	{
		*(x.m_data + i) = *(bet.m_data + i);
		for (k = i + 1; k < size; k++)
			*(x.m_data + i) -= *(alf.m_data + i * size + k) * *(x.m_data + k);
	}//end i
	return true;

}
/// <summary>
/// Формирование матриц alpha и gamma в компактной схеме исключения
/// Матрица alpha верхняя треугольная, gamma - нижняя треугольная
/// </summary>
/// <param name="alpha">формируемые матрицы</param>
/// <returns>значение определителя матрицы</returns>
template <typename T>
T MATRIX<T>::FormMatrixCompactScheme(MATRIX& alpha)
{
	if (this->m_columns != this->m_rows) return NAN;
	int n = this->m_rows;

	for (int i = 0; i < n; i++)
	{
		*(alpha.m_data + i * n) = *(this->m_data + i * n);
		if (i > 0) *(alpha.m_data + i) = *(this->m_data + i) / *this->m_data;
	}

	T sum = 0;
	int k = 1, i = 1;
	while (i < n)
	{
		if (k >= n)
		{
			k = 1; i++;
			if (i >= n) break;
		}
		if (i >= k)
		{
			sum = 0.0;
			for (int j = 0; j <= k - 1; j++)
			{
				sum += *(alpha.m_data + i * n + j) * *(alpha.m_data + j * n + k);
			}
			*(alpha.m_data + i * n + k) = *(this->m_data + i * n + k) - sum;
		}
		else
		{
			sum = 0.0;
			for (int j = 0; j <= i - 1; j++)
			{
				sum += *(alpha.m_data + i * n + j) * *(alpha.m_data + j * n + k);
			}
			if (abs(*(alpha.m_data + i * n + i)) < 1.0e-18)
			{
				return 0.0;
			}
			*(alpha.m_data + i * n + k) = (*(this->m_data + i * n + k) - sum) / *(alpha.m_data + i * n + i);
		}
		k++;
	}

	// вычисление определителя
	T det = 1;
	for (int i = 0; i < n; i++)
		det *= *(alpha.m_data + i * n + i);
	return det;
}

/// <summary>
/// Решение системы линейных уравнений компактной схемой исключения
/// </summary>
/// <param name="A">матрица системы уравнений</param>
/// <param name="b">вектор правой части системы уравнений</param>
/// <param name="x">решение системы уравнений</param>
template <typename T>
void CompactSchemeSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x)
{
	if (A.m_rows != A.m_columns) return;
	MATRIX<T> alpha(A.m_rows, A.m_columns);
	T det = A.FormMatrixCompactScheme(alpha);
	// решение только для неособенной матрицы
	if (det != 0.0)
		TriangleSolve(alpha, b, x);
}
/// <summary>
/// Решение системы линейных алгебраических уравнений A*x = b с верхней треугольной матрицей
/// </summary>
/// <param name="A">верхняя треугольная матрица системы</param>
/// <param name="b">вектор правой части системы</param>
/// <param name="x">решение системы</param>
template <typename T>
void TriangleSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x)
{
	int n = b.m_size;
	VECTOR<T> beta(n);
	T sum = 0;
	*beta.m_data = *b.m_data / *A.m_data;
	for (int i = 1; i < n; i++)
	{
		sum = 0.0;
		for (int j = 0; j <= i - 1; j++)
			sum += *(A.m_data + i * n + j) * *(beta.m_data + j);
		*(beta.m_data + i) = (*(b.m_data + i) - sum) / *(A.m_data + i * n + i);
	}

	// решение системы уравнений с труегольной матрицей		   
	*(x.m_data + n - 1) = *(beta.m_data + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		sum = 0.0;
		for (int j = n - 1; j > i; j--)
			sum += *(A.m_data + i * n + j) * *(x.m_data + j);
		*(x.m_data + i) = *(beta.m_data + i) - sum;
	}

}
/// <summary>
/// Решение СЛАУ с применением QR декомпозиции матрицы системы A
/// </summary>
/// <param name="A">матрица системы линейных уравнений</param>
/// <param name="b">вектор правой части системы линейных уравнений</param>
/// <param name="x">вектор решения системы линейныз уравнений</param>
template <typename T>
void QRDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x)
{
	if (A.m_rows != A.m_columns || A.m_rows != b.m_size) return;
	MATRIX<T> Q(A.m_rows, A.m_columns), R(A.m_rows, A.m_columns);
	int n = A.m_rows;
	A.QRDecomposition(Q, R);
	VECTOR <T>beta(b.m_size);
	beta = Q.Transpose() * b;

	// решение системы уравнений с верхней труегольной матрицей		   
	*(x.m_data + n - 1) = *(beta.m_data + n - 1) / *(R.m_data + (n - 1) * n + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		T sum = 0;
		for (int j = n - 1; j > i; j--)
			sum += *(R.m_data + i * n + j) * *(x.m_data + j);
		*(x.m_data + i) = (*(beta.m_data + i) - sum) / *(R.m_data + i * n + i);
	}


}
/// <summary>
/// Решение системы линейных алгебраических уравнений с применением метода LLT декомпозиции A = LLT
/// </summary>
/// <param name="A">матрица</param>
/// <param name="b">вектор правой части системы</param>
/// <param name="x">решение системы</param>
template <typename T>
void LLTDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x)
{
	if (A.m_rows != A.m_columns || A.m_rows != b.m_size) return;
	
	int n = A.m_rows;
	if (!A.IsSymmetric())
	{
		throw "Метод предназначен для решения системы с симметричной матрицей";
		return;
	}

	MATRIX<T> L(A.m_rows, A.m_columns);
	A.CholeskyDecomposition(L);
	// решение системы с верхней треугольной матрицей Ly = b
	VECTOR<T> y(n);
	for (int i = 0; i < n; i++)
	{
		*(y.m_data + i) = *(b.m_data + i);
		T sum = 0;
		for (int k = 0; k < i; k++)
		{
			sum += *(L.m_data + i * n + k) * *(y.m_data + k);
		}
		*(y.m_data + i) -= sum;
		*(y.m_data + i) /= *(L.m_data + i * n + i);

	}

	L = L.Transpose();
	// решение системы уравнений с нижней треугольной матрицей L^tx = y   
	*(x.m_data + n - 1) = *(y.m_data + n - 1) / *(L.m_data + (n - 1) * n + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		T sum = 0;
		for (int j = n - 1; j > i; j--)
			sum += *(L.m_data + i * n + j) * *(x.m_data + j);
		*(x.m_data + i) = (*(y.m_data + i) - sum) / *(L.m_data + i * n + i);
	}

}

/// <summary>
/// Декомпозиция матрицы A = LU
/// L - нижняя треугольная матрица, на главной диагонали которой расположены единицы
/// U - верхняя треугольная матрица
/// </summary>
/// <param name="alfa">Матрицы L и U, ниже гланой диагонали которой расположены внедиагональные элементы L, 
/// на главной диагонали и выше расположены элементы матрицы U</param>
/// <returns>Определитель матрицы</returns>
template <typename T>
T MATRIX<T>::LUDecomposition(MATRIX& alfa)
{
	T det = 0;
	if (m_rows != m_columns) return det;

	// декомпозиция матрицы
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_columns; j++)
		{
			T val = 0;
			if (i <= j)
			{
				for (int k = 0; k <= i - 1; k++)
					val += *(alfa.m_data + i * alfa.m_columns + k) * *(alfa.m_data + k * alfa.m_columns + j);
				*(alfa.m_data + i * alfa.m_columns + j) = *(m_data + i * m_columns + j) - val;
			}
			else
			{
				for (int k = 0; k <= j - 1; k++)
					val += *(alfa.m_data + i * alfa.m_columns + k) * *(alfa.m_data + k * alfa.m_columns + j);
				if (abs(*(alfa.m_data + j * alfa.m_columns + j)) < 1.0e-36) return 0;
				*(alfa.m_data + i * alfa.m_columns + j) = (*(m_data + i * m_columns + j) - val) / *(alfa.m_data + j * alfa.m_columns + j);

			}
		}
	det = 1;
	for(int i = 0; i < alfa.m_rows; i++)
		det *= *(alfa.m_data + i * alfa.m_columns + i);
	return det;
}

/// <summary>
/// Решение системы линейных алгебраических уравнений с применением метода LU декомпозиции A = LU
/// L - нижняя тругольная матрица
/// U - верхняя треугольная матрица 
/// </summary>
/// <param name="A">матрица</param>
/// <param name="b">вектор правой части системы</param>
/// <param name="x">решение системы</param>

template <typename T>
void LUDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x)
{
	if (A.m_rows != A.m_columns || A.m_rows != b.m_size) return;

	MATRIX<T> alfa(A.m_rows, A.m_columns);
	int n = A.m_rows;

	T det = A.LUDecomposition(alfa);
	if (abs(det) > 1.0e-36)
	{
		// решение СЛАУ
		VECTOR<T> y(n);
		for (int k = 0; k < n; k++)
		{
			T val = 0;
			for (int j = 0; j <= k - 1; j++)
				val += *(alfa.m_data + k * alfa.m_columns + j) * *(y.m_data + j);
			*(y.m_data + k) = *(b.m_data + k) - val;
		}

		for (int k = n - 1; k >= 0; k--)
		{
			T sum = 0;
			for (int j = n - 1; j > k; j--)
				sum += *(alfa.m_data + k * alfa.m_columns + j) * *(x.m_data + j);
			*(x.m_data + k) = (*(y.m_data + k) - sum) / *(alfa.m_data + k * alfa.m_columns + k);
		}
	}
}
/// <summary>
/// Вычисление минора квадаратной матрицы
/// </summary>
/// <param name="i">строка</param>
/// <param name="j">столбец</param>
/// <returns>значение минора, в случае ошибки NAN</returns>
template <typename T>
T MATRIX<T>::Minor(int i, int j)
{
	int n = m_rows;
	if (m_columns != m_rows)
	{
		throw "Миноры вычисляются только для квадратной матрицы";
		return NAN;
	}
	if (i < 0 || j < 0 || i > n - 1 || j > n - 1 || m_columns != m_rows)
	{
		throw "Выход индексов минора за допустимые пределы";
		return NAN;
	}
	MATRIX minor(n - 1, n - 1);

	// заполнение матрицы минора данными
	for (int k = 0; k < n; k++)
	{
		for (int m = 0; m < n; m++)
		{
			if (k < i && m < j)
				*(minor.m_data + k * (n - 1) + m) = *(m_data + k * n + m);
			else if (k < i && m > j)
				*(minor.m_data + k * (n - 1) + m - 1) = *(m_data + k * n + m);
			else if (k > i && m < j)
				*(minor.m_data + (k - 1) * (n - 1) + m) = *(m_data + k * n + m);
			else if (k > i && m > j)
				*(minor.m_data + (k - 1) * (n - 1) + m - 1) = *(m_data + k * n + m);
		}
	}

	return minor.Determinant();
}
/// <summary>
/// Возвращает обратную матрицу при помощи алгебраических дополнений
/// </summary>
/// <returns>обратную матрицу</returns>
template <typename T>
MATRIX<T> MATRIX<T>::Reverse()
{
	MATRIX A(m_rows, m_columns);
	if (m_rows != m_columns)
	{
		throw "Для вычисления обратной матрицы исходная матрица должна быть квадратной";
		return A;
	}
	T det = Determinant();

	if (abs(det) >= 1.0e-36)
	{
		auto Minors = [&](int row_begin, int row_end, int column_begin, int column_end)
			{
				for (int i = row_begin; i < row_end; i++)
					for (int j = column_begin; j < column_end; j++)
					{
						*(A.m_data + A.m_columns * i + j) = Minor(j, i) * ((i + j) % 2 > 0 ? -1.0 : 1.0) / det;
					}
			};
		int size = A.m_rows;
		if (size >= MIN_SIZE_FOR_THREAD)
		{
			int size2 = size / 2;
			thread t1(Minors, 0, size2, 0, size2);
			thread t2(Minors, size2, size, 0, size2);
			thread t3(Minors, 0, size2, size2, size);
			thread t4(Minors, size2, size, size2, size);

			t1.join();
			t2.join();
			t3.join();
			t4.join();
		}
		else
			Minors(0, size, 0, size);
	}


	return A;
}
/// <summary>
/// Обращение марицы при пощи решения систем уравнений A*X = E, колонки матрицы X - столбцы обратной матрицы
/// </summary>
/// <returns>обратную матрицу</returns>
template <typename T>
MATRIX<T> MATRIX<T>::Invert()
{
	MATRIX A(m_rows, m_columns);
	if (m_rows != m_columns)
	{
		throw "Для вычисления обратной матрицы исходная матрица должна быть квадратной";
		return A;
	}

	MATRIX alpha(m_rows, m_columns);
	T det = FormMatrixCompactScheme(alpha);

	if (abs(det) < 1.0e-36)
	{
		throw "Это вырожденная матрица";
		return A;
	}

	for (int i = 0; i < m_rows; i++)
	{
		VECTOR<T> e(m_rows), x(m_rows);
		e[i] = 1;
		TriangleSolve(alpha, e, x);
		A.CopyColumn(x, i);
	}
	return A;
}

/// <summary>
/// Копировать вектор в  j-ю колонку матрицы
/// </summary>
/// <param name="v">вектор</param>
/// <param name="j">номер колонки</param>
template <typename T>
void MATRIX<T>::CopyColumn(VECTOR<T>& v, int j)
{
	if (m_columns < j)
	{
		throw "Номер колонки превышает число колонок матрицы";
		return;
	}
	if (m_rows != v.size())
	{
		throw "Число строк (число элементов в колонке) матрицы не совпадает с размерностью вектора";
		return;
	}
	for (int i = 0; i < m_rows; i++)
		*(m_data + i * m_columns + j) = v[i];
}

/// <summary>
/// Копировать j-ю колонку матрицы в вектор
/// </summary>
/// <param name="j">номер колонки</param>
/// <returns>вектор содержащий колонку матрицы</returns>
template <typename T>
VECTOR<T> MATRIX<T>::CopyColumn2Vector(int j)
{
	VECTOR<T> v(m_rows);
	if (m_columns < j)
	{
		throw "Номер колонки превышает число колонок матрицы";
		return v;
	}
	for (int i = 0; i < m_rows; i++)
		v[i] = *(m_data + i * m_columns + j);
	return v;
}

// Характеристический полином x ^ n + koeff[1] * x ^ (n - 1) + ... + koeff[n-1] 
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
template <typename T>
complex<T> Polyfun(const complex<T>& z, VECTOR<T>& koeff)
{
	complex<T> poly(0, 0);
	int n = koeff.size();
	poly = -1.0 * koeff[n - 1] - koeff[n - 2] * z - koeff[n - 3] * z * z;
	T p = 0;
	for (int i = 3; i < n; i++)
	{
		p = (T)i;
		poly -= koeff[n - i - 1] * pow(z, p);

	}
	poly += pow(z, (T)n);
	return poly;
}

// Производная характеристического полинома x ^ n + koeff[1] * x ^ (n - 1) + ... + koeff[n-1]
template <typename T>
complex<T> PolyfunDerivative(const complex<T>& z, VECTOR<T>& koeff)
{
	complex<T> poly(0, 0);
	int n = koeff.size();
	poly = -1.0 * koeff[n - 2] - 2.0 * koeff[n - 3] * z - 3.0 * koeff[n - 4] * z * z;
	T p = 0 ;
	for (int i = 4; i < n; i++)
	{
		p = (T)(i - 1);
		poly -= ((T)i) * koeff[n - i - 1] * pow(z, p);

	}
	poly += (T)n * pow(z, (T)(n - 1));
	return poly;
}
/// <summary>
/// Генарация вещественных случайных чисел
/// </summary>
/// <typeparam name="T">тип</typeparam>
/// <param name="min">минимальное значение и нтервала</param>
/// <param name="max">максимальное значение интервала</param>
/// <returns>случайное число</returns>
template <typename T>
T rand_range(T min, T max) {
	return min + (max - min) * std::rand() / (RAND_MAX + 1);
}

template <typename T>
inline T randomDouble()
{
	return (T)(rand()) / (T)(rand());
}

// Заполнение массива начальных значений корней полинома init_values
// при помощи генератора случайных чисел
// n - степень полинома.
// Если полином нечётной степени, то первый элемент массива - вещественное число
template <typename T>
void FillInitialValues(complex<T>* init_values, int n)
{
	srand(1);
	int i = 0;
	while (i < n)
	{
		if (i == 0) init_values[i] = complex<T>(randomDouble<T>(), n % 2 != 0 ? 0 : randomDouble<T>());
		else
		{
			init_values[i] = complex<T>(randomDouble<T>(), randomDouble<T>());
			if (i == n - 1) break;
			i++;
			init_values[i] = n % 2 != 0 ? conj(init_values[i - 1]) : complex<T>(randomDouble<T>(), 0);
		}
		i++;
	}
}

// Поиск всех корней полинома степени n с вещественными коэффициентами 
// методом Аберта-Эрлиха (https://en.wikipedia.org/wiki/Aberth_method)
// Полином koeff[n]*x^n+koeff[n-1]*x^(n-1)+...+koeff[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n
template <typename T>
void Polyroots(VECTOR<T>& koeff, complex<T>* roots)
{
	int n = koeff.size();
	complex<T>* w = new complex<T>[n]; // массив начальных значений корней полинома, далее используется для чисел смещения
	FillInitialValues(roots, n);

	T err = 1.0;
	int iter = 0;
	while (err >= EPS && iter < MAX_ITER_NUMBER)
	{
		for (int i = 0; i < n; i++)
		{
			complex<T> pr(0, 0);
			for (int j = 0; j < n; j++)
			{
				if (i != j) pr += 1.0 / (roots[i] - roots[j]);
			}
			complex<T> pp = Polyfun(roots[i], koeff) / PolyfunDerivative(roots[i], koeff);
			w[i] = pp / (1.0 - pp * pr);
		}

		// проверка максимального модуля чисел смещений
		err = -1.0;
		T err_max = 0;

		for (int i = 0; i < n; i++)
		{
			err_max = abs(w[i]);
			if (err_max > err) err = err_max;
			roots[i] -= w[i];
		}
		iter++;
#ifdef _DEBUG
		cout << "Кол-во итераций: " << iter << ". Значение погрешности вычислений: " << err << endl;
#endif

	}
	delete[] w;
}

/// <summary>
/// Нахождение собственных значений и собственных векторов матрицы методом Крылова А.Н.
/// </summary>
/// <param name="lambda">массив собственных значений</param>
/// <param name="vect">массив собственных векторов</param>
template <typename T>
void MATRIX<T>::EigenvaluesAndVectorsKrylov(complex<T>* lambda, complex<T>** vect)
{
	if (m_rows != m_columns)
	{
		throw "Матрица не квадратная";
		return;
	}

	// поиск собственных значений
	int n = m_columns;
	VECTOR<T> mcolumn(n);
	MATRIX a(n, n);
	mcolumn[0] = 1.0;

	for (int i = n - 1; i >= 0; i--)
	{
		a.CopyColumn(mcolumn, i);
		mcolumn = *this * mcolumn;
	}

	VECTOR<T> p(n); // значения коэффициентов полинома
	QRDecompositionSolve(a, mcolumn, p);
#ifdef _DEBUG
	cout << "Матрица СЛАУ для определения коэффициентов полинома" << endl;
	cout << a << endl;

	cout << "Вектор правой части для определения коэффициентов характеристического полинома" << endl;
	cout << mcolumn << endl;

	cout << "Коэффициенты характеристического полинома" << endl;
	cout << p << endl;
#endif // _DEBUG

	Polyroots(p, lambda);

	// поиск собственных векторов 
	if (vect == nullptr) return; // только собственные значения
	for (int i = 0; i < n; i++)
	{
		mcolumn = a.CopyColumn2Vector(0);
		for (int j = 0; j < n; j++)
			vect[i][j] = mcolumn[j];


		complex<T> q(1, 0);
		for (int k = 1; k < n; k++)
		{
			mcolumn = a.CopyColumn2Vector(k);
			q = lambda[i] * q - p[k - 1];
			for (int j = 0; j < n; j++)
			{
				vect[i][j] += mcolumn[j] * q;
			}
		}
	}
}

/// <summary>
/// Нахождение собственных значений и собственных векторов матрицы методом Ле Веррье - Фаддеева Д.К.
/// </summary>
/// <param name="lambda">массив собственных значений</param>
/// <param name="vect">массив собственных векторов</param>
template <typename T>
void MATRIX<T>::EigenvaluesAndVectorsLeVerrierFaddeev(complex<T>* lambda, complex<T>** vect)
{
	if (m_rows != m_columns)
	{
		throw "Матрица не квадратная";
		return;
	}
	// поиск собственных значений
	int n = m_columns;
	VECTOR<T> p(n);
	vector<VECTOR<T>> b(n, VECTOR<T>(n));
	MATRIX A(n, n), B(n, n);

	for (int i = 0; i < n; i++)
	{
		A = i == 0 ? *this : *this * B;
		p[i] = A.Sp() / (i + 1);
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				if (j != k) B(j, k) = A(j, k);
				else B(j, k) = A(j, k) - p[i];

		b[i] = B.CopyColumn2Vector(0);
	}

#ifdef _DEBUG
	cout << "Коэффициенты характеристического полинома" << endl;
	cout << p << endl;
#endif // _DEBUG

	Polyroots(p, lambda);

	// поиск собственных векторов 
	if (vect == nullptr) return; // только собственные значения
	VECTOR<complex<T>> y(n);

	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < n; k++)
			y[k] = k == 0 ? complex<T>(1, 0) : complex<T>(0, 0);
		for (int j = 1; j < n; j++)
		{
			for (int k = 0; k < n; k++)
				y[k] = lambda[i] * y[k] + b[j - 1][k];

		}

		for (int k = 0; k < n; k++)
			vect[i][k] = y[k];
	}
	b.clear();
}

/// <summary>
/// Обращение матрицы методом Фаддеева Д.К.
/// </summary>
/// <returns>обратную матрицу</returns>
template <typename T>
MATRIX<T> MATRIX<T>::InvertFaddev()
{
	MATRIX Ainv(m_rows, m_columns);
	if (m_rows != m_columns)
	{
		throw "Матрица не квадратная";
		return Ainv;
	}
	// поиск собственных значений
	int n = m_columns;
	VECTOR<T> p(n);

	MATRIX A(n, n), B(n, n);

	for (int i = 0; i < n; i++)
	{
		A = i == 0 ? *this : *this * B;
		p[i] = A.Sp() / (i + 1);
		if (i == n - 1) Ainv = B;
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				if (j != k) B(j, k) = A(j, k);
				else B(j, k) = A(j, k) - p[i];
	}
	Ainv /= p[n - 1];
	return Ainv;
}

/// <summary>
/// Нахождение собственных значений и собственных векторов
/// </summary>
/// <param name="lambda">собственные значения</param>
/// <param name="vect">собственные векторы</param>
template <typename T>
void MATRIX<T>::Rotate(T* lambda, T** vect)
{
	if (!IsSymmetric())
	{
		throw "Метод вращений применяется только для симметричных матриц";
		return;
	}
	MATRIX A(m_rows, m_columns), U(m_rows, m_columns);
	A = *this;
	int n = A.m_columns;
	int iter = 1;
	while (true)
	{
		// поиск позиции максимального по модулю внедиагонального элемента
		int i0 = 0, j0 = 0;
		T val = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				T el = abs(*(A.m_data + i * A.m_columns + j));
				if (i != j && el > val)
				{
					val = el;
					i0 = i;
					j0 = j;
				}
			}
		if (val < EPS) break;
		// угол матрицы вращения
		T fi = 0.5 * atan(2.0 * *(A.m_data + i0 * A.m_columns + j0) / (*(A.m_data + i0 * A.m_columns + i0) - *(A.m_data + j0 * A.m_columns + j0)));
		T cs = cos(fi);
		T ss = sin(fi);

		// U0*U1*...Un
		MATRIX U0(n, n), B(n, n);

		if (iter > 1)
		{
			U0 = U;
			for (int i = 0; i < n; i++)
			{
				*(U.m_data + i * U.m_columns + i0) = *(U0.m_data + i * U0.m_columns + i0) * cs +
					*(U0.m_data + i * U0.m_columns + j0) * ss;
				*(U.m_data + i * U.m_columns + j0) = (-1.0) * *(U0.m_data + i * U0.m_columns + i0) * ss +
					*(U0.m_data + i * U0.m_columns + j0) * cs;
			}
		}
		else
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
				{
					if (i == j && i != i0 && j != j0) *(U0.m_data + i * U0.m_columns + j) = 1.0;
					else if (i == i0 && j == i0) *(U0.m_data + i * U0.m_columns + j) = cs;
					else if (i == i0 && j == j0) *(U0.m_data + i * U0.m_columns + j) = -ss;
					else if (i == j0 && j == j0) *(U0.m_data + i * U0.m_columns + j) = cs;
					else if (i == j0 && j == i0) *(U0.m_data + i * U0.m_columns + j) = ss;
					else *(U0.m_data + i * U0.m_columns + j) = 0.0;
				}

			U = U0;
		}

		// Ut*A*U
		B = A;
		for (int i = 0; i < n; i++)
		{
			*(B.m_data + i * B.m_columns + i0) = *(A.m_data + i * A.m_columns + i0) * cs +
				*(A.m_data + i * A.m_columns + j0) * ss;
			*(B.m_data + i * B.m_columns + j0) = (-1.0) * *(A.m_data + i * A.m_columns + i0) * ss +
				*(A.m_data + i * A.m_columns + j0) * cs;
		}
		A = B;
		for (int i = 0; i < n; i++)
		{
			*(A.m_data + i0 * A.m_columns + i) = *(B.m_data + i0 * B.m_columns + i) * cs +
				*(B.m_data + j0 * B.m_columns + i) * ss;
			*(A.m_data + j0 * A.m_columns + i) = (-1.0) * *(B.m_data + i0 * B.m_columns + i) * ss +
				*(B.m_data + j0 * B.m_columns + i) * cs;
		}

		iter++;
	}

	if (vect == nullptr) return;

	for (int i = 0; i < n; i++)
	{
		lambda[i] = A(i, i);
		for (int j = 0; j < n; j++)
			vect[i][j] = *(U.m_data + j * U.m_columns + i);
	}

}
