#pragma once
#include "VECTOR.h"

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
	memset(m_data, 0, sizeof m_data);

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
double VECTOR<T>::norm()
{
	double val = 0.0;
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
	if (this->m_size != src.m_size)
	{
		if (!this->m_data) delete[] this->m_data;
		this->m_data = new double[src.m_size];
		this->m_size = src.m_size;
	}
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
		prod += *(v1.m_data + i) *(v2.m_data + i);
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
		swap(*(m_data+start), *(m_data+end));
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
bool VECTOR<T>::writeToFile(const char* fileName, VECTOR<T> &vect)
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
