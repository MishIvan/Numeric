#include "VECTOR.h"
/// <summary>
///  онструктор вектора: выделение пам€ти под указатель и инициализаци€ компонент вектора значени€ми val
/// </summary>
/// <param name="n">размерность вектора</param>
VECTOR::VECTOR(int n, double val)
{
	m_data = 0; m_size = 0;
	if (n < 1) 
		 return;
	m_size = n;
	m_data = new double[m_size];
	for (int i = 0; i < m_size; i++)
		*(m_data + i) = val;

}
/// <summary>
///  онструктор копировани€ вектора src
/// </summary>
/// <param name="src">исходный вектор дл€ копировани€</param>
VECTOR::VECTOR(const VECTOR& src)
{
	if (!this->m_data) delete[] this->m_data;
	this->m_data = new double[src.m_size];
	this->m_size = src.m_size;
	for (int i = 0; i < this->m_size; i++)
		*(this->m_data + i) = *(src.m_data + i);
}

/// <summary>
/// ¬озвращает евклидову норму вектора
/// </summary>
/// <returns>yорма вектора</returns>
double VECTOR::norm()
{
	double val = 0.0;
	for (int i = 0; i < m_size; i++)
		val += *(m_data + i) * *(m_data + i);
	return sqrt(val);
}

/// <summary>
/// 
/// </summary>
/// <param name="src">исходный вестор дл€ присвоени€</param>
/// <returns></returns>
VECTOR& VECTOR::operator=(const VECTOR& src)
{
	if (!this->m_data) delete[] this->m_data;
	this->m_data = new double[src.m_size];
	this->m_size = src.m_size;
	for (int i = 0; i < this->m_size; i++)
		*(this->m_data + i) = *(src.m_data + i);
	return *this;
}
/// <summary>
/// ¬озвращает скал€рное произведение векторов v1 и v2
/// ¬ектора должны быть одной размерности
/// </summary>
/// <param name="v1">первый вектор</param>
/// <param name="v2">второй вектор</param>
/// <returns></returns>
double operator*(const VECTOR& v1, const VECTOR& v2)
{
	double prod = 0.0;
	if (v1.m_size != v2.m_size) return prod;
	for (int i = 0; i < v1.m_size; i++)
		prod += *(v1.m_data + i) * *(v2.m_data + i);
	return prod;
}
/// <summary>
/// ѕерегрузка операторов вычитани€ и сложени€ векторов v1 и v2
/// размерности векторов должны совпадать
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns>вектор - разность (сумму) векторов  </returns>
VECTOR operator-(const VECTOR& v1, const VECTOR& v2)
{
	if (v1.m_size != v2.m_size) return VECTOR(3);
	VECTOR res(v1.m_size);
	for (int i = 0; i < res.m_size; i++)
		*(res.m_data + i) = *(v1.m_data + i) - *(v2.m_data + i);
	return res;
}

VECTOR operator+(const VECTOR& v1, const VECTOR& v2)
{
	if (v1.m_size != v2.m_size) return VECTOR(3);
	VECTOR res(v1.m_size);
	for (int i = 0; i < res.m_size; i++)
		*(res.m_data + i) = *(v1.m_data + i) - *(v2.m_data + i);
	return res;
}

/// <summary>
/// ѕерегрузка оператора вывода на консоль или в файловый поток (файл нужно будет открыть дл€ чтени€)
/// </summary>
/// <param name="s"></param>
/// <param name="v">вектор</param>
/// <returns></returns>
ostream& operator<<(ostream& s, VECTOR& v)
{
	int n = v.m_size;
	for (int i = 0; i < n; i++)
		s << *(v.m_data + i) << ' ';
	return s;
}
/// <summary>
/// ѕерегрузка оператора ввода с консоли или в файловый поток (файл должен быть открыт дл€ записи) 
/// </summary>
/// <param name="s"></param>
/// <param name="v">вектор</param>
/// <returns></returns>
istream& operator>>(istream& s, VECTOR& v)
{
	int n = v.m_size;
	for (int i = 0; i < n; i++)
		s >> *(v.m_data + i); 
	return s;
}

/// <summary>
/// „тение вектора из текстового файла
/// </summary>
/// <param name="fileName">полное им€ файла</param>
/// <param name="vect">вектор, создаваемый по данным файла</param>
/// <returns>true - в случае успешного считывани€ данных, false - в случае ошибки </returns>
bool VECTOR::readFromFile(const char* fileName, VECTOR& vect)
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
/// «апись ветора в текстовый файл
/// </summary>
/// <param name="fileName">полное им€ текстового файла</param>
/// <returns>true - в случае успешного считывани€ данных, false - в случае ошибки</returns>
bool VECTOR::writeToFile(const char* fileName, VECTOR &vect)
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
/// ƒеструктор. ќсвобождение пам€ти
/// </summary>
VECTOR::~VECTOR()
{
	if (m_data) delete[] m_data;
}
