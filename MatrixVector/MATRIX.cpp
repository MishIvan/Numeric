#include "MATRIX.h"

/// <summary>
/// Генарация вещественных случайных чисел
/// </summary>
/// <typeparam name="T">тип</typeparam>
/// <param name="min">минимальное значение и нтервала</param>
/// <param name="max">максимальное значение интервала</param>
/// <returns>случайное число</returns>
double rand_range(double min, double max) {
	return min + (max - min) * std::rand() / (RAND_MAX + 1);
}


/// <summary>
/// Конструктор матрицы M x N, всем элементам которой присваивается значение val
/// </summary>
/// <param name="M">число строк матрицы</param>
/// <param name="N">число столбцов матрицы</param>
/// <param name="val"></param>
MATRIX::MATRIX(int M, int N, double val)
{
	m_columns = m_rows = 0;
	m_data = 0;
	if (M < 0 || N < 0)
	{
		throw invalid_argument("Задана отрицательная размерность матрицы");
		return;
	}
	m_rows = M; m_columns = N;
	this->m_data = new double[this->m_rows * this->m_columns];
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			*(m_data + i * m_columns + j) = val;
}
/// <summary>
/// Конструктор копирования матрицы из матрицы src 
/// </summary>
/// <param name="src">исходная матрица для копирования</param>
MATRIX::MATRIX(const MATRIX& src)
{
	if (!this->m_data) delete[] this->m_data;
	this->m_columns = src.m_columns;
	this->m_rows = src.m_rows;
	this->m_data = new double[this->m_rows * this->m_columns];
	for (int i = 0; i < this->m_rows; i++)
		for (int j = 0; j < this->m_columns; j++)
			*(this->m_data + i * this->m_columns + j) = *(src.m_data + i*src.m_columns +j);
}
/// <summary>
/// Перегрузка оператора присваивания
/// </summary>
/// <param name="src">матрица для присваивания</param>
/// <returns></returns>
MATRIX& MATRIX::operator=(const MATRIX& src)
{
	if (this->m_rows != src.m_rows || this->m_columns != src.m_columns)
	{
		if (!this->m_data) delete[] this->m_data;
		this->m_columns = src.m_columns;
		this->m_rows = src.m_rows;
		this->m_data = new double[this->m_rows * this->m_columns];
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
MATRIX operator*(const MATRIX& matr1, const MATRIX& matr2)
{
	MATRIX pmatr(matr1.m_rows, matr2.m_columns);
	if (matr1.m_columns != matr2.m_rows)
	{
		throw("Несоответствие размерности матриц");
		return pmatr;
	}
	int i, j;
	for(i = 0; i < matr1.m_rows; i++)
		for (j = 0; j < matr2.m_columns; j++)
		{
			double prod = 0.0;
			for (int k = 0; k < matr1.m_columns; k++)
			{
				if (*(matr1.m_data + i * matr1.m_columns + k) == 0.0 || *(matr2.m_data + k * matr2.m_columns + j) == 0.0)
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
ostream& operator<<(ostream& s, MATRIX& matr)
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
istream& operator>>(istream& s, MATRIX& matr)
{
	int m = matr.m_rows;
	int n = matr.m_columns;
	double* buff = new double[n];
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
double MATRIX::Determinant()
{
	if(m_rows != m_columns) return NAN;
	MATRIX alpha(m_rows, m_columns);
	return FormMatrixCompactScheme(alpha);
}

/// <summary>
/// Вычисляет первую норму матрицы
/// </summary>
/// <returns>значение первой нормы матрицы</returns>
double MATRIX::normI()
{
	double norm = 0;
	for (int i = 0; i < m_rows; i++)
	{
		double sum = 0;
		for (int j = 0; j < m_columns; j++)
			sum += abs(*(m_data + i * m_columns + j));
		if (sum > norm) norm = sum;
	}
	return norm;
}

/// <summary>
/// Считывание данных их файла и загрузка их в матрицу matr
/// </summary>
/// <param name="fileName">полное имя файла</param>
/// <param name="matr">матрица</param>
/// <returns>true - в случае успешного считывания данных, false - в случае ошибки</returns>
bool MATRIX::readFromFile(const char* fileName, MATRIX& matr)
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
bool MATRIX::writeToFile(const char* fileName, MATRIX& matr)
{
	ofstream fs;
	fs.open(fileName);
	if(fs.is_open()){
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
 /// <summary>
/// Вычислить след матрицы, сумму её колонок по диагонали
/// </summary>
/// <returns>значение следа матрицы</returns>
double MATRIX::Sp()
{
	double val = 0.0;
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
 MATRIX operator+(const MATRIX& matr1, const MATRIX& matr2)
{
	MATRIX prod(matr1.m_rows, matr1.m_columns);
	if (matr1.m_columns != matr2.m_columns && matr1.m_rows != matr2.m_rows)
	{
		throw("Размерности матриц не совпадают");
		return prod;
	}
	for (int i = 0; i < matr1.m_rows; i++)
		for (int j = 0; j < matr1.m_columns; j++)
			*(prod.m_data + i*prod.m_columns + j) = *(matr1.m_data + i * matr1.m_columns + j) + *(matr2.m_data + i * matr2.m_columns + j);
	return prod;
}

 /// <summary>
/// Вычитание марицы matr2 из матрицы matr1
/// </summary>
/// <param name="matr1">матрица</param>
/// <param name="matr2">матрица</param>
/// <returns>матрица - результат вычитания</returns>
 MATRIX operator-(const MATRIX& matr1, const MATRIX& matr2)
 {
	 MATRIX prod(matr1.m_rows, matr1.m_columns);
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
VECTOR operator*(const MATRIX& matr, const VECTOR& v)
{
	VECTOR prod(matr.m_rows);
	if (matr.m_columns != v.m_size) return prod;
	for (int i = 0; i < matr.m_rows; i++)
	{
		double val = 0.0;
		for (int j = 0; j < matr.m_columns; j++)
		{
			if ( *(matr.m_data + i * matr.m_columns + j) == 0.0 || *(v.m_data + j) == 0.0) continue;
			val += *(matr.m_data + i * matr.m_columns + j) * *(v.m_data + j);
		}
		*(prod.m_data + i) = val;
	}
	return prod;
}

/// <summary>
/// Умножение матрицы на матрицу matr справа. Результат записывается в матрицу
/// </summary>
/// <param name="matr">матрица, на которую умножают справа</param>
/// <returns>эта матрица как результат умножения</returns>
MATRIX& MATRIX::operator*=(const MATRIX& matr)
{
	if (m_columns != matr.m_rows)
	{
		throw("Несоответсвие размерности матриц");
		return *this;
	}
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_columns; j++)
		{
			double val = 0.0;
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
MATRIX& MATRIX::operator/=(double val)
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
MATRIX operator*(double alf, const MATRIX& matr)
{
	MATRIX prod(matr.m_rows, matr.m_columns );
	for (int i = 0; i < matr.m_rows; i++)	
		for (int j = 0; j < matr.m_columns; j++)
			*(prod.m_data +i * prod.m_columns +j) += *(matr.m_data + i * matr.m_columns + j) * alf;	
	return prod;
}

/// <summary>
/// QR разложение квадратной матрицы
/// </summary>
/// <param name="Q">матрица Q</param>
/// <param name="R">матрица R</param>
bool MATRIX::QRDecomposition(MATRIX& Q, MATRIX& R)
{
	if (m_columns != m_rows) return false;
	int n = m_rows;
	double sum = 0.0;

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
				sum += *(Q.m_data + k* n + i) * *(m_data + k * n + j);
			*(R.m_data + i * n + j) = sum;

			// r(i,j)*q(i)
			for(int k=0; k < n; k++)
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
bool MATRIX::CholeskyDecomposition(MATRIX& L)
{
	if (m_rows != m_columns) return false;
	int n = m_rows;
	for (int i = 0; i < n; i++) 
		for (int j = 0; j <= i; j++) {

			double sum = 0;
			for (int k = 0; k < j; k++) 
				sum += *(L.m_data + i*n + k) * *(L.m_data + j*n + k);

			if (i == j)
			{
				*(L.m_data + i * n + j) = sqrt(*(m_data + i * n + i) - sum);
				if ( *(L.m_data + i * n + j) <= 0.0 || isnan(*(L.m_data + i * n + j)) ) return false;
			}
			else 
				*(L.m_data + i * n + j) = (*(m_data + i * n + j) - sum)/(*(L.m_data + j * n + j));
		}
	
	return true;
}
/// <summary>
/// Получение транспонированной матрицы 
/// </summary>
/// <returns>транспонированную матрицу</returns>
MATRIX MATRIX::Transpose()
{
	MATRIX tr(m_columns, m_rows);
	for (int i = 0; i < m_columns; i++)
		for (int j = 0; j < m_rows; j++)
			*(tr.m_data + i * tr.m_columns + j) = *(m_data + j * m_columns + i);
	return tr;
}

/// <summary>
/// Копировать вектор в  j-ю колонку матрицы
/// </summary>
/// <param name="v">вектор</param>
/// <param name="j">номер колонки</param>
void MATRIX::CopyColumn(VECTOR& v, int j)
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
/// Деструктор: освобождение памяти, занимаемой под матрицу
/// </summary>
MATRIX::~MATRIX()
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
bool Gauss(MATRIX &a, VECTOR &b, VECTOR &x)
{
	int i, k, m;
	long double amm, aim;

	// матрица должна быть квадратной и размерность вектора должна совпадать 
	// с размерностью матрицы
	if (a.columns() != b.size() || a.columns() != a.rows()) return false;
	int size = a.rows();

	// сведение исходной системы к системе с верхней треугольной матрицей
	MATRIX alf(size, size);
	VECTOR bet(size);
	alf = a;
	bet = b;
	for (m = 0; m <= size - 2; m++)
	{
		amm =  alf(m, m);
		for (k = m; k <= size - 1; k++)
			alf(m,k) /= amm;
		bet[m]/= amm;
		for (i = m + 1; i <= size - 1; i++)
		{
			aim = alf(i,m);
			for (k = m; k <= size - 1; k++)
				alf(i,k) -= alf(m,k) * aim;
			bet[i] -= bet[m] * aim;
		}//end i 
	}//end m 

	// нахождение решения СЛАУ с верхней треугольной матрицей
	x[size - 1] = bet[size - 1] / alf(size - 1, size - 1);
	for (i = size - 2; i >= 0; i--)
	{
		x[i] = bet[i];
		for (k = i + 1; k < size; k++)
			x[i] -= alf(i,k)* x[k];
	}//end i
	return true;
	
}
/// <summary>
/// Формирование матриц alpha и gamma в компактной схеме исключения
/// Матрица alpha верхняя треугольная, gamma - нижняя треугольная
/// </summary>
/// <param name="alpha">формируемые матрицы</param>
/// <returns>значение определителя матрицы</returns>
double MATRIX::FormMatrixCompactScheme(MATRIX& alpha)
{
	if (this->m_columns != this->m_rows) return NAN;
	int n = this->m_rows;

	for (int i = 0; i < n; i++)
	{
		*(alpha.m_data +i * n) = *(this->m_data + i * n);
		if (i > 0) *(alpha.m_data + i) = *(this->m_data +i) / *this->m_data;
	}

	double sum = 0.0;
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
				sum += *(alpha.m_data +i * n + j) * *(alpha.m_data +j * n + k);
			}
			*(alpha.m_data +i * n + k) = *(this->m_data +i * n + k) - sum;
		}
		else
		{
			sum = 0.0;
			for (int j = 0; j <= i - 1; j++)
			{
				sum += *(alpha.m_data + i * n + j) * *(alpha.m_data +j * n + k);
			}
			if (abs(*(alpha.m_data + i * n + i)) < 1.0e-18)
			{
				return 0.0;
			}
			*(alpha.m_data + i * n + k) = (*(this->m_data +i * n + k) - sum) / *(alpha.m_data +i * n + i);
		}
		k++;
	}

	// вычисление определителя
	double det = 1.0;
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
void CompactSchemeSolve(MATRIX &A, VECTOR& b, VECTOR& x)
{
	if (A.rows() != A.columns()) return;
	int n = A.rows();
	MATRIX alpha(n , n);
	double det = A.FormMatrixCompactScheme(alpha);
	// решение только для неособенной матрицы
	if (abs(det) < 1.0e-36)
		TriangleSolve(alpha, b, x);
	
}
/// <summary>
/// Решение системы линейных алгебраических уравнений A*x = b с верхней треугольной матрицей
/// </summary>
/// <param name="A">верхняя треугольная матрица системы</param>
/// <param name="b">вектор правой части системы</param>
/// <param name="x">решение системы</param>
void TriangleSolve(MATRIX& A, VECTOR& b, VECTOR& x)
{
	int n = b.size();
	VECTOR beta(n);
	double sum = 0.0;
	beta[0] /= A(0,0);
	for (int i = 1; i < n; i++)
	{
		sum = 0.0;
		for (int j = 0; j <= i - 1; j++)
			sum += A(i,j) * beta[j];
		beta[i] = (b[i] - sum) / A(i,i);
	}

	// решение системы уравнений с труегольной матрицей		   
	x[n - 1] = beta[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		sum = 0.0;
		for (int j = n - 1; j > i; j--)
			sum += A(i,j) * x[j];
		x[i] = beta[i] - sum;
	}

 }
/// <summary>
/// Решение СЛАУ с применением QR декомпозиции матрицы системы A
/// </summary>
/// <param name="A">матрица системы линейных уравнений</param>
/// <param name="b">вектор правой части системы линейных уравнений</param>
/// <param name="x">вектор решения системы линейныз уравнений</param>
void QRDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x)
{ 
	if (A.rows() != A.columns() || A.rows() != b.size()) return;
	MATRIX Q(A.rows(), A.columns()), R(A.rows(), A.columns());
	int n = A.rows();
	A.QRDecomposition(Q, R);
	VECTOR beta(Q.Transpose() * b);

	// решение системы уравнений с верхней труегольной матрицей		   
	x[n - 1] = beta[n - 1]/ R(n - 1,n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = n - 1; j > i; j--)
			sum += R(i, j) * x[j];
		x[i] = (beta[i] - sum) / R(i, i);
	}


}
/// <summary>
/// Решение системы уравнений с применением метода LLT декомпозиции A = LLT
/// </summary>
/// <param name="A">матрица</param>
/// <param name="b"></param>
/// <param name="x"></param>
void LLTDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x)
{
	if (A.rows() != A.columns() || A.rows() != b.size()) return;

	MATRIX L(A.rows(), A.columns()), Anorm(A.rows(), A.columns());
	VECTOR bet(A.rows());

	int n = A.rows();
	if (!A.IsSymmetric())
	{
		MATRIX At(A.Transpose());
		Anorm =  At * A; // нормализация матрицы
		bet = At * b;		
	}
	else
	{
		Anorm = A; bet = b;
	}
	Anorm.CholeskyDecomposition(L);
	// решение системы с верхней треугольной матрицей Ly = b
	VECTOR y(n);
	for (int i = 0; i < n; i++)
	{
		y[i] = bet[i];
		for (int k = 0; k < i; k++)
			y[i] -= L(i, k) * y[k];
		y[i] /= L(i,i);
		
	}

	L = L.Transpose();
	// решение системы уравнений с верхней труегольной матрицей L^tx = y   
	x[n - 1] = y[n - 1] / L(n - 1, n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = y[i];
		for (int j = n - 1; j > i; j--)
			x[i] -= L(i, j) * x[j];
		x[i] /= L(i,i);
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
double MATRIX::LUDecomposition(MATRIX& alfa)
{
	double det = 0;
	if (m_rows != m_columns) return det;

	// декомпозиция матрицы
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_columns; j++)
		{
			double val = 0;
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
	for (int i = 0; i < alfa.m_rows; i++)
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

void LUDecompositionSolve(MATRIX& A, VECTOR& b, VECTOR& x)
{
	if (A.rows() != A.columns() || A.rows() != b.size()) return;

	MATRIX alfa(A.rows(), A.columns());
	int n = A.rows();

	double det = A.LUDecomposition(alfa);
	if (abs(det) < 1.0e-36) return;
	// решение СЛАУ
	VECTOR y(n);
	for (int k = 0; k < n; k++)
	{
		y[k] = b[k];
		for (int j = 0; j <= k - 1; j++)
			y[k]-= alfa(k, j) * y[j];
	}

	for (int k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int j = n - 1; j > k; j--)
			x[k] -= alfa(k, j) * x[j];
		x[k] /= alfa(k, k);
	}
	
}

/// <summary>
/// Вычисление минора квадаратной матрицы
/// </summary>
/// <param name="i">строка</param>
/// <param name="j">столбец</param>
/// <returns>значение минора, в случае ошибки NAN</returns>
double MATRIX::Minor(int i, int j)
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
	MATRIX minor(n - 1,n - 1);

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
MATRIX MATRIX::Reverse()
{
	MATRIX A(m_rows, m_columns);
	if (m_rows != m_columns)
	{
		throw "Для вычисления обратной матрицы исходная матрица должна быть квадратной";
		return A;
	}
	double det = Determinant();
	
	if (abs(det) >= 1.0e-36)
	{
		auto Minors = [&](int row_begin, int row_end, int column_begin, int column_end)
			{
				for (int i = row_begin; i < row_end; i++)
					for (int j = column_begin; j < column_end; j++)
					{
						*(A.m_data + A.m_columns * i + j) = Minor(j, i) * ( (i + j) % 2 > 0 ? -1.0 : 1.0)  / det;
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
MATRIX MATRIX::Invert()
{
	MATRIX A(m_rows, m_columns);
	if (m_rows != m_columns)
	{
		throw "Для вычисления обратной матрицы исходная матрица должна быть квадратной";
		return A;
	}

	MATRIX alpha(m_rows, m_columns);
	double det = FormMatrixCompactScheme(alpha);

	if (abs(det) < 1.0e-36)
	{
		throw "Это вырожденная матрица";
		return A;
	}
	
	for (int i = 0; i < m_columns; i++)
	{
		VECTOR e(m_columns), x(m_columns);
		for (int j = 0; j < m_columns; j++)
		{
			if(i == j) e[j] = 1.0;
		}
		TriangleSolve(alpha, e, x);
		A.CopyColumn(x, i);
	}
	
	// итерационное уточнение обратной матрицы
	if (m_rows >= MIN_SIZE_INVERSION)
	{
		MATRIX E(m_rows, m_columns);
		for (int i = 0; i < m_rows; i++)
			*(E.m_data + i * m_columns + i) = 1;
		int iter = 0;
		double norm;
		do
		{
			MATRIX R(E - *this * A);
			A = A * (E + R);
			norm = R.normI();
			if (++iter > MAX_ITER_NUMBER) break;
		} while (norm < 1.0e-17);
	}
	return A;
}
/// <summary>
/// Обращение матриы с пощью LU разложения
/// </summary>
/// <returns>обратную матрицу</returns>
MATRIX MATRIX::InvertLU()
{
	MATRIX D(m_rows, m_columns), alpha(m_rows, m_columns);
	if (m_rows != m_columns) // матрица должна быть квадратной
		return D;
	double det = LUDecomposition(alpha); // декомпозиция A = LU
	if (abs(det) < 1.0e-36)
		return D;
	// формирование элементов обратной матрицы
	int n = m_rows;
	for (int i = n - 1; i >= 0; i--)
		for (int j = n - 1; j >= 0; j--)
		{
			double sum = 0;
			if (i < j)
			{
				for (int k = i + 1; k < n; k++)
					sum -= *(alpha.m_data + i * alpha.m_columns + k) * *(D.m_data + k * D.m_columns + j);
				*(D.m_data + i * D.m_columns + j) = sum / *(alpha.m_data + i * alpha.m_columns + i);

			}
			else if (i == j)
			{
				for (int k = j + 1; k < n; k++)
					sum += *(alpha.m_data + j * alpha.m_columns + k) * *(D.m_data + k * D.m_columns + j);
				*(D.m_data + j * D.m_columns + i) = (1 - sum) / *(alpha.m_data + j * alpha.m_columns + j);

			}
			else
			{
				for (int k = j + 1; k < n; k++)
					sum -= *(alpha.m_data + k * alpha.m_columns + j) * *(D.m_data + i * D.m_columns + k);
				*(D.m_data + i * D.m_columns + j) = sum;

			}
		}
	// итерационное уточнение обратной матрицы
	if (m_rows >= MIN_SIZE_INVERSION)
	{
		MATRIX E(m_rows, m_columns);
		for (int i = 0; i < m_rows; i++)
			*(E.m_data + i * m_columns + i) = 1;
		int iter = 0;
		double norm;
		do
		{
			MATRIX R(E - *this * D);
			D = D * (E + R);
			norm = R.normI();
			if (++iter > MAX_ITER_NUMBER) break;
		} while (norm < 1.0e-17);
	}

	return D;
}

/// <summary>
/// Копировать j-ю колонку матрицы в вектор
/// </summary>
/// <param name="j">номер колонки</param>
/// <returns>вектор содержащий колонку матрицы</returns>
VECTOR MATRIX::CopyColumn2Vector(int j)
{
	VECTOR v(m_rows);
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

complex<double> Polyfun(const complex<double>& z, VECTOR& koeff)
{
	complex<double> poly(0, 0);
	int n = koeff.size();
	poly = -1.0*koeff[n-1] - koeff[n-2] * z - koeff[n-3] * z * z;
	double p = 0.0;
	for (int i = 3; i < n; i++)
	{
		p = (double)i;
		poly -= koeff[n-i-1] * pow(z, p);

	}
	poly += pow(z, (double)n);
	return poly;
}

// Производная характеристического полинома x ^ n + koeff[1] * x ^ (n - 1) + ... + koeff[n-1]
complex<double> PolyfunDerivative(const complex<double>& z, VECTOR &koeff)
{
	complex<double> poly(0, 0);
	int n = koeff.size();
	poly = -1.0*koeff[n-2] - 2.0 * koeff[n-3] * z - 3.0 * koeff[n-4] * z * z;
	double p = 0.0;
	for (int i = 4; i < n; i++)
	{
		p = (double)(i - 1);
		poly -= ((double)i) * koeff[n-i-1] * pow(z, p);

	}
	poly += (double)n * pow(z, (double)(n - 1));
	return poly;
}

#define EPS 1.0e-12
#define MAX_ITER_NUMBER 30000
inline double randomDouble()
{
	return (double)(rand()) / (double)(rand());
}

// Заполнение массива начальных значений корней полинома init_values
// при помощи генератора случайных чисел
// n - степень полинома.
// Если полином нечётной степени, то первый элемент массива - вещественное число

void FillInitialValues(complex<double>* init_values, int n)
{
	srand(1);
	int i = 0;
	while (i < n)
	{
		if (i == 0) init_values[i] = complex<double>(randomDouble(), n % 2 != 0 ? 0 : randomDouble());
		else
		{
			init_values[i] = complex<double>(randomDouble(), randomDouble());
			if (i == n - 1) break;
			i++;
			init_values[i] = n % 2 != 0 ? conj(init_values[i - 1]) : complex<double>(randomDouble(), 0);
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

void Polyroots(VECTOR &koeff, complex<double>* roots)
{
	int n = koeff.size();
	complex<double>* w = new complex<double>[n]; // массив начальных значений корней полинома, далее используется для чисел смещения
	FillInitialValues(roots, n);

	double err = 1.0;
	int iter = 0;
	while (err >= EPS && iter < MAX_ITER_NUMBER)
	{
		for (int i = 0; i < n; i++)
		{
			complex<double> pr(0, 0);
			for (int j = 0; j < n; j++)
			{
				if (i != j) pr += 1.0 / (roots[i] - roots[j]);
			}
			complex<double> pp = Polyfun(roots[i], koeff) / PolyfunDerivative(roots[i], koeff);
			w[i] = pp / (1.0 - pp * pr);
		}

		// проверка максимального модуля чисел смещений
		err = -1.0;
		double err_max = 0.0;

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
void MATRIX::EigenvaluesAndVectorsKrylov(complex<double>* lambda, complex<double> **vect)
{
	if (m_rows != m_columns)
	{
		throw "Матрица не квадратная";
		return;
	}

	// поиск собственных значений
	int n = m_columns;
	VECTOR mcolumn(n);
	MATRIX a(n, n);
	mcolumn[0] = 1.0;	

	for (int i = n-1; i >=0; i--)
	{
		a.CopyColumn(mcolumn, i);
		mcolumn = *this * mcolumn;
	}

	VECTOR p(n); // значения коэффициентов полинома
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


		complex<double> q(1, 0);
		for (int k = 1; k < n; k++)
		{
			mcolumn = a.CopyColumn2Vector(k);
			q = lambda[i] * q - p[k-1];
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
void MATRIX::EigenvaluesAndVectorsLeVerrierFaddeev(complex<double>* lambda, complex<double>** vect)
{
	if (m_rows != m_columns)
	{
		throw "Матрица не квадратная";
		return;
	}
	// поиск собственных значений
	int n = m_columns;
	VECTOR p(n);
	vector<VECTOR> b(n, VECTOR(n));
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
	complex<double>* y = new complex<double>[n];
	
	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < n; k++)
			y[k] = k == 0 ? complex<double>(1, 0) : complex<double>(0, 0);
		for (int j = 1; j < n; j++)
		{
			for (int k = 0; k < n; k++)
				y[k] = lambda[i] * y[k] + b[j-1][k];

		}

		for (int k = 0; k < n; k++)
			vect[i][k] = y[k];
	}
	delete[] y;
}

/// <summary>
/// Обращение матрицы методом Фаддеева Д.К.
/// </summary>
/// <returns>обратную матрицу</returns>
MATRIX MATRIX::InvertFaddev()
{
	MATRIX Ainv(m_rows, m_columns);
	if (m_rows != m_columns)
	{
		throw "Матрица не квадратная";
		return Ainv;
	}
	// поиск собственных значений
	int n = m_columns;
	VECTOR p(n);

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
void MATRIX::Rotate(double* lambda, double** vect)
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
		double val = 0.0;
		for(int i=0;i<n;i++)
			for (int j = 0; j < n; j++)
			{
				double el = abs(*(A.m_data + i * A.m_columns + j));
				if (i != j && el > val)
				{
					val = el;
					i0 = i;
					j0 = j;
				}
			}
		if (val < EPS) break;
		// угол матрицы вращения
		double fi = 0.5 * atan(2.0 * *(A.m_data + i0 * A.m_columns + j0) / ( *(A.m_data + i0 * A.m_columns + i0) - *(A.m_data + j0 * A.m_columns + j0)));
		double cs = cos(fi);
		double ss = sin(fi);

		// U0*U1*...Un
		MATRIX U0(n, n), B(n,n);
		
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

	for (int i = 0; i < n; i++)
	{
		lambda[i] = A(i, i);
		for (int j = 0; j < n; j++)
			vect[i][j] = *(U.m_data + j * U.m_columns + i);
	}

}