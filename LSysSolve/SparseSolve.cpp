#include "SparseSolve.h"
#include "SysSolve.h"
using namespace std;

/// <summary>
/// Поиск значения элемента на строке i и столбце j
/// </summary>
/// <param name="matrix">разреженная матрица</param>
/// <param name="row">строка</param>
/// <param name="column">столбец</param>
/// <returns>значение элемента</returns>
double FindElement(const vector<SparseElement>& matrix,int _row, int _column = 1)
{
	auto iter = find_if(matrix.begin(), matrix.end(),
		[_row,_column](const SparseElement &el) { return el.row == _row && el.column == _column; }
	);
	if (iter != matrix.end())
		return iter->value;
	else
		return NAN;
}

/// <summary>
/// Вычисление размера разреженной матрицы
/// </summary>
/// <param name="matrix">разреженная матрица</param>
/// <param name="rows">число строк</param>
/// <param name="columns">число столбцов</param>
void MatrixSize(const vector<SparseElement>& matrix, int& rows, int& columns)
{
	// максимальное значение строки
	vector<SparseElement>::const_iterator iter = max_element(matrix.begin(), matrix.end(),
		[](const SparseElement& one, const SparseElement& two)
		{
				return one.row < two.row;
		});
	if (iter == matrix.end()) { rows = -1; columns = -1; return; }
	rows = iter->row;

	// максимальное значение колонки
	iter = max_element(matrix.begin(), matrix.end(),
		[](const SparseElement& one, const SparseElement& two)
		{
			return one.column < two.column;
		});
	if (iter == matrix.end()) { rows = -1; columns = -1; return; }
	columns = iter->column;

}
/// <summary>
/// Степень заполнения матрицы в процентах
/// </summary>
/// <param name="matrix">матрица</param>
/// <returns>процент заполнения матрицы</returns>
double FullnessDegree(const vector<SparseElement> matrix)
{
	int m = 0, n = 0;
	MatrixSize(matrix, m, n);
	if (m <= 0 || n <= 0) 
		return NAN;
	else 
		(double)matrix.size() * 100.0 / (double)(m * n);
}

/// <summary>
/// Евклидова (сферическая норма) вектора
/// </summary>
/// <param name="vect">разреженный вектор</param>
/// <returns>значение нормы</returns>
double SparseNormv(const vector<SparseElement>& vect) 
{
	if (vect.empty()) return 0.0;
	double norm = 1.0;
	for (const auto& elem : vect)
		if(elem.column == 1)
			norm *= elem.value * elem.value;
	return sqrt(norm);
}
/// <summary>
/// Погрешность решения - норма невязки A*x - b решения СЛАУ
/// </summary>
/// <param name="A">разреженная матрица СЛАУ</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">полученное решение СЛАУ</param>
/// <returns>погре</returns>
double ErrorMeasure(const vector<SparseElement>& A, const vector<SparseElement>& b, const vector<SparseElement>& x)
{
	double err = NAN;
	vector <SparseElement> ax, v_err;
	
	// A*x
	double sum = 0;
	for (const auto& elem_A : A)
	{
		for (const auto& elem_x : x)
			if (elem_A.column == elem_x.row)
				sum += elem_A.value * elem_x.value;
		SparseElement el_ax{ elem_A.row,1,sum };
		ax.push_back(el_ax);
	} 
		
	
	//  A*x - b доработать!!! 
	SparseElement elem_err;
	for (const auto& elem_ax : ax)
	{
		vector<SparseElement>::const_iterator iter;
		for (iter = b.begin(); iter != b.end(); ++iter)
			if (elem_ax.row == iter->row && elem_ax.column == 1 && iter->column == 1)
			{
				elem_err.row = elem_ax.row,
				elem_err.column = 1;
				elem_err.value = elem_ax.value - iter->value;
				v_err.push_back(elem_err);
				break;
			}
			
		if (iter == b.end()) 
		{
			elem_err.row = elem_ax.row,
			elem_err.column = 1;
			elem_err.value = elem_ax.value;
			v_err.push_back(elem_err);
		}
	}

	return SparseNormv(v_err);
}

/// <summary>
/// Решение СЛАУ методом вращения c разреженной матрицей, решение сводится
/// к решению СЛАУ с верхней треугольной матрицей
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <returns>-1 - элемент не найден,
/// -2 - число элементов в строке матрицы не совпадает с размерностью вектора
/// или матрица не квадратная
/// 0 - превышено максимальное число итераций, 
/// 1 - решение завершено успешно</returns>
int SparseRotationSolve(const vector<SparseElement>& A, 
						const vector<SparseElement>& b, 
						vector<SparseElement>& x)
{
	int m, n, p, q; // размерность матрицы и вектора
	// проверка совпадения числа элементов в строке 
	// матрицы A с размерностью вектора b 
	MatrixSize(A, m, n);
	if (m <= 0 || n <= 0 || m != n) return -2;
	MatrixSize(b, p, q);
	if (p <= 0 || q !=1 || p != n) return -2;

	// верхняя тругольная матрица
	vector<SparseElement> T = A; // на текущей итерации

	vector<SparseElement> T0 = A; // на предыдущей итерации

	// преобразованный вектор правой части
	vector<SparseElement> bet0 = b; // на предыдущей итерации

	vector<SparseElement> bet = b; // на текущей итерации

	double val = 0;
	int i0 = 0, j0 = 0;
	int iter_num = 1; // число итерации, ограничивается максимальным числом итераций
	while (true)
	{
		// поиск позиции максимального по модулю элемента ниже главной диагонали
		vector<SparseElement>::iterator iter = max_element(T.begin(), T.end(),
			[](const SparseElement& one, const SparseElement& two)
			{
				if (one.row < one.column)
				{
					return abs(one.value) < abs(two.value);
				}
				else
					return false;
			}
		); 
		if (iter == T.end()) return -1;
		
		j0 = iter->row;
		i0 = iter->column;
		val = iter->value;
		
		if (val < 1.0e-17) break;

		// угол матрицы вращения, находится по условию T(j0, i0) = 0
		double t1 = FindElement(T, j0, i0);
		if (isnan(t1)) return -1;
		double t2 = FindElement(T, i0, i0);
		if (isnan(t2)) return -1;

		double fi = atan(t1 / t2);
		double cs = cos(fi);
		double ss = sin(fi);

		// T(k) = U*T(k-1), U - матрица вращения
		// bet(k) = U*bet(k-1)
		T = T0;
		bet = bet0;

		// заполнение строк T(i0,k) и T(j0,k)
		for (auto iter_t = T.begin(); iter_t != T.end(); ++iter_t)
		{
			double t_i0 = NAN, t_j0 = NAN;
			int k = iter_t->column;

			if (iter_t->row == i0)
			{
				// T0(i0,i)
				t_i0 = FindElement(T0, i0, k);
				if (isnan(t_i0)) return -1;		

				// T0(j0,i)
				t_j0 = FindElement(T0, j0, k);
				if (isnan(t_i0)) return -1;

				// T(i0,k) = T0(i0,k)*cos(phi) + T0(j,k)*sin(phi)
				iter_t->value = t_i0 * cs + t_j0 * ss;
			}

			if (iter_t->row == j0)
			{
				// T0(i0,i)
				t_i0 = FindElement(T0, i0, k);
				if (isnan(t_i0)) return -1;

				// T0(j0,i)
				t_j0 = FindElement(T0, j0, k);
				if (isnan(t_i0)) return -1;

				// T(j0,k) = T0(i0,k)*cos(phi) + T0*sin(phi)
				iter_t->value = -t_i0 * ss + t_j0 * cs;
			}
		}
		
		// b0(i0)
		double b_i0 = FindElement(bet0, i0, 1);
		if (isnan(b_i0)) return -1;

		// b0(j0)
		double b_j0 = FindElement(bet0, j0, 1);
		if (isnan(b_j0)) return -1;

		for (auto iter_b = bet.begin(); iter_b != bet.end(); ++iter_b)
		{
			if (iter_b->row == i0)
				iter->value = b_i0 * cs + b_j0 * ss;
			if (iter_b->row == j0)
				iter->value = -b_i0 * ss + b_j0 * cs;
		}
		// результаты для следующей итерации
		T0 = T;
		bet0 = bet;
		if (++iter_num > MAX_ITERATION_NUMBER) return 0;
	}
	T0.clear();
	bet0.clear();

	// вычисление вектора решения
	// решение системы уравнений с верхней треугольной матрицей T(n)*x = bet(n)   

	if (!x.empty()) x.clear();
	val = FindElement(T, n, n);
	if (isnan(val)) return -1;

	SparseElement el{ n ,1,val };
	x.push_back(el);
	for (int i = n - 1; i >= 1; i--)
	{ 

		el.row = i;
		el.column = 1;
		val = FindElement(bet, i);
		if (isnan(val)) return -1;
	
		double tval = NAN, xval = NAN;
		for (int j = n; j > i; j--)
		{
			tval = FindElement(T, i, j);
			if (isnan(tval)) return -1;

			xval = FindElement(x, j);
			if (isnan(xval)) return -1;

			val -= tval * xval;
		}
		tval = FindElement(T, i, i);
		if (isnan(tval)) return -1;

		val /= tval;
		x.push_back(el);
	}
	
	reverse(x.begin(), x.end());
	T.clear();
	bet.clear();	
	return 1;
}