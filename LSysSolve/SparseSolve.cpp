#include "SparseSolve.h"
#include "SysSolve.h"
using namespace std;

/// <summary>
/// Поиск значения элемента в строке i и столбце j
/// </summary>
/// <param name="matrix">разреженная матрица</param>
/// <param name="_row">строка</param>
/// <param name="_column">столбец</param>
/// <returns>ненулевое значение элемента, 0 - если элемент не найден</returns>
double FindElement(const vector<SparseElement>& matrix,int _row, int _column = 1)
{
	auto iter = find_if(matrix.begin(), matrix.end(),
		[_row,_column](const SparseElement &el) { return el.row == _row && el.column == _column; }
	);
	if (iter != matrix.end())
		return iter->value;
	else
		return 0.0;
}

/// <summary>
/// Вычисление размера разреженной матрицы
/// предполагается, что все элементы на гланой диагонали ненулевые
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
		return (double)matrix.size() * 100.0 / (double)(m * n);
}
/// <summary>
/// Вывод на консоль матрицы или вектора
/// </summary>
/// <param name="matrix">разреженная матрица</param>
void PrintMatrix(vector<SparseElement> matrix)
{
	int m = 0, n = 0;
	MatrixSize(matrix, m, n); // определение размера матрицы
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			double val = FindElement(matrix, i, j);
			cout << val << '\t';
		}
		cout << endl;
	}
}
/// <summary>
/// Евклидова (сферическая норма) вектора
/// </summary>
/// <param name="vect">разреженный вектор</param>
/// <returns>значение нормы</returns>
double SparseNormv(const vector<SparseElement>& vect) 
{
	double norm = 0.0;
	for (const auto& elem : vect)
		if(elem.column == 1)
			norm += elem.value * elem.value;
	return sqrt(norm);
}
/// <summary>
/// Погрешность решения - норма невязки A*x - b решения СЛАУ
/// </summary>
/// <param name="A">разреженная матрица СЛАУ</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">полученное решение СЛАУ</param>
/// <returns>погрешность</returns>
double ErrorMeasure(const vector<SparseElement>& A, const vector<SparseElement>& b, const vector<SparseElement>& x)
{
	vector <SparseElement> ax, v_err;
	
	// A*x
	double sum = 0;
	for (const auto& elem_A : A)
	{
		sum = 0.0;
		for (const auto& elem_x : x)
			if (elem_A.column == elem_x.row)
				sum += elem_A.value * elem_x.value;
		SparseElement el_ax{ elem_A.row,1,sum };
		ax.push_back(el_ax);
	} 
	PrintMatrix(ax);
	cout << "***" << endl;
	PrintMatrix(b);
	//  A*x - b 
	SparseElement elem_err;
	for (const auto& elem_ax : ax)
	{
		for (auto iter = b.begin(); iter != b.end(); ++iter)
			if (elem_ax.row == iter->row && elem_ax.column == 1 && iter->column == 1)
			{
				elem_err.row = elem_ax.row,
				elem_err.column = 1;
				elem_err.value = elem_ax.value - iter->value;
				v_err.push_back(elem_err);
			}
	}
			
	for (const auto& elem_ax : ax)
	{
		if (FindElement(v_err, elem_ax.row) == 0)
		{
			SparseElement el_err{ elem_ax.row,1,elem_ax.value };
			v_err.push_back(el_err);
		}
	}

	for (const auto& elem_b : b)
	{
		if (FindElement(v_err, elem_b.row) == 0)
		{
			SparseElement el_err{ elem_b.row,1,-elem_b.value };
			v_err.push_back(el_err);
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
/// <returns>-1 - деление на ноль, нулевой элемент на главной диагонали
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
	int iter_num = 0; // число итерации, ограничивается максимальным числом итераций
	while (true)
	{
		vector<SparseElement>::iterator iter;
		// любой элемент ниже главной диагонали
		for (iter = T.begin(); iter != T.end(); ++iter)
			if (iter->column < iter->row) break;
		// поиск позиции максимального по модулю элемента ниже главной диагонали
		 iter = max_element(iter, T.end(),
			[](const SparseElement& one, const SparseElement& two)
			{
				if (one.column < one.row && two.column < two.row)
					return abs(one.value) < abs(two.value);
				else
					return false;
			}
		); 
		
		j0 = iter->row;
		i0 = iter->column;
		val = iter->value;
		
		if (abs(val) < DBL_EPSILON) break;

		// угол матрицы вращения, находится по условию T(j0, i0) = 0
		double t2 = FindElement(T, i0, i0);
		if (t2 == 0.0) return -1;

		double fi = atan(val / t2);
		double cs = cos(fi);
		double ss = sin(fi);

		// T(k) = U*T(k-1), U - матрица вращения
		// bet(k) = U*bet(k-1)
		T = T0;
		bet = bet0;

		// заполнение строк T(i0,k) и T(j0,k)
		double t_i0 = 0.0, t_j0 = 0.0;
		int k = 0;
		for (auto iter_t = T.begin(); iter_t != T.end(); ++iter_t)
		{			
			k = iter_t->column;

			if (iter_t->row == i0)
			{
				// T0(i0,k)
				t_i0 = FindElement(T0, i0, k);

				// T0(j0,k)
				t_j0 = FindElement(T0, j0, k);

				// T(i0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
				iter_t->value = t_i0 * cs + t_j0 * ss;
			}

			if (iter_t->row == j0)
			{
				// T0(i0,k)
				t_i0 = FindElement(T0, i0, k);

				// T0(j0,k)
				t_j0 = FindElement(T0, j0, k);

				// T(j0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
				iter_t->value = -t_i0 * ss + t_j0 * cs;
			}
		}
		
		// b0(i0)
		double b_i0 = FindElement(bet0, i0);

		// b0(j0)
		double b_j0 = FindElement(bet0, j0);

		for (auto iter_b = bet.begin(); iter_b != bet.end(); ++iter_b)
		{
			if (iter_b->row == i0)
				iter_b->value = b_i0 * cs + b_j0 * ss;
			if (iter_b->row == j0)
				iter_b->value = -b_i0 * ss + b_j0 * cs;
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

	SparseElement el;
	for (int i = n; i >= 1; i--)
	{ 
		val = FindElement(bet, i);
	
		double tval = 0.0, xval = 0.0;
		for (int j = i + 1; j <= n; j++)
		{
			tval = FindElement(T, i, j);
			xval = FindElement(x, j);
			val-= tval * xval;
		}

		tval = FindElement(T, i, i);
		if (tval == 0.0) return -1;

		if (abs(val) <= DBL_MIN) continue;

		el.row = i;
		el.column = 1;
		el.value = val/tval;
		x.push_back(el);
	}
	
	reverse(x.begin(), x.end());
	T.clear();
	bet.clear();	
	return 1;
}