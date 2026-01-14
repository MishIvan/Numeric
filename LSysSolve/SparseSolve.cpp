#include "SparseSolve.h"
#include "SysSolve.h"
using namespace std;

/// <summary>
/// Поиск значения элемента на строке i и столбце j
/// </summary>
/// <param name="matrix">разреженная матрица</param>
/// <param name="i">строка</param>
/// <param name="j">столбец</param>
/// <returns>значение элемента</returns>
double FindElement(const vector<SparseElement>& matrix,int i, int j = 1)
{
	auto iter = find_if(matrix.begin(), 
		matrix.end(),
		[i,j](const SparseElement &el) { return el.row == i && el.column == j; }
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
/// Решение СЛАУ методом вращения c разреженной матрицей, решение сводится
/// к решению СЛАУ с верхней треугольной матрицей
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <returns>-1 - элемент не найден,
/// 0 - превышено максимальное число итераций, 
/// 1 - решение завершено успешно</returns>
int SparseRotationSolve(const vector<SparseElement>& A, 
						const vector<SparseElement>& b, 
						vector<SparseElement>& x)
{

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
		/*for (int i = 0; i < n; i++)
		{
			*(T + i0 * n + i) = *(T0 + i0 * n + i) * cs +
				*(T0 + j0 * n + i) * ss;
			*(T + j0 * n + i) = (-1.0) * *(T0 + i0 * n + i) * ss +
				*(T0 + j0 * n + i) * cs;

			*(bet + i0) = *(bet0 + i0) * cs +
				*(bet0 + j0) * ss;
			*(bet + j0) = (-1.0) * *(bet0 + i0) * ss +
				*(bet0 + j0) * cs;
		} */
		// T(i0,i)
		iter = find_if(T.begin(), T.end(),
			[i0](const SparseElement& el) { return el.row == i0; }
			);
		if (iter == T.end()) return -1;
		
		// T0(i0,i)
		auto iter_i0 = find_if(T0.begin(), T0.end(),
			[i0](const SparseElement& el) { return el.row == i0; }
		);
		if (iter_i0 == T.end()) return -1;

		// T0(j0,i)
		auto iter_j0 = find_if(T0.begin(), T0.end(),
			[j0](const SparseElement& el) { return el.row == j0; }
		);
		if (iter_j0 == T.end()) return -1;
		 
		// T(i0,k) = T0(i0,k)*cos(phi) + T0*sin(phi)
		iter->value = iter_i0->value * cs + iter_j0->value * ss;

		// результаты для следующей итерации
		T0 = T;
		bet0 = bet;
		if (++iter_num > MAX_ITERATION_NUMBER) return 0;
	}
	T0.clear();
	bet0.clear();

	// вычисление вектора решения
	// решение системы уравнений с верхней треугольной матрицей T(n)*x = bet(n)   
	/* *(x + n - 1) = *(bet + n - 1) / *(T + (n - 1) * n + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		*(x + i) = *(bet + i);
		for (int j = n - 1; j > i; j--)
			*(x + i) -= *(T + i * n + j) * *(x + j);
		*(x + i) /= *(T + i * n + i);
	}*/
	return 1;
}