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
/// Умнoжение матриц first и second
/// </summary>
/// <param name="first">матрица слева</param>
/// <param name="second">матрица справа</param>
/// <returns>произведение матриц</returns>
vector <SparseElement> SparseMultiply(const vector<SparseElement>& first,
	const vector<SparseElement>& second)
{
	map<int, map<int, double>> matrix;
	vector <SparseElement> result;
	for (const auto& el_first : first)
	{
		map<int, double> mrow;
		if (matrix.count(el_first.row) == 0)
			matrix[el_first.row] = mrow;
		for(const auto& el_second : second)
			if (el_first.column == el_second.row)
			{
				if (matrix[el_first.row].count(el_second.column) == 0)
					matrix[el_first.row][el_second.column] = el_first.value * el_second.value;
				else
					matrix[el_first.row][el_second.column] += el_first.value * el_second.value;
			}
	}

	for (const auto matrix_el : matrix)
	{
		map<int, double> mrow = matrix_el.second;
		int row = matrix_el.first;
		for (const auto& mrow_el : mrow)
		{
			SparseElement res_el{ row, mrow_el.first, mrow_el.second };
			result.push_back(res_el);
		}
	}

	return result;
}
/// <summary>
/// Вывод на консоль матрицы или вектора
/// </summary>
/// <param name="matrix">разреженная матрица</param>
/// <param name="n">порядок матрицы</param>
void PrintMatrix(vector<SparseElement> matrix, int n)
{
	for (int i = 1; i <= n; i++)
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
/// для сортировки матриц и векторов
/// </summary>
/// <param name="one"></param>
/// <param name="two"></param>
/// <returns></returns>
bool MatrixLexiograhicCompare(const SparseElement& one,
	const SparseElement& two)
{
	return one.row < two.row ||
		(one.row == two.row && one.column < two.column);
}
bool VectorCompare(const SparseElement& one, const SparseElement& two)
{
	return one.row < two.row;
}

/// <summary>
/// Погрешность решения - норма невязки A*x - b решения СЛАУ
/// </summary>
/// <param name="A">разреженная матрица СЛАУ</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">полученное решение СЛАУ</param>
/// <param name="n">порядок матрицы и размерность вектора правой части</param>
/// <returns>погрешность</returns>
double ErrorMeasure(vector<SparseElement>& A, 
					vector<SparseElement>& b, 
					vector<SparseElement>& x,
					int n)
{
	vector <SparseElement> ax, v_err;
	if(!is_sorted(A.begin(), A.end(), 
		MatrixLexiograhicCompare))
	sort(A.begin(), A.end(), MatrixLexiograhicCompare);

	if(!is_sorted(x.begin(), x.end(), VectorCompare))
		sort(x.begin(), x.end(), VectorCompare);
	if (!is_sorted(b.begin(), b.end(), VectorCompare))
		sort(b.begin(), b.end(), VectorCompare);

	// A*x
	SparseElement el;
	ax = SparseMultiply(A, x);

	for (const auto& elem_ax : ax)
		cout << elem_ax.row  << '\t' << elem_ax.value << endl;
	cout << "***" << endl;
	for (const auto& elem_b : b)
		cout << elem_b.row << '\t' << elem_b.value << endl;
	cout << "***" << endl;
	for (const auto& elem_x : x)
		cout << elem_x.row << '\t' << elem_x.value << endl;


	// A*x - b
	for (int i = 1; i <= n; i++)
	{
		double val = FindElement(ax, i) - FindElement(b, i);
		if (abs(val) <= DBL_MIN) continue;
		el.row = i;
		el.column = 1;
		el.value = val;
		v_err.push_back(el);
	}

	//cout << "***" << endl;
	//for (const auto& elem_err : v_err)
	//	cout << elem_err.row << '\t' << elem_err.column << '\t' << elem_err.value << endl;

	return SparseNormv(v_err);
}

/// <summary>
/// Решение СЛАУ методом вращения c разреженной матрицей, решение сводится
/// к решению СЛАУ с верхней треугольной матрицей
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <param name="n">порядок матрицы и размерность вектора правой части</param>
/// <returns>-1 - деление на ноль, нулевой элемент на главной диагонали
/// -2 - число элементов в строке матрицы не совпадает с размерностью вектора
/// или матрица не квадратная
/// 0 - превышено максимальное число итераций, 
/// 1 - решение завершено успешно</returns>
int SparseRotationSolve(const vector<SparseElement>& A, 
						const vector<SparseElement>& b, 
						vector<SparseElement>& x,
						int n)
{
	// верхняя тругольная матрица
	vector<SparseElement> T; // на текущей итерации

	vector<SparseElement> T0 = A; // на предыдущей итерации

	// преобразованный вектор правой части
	vector<SparseElement> bet0 = b; // на предыдущей итерации

	vector<SparseElement> bet; // на текущей итерации

	double val = 0;
	int i0 = 0, j0 = 0;
	int iter_num = 0; // число итерации, ограничивается максимальным числом итераций
	while (true)
	{
		vector<SparseElement>::iterator iter;
		// любой элемент ниже главной диагонали
		for (iter = T0.begin(); iter != T0.end(); ++iter)
			if (iter->column < iter->row) break;
		// поиск позиции максимального по модулю элемента ниже главной диагонали
		 iter = max_element(iter, T0.end(),
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

		// угол матрицы вращения, находится по условию T(k-1)(j0, i0) = 0
		double t2 = FindElement(T0, i0, i0);
		if (t2 == 0.0) return -1;

		double fi = atan(val / t2);
		double cs = cos(fi);
		double ss = sin(fi);

		// T(k) = U*T(k-1), U - матрица вращения
		T = T0;
		// заполнение строк T(i0,k) и T(j0,k)
		double t_i0 = 0.0, t_j0 = 0.0;
		int k = 0;
		for (auto iter_t = T.begin(); iter_t != T.end(); ++iter_t)
		{			
			k = iter_t->column;

			if (iter_t->row == i0 || iter_t->row == j0)
			{
				// T0(i0,k)
				t_i0 = FindElement(T0, i0, k);

				// T0(j0,k)
				t_j0 = FindElement(T0, j0, k);
			}

			// T(i0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
			if (iter_t->row == i0)
				iter_t->value = t_i0 * cs + t_j0 * ss;

			// T(j0,k) = -T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
			if (iter_t->row == j0)
				iter_t->value = -t_i0 * ss + t_j0 * cs;
		}

		// bet(k) = U*bet(k-1)
		bet = bet0;
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
		if (abs(tval) <= DBL_MIN) return -1;

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