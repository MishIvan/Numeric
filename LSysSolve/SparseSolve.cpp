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
double FindElement(const vector<SparseElement>& matrix,int _row, int _column)
{
	auto iter = find_if( matrix.begin(), matrix.end(),
		[_row,_column](const SparseElement &el) { return el.row == _row && el.column == _column; }
	);
	if (iter != matrix.end())
		return iter->value;
	else
		return 0.0;
} 

/// <summary>
/// Присвоить значение элемента матрицы: если найден - присвоить значение,
/// если нет и значение ненулевое - создать элемент
/// </summary>
/// <param name="matrix">матрица или вектор</param>
/// <param name="_row">строка</param>
/// <param name="_column">столбец</param>
/// <param name="_value">значение</param>
void SetValue(vector<SparseElement>& matrix, int _row, int _column, double _value)
{
	auto iter = find_if(matrix.begin(), matrix.end(),
		[_row, _column](const SparseElement& el) { return el.row == _row && el.column == _column; }
	);
	if (iter != matrix.end())
		iter->value = _value;
	else
	{
		if (abs(_value) <= DBL_MIN) return;
		SparseElement el{ _row,_column, _value };
		matrix.push_back(el);
	}
}

/// <summary>
/// Транспонированную матрицу
/// </summary>
/// <param name="matrix">исходная матрица</param>
/// <param name="n">порядок матрицы</param>
/// <returns>транспонированну  matrix матрицу</returns>
vector<SparseElement> SparseTranspose(const vector<SparseElement>& matrix)
{
	vector <SparseElement> result;
	for (const auto& elem : matrix)
	{
		SetValue(result, elem.column, elem.row, elem.value);		
	}

	return result;
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
		for(const auto& el_second : second)
			if (el_first.column == el_second.row)
			{
				matrix[el_first.row][el_second.column] += el_first.value * el_second.value;
			}
	}

	for (const auto& matrix_el : matrix)
	{
		int row = matrix_el.first;
		for (const auto& mrow_el : matrix_el.second)
		{
			SparseElement res_el{ row, mrow_el.first, mrow_el.second };
			result.push_back(res_el);
		}
	}

	return result;
}

/// <summary>
/// Разность двух матриц first и second
/// </summary>
/// <param name="first">разреженная матрица</param>
/// <param name="second">разреженная матрица</param>
/// <returns>разреженную матрицу first - second</returns>
vector<SparseElement> SparseDifference(const vector<SparseElement>& first,
	const vector<SparseElement>& second)
{
	map<int, map<int, double>> buffer;
	
	for (const auto& elem_first : first)
		buffer[elem_first.row][elem_first.column]  = elem_first.value;
	for (const auto& elem_second : second)
		buffer[elem_second.row][elem_second.column] -= elem_second.value;

	vector<SparseElement> result;
	
	for (const auto& buffer_el : buffer)
	{
		int row = buffer_el.first;
		for (const auto& mrow_el : buffer_el.second)
		{
			if (abs(mrow_el.second) <= DBL_MIN) continue;
			result.push_back({ row, mrow_el.first, mrow_el.second });
		}
	}

	return result;
}

map<int, map<int, double>> Matrix2Map(const vector<SparseElement> matrix)
{
	map<int, map<int, double>> result;
	for (const auto& elem : matrix)
		result[elem.row][elem.column] = elem.value;
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

double SparseScalarProduct(const vector<SparseElement>& vect1,
	const vector<SparseElement>& vect2)
{
	vector<SparseElement> value = SparseMultiply(SparseTranspose(vect1), vect2);
	return value[0].value;
}
/// <summary>
/// Евклидова (сферическая норма) вектора
/// </summary>
/// <param name="vect">разреженный вектор</param>
/// <returns>значение нормы</returns>
double SparseNormv(const vector<SparseElement>& vect) 
{
	return sqrt(SparseScalarProduct(vect, vect));
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
/// <returns>погрешность</returns>
double ErrorMeasure(const vector<SparseElement>& A, 
					const vector<SparseElement>& b, 
					const vector<SparseElement>& x)
{
	// A*x
	vector <SparseElement> ax = SparseMultiply(A, x);

	return SparseNormv(SparseDifference(ax, b));
}

/// <summary>
/// Решение СЛАУ методом вращения c разреженной матрицей, решение сводится
/// к решению СЛАУ с верхней треугольной матрицей
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <param name="n">порядок матрицы и размерность вектора правой части</param>
/// <returns>
/// -1 - деление на ноль, нулевой элемент на главной диагонали
/// -2 - число элементов в строке матрицы не совпадает с размерностью вектора
/// или матрица не квадратная
/// -3 - элемент не найден
///   0 - превышено максимальное число итераций, 
/// 1 - решение завершено успешно
/// </returns>
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
		
		j0 = iter->column;
		i0 = iter->row;
		val = iter->value;

 		if (abs(val) < DBL_EPSILON) break;
		// угол матрицы вращения, находится по условию T(k-1)(i0, j0) = 0
		double t2 = FindElement(T0, j0, j0);
		if (t2 == 0.0) return -1;

		double fi = atan(-val / t2);
		double cs = cos(fi);
		double ss = sin(fi);

		// T(k) = U*T(k-1), U - матрица вращения
		T = T0;
		// заполнение строк T(i0,k) и T(j0,k)	
		double t_i0 = 0.0, t_j0 = 0.0;
		for (int k = 1; k <= n; k++)
		{
			// T0(i0,k)
			t_i0 = FindElement(T0, i0, k);

			// T0(j0,k)
			t_j0 = FindElement(T0, j0, k);

			// T(i0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
			val = t_i0 * cs + t_j0 * ss;
			SetValue(T, i0, k, val);

			// T(j0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
			val = -t_i0 * ss + t_j0 * cs;
			SetValue(T, j0, k, val);
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
			if (tval == 0) continue;
			xval = FindElement(x, j);
			if (xval == 0) continue;
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

/// <summary>
/// Решение системы уравнений A*x = b методом  релаксации
/// для симметричной положительно определённой матрицы A
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <param name="omega">релаксационный множитель</param>
/// <returns>
///  1 - если сходимость была достигнута,
/// -1 - элемент матрицы на главной диагонали равен нулю,
///  0 - превышено максимальное число итераций
/// </returns>
int SparseRelaxation(vector<SparseElement>& A,
					vector<SparseElement>& b, 
					vector<SparseElement>& x, 
					double omega)
{
	vector<SparseElement> x0;
	x0.push_back({ 1, 1, 1.0 });
	map<int, map<int, double>> Amap;
	Amap = Matrix2Map(A);

	int iter = 0; // число итераций
	double err_norm = 0.0; // погрешность
	do
	{		
		double sum = 0.0, val_d = 0.0;

		for (const auto& elem_a : Amap)
		{	
			map<int,double> Arow = elem_a.second;
			double sum = 0.0;
			int row = elem_a.first;
			
			for (const auto& elem_x : x)
			{
				if (elem_x.row <= row - 1 && Arow.count(elem_x.row) > 0)
					sum += Arow[elem_x.row] * elem_x.value;
			}

			for (const auto& elem_x : x0)
			{ 
				if (elem_x.row >= row && Arow.count(elem_x.row) > 0)
					sum += Arow[elem_x.row] * elem_x.value;
			}

			if (Arow.count(row) == 0) return -1;
			double val_d = Amap[row][row];	
			double val = FindElement(x0, row) - omega * (sum - FindElement(b, row)) / val_d;
			SetValue(x, row, 1, val);
		}


		// норма разности между значением на текущей и предыдущей итерации
		err_norm = SparseNormv(SparseDifference(x, x0));
		if (++iter > MAX_ITERATION_NUMBER) return 0;

		if (!is_sorted(x.begin(), x.end(), VectorCompare))
			sort(x.begin(), x.end(), VectorCompare);
		x0 = x;

	} while (err_norm > DBL_EPSILON*100);

	return 1;
}

/// <summary>
/// Решение системы уравнений A*x = b методом градиентного спуска
/// для симметричной положительно определённой матрицы A
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <returns>
///  1 - если сходимость была достигнута,
///  0 - превышено максимальное число итераций
/// </returns>
int SparseGradientDescent(const vector<SparseElement>& A, 
	const vector<SparseElement>& b, 
	vector<SparseElement>& x)
{

	// градиент: 2*(Ax0 - b)
	auto Gradient = [A, b](const vector<SparseElement>& x)
		{		
			vector<SparseElement> result = SparseDifference(SparseMultiply(A, x), b);
			for (auto& elem : result)
				elem.value *= 2.0;
			return result;
		};

	// A*(Ax0 - b)
	auto VectorRight = [A, b](const vector<SparseElement>& x)
		{
			vector<SparseElement> ax = SparseDifference(SparseMultiply(A, x), b);
			vector<SparseElement> result = SparseMultiply(A, ax);
			return result;
		};

	vector<SparseElement> x0;
	x0.push_back({ 1, 1, 1.0 });
	vector<SparseElement> grad;

	int iter = 0;
	double err_norm = 0.0, t = 0.0;
	do
	{
		grad = Gradient(x0);
		double sc1 = SparseScalarProduct(grad, grad) * 0.25;

		double sc2 = SparseScalarProduct(grad, VectorRight(x0));
		t = sc1 / sc2;

		for (auto& elem : grad)
			elem.value *= t;

		x = SparseDifference(x0, grad);

		// норма разности между значением на текущей и предыдущей итерации
		err_norm = SparseNormv(SparseDifference(x, x0));
		if (++iter > MAX_ITERATION_NUMBER) return 0;
		x0 = x;
	} while (err_norm > DBL_EPSILON);
	return 1;
}