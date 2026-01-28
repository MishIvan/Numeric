#include "SparseMapSolve.h"

using namespace std;

/// <summary>
/// Транспонированную матрицу
/// </summary>
/// <param name="matrix">исходная матрица</param>
/// <returns>транспонированну  matrix матрицу</returns>
map<int,map<int,double>> SparseTranspose(const map<int, map<int, double>>& matrix)
{
	map<int, map<int, double>> result;
	for (const auto& elem_row : matrix)
	{
		map<int, double> mrow = elem_row.second;
		for (const auto& elem_column : mrow)
			result[elem_column.first][elem_row.first] = elem_column.second;
	}

	return result;
}

/// <summary>
/// Умнoжение матрицы first и вектора second
/// </summary>
/// <param name="first">матрица</param>
/// <param name="second">вектор</param>
/// <returns>вектор - результат умножения</returns>
map<int, double> SparseMultiply(const map<int, map<int, double>>& first,
	const map<int, double>& second)
{
	map<int, double> result;

	for (const auto& elem_row_first : first)
	{
		map<int, double> row_first = elem_row_first.second;
		for (const auto& elem_column_first : row_first)

			for (const auto& elem_second : second)
				if (elem_second.first == elem_column_first.first)
					result[elem_row_first.first] += elem_column_first.second * elem_second.second;
	}

	return result;
}

/// <summary>
/// Умнoжение матриц first и second
/// </summary>
/// <param name="first">матрица слева</param>
/// <param name="second">матрица справа</param>
/// <returns>матрица - результат умножения</returns>
map<int,map<int, double>> SparseMultiply(const map<int, map<int, double>>& first,
	const map<int, map<int, double>>& second)
{
	map<int,map<int, double>> result;

	for (const auto& elem_row_first : first)
	{
		map<int, double> row_first = elem_row_first.second;
		for (const auto& elem_column_first : row_first)
		{
			for (const auto& elem_row_second : second)
			{
				map<int, double> row_second = elem_row_first.second;
				for(const auto& elem_column_second : row_second)
					if (elem_column_first.first == elem_row_second.first)
						result[elem_row_first.first][elem_column_second.first] += 
							elem_column_first.second * elem_column_second.second;
			}
		}

	}

	return result;
}

/// <summary>
/// Разность двух векторов first и second
/// </summary>
/// <param name="first">вектор</param>
/// <param name="second">вектор</param>
/// <returns>разность векторов first и second</returns>
map<int, double> SparseDifference(const map<int, double>& first,
	const map<int, double>& second)
{
	map<int, double> result;

	for (const auto& row_first : first)
		result[row_first.first] = row_first.second;
	for (const auto& row_second : second)
		result[row_second.first] -= row_second.second;

	return result;
}

/// <summary>
/// Скалярное произведение vect1 и vect2
/// </summary>
/// <param name="vect1">вектор</param>
/// <param name="vect2">вектор</param>
/// <returns>значение скалярного произведения vect1 и vect2</returns>
double SparseScalarProduct(const map<int, double>& vect1,
	const map<int, double>& vect2)
{
	double result = 0.0;

	for (const auto& row_vect1 : vect1)
	{
		for (const auto& row_vect2 : vect2)
			if (row_vect1.first == row_vect2.first)
				result += row_vect1.second * row_vect2.second;
	}
	
	return result;
}
/// <summary>
/// Евклидова (сферическая норма) вектора
/// </summary>
/// <param name="vect">разреженный вектор</param>
/// <returns>значение нормы</returns>
double SparseNormv(const map<int, double>& vect)
{
	return sqrt(SparseScalarProduct(vect, vect));
}

/// <summary>
/// Погрешность решения - норма невязки A*x - b решения СЛАУ
/// </summary>
/// <param name="A">разреженная матрица СЛАУ</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">полученное решение СЛАУ</param>
/// <returns>значение погрешность</returns>
double ErrorMeasure(const map<int,map<int, double>>& A,
	const map<int, double>& b,
	const map<int, double>& x)
{
	// A*x
	map<int, double> ax = SparseMultiply(A, x);

	return SparseNormv(SparseDifference(ax, b));
}

/// <summary>
/// Найден ли элемент матрицы (row,column)
/// </summary>
/// <param name="matrix"></param>
/// <param name="row">строка</param>
/// <param name="column">столбец</param>
/// <returns>true - элемент найден, false - иначе</returns>
bool Find(const map<int, map<int, double>>& matrix, int row, int column)
{
	map<int, map<int, double>>::const_iterator iter = matrix.find(row);
	if (iter != matrix.end())
	{
		map<int, double> matrix_row = iter->second;
		return matrix_row.count(column) > 0;
	}
	else
		return false;
		
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
int SparseRotationSolve(const map<int, map<int, double>>& A,
	const map<int, double>& b,
	map<int, double>& x,
	int n)
{
	// верхняя тругольная матрица
	map<int, map<int, double>> T; // на текущей итерации

	map<int, map<int, double>> T0 = A; // на предыдущей итерации

	// преобразованный вектор правой части
	map<int, double> bet0 = b; // на предыдущей итерации

	map<int, double> bet; // на текущей итерации

	double val = 0;
	int i0 = 0, j0 = 0;
	int iter_num = 0; // число итерации, ограничивается максимальным числом итераций
	while (true)
	{
		// поиск максимального по модулю элемента матрицы T0
		i0 = j0 = 0;
		val = 0.0;
		for (const auto& elem_row : T0)
		{
			map<int, double> vect_column = elem_row.second;
			for (const auto& elem_column : vect_column)
			{
				double el = abs(elem_column.second);
				if (elem_row.first > elem_column.first && el > val)
				{
					i0 = elem_row.first;
					j0 = elem_column.first;
					val = el;
				}

			}
		}
		
		val = T0[i0][j0];

		if (abs(val) < DBL_EPSILON) break;
		// угол матрицы вращения, находится по условию T(k-1)(i0, j0) = 0
		double t2 = Find(T0, j0, j0) ? T0[j0][j0] : 0.0;
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
			t_i0 = Find(T0, i0, k) ? T0[i0][k] : 0.0;

			// T0(j0,k)
			t_j0 = Find(T0, j0, k) ? T0[j0][k] : 0.0;

			// T(i0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
			val = t_i0 * cs + t_j0 * ss;
			if(val != 0.0) T[i0][k] = val;

			// T(j0,k) = T0(i0,k)*cos(phi) + T0(j0,k)*sin(phi)
			val = -t_i0 * ss + t_j0 * cs;
			if (val != 0.0) T[j0][k] = val;
		}

		// bet(k) = U*bet(k-1)
		bet = bet0;
		// b0(i0)
		double b_i0 = bet0.count(i0) > 0 ? bet0[i0] : 0.0;

		// b0(j0)
		double b_j0 = bet0.count(j0) > 0 ? bet0[j0] : 0.0;

		val = b_i0 * cs + b_j0 * ss;
		if (val != 0.0)
			bet[i0] = val;
		val = -b_i0 * ss + b_j0 * cs;
		if (val != 0.0)
			bet[j0] = val;

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

	for (int i = n; i >= 1; i--)
	{
		val = bet.count(i) > 0 ? bet[i] : 0;

		double tval = 0.0, xval = 0.0;
		for (int j = i + 1; j <= n; j++)
		{
			tval = Find(T, i, j) ? T[i][j] : 0;
			if (tval == 0) continue;
			xval = x.count(j) > 0 ? x[j] : 0;
			if (xval == 0) continue;
			val -= tval * xval;
		}

		tval = Find(T, i, i);
		if (abs(tval) <= DBL_MIN) return -1;

		if (abs(val) <= DBL_MIN) continue;

		x[i] = val / tval;
	}

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
int SparseRelaxation(const map<int, map<int, double>>& A,
	const map<int, double>& b,
	map<int, double>& x,
	double omega)
{
	map<int, double> x0;
	x0[1] = 1.0;

	int iter = 0; // число итераций
	double err_norm = 0.0; // погрешность
	do
	{
		double sum = 0.0, val_d = 0.0;

		for (const auto& elem_row_a : A)
		{
			map<int, double> Arow = elem_row_a.second;
			double sum = 0.0;
			int row = elem_row_a.first;
			if (Arow.count(row) == 0) return -1; // нулевой диагональный элемент

			for (const auto& elem_x : x)
			{
				if (elem_x.first <= row - 1 && Arow.count(elem_x.first) > 0)
					sum += Arow[elem_x.first] * elem_x.second;
			}

			for (const auto& elem_x : x0)
			{
				if (elem_x.first >= row && Arow.count(elem_x.first) > 0)
					sum += Arow[elem_x.first] * elem_x.second;
			}			

			// b[row]
			auto iter_b = find_if(b.begin(), b.end(),
				[row](const auto& pair) {return pair.first == row; }
			);
			double val_b = iter_b != b.end() ? iter_b->second : 0.0;

			// x0[row]
			auto iter_x = find_if(x0.begin(), x0.end(),
				[row](const auto& pair) {return pair.first == row; }
			);
			double val_x = iter_x != x0.end() ? iter_x->second : 0.0;

			// A[row][row]
			auto iter_a = find_if(Arow.begin(), Arow.end(),
				[row](const auto& pair) {return pair.first == row; }
			);
			double val_a = iter_a != Arow.end() ? iter_a->second : 0.0;

			double val = val_x - omega * (sum - val_b) / val_a;
			if (abs(val) <= DBL_MIN) continue;
			x[row] = val;
		}
		err_norm = SparseNormv(SparseDifference(x, x0));
		x = x0;
		if (++iter > MAX_ITERATION_NUMBER) return 0;
	} while (err_norm > DBL_EPSILON);
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
int SparseGradientDescent(const map<int,map<int,double>>& A,
	const map<int, double>& b,
	map<int, double>& x)
{

	// градиент: 2*(Ax0 - b)
	auto Gradient = [A, b](const map<int, double>& x)
		{
			map<int, double> result = SparseDifference(SparseMultiply(A, x), b);
			for (auto& elem : result)
				elem.second *= 2.0;
			return result;
		};

	// A*(Ax0 - b)
	auto VectorRight = [A, b](const map<int, double>& x)
		{
			map<int, double> ax = SparseDifference(SparseMultiply(A, x), b);
			map<int, double> result = SparseMultiply(A, ax);
			return result;
		};

	map<int, double> x0;
	x0[1] = 1.0;
	map<int, double> grad;

	int iter = 0;
	double err_norm = 0.0, t = 0.0;
	do
	{
		grad = Gradient(x0);
		double sc1 = SparseScalarProduct(grad, grad) * 0.25;

		double sc2 = SparseScalarProduct(grad, VectorRight(x0));
		t = sc1 / sc2;

		for (auto& elem : grad)
			elem.second *= t;

		x = SparseDifference(x0, grad);

		// норма разности между значением на текущей и предыдущей итерации
		err_norm = SparseNormv(SparseDifference(x, x0));
		if (++iter > MAX_ITERATION_NUMBER) return 0;
		x0 = x;
	} while (err_norm > DBL_EPSILON);
	return 1;
}


