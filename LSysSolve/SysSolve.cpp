#include "SysSolve.h"
/// <summary>
/// Генарация вещественных случайных чисел
/// </summary>
/// <param name="min">минимальное значение и нтервала</param>
/// <param name="max">максимальное значение интервала</param>
/// <returns>случайное число в указанном интервале</returns>
double rand_range(double min, double max) 
{
	return min + (max - min) * std::rand() / (RAND_MAX + 1);
}

/// <summary>
/// Вычисление сферической нормы n-мерного вектора v
/// </summary>
/// <param name="v">вектор</param>
/// <param name="n">размерность вектора</param>
/// <returns>значение нормы</returns>
double norm(const double* v, int n)
{
	double nrm = 0.0;
	for (int i = 0; i < n; i++)
		nrm += *(v + i) * *(v + i);
	return sqrt(nrm);
}

/// <summary>
/// Является ли квадратная матрица симметричной
/// </summary>
/// <param name="A">матрица</param>
/// <param name="n">порядок матрицы</param>
/// <returns>true - матрица симметричная, false - матрица не симметричная</returns>
bool IsSymmetric(const double* A, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i != j && *(A + i * n + j) != *(A + j * n + i)) return false;
		}
	return true;
}

/// <summary>
/// Переобразование СЛАУ A*x = b к виду At*A*x = At*b
/// </summary>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="Anorm">матрица At*A</param>
/// <param name="bet">вектор At*b</param>
/// <param name="n">порядок матрицы A</param>
void TransformLinearSystem(const double* A, const double* b, double* Anorm, double* bet, int n)
{
	// Anorm = At * A; 
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			*(Anorm + i * n + j) = 0.0;
			for (int k = 0; k < n; k++)
				*(Anorm + i * n + j) += *(A + k * n + i) * *(A + k * n + j);
		}
	// bet = At * b;
	for (int i = 0; i < n; i++)
	{
		*(bet + i) = 0.0;
		for (int k = 0; k < n; k++)
			*(bet + i) += *(A + k * n + i) * *(b + k);
	}

}
/// <summary>
/// Вывод квадратной матрицы A порядка n на консоль
/// </summary>
/// <param name="A"></param>
/// <param name="n"></param>
void Matrix2Console(const double* A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) // вывод строки в поток
			std::cout << *(A + i * n + j) << ' ';
		std::cout << std::endl;
	}

}

/// <summary>
/// Решение системы линейных алгебраических уравнений (СЛАУ) методом Гаусса 
/// с вычислением определителя матрицы системы 
/// </summary>
/// <param name="a">матрица коэффициентов СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">решение СЛАУ</param>
/// <param name="size">порядок матрицы СЛАУ</param>
/// <returns>определитель матрицы a</returns>
double Gauss(const double *a, const double* b, double* x, int size)
{
	int i, k, m;
	long double amm, aim;

	// сведение исходной системы к системе с верхней треугольной матрицей
	double* alf = new double[size * size*sizeof(double)];
	double* bet = new double[size * sizeof(double)];
	memcpy(alf, a, size*size*sizeof(double));
	memcpy(bet, b, size * sizeof(double));
	for (m = 0; m <= size - 2; m++)
	{
		amm = *(alf + m*size + m);
		if (std::abs(amm) < 1.0e-36)
		{
			delete[] alf;
			delete[] bet;
			return 0.0;
		}
		for (k = m; k <= size - 1; k++)
			*(alf + m*size + k) /= amm;
		*(bet +m) /= amm;
		for (i = m + 1; i <= size - 1; i++)
		{
			aim = *(alf + i*size  + m);
			for (k = m; k <= size - 1; k++)
				*(alf + i*size + k) -= *(alf +m*size + k) * aim;
			*(bet +i) -= *(bet +m) * aim;
		}//end i 
	}//end m 

	// нахождение решения СЛАУ с верхней треугольной матрицей
	// на главной диагонали которой находятся единицы
	*(x + size - 1) = *(bet + size - 1) / *(alf + (size - 1)*size + size - 1);
	for (int i = size - 2; i >= 0; i--)
	{
		*(x + i) = *(bet +i);
		for (k = i + 1; k < size; k++)
			*(x + i) -= *(alf + i*size + k) * *(x + k);
	}//end i

	// вычисление определителя
	double det = 1.0;
	for (i = 0; i < size; i++)
		det *= *(alf + i * size + i);

	delete[] alf;
	delete[] bet;
	return det;
}

/// <summary>
/// Формирование матриц alpha и gamma в компактной схеме исключения
/// Матрица alpha верхняя треугольная, gamma - нижняя треугольная
/// </summary>
/// <param name="alpha">формируемые матрицы</param>
/// <param name="n">порядок матриц</param>
/// <returns>значение определителя матрицы</returns>
double FormMatrixCompactScheme(const double *A, double *alpha, int n)
{
	for (int i = 0; i < n; i++)
	{
		*(alpha + i * n) = *(A + i * n);
		if (i > 0) *(alpha + i) = *(A + i) / *A;
	}

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
			*(alpha + i * n + k) = *(A + i * n + k);
			for (int j = 0; j <= k - 1; j++)
				*(alpha + i * n + k) -= *(alpha + i * n + j) * *(alpha + j * n + k);
		}
		else
		{
			if (std::abs(*(alpha + i * n + i)) < 1.0e-36) return 0.0;
			*(alpha + i * n + k) = *(A + i * n + k);
			for (int j = 0; j <= i - 1; j++)
			{
				*(alpha + i * n + k) -= *(alpha + i * n + j) * *(alpha + j * n + k);
			}
			*(alpha + i * n + k) /=  *(alpha + i * n + i);
		}
		k++;
	}

	// вычисление определителя
	double det = 1.0;
	for (int i = 0; i < n; i++)
		det *= *(alpha + i * n + i);
	return det;
}

/// <summary>
/// Решение системы линейных уравнений компактной схемой исключения
/// </summary>
/// <param name="A">матрица системы уравнений</param>
/// <param name="b">вектор правой части системы уравнений</param>
/// <param name="x">решение системы уравнений</param>
/// <param name="n">порядок матрицы системы</param>
void CompactSchemeSolve(const double* A, const double* b, double* x, int n)
{
	double* alpha = new double[n*n*sizeof(double)];
	double det = FormMatrixCompactScheme(A, alpha, n);
	// решение только для неособенной матрицы
	if (abs(det) < 1.0e-36)
	{
		delete[] alpha;
		return;
	}
	double *beta = new double[n*sizeof(double)];
	*beta = *b / *A;
	for (int i = 1; i < n; i++)
	{
		*(beta + i) = *(b + i);
		for (int j = 0; j <= i - 1; j++)
			*(beta + i) -= *(alpha + i*n + j) * *(beta +j);
		*(beta +i) /= *(alpha + i*n + i);
	}

	// решение системы уравнений с труегольной матрицей		   
	*(x + n - 1) = *(beta + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		*(x + i) = *(beta + i);
		for (int j = n - 1; j > i; j--)
			*(x + i) -= *(alpha + i*n + j) * *(x + j);
	}

	delete[] alpha;
	delete[] beta;
}

/// <summary>
/// QR разложение квадратной матрицы A, A = QR
/// где Q - ортогональная, а R - верхняя треугольная матрица
/// </summary>
/// <param name="Q">ортогональная матрица Q</param>
/// <param name="R">верхняя треугольная матрица R</param>
/// <param name="A">исходная матрица</param>
/// <param name="n">порядок матрицы A</param>
/// <returns>false - если на главной диагонали матрицы R 
/// хотя бы один элемент равен нулю (A - вырожденная матрица), иначе - true</returns>
bool QRDecomposition(const double* A, double* Q, double* R, int n)
{
	double sum = 0.0;

	for (int j = 0; j < n; j++)
	{
		// q(j) = a(j)
		for (int k = 0; k < n; k++)
			*(Q + k * n + j) = *(A + k * n + j);

		for (int i = 0; i <= j - 1; i++)
		{
			// rij = q(i)^T*a(j)
			sum = 0.0;
			for (int k = 0; k < n; k++)
				sum += *(Q + k * n + i) * *(A + k * n + j);
			*(R + i * n + j) = sum;

			// r(i,j)*q(i)
			for (int k = 0; k < n; k++)
				*(Q + k * n + j) -= *(R + i * n + j) * *(Q + k * n + i);
		}

		// r(j,j) = || q(j)) ||2
		sum = 0.0;
		for (int k = 0; k < n; k++)
			sum += *(Q + k * n + j) * *(Q + k * n + j);
		*(R + j * n + j) = sqrt(sum);

		if (abs(*(R + j * n + j)) < 1.0e-36) return false;

		for (int k = 0; k < n; k++)
			*(Q + k * n + j) /= *(R + j * n + j);

	}

	return true;
}

/// <summary>
/// Решение СЛАУ с применением QR декомпозиции матрицы системы A
/// </summary>
/// <param name="A">матрица системы линейных уравнений</param>
/// <param name="b">вектор правой части системы линейных уравнений</param>
/// <param name="x">вектор решения системы линейных уравнений</param>
/// <param name="n">порядок матрицы A</param>
void QRDecompositionSolve(const double* A, const double* b, double* x, int n)
{
	double* Q = new double[n * n * sizeof(double)];
	double* R = new double[n * n * sizeof(double)]; 
	if (!QRDecomposition(A, Q, R, n))
	{
		delete[] Q;
		delete[] R;
		return;
	}
	double * beta = new double[n * sizeof(double)];
	//  beta = Qt*b
	for (int i = 0; i < n; i++)
	{
		*(beta + i) = 0.0;
		for (int k = 0; k < n; k++)
			*(beta + i) += *(Q + k * n + i) * *(b + k);

	}
	delete[] Q;
	// решение системы уравнений с верхней труегольной матрицей
	// R*x = Qt*beta		   
	*(x + n - 1) = *(beta + n - 1) / *(R + (n - 1)*n + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		*(x + i) = *(beta + i);
		for (int j = n - 1; j > i; j--)
			*(x + i) -= *(R + i * n + j) * *(x + j);
		*(x + i) /= *(R + i * n + i);
	}
	delete [] R;
	delete[] beta;
}

/// <summary>
/// Декомпозиция матрицы A = LU
/// L - нижняя треугольная матрица, на главной диагонали которой расположены единицы
/// U - верхняя треугольная матрица
/// </summary>
/// <param name="alfa">Матрицы L и U, ниже гланой диагонали которой расположены внедиагональные элементы L, 
/// на главной диагонали и выше расположены элементы матрицы U</param>
/// <param name="A">исходная матрица</param>
/// <param name="n">порядок матрицы</param>
/// <returns>Определитель матрицы</returns>
double LUDecomposition(const double* A, double* alfa, int n)
{

	// декомпозиция матрицы
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i <= j)
			{
				*(alfa + i * n + j) = *(A + i * n + j);
				for (int k = 0; k <= i - 1; k++)
					*(alfa + i * n + j) -= *(alfa + i * n + k) * *(alfa + k * n + j);
			}
			else
			{
				if (abs(*(alfa + j * n + j)) < 1.0e-36) return 0;
				*(alfa + i * n + j) = *(A + i * n + j);
				for (int k = 0; k <= j - 1; k++)
					*(alfa + i * n + j) -= *(alfa + i * n + k) * *(alfa + k * n + j);
				*(alfa + i * n + j) /= *(alfa + j * n + j);

			}
		}
	double det = 1;
	for (int i = 0; i < n; i++)
		det *= *(alfa + i * n + i);
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
void LUDecompositionSolve(const double* A, const double* b, double* x, int n)
{

	double *alfa = new double[n*n*sizeof(double)];
	
	double det = LUDecomposition(A, alfa, n);
	if (abs(det) < 1.0e-36) return;
	// решение СЛАУ
	double* y = new double[n * sizeof(double)];
	for (int k = 0; k < n; k++)
	{
		*(y + k) = *(b + k);
		for (int j = 0; j <= k - 1; j++)
			y[k] -= *(alfa + k*n +j) * *(y + j);
	}

	for (int k = n - 1; k >= 0; k--)
	{
		*(x + k) = *(y + k);
		for (int j = n - 1; j > k; j--)
			*(x + k) -= *(alfa + k*n + j) * *(x + j);
		*(x + k) /= *(alfa + k*n + k);
	}

	delete[] alfa;
}

/// <summary>
///  Разложение Холецкого симметричной положительно определённой матрицы A = L*L^t
/// </summary>
/// <param name="L">матрица в разложении Холецкого</param>
/// <param name="A">исходная матрица</param>
/// <param name="n">порядок матрицы A</param>
/// <returns></returns>
bool CholeskyDecomposition(const double*A, double* L, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++) {

			double sum = 0;
			for (int k = 0; k < j; k++)
				sum += *(L + i * n + k) * *(L + j * n + k);

			if (i == j)
			{
				*(L + i * n + j) = sqrt(*(A + i * n + i) - sum);
				if ( isnan(*(L + i * n + j)) ) return false;
			}
			else
				*(L + i * n + j) = (*(A + i * n + j) - sum) / (*(L + j * n + j));
		}

	return true;
}

/// <summary>
/// Решение системы уравнений с применением метода LLT декомпозиции A = LLT
/// Метод применим только для положительно определённых симметричных матриц
/// </summary>
/// <param name="A">матрица</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">вектор решения СЛАУ</param>
void LLTDecompositionSolve(const double* A, const double* b, double* x, int n)
{	
	
	double* L = new double[n * n * sizeof(double)];
	memset(L, 0, n * n * sizeof(double));
	CholeskyDecomposition(A, L, n);

	// решение системы с нижней треугольной матрицей Ly = b
	double* y = new double[n * sizeof(double)];
	for (int i = 0; i < n; i++)
	{
		*(y + i) = *(b + i);
		for (int k = 0; k < i; k++)
			*(y + i) -= *(L +i*n + k) * *(y + k);
		*(y + i) /= *(L + i * n + i);

	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i < j) *(L + i * n + j) = *(L + j * n + i);

	// решение системы уравнений с верхней труегольной матрицей L^t*x = y   
	*(x + n - 1) = *(y + n - 1) / *(L + (n - 1)*n + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		*(x + i) = *(y + i);
		for (int j = n - 1; j > i; j--)
			*(x + i) -= *(L + i * n + j) * *(x + j);
		*(x + i) /= *(L + i * n + i);
	}

	delete[] y;
	delete[] L;
}

/// <summary>
/// Решение системы уравнений A*x = b методом вехней релаксации
/// для симметричной положительно определённой матрицы A
/// </summary>
/// <param name="A">матрица</param>
/// <param name="b">вектор правой части</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <returns>true - если сходимость была достигнута, false - иначе</returns>
bool UpperRelaxation(const double* A, const double* b, double* x, int n)
{
	double* x0 = new double[n * sizeof(double)];
	double *errv = new double[n * sizeof(double)];
	memcpy(x0, b, n * sizeof(double));

	int iter = 0; // число итераций
	double err_norm = 0.0; // погрешность
	double omega = 0.99; // коэффициент релаксации
	do
	{		
		for (int i = 0; i < n; i++)
		{
			double sum = 0.0;
			for (int k = 0; k <= i - 1; k++)
				sum += *(A + i * n + k) * *(x + k);

			for (int k = i; k < n; k++)
				sum += *(A + i * n + k) * *(x0 + k);

			*(x + i) = *(x0 + i) - omega * (sum - *(b + i)) / *(A + i * n + i);						
		}
		// норма разности между значением на текущей и предыдущей итерации
		for(int k =0; k < n; k++)
			*(errv + k) = *(x + k) - *(x0 + k);
		err_norm = norm(errv, n);
		if (++iter > MAX_ITERATION_NUMBER)
		{
			delete[] errv;
			delete[] x0;
			return false;
		}
		memcpy(x0, x, n * sizeof(double));
	} while (err_norm > 1.0e-18);
	
	delete[] errv;
	delete[] x0;
	return true;
}

/// <summary>
/// Решение СЛАУ методом вращения, решение сводится
/// к решению СЛАУ с верхней треугольной матрицей
/// </summary>
/// <param name="Asrc">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
/// <param name="x">вектор решения СЛАУ</param>
/// <param name="n">порядок матрицы</param>
void RotationSolve(const double *Asrc, const double *b, double *x, int n)
{

	// верхняя тругольная матрица
	double* A = new double[n * n * sizeof(double)]; // на текущей итерации
	memcpy(A, Asrc, n * n * sizeof(double));

	double* A0 = new double[n * n * sizeof(double)]; // на предыдущей итерации
	memcpy(A0, Asrc, n * n * sizeof(double));

	// преобразованный вектор правой части
	double* bet0 = new double[n * sizeof(double)]; // на предыдущей итерации
	memcpy(bet0, b, n * sizeof(double));

	double* bet = new double[n * sizeof(double)]; // на текущей итерации
	memcpy(bet, b, n * sizeof(double));

	int iter = 1;
	while (true)
	{
		// поиск позиции максимального по модулю элемента ниже главной диагонали
		int i0 = 0, j0 = 0;
		double val = 0.0, el = 0.0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
			{
				if (i < j)
				{
					el = abs(*(A + j * n + i));
					if (el > val)
					{
						val = el;
						i0 = i;
						j0 = j;
					}
				}
			}
		if (val < 1.0e-17) break;
		
		// угол матрицы вращения, находится по условию A(j0, i0) = 0
		double fi = atan( *(A + j0 * n + i0) /  *(A + i0 * n + i0) );
		double cs = cos(fi);
		double ss = sin(fi);

		// A(n) = U*A(n-1), U - матрица вращения
		// bet(n) = U*bet(n-1)
		memcpy(A, A0, n * n * sizeof(double));
		memcpy(bet, bet0, n * sizeof(double));
		for (int i = 0; i < n; i++)
		{
			*(A + i0 * n + i) = *(A0 + i0 * n + i) * cs +
				*(A + j0 * n + i) * ss;
			*(A + j0 * n + i) = (-1.0) * *(A0 + i0 * n + i) * ss +
				*(A0 + j0 * n + i) * cs;

			*(bet + i0) = *(bet0 + i0) * cs +
				*(bet0 + j0) * ss;
			*(bet + j0) = (-1.0) * *(bet0 + i0) * ss +
				*(bet0 + j0) * cs;
		}		
		iter++;
		// результаты для следующей итерации
		memcpy(A0, A, n * n * sizeof(double));
		memcpy(bet0, bet, n * sizeof(double));
	}

	delete[] A0;
	delete[] bet0;

	// вычисление вектора решения
	// решение системы уравнений с верхней труегольной матрицей A(n)*x = bet   
	*(x + n - 1) = *(bet + n - 1) / *(A + (n - 1) * n + n - 1);
	for (int i = n - 2; i >= 0; i--)
	{
		*(x + i) = *(bet + i);
		for (int j = n - 1; j > i; j--)
			*(x + i) -= *(A + i * n + j) * *(x + j);
		*(x + i) /= *(A + i * n + i);
	}

	delete[] A;
	delete[] bet;
}