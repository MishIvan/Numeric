#include "SysSolve.h"

const char* sys_methods[7]
{
    "Метод Гаусса","LU декомпозиция",
    "Компактная схема\nисключения",
    "QR декомпозиция", 
    "LLT\nдекомпозиция",
    "Верхняя\nрелаксация",
    "Метод\nвращений"
};

// Тестирование решения СЛАУ различного порядка
void TestLinearSystemSolve2()
{
    int  n = 6;
    srand(10);
    // заполнение матрицы коэффициентов СЛАУ и вектора правой части с помощью генерации случайных чисел
    double* A = new double[n*n*sizeof(double)];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
             *(A + i * n + j) = rand_range(-1.0e3, 1.0e3);

    double *v = new double[n * sizeof(double)]; 
    double *x = new double[n * sizeof(double)];
    for (int i = 0; i < n; i++)
        *(v +i) = rand_range(-100.0, 200);

    double* Anorm = nullptr;
    double* bet = nullptr;

    double* errv = new double[n * sizeof(double)];
    std::cout << "Порядок матрицы n = " << n << std::endl;
    std::cout << "Метод" << "\t\t" << "Время, сек" << '\t' << "Норма невязки" << std::endl;

    for (int i = 0; i < 7; i++)
    {
        auto start = std::chrono::steady_clock::now();

        switch (i)
        {
        case 0:
            Gauss(A, v, x, n); break;
        case 1:
            LUDecompositionSolve(A, v, x, n); break;
        case 2:
            CompactSchemeSolve(A, v, x, n); break;
        case 3:
            QRDecompositionSolve(A, v, x, n); break;
        case 4:
            Anorm = new double[n * n * sizeof(double)];
            bet = new double[n * sizeof(double)];
            TransformLinearSystem(A, v, Anorm, bet, n);
            LLTDecompositionSolve(Anorm, bet, x, n); 
            delete[] Anorm;
            delete[] bet;
            break;
        case 5:
            Anorm = new double[n * n * sizeof(double)];
            bet = new double[n * sizeof(double)];
            TransformLinearSystem(A, v, Anorm, bet, n);
            UpperRelaxation(Anorm, bet, x, n);
            delete[] Anorm;
            delete[] bet;
            break;
        case 6:
            bool is_symmetric = IsSymmetric(A, n);
            if (!is_symmetric)
            {
                Anorm = new double[n * n * sizeof(double)];
                bet = new double[n * sizeof(double)];
                TransformLinearSystem(A, v, Anorm, bet, n);
            }
            else
            {
                Anorm = A; bet = v;
            }
            RotationSolve(Anorm, bet, x, n);
            if (!is_symmetric)
            {
                delete[] Anorm;
                delete[] bet;
            }
            break;

        }

        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        std::chrono::duration<double> secs = duration;

        // вычисление нормы вектора Ax - b
        for (int i = 0; i < n; i++)
        {
            *(errv + i) = -1.0 * *(v + i);
            for (int j = 0; j < n; j++)
                *(errv + i) += *(A + i * n + j) * *(x + j);
        }

        std::cout << sys_methods[i] << '\t' << secs.count() << '\t' << norm(errv, n) << std::endl;
    }

    delete[] A;
    delete[] v;
    delete[] x;
    delete[] errv;
}

int main()
{
    setlocale(LC_ALL, ""); // для от ображения кириллицы
    TestLinearSystemSolve2();
}

