#include "SysSolve.h"
#include "resource.h"
#define _CRT_SECURE_NO_WARNINGS

const char* sys_methods[7]
{
    "Метод Гаусса",
    "LU декомпозиция",
    "Компактная схема\nисключения",
    "QR декомпозиция", 
    "LLT\nдекомпозиция",
    "Верхняя\nрелаксация",
    "Метод\nвращений"
};

std::string LPWSTRToStdString(LPWSTR lpws) {
    if (!lpws) return "";

    size_t len = wcslen(lpws);
    char* buff = new char[len * 4];

    size_t converted = 0;
    std::string result = "";
    wcstombs_s(&converted, buff, len * 4, lpws, len);
    if (converted < len - 1) {
        delete[] buff;  
        return result;
    } // ошибка

    result = buff;
    delete[] buff;
    return result;
}
std::string GetMethodName(int id_res)
{
    LPWSTR buff = new WCHAR[1024];
    std::string result;
    int count = ::LoadString(::GetModuleHandle(NULL), id_res, buff, 1024*sizeof(WCHAR));
    if (count > 0)
        result = LPWSTRToStdString(buff);
    return result;
}
// Тестирование решения СЛАУ различного порядка
void TestLinearSystemSolve2()
{
    int n = 0;
    std::cout << "Введите число уравнений СЛАУ: ";
    std::cin >> n;
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
    bool conv = false;

    double* errv = new double[n * sizeof(double)];
    std::cout << "Порядок матрицы n = " << n << std::endl;
    std::cout << "Метод" << "\t\t" << "Время, сек" << '\t' << "Норма невязки" << std::endl;
    int res_id = 0;
    std::string method_name;

    for (int i = 0; i < 7; i++)
    {
        if ((n >= 30 && i == 5) || (n >=300 && i == 6)) continue;
        auto start = std::chrono::steady_clock::now();

        switch (i)
        {
        case 0:
            Gauss(A, v, x, n); 
            res_id = IDS_GAUSS;
            break;
        case 1:
            LUDecompositionSolve(A, v, x, n);
            res_id = IDS_LU_DECOMPOSITION;
            break;
        case 2:
            CompactSchemeSolve(A, v, x, n); 
            res_id = IDS_COMPACT_SCHEME;
            break;
        case 3:
            QRDecompositionSolve(A, v, x, n); 
            res_id = IDS_QR_DECOMPOSITION;
            break;
        case 4:
            Anorm = new double[n * n * sizeof(double)];
            bet = new double[n * sizeof(double)];
            TransformLinearSystem(A, v, Anorm, bet, n);
            LLTDecompositionSolve(Anorm, bet, x, n); 
            delete[] Anorm;
            delete[] bet;
            res_id = IDS_LLT_DECOMPOSITION;
            break;
        case 5:
            Anorm = new double[n * n * sizeof(double)];
            bet = new double[n * sizeof(double)];
            TransformLinearSystem(A, v, Anorm, bet, n);
            conv = Relaxation(Anorm, bet, x, n, 0.8);
            delete[] Anorm;
            delete[] bet;
            res_id = IDS_RELAXATION;
            if (!conv)
            {
                method_name = GetMethodName(res_id);
                std::cout << method_name << ": превышено максимальное число итерации\nдля метода верхней релаксации" << std::endl;
                std::cout << "Метод не сходится" << std::endl;
                continue;
            }
            break;
        case 6:
            RotationSolve(A, v, x, n);
            res_id = IDS_ROTATION;
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
        method_name = GetMethodName(res_id);
        std::cout << method_name << '\t' << secs.count() << '\t' << norm(errv, n) << std::endl;
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

