#include "SysSolve.h"
#include "SparseSolve.h"
//#include "SparseMapSolve.h"
#include "resource.h"
#include "LSysSolve.h"

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
std::string GetResourceString(int id_res)
{
    LPWSTR buff = new WCHAR[1024];
    std::string result = "";
    int count = ::LoadString(::GetModuleHandle(NULL), id_res, buff, 1024*sizeof(WCHAR));
    if (count > 0)
        result = LPWSTRToStdString(buff);
    delete[] buff;
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

    for (int i = 0; i < 8; i++)
    {
        //if ((n >= 30 && i == 5) || (n >=300 && i == 6)) continue;
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
            conv = Relaxation(Anorm, bet, x, n, 0.9);
            delete[] Anorm;
            delete[] bet;
            res_id = IDS_RELAXATION;
            if (!conv)
            {
                method_name = GetResourceString(res_id);
                std::string err_message = GetResourceString(IDS_EXCEED_MAX_NUMBER_ITERATION);
                std::cout << method_name << ". " << err_message << std::endl;
                continue;
            }
            break;
        case 6:
            conv = RotationSolve(A, v, x, n);
            res_id = IDS_ROTATION;
            if (!conv)
            {
                method_name = GetResourceString(res_id);
                std::string err_message = GetResourceString(IDS_EXCEED_MAX_NUMBER_ITERATION);
                std::cout << method_name << ". " << err_message << std::endl;
                continue;
            }

            break;
        case 7:
            Anorm = new double[n * n * sizeof(double)];
            bet = new double[n * sizeof(double)];
            TransformLinearSystem(A, v, Anorm, bet, n);
            conv = GradientDescent(Anorm, bet, x, n);
            delete[] Anorm;
            delete[] bet;
            res_id = IDS_GRADIENT_DESCENT;
            if (!conv)
            {
                method_name = GetResourceString(res_id);
                std::string err_message = GetResourceString(IDS_EXCEED_MAX_NUMBER_ITERATION);
                std::cout << method_name << ". " << err_message << std::endl;
                continue;
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
        method_name = GetResourceString(res_id);
        std::cout << method_name << '\t' << secs.count() << '\t' << norm(errv, n) << std::endl;
    }

    delete[] A;
    delete[] v;
    delete[] x;
    delete[] errv;
}

std::string GetErrorMessage(int ret_code)
{
    int res_id = 0;
    switch (ret_code)
    {
        case 0:
            res_id = IDS_EXCEED_MAX_NUMBER_ITERATION;
            break;
        case -2:
            res_id = IDS_NOT_MATCHED;
            break;
        case -1:
            res_id = IDS_ZERO;
            break;
        default:
            return "";
    }
    
    return GetResourceString(res_id);
}

/// <summary>
/// Заполнение матрицы СЛАУ и вектора правой части эдементами
/// </summary>
/// <param name="n">число уравнений</param>
/// <param name="A">матрица СЛАУ</param>
/// <param name="b">вектор правой части СЛАУ</param>
void FillData1(std::vector<SparseElement>& A, std::vector<SparseElement>& b, int n)
{
    int p = 0;
    std::cout << "Число внедиагональных элементов (1 - " << n - 1 << "): ";
    std::cin >> p;

    srand(10);
  
    for (int i = 1; i <= n; i++)
    {
        A.push_back({ i , i , rand_range(-1.0e3, 1.0e3) });
        int k = 1;
        for (int j = 1; j <= (p >> 1); j++)
        {
            if (i - k > 0 && k <= (p << 1))
                A.push_back({ i , i - k , rand_range(-1.0e3, 1.0e3) });

            if (k <= (p >> 1) && i + k <= n)
                A.push_back({ i , i + k , rand_range(-1.0e3, 1.0e3) });

            k++;
        }

    }

    for (int i = 1; i <= n; i++)
        b.push_back({ i , 1 , rand_range(-100.0, 200) });

}

void FillData2(std::vector<SparseElement>& A, std::vector<SparseElement>& b, int n)
{
    int p = 0;
    std::cout << "Число элементов выше главной диагонали (1 - " << n - 1 << "): ";
    std::cin >> p;

    srand(10);

    for (int i = 1; i <= n; i++)
    {
        A.push_back({ i , i , rand_range(1.0e1, 1.0e3) });
        int k = 1;
        while(k <=p)
        {
            if (i + k > n) break;
            double val = rand_range(-100.0, 100.0);
            A.push_back({ i , i + k , val / n });
            A.push_back({ i + k, i , val / n });
            k++;
        }

    }

    for (int i = 1; i <= n; i++)
        b.push_back({ i , 1 , rand_range(-100.0, 200) });

}

void FillData2(std::map<int,std::map<int,double>>& A, std::map<int, double>& b, int n)
{
    int p = 0;
    std::cout << "Число элементов выше главной диагонали (1 - " << n - 1 << "): ";
    std::cin >> p;

    srand(10);

    for (int i = 1; i <= n; i++)
    {
        A[i][i] = rand_range(1.0e1, 1.0e3);
        int k = 1;
        while (k <= p)
        {
            if (i + k > n) break;
            double val = rand_range(-100.0, 100.0);
            A[i][i + k] = val / n;
            A[i + k][i] = val / n;
            k++;
        }

    }

    for (int i = 1; i <= n; i++)
        b[i] = rand_range(-100.0, 200);

}


void SparseTestSolve()
{
    int n = 0;
    std::cout << " Число уравнений СЛАУ: ";
    std::cin >> n;

    // заполнение матрицы коэффициентов СЛАУ и вектора правой части с помощью генерации случайных чисел
    // для разреженных матриц и векторов нумерация элементов начинается с 1
    std::vector<SparseElement> A, v, x, At;
    FillData2(A, v, n);

    //for (const auto elem : A)
    //   std::cout << elem.row << '\t' << elem.column << '\t' << elem.value << std::endl;
    //PrintMatrix(A, n);
    std::cout << "Число уравнений: " << n << std::endl;
    double fullness_degree = (double)A.size() * 100.0 / (double)(n * n);
    std::cout << "Степень заполнения матрицы: " << fullness_degree << " %" << std::endl;
    
    int conv = 0, res_id = 0;
    std::string method_name;
    std::vector<SparseElement> Anorm, bet;

    for (int i = 0; i < 3; i++)
    {
        if (i != 1 && n >= 200) continue;
        auto start = std::chrono::steady_clock::now();
        switch (i)
        {
        case 0:
            
            conv = SparseRotationSolve(A, v, x, n);
            res_id = IDS_ROTATION;
            if (conv < 1)
            {
                method_name = GetResourceString(res_id);
                std::string err_message = GetErrorMessage(conv);
                std::cout << method_name << ". " << err_message << std::endl;
                continue;
            }

            break;
        case 1:
            //At = SparseTranspose(A);
            //Anorm = SparseMultiply(At, A);
            //bet = SparseMultiply(At, v);
            res_id = IDS_RELAXATION;
            conv = SparseRelaxation(A, v, x, 0.9);
            if (conv < 1)
            {
                method_name = GetResourceString(res_id);
                std::string err_message = GetErrorMessage(conv);
                std::cout << method_name << ". " << err_message << std::endl;
                continue;
            }
            break;
        case 2:
            //At = SparseTranspose(A);
            //Anorm = SparseMultiply(At, A);
            //bet = SparseMultiply(At, v);
            res_id = IDS_GRADIENT_DESCENT;
            conv = SparseGradientDescent(A, v, x);
            if (conv < 1)
            {
                method_name = GetResourceString(res_id);
                std::string err_message = GetErrorMessage(conv);
                std::cout << method_name << ". " << err_message << std::endl;
                continue;
            }
            break;

        }

        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        std::chrono::duration<double> secs = duration;

        method_name = GetResourceString(res_id);
        std::cout << method_name << '\t' << secs.count() << '\t' << ErrorMeasure(A, v, x) << std::endl;


    }
}

//void SparseTestSolveM()
//{
//    int n = 0;
//    std::cout << " Число уравнений СЛАУ: ";
//    std::cin >> n;
//
//    // заполнение матрицы коэффициентов СЛАУ и вектора правой части с помощью генерации случайных чисел
//    // для разреженных матриц и векторов нумерация элементов начинается с 1
//    std::map<int,std::map<int,double>> A, At, Anorm;
//    std::map<int, double> v, x, bet;
//    FillData2(A, v, n);
//
//    //for (const auto elem : A)
//    //   std::cout << elem.row << '\t' << elem.column << '\t' << elem.value << std::endl;
//    //PrintMatrix(A, n);
//    std::cout << "Число уравнений: " << n << std::endl;
//    double fullness_degree = (double)A.size() * 100.0 / (double)(n * n);
//    std::cout << "Степень заполнения матрицы: " << fullness_degree << " %" << std::endl;
//
//    int conv = 0, res_id = 0;
//    std::string method_name;
//
//    for (int i = 0; i < 3; i++)
//    {
//        auto start = std::chrono::steady_clock::now();
//        switch (i)
//        {
//        case 0:
//
//            conv = SparseRotationSolve(A, v, x, n);
//            res_id = IDS_ROTATION;
//            if (conv < 1)
//            {
//                method_name = GetResourceString(res_id);
//                std::string err_message = GetErrorMessage(conv);
//                std::cout << method_name << ". " << err_message << std::endl;
//                continue;
//            }
//
//            break;
//        case 1:
//            //At = SparseTranspose(A);
//            //Anorm = SparseMultiply(At, A);
//            //bet = SparseMultiply(At, v);
//            res_id = IDS_RELAXATION;
//            conv = SparseRelaxation(A, v, x, 1.9);
//            if (conv < 1)
//            {
//                method_name = GetResourceString(res_id);
//                std::string err_message = GetErrorMessage(conv);
//                std::cout << method_name << ". " << err_message << std::endl;
//                continue;
//            }
//            break;
//        case 2:
//            //At = SparseTranspose(A);
//            //Anorm = SparseMultiply(At, A);
//            //bet = SparseMultiply(At, v);
//            res_id = IDS_GRADIENT_DESCENT;
//            conv = SparseGradientDescent(A, v, x);
//            if (conv < 1)
//            {
//                method_name = GetResourceString(res_id);
//                std::string err_message = GetErrorMessage(conv);
//                std::cout << method_name << ". " << err_message << std::endl;
//                continue;
//            }
//            break;
//
//        }
//
//        auto end = std::chrono::steady_clock::now();
//        auto duration = end - start;
//        std::chrono::duration<double> secs = duration;
//
//        method_name = GetResourceString(res_id);
//        std::cout << method_name << '\t' << secs.count() << '\t' << ErrorMeasure(A, v, x) << std::endl;
//
//
//    }
//}

int main()
{
    setlocale(LC_ALL, ""); // для от ображения кириллицы
    //TestLinearSystemSolve2();
    SparseTestSolve();
}

