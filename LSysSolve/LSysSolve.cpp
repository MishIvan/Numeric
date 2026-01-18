#include "SysSolve.h"
#include "SparseSolve.h"
#include "resource.h"

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

    for (int i = 0; i < 7; i++)
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

void SparseTestSolve()
{
    int n = 0, p = 0;
    std::cout << "Введите число уравнений СЛАУ: ";
    std::cin >> n;
    srand(10);
    // заполнение матрицы коэффициентов СЛАУ и вектора правой части с помощью генерации случайных чисел
    // для разреженных матриц и векторов нумерация элементов начинается с 1
    std::vector<SparseElement> A, v, x;
    p = n < 20 ? n : n / 10;
    for (int i = 1; i <= n; i++)
    {
        SparseElement sp_el{ i , i , rand_range(-1.0e3, 1.0e3) };
        A.push_back(sp_el);        
        for (int j = 1; j <= p; j++)
        {
            if (i != j)
            {
                sp_el.row = i;
                sp_el.column = n >= 20 ? irand_range(1, n) : j;
                sp_el.value = rand_range(-1.0e3, 1.0e3);
                A.push_back(sp_el);
            }
        }
    }

    PrintMatrix(A);
    std::cout << "Степень заполнения матрицы: " << FullnessDegree(A) << " %" << std::endl;
    
    int k = 0;
    p = n < 20 ? n : n / 3;
    for (int i = 0; i < p; i++)
    {
        SparseElement sp_el{ ++k , 1 , rand_range(-100.0, 200) };
        v.push_back(sp_el);
    }

    int conv = 0, res_id = 0;
    std::string method_name;
    std::vector<SparseElement> Anorm, bet;

    for (int i = 0; i < 1; i++)
    {
        auto start = std::chrono::steady_clock::now();
        switch (i)
        {
        case 0:
            
            conv = SparseRotationSolve(A, v, x);
            res_id = IDS_ROTATION;
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
int main()
{
    setlocale(LC_ALL, ""); // для от ображения кириллицы
    //TestLinearSystemSolve2();
    SparseTestSolve();
}

