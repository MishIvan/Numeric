// MatrixVector.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <ctime>
#include <chrono>
#include "MATRIX.h"
#include "MatrixVector.h"
const char* sys_methods[4]{ "Метод Гаусса","LU декомпозиция","Компактная схема исключения","QR декомпозиция" };

/// <summary>
/// Получить полный путь файла в папке, из которой запускается исполняемый файл программы
/// </summary>
/// <param name="fullExePath">полный путь запуска программы</param>
/// <param name="dataFileName">имя файла данных в каталоге, где расположен исполняемый файл</param>
/// <param name="fullFileName">имя файла данных в каталоге, где расположен исполняемый файл</param>
/// <returns></returns>
void GetFullPathInWD(char* fullExePath, const char* dataFileName, char* fullFileName)
{
    strcpy_s(fullFileName, 1024, fullExePath);
    std::string s1 = fullFileName;
    int k = s1.find_last_of('\\');
    if (k == std::string::npos)
        strcpy_s(fullFileName, 1024, dataFileName);
    else
    {
        s1.replace(k + 1, s1.size() - 1, dataFileName);
        strcpy_s(fullFileName, 1024, s1.data());
    }
}
/// <summary>
/// Тест перегрузки операции умножения матриц
/// </summary>
void TestMatrixMultiplication(char *appPath, char * path)
{
    GetFullPathInWD(appPath, "MatrixMult_in.txt", path);
    // считывание данных
    ifstream fs;
    fs.open(path);
    if (fs.is_open())
    {
        int m, n;
        fs >> m >> n;
        MATRIX A(m, n);
        fs >> A;
        std::cout << "Матрица A" << endl;
        std::cout << A;

        fs >> m >> n;
        MATRIX B(m, n);
        fs >> B;
        std::cout << endl << "Матрица B" << endl;
        std::cout << B;

        MATRIX C(A.rows(), B.columns());
        C = A * B;
        std::cout << endl << "Матрица C = A*B" << endl;
        std::cout << C;

        fs.close();
    }
}

/// <summary>
/// Тест перегрузки операции умножения матрицы на вектор 
/// </summary>
void TestMatrixVectorMultiplication(char* appPath, char* path)
{
    GetFullPathInWD(appPath, "MatrixVectMult_in.txt", path);
    // считывание данных
    ifstream fs;
    fs.open(path);
    if (fs.is_open())
    {
        int m, n;
        fs >> m >> n;
        MATRIX A(m, n);
        fs >> A;
        std::cout << "Матрица A" << endl;
        std::cout << A;

        fs >> n;
        VECTOR v(n);
        fs >> v;
        std::cout << endl << "Вектор v" << endl;
        std::cout << v;

        VECTOR vp(A.rows());
        vp = A* v;
        std::cout << endl << "Вектор  vp= A*v" << endl;
        std::cout << vp;

        fs.close();
    }
}

/// <summary>
/// Тест функционала матрицы и вектора
/// </summary>
/// <param name="appPath">полный путь к программе</param>
/// <param name="path">указатель на полные путь файла</param>
void TestMatrixNVector(char*appPath, char*  path)
{
    /// Вектор. Считывание данных
    GetFullPathInWD(appPath, "Vector_in.txt", path);
    VECTOR v(1);
    bool success = VECTOR::readFromFile(path, v);
    if (success)
    {
        std::cout << "Вектор размерностью " << v.size() << " считан из файла" << endl;
        std::cout << v << endl;
    }
    else
        std::cout << "Неудачная попытка чтения вектора из файла" << endl;

    /// Вектор. Считывание данных с консоли и запись в файл
    VECTOR v_in(5);
    std::cout << "Введите вектор размерностью 5 с консоли" << endl;
    cin >> v_in;
    std::cout << "Вектор разменостью " << v_in.size() << endl;
    std::cout << v_in << endl;
    GetFullPathInWD(appPath, "Vector_out.txt", path);
    success = VECTOR::writeToFile(path, v_in);
    if (success)
        std::cout << "Вектор, введённый с консоли, записан в файл" << endl;
    else
        std::cout << "Неудачная попытка записи вектора, введённого с консоли" << endl;


    /// Матрица. Считывание из файла
    GetFullPathInWD(appPath, "Matrix_in.txt", path);
    MATRIX matr(3, 3);
    success = MATRIX::readFromFile(path, matr);
    if (success)
    {
        std::cout << "Считаны данные матрицы " << matr.rows() << " x " << matr.columns() << endl;
        std::cout << matr << endl;
    }
    else
        std::cout << "Неудачная попытка считывания данных матрицы" << endl;

    ///  Матрица. Ввод с консоли и запись в файл
    MATRIX matr_in(3, 5);
    std::cout << "Введите матрицу 3 х 5 с консоли" << endl;
    std::cin >> matr_in;
    //cout << "Введена матрица" << endl;
    //cout << matr_in << endl;
    GetFullPathInWD(appPath, "Matrix_out.txt", path);
    success = MATRIX::writeToFile(path, matr_in);
    if (success)
        std::cout << "Матрица, введённая с консоли, записана в файл" << endl;
    else
        std::cout << "Неудачная попытка записи матрицы, введённой с консоли" << endl;

    /// Матрица. Умножение матриц
    TestMatrixMultiplication(appPath, path);

    /// Матрица. Умножение вестора на матрицу
    TestMatrixVectorMultiplication(appPath, path);
}
/// <summary>
/// QR разложение
/// </summary>
/// <param name="A"></param>
/// <param name="res"></param>
void TestQRDecomposition(MATRIX& A)
{
    // 
    MATRIX Q(A.rows(), A.columns()), R(A.rows(), A.columns());
    if (!A.QRDecomposition(Q, R))
    {
        std::cout << "Матрица R вырождена" << endl;
        return;
    }


    std::cout << endl << "Matrix Q" << endl;
    std::cout << Q << endl;

    std::cout << endl << "Matrix R" << endl;
    std::cout << R << endl;

    MATRIX rt(R.columns(), R.rows());
    rt = Q * R;

    std::cout << endl << "Matrix A = Q*R" << endl;
    std::cout << rt << endl;

    std::cout << endl << "Matrix E = Q^t*Q" << endl;
    rt = Q.Transpose() * Q;
    std::cout << rt << endl;
    
}
/// <summary>
/// Тест разложения Холецкого.  Разложение Холецкого работает для симметричнйо и положительно определённой матрицы
/// </summary>
/// <param name="A"></param>
void TestCholeskyDecomposuition(MATRIX& A)
{
    MATRIX L(A.rows(), A.columns()), Anorm(A.rows(), A.columns());
    Anorm = A.Transpose() * A;

    std::cout << endl << "Matrix A^t*A" << endl;
    std::cout << Anorm << endl;

    if (!Anorm.CholeskyDecomposition(L))
    {
        std::cout << "Разложение Холецкого неприменимо для матрицы" << endl;
        return;
    }

    std::cout << endl << "Matrix L" << endl;
    std::cout << L << endl;

    MATRIX rt(A.columns(), A.rows());
    rt = L * L.Transpose();

    std::cout << endl << "Matrix A^t*A = L*L^t" << endl;
    std::cout << rt << endl;

}

void TestMatrixReverse()
{
    int size = 6;
    MATRIX A(size, size), A1(size, size);
    srand(10);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            A(i, j) = rand_range(-1000,1000);
 
    if (size < 20)
    {
        cout << endl << "Matrix A" << endl;
        cout << A << endl;
    }
    for (int i = 0; i < 2; i++)
    {
        auto start = std::chrono::steady_clock::now();
        A1 = i== 0 ? A.Invert() : A.InvertLU();    // обращение матрицы
        cout << (i == 0 ? "Решение AX = E" : "LU декомпозиция") << endl;
        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        std::chrono::duration<double> secs = duration;
        cout << "Время решения " << secs.count() << " сек" << endl;

        if (size < 20)
        {
            cout << endl << "Matrix A(-1)" << endl;
            cout << A1 << endl;
        }

        MATRIX E(A1 * A);
        // максимальный внедиагональный элемент по модулю матрицы E
        int n = E.rows();

        if (size < 20)
        {
            cout << endl << "Matrix A(-1)*A" << endl;
            cout << E << endl;
        }
        for (int k = 0; k < n; k++)
            E(k, k) -= 1.0;

        cout << "Погрешность вычисления A(-1)*A - E: " << E.normI() << endl;
    }

}

void TestLinearSystemSolve(char* appPath, char* path)
{
    GetFullPathInWD(appPath, "MatrixVectMult_in.txt", path);
    // считывание данных
    ifstream fs;
    fs.open(path);
    if (fs.is_open())
    {
        int m, n;
        fs >> m >> n;
        MATRIX A(m, n);
        fs >> A;
        cout << "Матрица СЛАУ A" << endl;
        cout << A;

        fs >> n;
        VECTOR v(n);
        fs >> v;
        cout << endl << "Вектор правой части СЛАУ v" << endl;
        cout << v << endl;

        fs.close();

        VECTOR x(A.rows());
        double det = A.Determinant();
        bool res = true; 
        auto start = std::chrono::steady_clock::now();
        //res = Gauss(A, v, x);
        //CompactSchemeSolve(A, v, x);
        QRDecompositionSolve(A, v, x);
        //LLTDecompositionSolve(A, v, x);
        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        std::chrono::duration<double> secs = duration;

        if ((det != 0.0 || !isnan(det)) && res)
        {
            cout << "Время решения " << secs.count() << " сек" << endl;
            cout << endl << "Вектор  решения CЛАУ A*x = v" << endl;
            cout << x;

            cout << endl << "Определитель марицы A равен " << det << endl;

            //VECTOR vn(x.size());
            // vn = A*x - v;
            double norm = (A * x - v).norm();
            cout << endl << "Норма невязки " << norm << endl;

        }

        
    }
}
// Тестирование решения СЛАУ различного порядка
void TestLinearSystemSolve2()
{
    int  n = 1000;
    srand(10);
    // заполнение матрицы коэффициентов СЛАУ и вектора правой части с помощью генерации случайных чисел
    MATRIX A(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A(i, j) = rand_range(-1.0e3, 1.0e3);

    double det = A.Determinant();
    if (abs(det) < 1.0e-36)
    {
        cout << "Матрица СЛАУ вырождена" << endl;
        return;
    }

    VECTOR v(n), x(n);
    for (int i = 0; i < n; i++)
        v[i] = rand_range(-10.0, 200);

    cout << "Метод" << '\t' << "Время, сек" << '\t' << "Норма невязки" << endl;

    for (int i = 0; i < 4; i++)
    {
        auto start = std::chrono::steady_clock::now();
        switch (i)
        {
        case 0:
            Gauss(A, v, x); break;
        case 1:
            LUDecompositionSolve(A, v, x); break;
        case 2:
            CompactSchemeSolve(A, v, x); break;
        case 3:
            QRDecompositionSolve(A, v, x); break;
        }

        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        std::chrono::duration<double> secs = duration;

        cout << sys_methods[i] << '\t' << secs.count() << '\t' << (A * x - v).norm() << endl;
    }

}
void TestEigeValuesAndVectors(char* appPath, char* path)
{
    GetFullPathInWD(appPath, "Matrix_in.txt", path);
    // считывание данных
    ifstream fs;
    fs.open(path);
    if (fs.is_open())
    {
        int m, n;
        fs >> m >> n;
        MATRIX A(m, n);
        fs >> A;
        cout << "Матрица для поиска собственных значений" << endl;
        cout << A << endl;
        fs.close();

        complex<double>* lambda = new complex<double>[n];
        complex<double>** vects;
        vects = new complex<double>*[n];
        for (int i = 0; i < n; i++)
            vects[i] = new complex<double>[n];

        A.EigenvaluesAndVectorsLeVerrierFaddeev(lambda, vects);

        cout << "Вектор собственных значений" << endl;
        for (int i = 0; i < n; i++)
            cout << lambda[i] << '\t';
        cout << endl;

        cout << "Собственные вектора матрицы" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                complex<double> val = vects[i][j];
                cout << val << '\t';
            }
            cout << endl;

        }
        delete[] lambda;
        for (int i = 0; i < n; i++)
            delete[] vects[i];
        delete[] vects;

    }
}

    void TestRotations(char* appPath, char* path)
    {
        GetFullPathInWD(appPath, "Matrix_in.test", path);
        // считывание данных
        ifstream fs;
        fs.open(path);
        if (fs.is_open())
        {
            int m, n;
            fs >> m >> n;
            MATRIX A(m, n);
            fs >> A;
            cout << "Матрица для поиска собственных значений" << endl;
            cout << A << endl;
            fs.close();

            double* lambda = new double[n];
            double** vects = new double*[n];
            for (int i = 0; i < n; i++)
                vects[i] = new double[n];

            A.Rotate(lambda, vects);

            cout << "Вектор собственных значений" << endl;
            for (int i = 0; i < n; i++)
                cout << lambda[i] << '\t';
            cout << endl;

            cout << "Собственные вектора матрицы" << endl;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double val = vects[i][j];
                    cout << val << '\t';
                }
                cout << endl;

            }
            delete[] lambda;
            for (int i = 0; i < n; i++)
                delete[] vects[i];
            delete[] vects;

        }

}

int main(int argc, char *argv[])
{
    setlocale(LC_ALL, ""); // для от ображения кириллицы
    char path[1024]; // буфер пути файла данных
    //TestMatrixReverse();
    TestLinearSystemSolve2();
    //TestMatrixNVector(argv[0], path);
    //TestLinearSystemSolve(argv[0], path);
    //TestEigeValuesAndVectors(argv[0], path);
    //TestRotations(argv[0], path);
    
}