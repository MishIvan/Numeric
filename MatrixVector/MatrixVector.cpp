// MatrixVector.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <ctime>
#include "MATRIX.h"
#include "MatrixVector.h"


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
        cout << "Матрица A" << endl;
        cout << A;

        fs >> m >> n;
        MATRIX B(m, n);
        fs >> B;
        cout << endl << "Матрица B" << endl;
        cout << B;

        MATRIX C(A.rows(), B.columns());
        C = A * B;
        cout << endl << "Матрица C = A*B" << endl;
        cout << C;

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
        cout << "Матрица A" << endl;
        cout << A;

        fs >> n;
        VECTOR v(n);
        fs >> v;
        cout << endl << "Вектор v" << endl;
        cout << v;

        VECTOR vp(A.rows());
        vp = A* v;
        cout << endl << "Вектор  vp= A*v" << endl;
        cout << vp;

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
        cout << "Вектор размерностью " << v.size() << " считан из файла" << endl;
        cout << v << endl;
    }
    else
        cout << "Неудачная попытка чтения вектора из файла" << endl;

    /// Вектор. Считывание данных с консоли и запись в файл
    VECTOR v_in(5);
    cout << "Введите вектор размерностью 5 с консоли" << endl;
    cin >> v_in;
    cout << "Вектор разменостью " << v_in.size() << endl;
    cout << v_in << endl;
    GetFullPathInWD(appPath, "Vector_out.txt", path);
    success = VECTOR::writeToFile(path, v_in);
    if (success)
        cout << "Вектор, введённый с консоли, записан в файл" << endl;
    else
        cout << "Неудачная попытка записи вектора, введённого с консоли" << endl;


    /// Матрица. Считывание из файла
    GetFullPathInWD(appPath, "Matrix_in.txt", path);
    MATRIX matr(3, 3);
    success = MATRIX::readFromFile(path, matr);
    if (success)
    {
        cout << "Считаны данные матрицы " << matr.rows() << " x " << matr.columns() << endl;
        cout << matr << endl;
    }
    else
        cout << "Неудачная попытка считывания данных матрицы" << endl;

    ///  Матрица. Ввод с консоли и запись в файл
    MATRIX matr_in(3, 5);
    cout << "Введите матрицу 3 х 5 с консоли" << endl;
    cin >> matr_in;
    //cout << "Введена матрица" << endl;
    //cout << matr_in << endl;
    GetFullPathInWD(appPath, "Matrix_out.txt", path);
    success = MATRIX::writeToFile(path, matr_in);
    if (success)
        cout << "Матрица, введённая с консоли, записана в файл" << endl;
    else
        cout << "Неудачная попытка записи матрицы, введённой с консоли" << endl;

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
        cout << "Матрица R вырождена" << endl;
        return;
    }


    cout << endl << "Matrix Q" << endl;
    cout << Q << endl;

    cout << endl << "Matrix R" << endl;
    cout << R << endl;

    MATRIX rt(R.columns(), R.rows());
    rt = Q * R;

    cout << endl << "Matrix A = Q*R" << endl;
    cout << rt << endl;

    cout << endl << "Matrix E = Q^t*Q" << endl;
    rt = Q.Transpose() * Q;
    cout << rt << endl;
    
}
/// <summary>
/// Тест разложения Холецкого.  Разложение Холецкого работает для симметричнйо и положительно определённой матрицы
/// </summary>
/// <param name="A"></param>
void TestCholeskyDecomposuition(MATRIX& A)
{
    MATRIX L(A.rows(), A.columns()), Anorm(A.rows(), A.columns());
    Anorm = A.Transpose() * A;

    cout << endl << "Matrix A^t*A" << endl;
    cout << Anorm << endl;

    if (!Anorm.CholeskyDecomposition(L))
    {
        cout << "Разложение Холецкого неприменимо для матрицы" << endl;
        return;
    }

    cout << endl << "Matrix L" << endl;
    cout << L << endl;

    MATRIX rt(A.columns(), A.rows());
    rt = L * L.Transpose();

    cout << endl << "Matrix A^t*A = L*L^t" << endl;
    cout << rt << endl;

}

void TestMatrixReverse()
{
    int size = 10;
    MATRIX A(size, size), A1(size, size);
    srand(1);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            // A(i, j) = fmod((double)(i + 1) ,(double)(j + 2)) +0.25;
            A(i, j) = (double)rand() / (double)rand();
 
    if (size < 20)
    {
        cout << endl << "Matrix A" << endl;
        cout << A << endl;
    }
    
    clock_t  time_end;
    //A1 = A.Reverse(); // вычислением алгебраических дополнений
    //A1 = A.Invert();    // решением системы
    A1 = A.InvertFaddev();
    time_end = clock();
    double secs = (double)time_end / CLOCKS_PER_SEC;
    cout << "Время решения " << secs << " сек" << endl;

    if (size < 20)
    {
        cout << endl << "Matrix A(-1)" << endl;
        cout << A1 << endl;
    }

    MATRIX E(A.rows(), A.columns());
    E = A1 * A;
     
    // максимальный внедиагональный элемент по модулю матрицы E
    int n = E.rows();
    int m = E.columns();
    double emax = 0.0;
    for(int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
        {
            if (i != j && fabs(E(i, j)) > emax) emax = fabs(E(i, j));
        }
   
    if (size < 20)
    {
        cout << endl << "Matrix A(-1)*A" << endl;
        cout << E << endl;
    }

    cout << "Максимальный по модулю внедиагональный элемент матрицы A(-1)*A: " << emax << endl;


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
        clock_t  time_end;
        //res = Gauss(A, v, x);
        //CompactSchemeSolve(A, v, x);
        QRDecompositionSolve(A, v, x);
        //LLTDecompositionSolve(A, v, x);

        time_end = clock();
        double secs = (double)time_end/CLOCKS_PER_SEC;
        if ((det != 0.0 || !isnan(det)) && res)
        {
            cout << "Время решения " << secs << " сек" << endl;
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
        delete [] lambda;
        for (int i = 0; i < n; i++)
            delete [] vects[i];
        delete[] vects;

        
    }

}

int main(int argc, char *argv[])
{
    setlocale(LC_ALL, ""); // для от ображения кириллицы
    char path[1024]; // буфер пути файла данных
    TestMatrixReverse();
    //TestMatrixNVector(argv[0], path);
    //TestLinearSystemSolve(argv[0], path);
    //TestEigeValuesAndVectors(argv[0], path);
    
}