# Шаблоны классов MATRIX\<T\>("Матрица") и VECTOR\<T\> ("Вектор") и их приложение для решения системы линейных алгебраических уравнений, обращения матрицы и поиска собственных значений и собственных векторов матрицы
Реализация шаблонов классов MATRIX\<T\>("Матрица") и VECTOR\<T\>("Вектор")в проекте _MatrixVector_ в файле _MatrixVectorTemplate.hpp_. Для того, чтобы воспользоваться функциями классов в своём проекте, следует вставить строку в файл своего исходного кода _#include "MatrixVectorTemplate.hpp"_.

Элементы вектора и матрицы хранятся в сплошных массивах _m_data_, указатель на тип T, скрытый член класса. i-ый элемент вектора: _\*(m_data + i)_, элементы матрицы $A_{ij}$ распологаются 
в массиве построчно: _\*(m_data + i\*m_columns+j)_, где _m_columns_ - число столбцов матрицы. Шаблон класса содержит перегрузку операторов _=, \*, +, -, <<, >>, /=_:
```CPP
	template <typename T>friend MATRIX<T> operator*(const MATRIX<T>& , const MATRIX<T>& );
	template <typename T> friend ostream& operator<<(ostream& s, MATRIX<T>& );
	template <typename T> friend istream& operator>>(istream& s, MATRIX<T>& );
	template <typename T> friend VECTOR<T> operator*(const MATRIX<T>& , VECTOR<T>& );
	template <typename T> friend MATRIX<T> operator*(T , const MATRIX<T>& );
	template <typename T> friend MATRIX<T> operator+(const MATRIX<T>& , const MATRIX<T>& );
	template <typename T> friend MATRIX<T> operator-(const MATRIX<T>& , const MATRIX<T>& );
	MATRIX& operator/=(T );
	MATRIX& operator*=(const MATRIX& matr);
```
Таким образом, реализованы сложение, вычитание матриц и умножение матриц, умножение матрицы на вектор, умножение матрицы на скаляр и деление элементов матрицы на скаляр. Также реализован ввод-вывод 
элементов матрицы в поток вводв-вывода. Элементы марицы вводятся и выводятся построчно. Для вектора перегружены операции _=, \*, +, -, <<, >>_:
```CPP
VECTOR& operator=(const VECTOR& );
template <typename T> friend T operator*(const VECTOR<T>& , const VECTOR<T>& );
template <typename T> friend VECTOR<T> operator*(T a, const VECTOR<T>& );
template <typename T> friend VECTOR<T> operator+(const VECTOR<T>& , const VECTOR<T>& );
template <typename T> friend VECTOR<T> operator-(const VECTOR<T>& , const VECTOR<T>& );
template <typename T> friend ostream& operator<<(ostream& s, VECTOR<T>& );
template <typename T> friend istream& operator>>(istream& s, VECTOR<T>& );
```
Умножение вектора на вектор возвращает скалярное произведение векторов. Ввод-вывод в поток осушествляется в одной строке. 
# РЕШЕНИЕ системы линейных алгебраических уравнений (СЛАУ)
С помощью классов MATRIX\<T\> и VECTOR\<T\> реализовано решение СЛАУ несколькими методами: методом Гаусса, методом LU декомпозиции, LL<sup>T</sup> для симметричных положительно определённых матриц, QR декомпозиции и методом компактной схемы 
исключения [^1] в виде дружественных функций классов MATRIX<T> и VECTOR<T>:
```CPP
// Решение СЛАУ Ax = b
// A - матрица коэффициентов системы уравнений
// b - вектор правой части системы уравнений
// x - вуктор решений системы уравнений
// Матрица СЛАУ должна быть квадратной, а число элементов в векторе правой части и векторе решения
// должно совпадать с числом строк матрицы

// Решение СЛАУ методом Гаусса
template <typename T> friend bool Gauss(const MATRIX<T>& A, const VECTOR<T>& b, VECTOR<T>& x);

// Решение СЛАУ компактной схемой исключения
template <typename T> friend void CompactSchemeSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

// Решение СЛАУ методом декомпозиции A = QR, Q - ортогональная матрица, R - верхняя треугольная матрица
template <typename T> friend void QRDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

// Решение СЛАУ методом декомпозиции симметриченой положительно определённой матрицы A = LLT,
// L - нижняя треугольная матрица, LT - транспонированная к матрице L
template <typename T> friend void LLTDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

// Решение СЛАУ методом декомпозиции A = LU, где L - нижняя треугольная матрица, U - верхняя треугольная матрица
template <typename T> friend void LUDecompositionSolve(MATRIX<T>& A, VECTOR<T>& b, VECTOR<T>& x);

```
# Поиск обратной матрицы
Поиск обратной матрицы осуществлён в двух функциях шаблона класса MATRIX<T>:
```CPP
	MATRIX Reverse(); // обращение матрицы вычислением алгебраических дополнений
	MATRIX Invert(); // Решением система AX = I, I - единичная матрица, X - обратная матрица

```
# Поиск собственных значений матрицы 
Поиск собственных значений осуществлён в трёх функциях шаблона класса MATRIX<T>: методом Крылова А.Н., Леверье-Фаддева и методом вращений для симметричной матрицы [^2]\:
```CPP
	// lambda - массив собственных значений матрицы
	// vect - массив собственных векторов, соответствующих собственным значениям матрицы
	// Поиск собственных значений и собственных векторов матрицы методом Крылова А.Н.
	void EigenvaluesAndVectorsKrylov(complex<T>* lambda, complex<T>** vect);
	// Поиск собственных значений и собственных векторов матрицы методом Леверье-Фаддеева
	void EigenvaluesAndVectorsLeVerrierFaddeev(complex<T>* lambda, complex<T>** vect);
	// Поиск собственных значений и собственных векторов симметричной матрицы методом Леверье-Фаддеева
	void Rotate(T* lambda, T** vect);

```
В общем случае, собственные значения матрицы могут быть как комплексными, так и вещественными числами.
## Примеры
### Решение СЛАУ
```CPP
// Файл "MatrixVectMult_in.txt" для считывания элементо матрицы и вектора правой части
// Матрица считывается построчно, вектор правой части системы - отдельной строкой
/*
6 6
1789.0 17.0 45.0 1205.0 13.0 23.09
16.0 22.0 48.0 34.0 56088.0001 3456.0
99.0 21.0 14.0 7.0 2.0 782.0
112.0 117.0 29.0 456.0234 22.0 435.0
0.345 0.0089 1.004 2.67 16.42 0.009
0.009 0.003 240.34 0.998 0.87 0.22
6
19.0 77.0 66.0 11.0 22.0 234.0
*/

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

    GetFullPathInWD(appPath, "MatrixVectMult_in.txt", path);
    // считывание данных
    ifstream fs;
    fs.open(path);
    if (fs.is_open())
    {
        int m, n;
        fs >> m >> n;
        MATRIX<double> A(m, n);
        fs >> A;
        cout << "Матрица СЛАУ A" << endl;
        cout << A;

        fs >> n;
        VECTOR<double> v(n);
        fs >> v;
        cout << endl << "Вектор правой части СЛАУ v" << endl;
        cout << v << endl;

        fs.close();

        VECTOR<double> x(A.rows());
        double det = A.Determinant();
        bool res = true; 
        auto start = std::chrono::steady_clock::now();

        //res = Gauss(A, v, x); // метод Гаусса
        CompactSchemeSolve(A, v, x); // компактная схема исключения
        //QRDecompositionSolve(A, v, x); // QR декомпозиция
        LUDecompositionSolve(A, v, x); // LU декомпозиция
        // LLT декомпозиция
        /*MATRIX<double> At(A.Transpose());
        MATRIX<double> a(A.rows(), A.columns());
        a = At * A;
        VECTOR<double> b(v.size());
        b = At * v;
        LLTDecompositionSolve(a, b, x);*/

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
            double norm = (A * x - v).norm(); // евклидова норма вектора невязки
            cout << endl << "Норма невязки " << norm << endl;

        }

```
### Обращение матрицы
```CPP
    int size = 100;
    MATRIX<double> A(size, size), A1(size, size);
    srand(10);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            A(i, j) = (double)rand() / (double)rand();
 
    if (size < 20)
    {
        cout << endl << "Matrix A" << endl;
        cout << A << endl;
    }
    
    //clock_t  time_end;
    auto start = std::chrono::steady_clock::now();
    //A1 = A.Reverse(); // вычислением алгебраических дополнений
    //A1 = A.Invert();    // решением системы
    auto end = std::chrono::steady_clock::now();
    //time_end = clock();
    //double secs = (double)time_end / CLOCKS_PER_SEC;
    auto duration = end - start;
    std::chrono::duration<double> secs = duration;
    cout << "Время решения " << secs.count() << " сек" << endl;

    if (size < 20)
    {
        cout << endl << "Matrix A(-1)" << endl;
        cout << A1 << endl;
    }

    MATRIX<double> E(A.rows(), A.columns());
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

    cout << "Максимальный по модулю внедиагональный элемент A(-1)*A: " << emax << endl;


```
# Поиск собственных значений и собственных векторов матрицы
```CPP
// *** Поиск собственных значений квадраьной матрицы общего вида ***
// Файл "Matrix_in.test"
/*
4 4
2.2 1.0 0.5 2.0
1.0 1.3 2.0 1.0
0.5 2.0 0.5 1.6
2.0 1.0 1.6 2.0
*/
    GetFullPathInWD(appPath, "Matrix_in.test", path);
    // считывание данных
    ifstream fs;
    fs.open(path);
    if (fs.is_open())
    {
        int m, n;
        fs >> m >> n;
        MATRIX<double> A(m, n);
        fs >> A;
        cout << "Матрица для поиска собственных значений" << endl;
        cout << A << endl;
        fs.close();

        complex<double>* lambda = new complex<double>[n];
        complex<double>** vects = new complex<double>*[n];
        for (int i = 0; i < n; i++)
            vects[i] = new complex<double>[n];

        A.EigenvaluesAndVectorsLeVerrierFaddeev(lambda, vects);
        //A.EigenvaluesAndVectorsKrylov(lambda, vects);

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
            delete vects[i];
        delete[] vects;

    }

// *** Поиск собственных значений симметричной матрицы методом вращений ***
        GetFullPathInWD(appPath, "Matrix_in.test", path);
        // считывание данных
        ifstream fs;
        fs.open(path);
        if (fs.is_open())
        {
            int m, n;
            fs >> m >> n;
            MATRIX<double> A(m, n);
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

```
[^1]: Корн, Г. Справочник по математике для научных работников и инженеров : определения, теоремы, 
формулы / Г. Корн, Т. Корн ; пер. с англ. под общ. ред. И. Г. Арамановича. — 3‑е изд. — Москва : Наука, 1973. — 832  с.
[^2]: Долгополоа Д. В. Метод нахождения собственных значений и собственных векторов матриц / Долгополоа Д. В.; 
Методические указания - СПб: Санкт-Петербургский государственный технологический институт
(Технический университет), ИК «Синтез», Кафедра прикладной математики, 2005.- 40 с.
URL: https://sa.spbti.ru/files/chumakov/Sobstvennye%20znachenija.pdf (дата обращения 8.12.2025).
