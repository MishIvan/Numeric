# Шаблоны классов MATRIX<T>("Матрица") и VECTOR<T> ("Вектор"). Их приложение для решения системы линейных алгебраических уравнений, обращения матрицы и поиска собственных значений и собственных векторов матрицы
Реализация шаблонов классов MATRIX<T>("Матрица") и VECTOR<T> ("Вектор"). Элементы вектора и матрицы хранятся в сплошных массивах m_data. i-ый элемент вектора: _\*(m_data + i)_, элементы матрицы $A_{ij}$ распологаются 
в массиве построчно: _\*(m_data + i\*m_columns+j)_, где _m_columns_ - число столбцов матрицы. Шаблон класса содержит перегрузку операторов _=, \*, +, -, <<, >>, /=_:
```
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
```
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
С помощью классов MATRIX<T> и VECTOR<T> реализовано решение СЛАУ несколькими методами: методом Гаусса, методом LU декомпозиции, LL<sup>T</sup> для симметричных положительно определённых матриц,
QR декомпозиции и методом компактной схемы исключения [^1]\ в виде дружественных функций классов MATRIX<T> и VECTOR<T>:
```
// Решение СЛАУ методом Гаусса
template <typename T> friend bool Gauss(const MATRIX<T>& a, const VECTOR<T>& b, VECTOR<T>& x);

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
```
	MATRIX Reverse(); // обращение матрицы вычислением алгебраических дополнений
	MATRIX Invert(); // Решением система AX = I, I - единичная матрица, X - обратная матрица

```
## Примеры
### Решение СЛАУ
```
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
```

```
[^1]: Корн, Г. Справочник по математике для научных работников и инженеров : определения, теоремы, 
формулы / Г. Корн, Т. Корн ; пер. с англ. под общ. ред. И. Г. Арамановича. — 5‑е изд. — Москва : Наука, 1984. — 831 с.
