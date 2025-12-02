#include <iostream>
#include <complex>
using namespace std;
#define EPS 1.0e-12
#define MAX_ITER_NUMBER 30000

// генератор вещественных случайных чисел
double randomDouble()
{
    return (double)(rand()) / (double)(rand());
}

// Заполнение массива начальных значений корней полинома init_values
// при помощи генератора случайных чисел
// n - степень полинома.
// Если полином нечётной степени, то первый элемент массива - вещественное число

void FillInitialValues(complex<double>* init_values, int n)
{
    srand(time(0));
    int i = 0;
    while (i < n)
    {
        if (i == 0) init_values[i] = complex<double>(randomDouble(), n % 2 != 0 ? 0 : randomDouble());
        else
        {
            init_values[i] = complex<double>(randomDouble(), randomDouble());
            i++;
            init_values[i] = conj(init_values[i - 1]);
        }
        i++;
    }
}


// Полином koeff[n]*x ^ n + koeff[n - 1] * x ^ (n - 1) + ... + koeff[0] 
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1

complex<double> Polyfun(const complex<double> &z, const double* koeff, int n)
{
    complex<double> poly(0, 0);
    poly += koeff[0] + koeff[1] * z + koeff[2] * z * z;
    double p = 0.0;
    for (int i = 3; i <= n; i++)
    {
        p = (double)i;
        poly += koeff[i] * pow(z, p);

    }
    return poly;
}

// Производная полинома n-ой степени koeff[n]*x ^ n + koeff[n - 1] * x ^ (n - 1) + ... + koeff[0]

complex<double> PolyfunDerivative(const complex<double> &z, const double* koeff, int n)
{
    complex<double> poly(0, 0);
    poly += koeff[1] + 2.0 * koeff[2] * z + 3.0 * koeff[3] * z * z;
    double p = 0.0;
    for (int i = 4; i <= n; i++)
    {
        p = (double)(i-1);
        poly += ((double)i)*koeff[i] * pow(z, p);

    }
    return poly;
}


// Поиск всех корней полинома степени n с вещественными коэффициентами 
// методом Дюрана-Кернера (https://en.wikipedia.org/wiki/Durand–Kerner_method)
// Полином koeff[n]*x^n+koeff[n-1]*x^(n-1)+...+koeff[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n

void PolyrootsDC(const double* koeff, int n, complex<double>* roots)
{
    complex<double>* prev_iter = new complex<double>[n]; // массив начальных значений корней полинома
    FillInitialValues(prev_iter, n);
    for (int i = 0; i < n; i++)
        roots[i] = prev_iter[i];

    double err = 1.0;
    int iter = 0;
    while (err >= EPS && iter < MAX_ITER_NUMBER)
    {
        for (int i = 0; i < n; i++)
        {            
            complex<double> pr(1, 0);
            for (int j = 0; j < n; j++)
            {
                if (i != j) pr *= roots[i] - roots[j];
            }
            roots[i] = prev_iter[i] - (Polyfun(prev_iter[i], koeff, n) / koeff[n]) / pr;            
        }
        
         err = 0;

        for (int i = 0; i < n; i++)
        {
            err += pow(abs(roots[i] - prev_iter[i]), 2.0);
            prev_iter[i] = roots[i];
        }
        err = sqrt(err);
        iter++;
#ifdef _DEBUG
        cout << "Кол-во итераций: " << iter << ". Значение погрешности вычислений: " << err << endl;
#endif

    } 
    delete[] prev_iter;
}

// Поиск всех корней полинома степени n с вещественными коэффициентами 
// методом Аберта-Эрлиха (https://en.wikipedia.org/wiki/Aberth_method)
// Полином koeff[n]*x^n+koeff[n-1]*x^(n-1)+...+koeff[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n

void PolyrootsAE(const double* koeff, int n, complex<double>* roots)
{
    complex<double> *w = new complex<double>[n]; // массив начальных значений корней полинома, далее используется для чисел смещения
    FillInitialValues(w, n);
    for (int i = 0; i < n; i++)
        roots[i] = w[i];

    double err = 1.0;
    int iter = 0;
    while (err >= EPS && iter < MAX_ITER_NUMBER)
    {
        for (int i = 0; i < n; i++)
        {
            complex<double> pr(0, 0), u(1,0);
            for (int j = 0; j < n; j++)
            {
                if (i != j) pr += u / (roots[i] - roots[j]);
            }
            complex<double> pp = Polyfun(roots[i], koeff, n) / PolyfunDerivative(roots[i], koeff, n);
            w[i] = pp / (u - pp*pr);            
        }

        err = 0;

        for (int i = 0; i < n; i++)
        {
            err += pow(abs(w[i]), 2.0);
            roots[i] -= w[i];
        }
        err = sqrt(err);
        iter++;
#ifdef _DEBUG
        cout << "Кол-во итераций: " << iter << ". Значение погрешности вычислений: " << err << endl;
#endif

    }
    delete[] w;
}

int main()
{
    setlocale(LC_ALL, "");
    double a[8]{ 201,1,-16,-12,0.031,121,9,8 }; // массив вещественных коэффициентов полинома 
    int n = 7; // степень полинома
    complex<double> *roots = new complex<double>[n]; // массив корней полинома

    PolyrootsDC(a, n, roots); // поиск корней

    cout << "*** Корни полинома ***" << endl;
    for (int i=0; i < n; i++)
    {
        cout << roots[i] << "\t";
    }
    cout << endl;

    cout << "*** Ошибка вычисления ***" << endl;
    for (int i=0; i < n; i++)
    {
        cout << abs(Polyfun(roots[i], a, n)) << "\t";
    }
    cout << endl;


    delete[] roots;
}


