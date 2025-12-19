#include <iostream>
#include <complex>
using namespace std;
#define EPS 1.0e-12
#define MAX_ITER_NUMBER 30000

template <typename T>
inline T randomDouble()
{
    return (T)(rand()) / (T)(rand());
}

// Заполнение массива начальных значений корней полинома init_values
// при помощи генератора случайных чисел
// n - степень полинома.
// Если полином нечётной степени, то первый элемент массива - вещественное число
template <typename T>
void FillInitialValues(complex<T>* init_values, int n)
{
    srand(1);
    int i = 0;
    while (i < n)
    {
        if (i == 0) init_values[i] = complex<T>(randomDouble<T>(), n % 2 != 0 ? 0 : randomDouble<T>());
        else
        {
            init_values[i] = complex<T>(randomDouble<T>(), randomDouble<T>());
            if (i == n - 1) break;
            i++;
            init_values[i] = n % 2 != 0 ? conj(init_values[i - 1]) : complex<T>(randomDouble<T>(), 0);
        }
        i++;
    }
}


// Полином koeff[n]*x ^ n + koeff[n - 1] * x ^ (n - 1) + ... + koeff[0] 
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
template <typename T>
complex<T> Polyfun(const complex<T>& z, const T* koeff, int n)
{
    complex<T> poly(0, 0);
    poly += koeff[0] + koeff[1] * z + koeff[2] * z * z;
    T p = 0;
    for (int i = 3; i <= n; i++)
    {
        p = (T)i;
        poly += koeff[i] * pow(z, p);

    }
    return poly;
}

// Производная полинома n-ой степени koeff[n]*x ^ n + koeff[n - 1] * x ^ (n - 1) + ... + koeff[0]
template <typename T>
complex<T> PolyfunDerivative(const complex<T>& z, const T* koeff, int n)
{
    complex<T> poly(0, 0);
    poly += koeff[1] + 2.0 * koeff[2] * z + 3.0 * koeff[3] * z * z;
    T p = 0;
    for (int i = 4; i <= n; i++)
    {
        p = (T)(i - 1);
        poly += ((T)i) * koeff[i] * pow(z, p);

    }
    return poly;
}


// Поиск всех корней полинома степени n с вещественными коэффициентами 
// методом Дюрана-Кернера (https://en.wikipedia.org/wiki/Durand–Kerner_method)
// Полином koeff[n]*x^n+koeff[n-1]*x^(n-1)+...+koeff[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n
template <typename T>
void PolyrootsDC(const T* koeff, int n, complex<T>* roots)
{
    complex<T>* prev_iter = new complex<T>[n]; // массив начальных значений корней полинома
    FillInitialValues<T>(prev_iter, n);
    for (int i = 0; i < n; i++)
        roots[i] = prev_iter[i];

    T err = 1;
    int iter = 0;
    while (err >= EPS && iter < MAX_ITER_NUMBER)
    {
        for (int i = 0; i < n; i++)
        {
            complex<T> pr(1, 0);
            for (int j = 0; j < n; j++)
            {
                if (i != j) pr *= roots[i] - roots[j];
            }
            roots[i] = prev_iter[i] - (Polyfun<T>(prev_iter[i], koeff, n) / koeff[n]) / pr;
        }

        // Проверка максимального модуля разности текущей и предыдущей итерации корней
        err = -1.0;
        T err_max = 0;

        for (int i = 0; i < n; i++)
        {
            err_max = abs(roots[i] - prev_iter[i]);
            if (err_max > err) err = err_max;
            prev_iter[i] = roots[i];
        }
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
template <typename T>
void PolyrootsAE(const T* koeff, int n, complex<T>* roots)
{
    complex<T>* w = new complex<T>[n]; // массив начальных значений корней полинома, далее используется для чисел смещения
    FillInitialValues(w, n);
    for (int i = 0; i < n; i++)
        roots[i] = w[i];

    T err = 1.0;
    int iter = 0;
    while (err >= EPS && iter < MAX_ITER_NUMBER)
    {
        for (int i = 0; i < n; i++)
        {
            complex<T> pr(0, 0);
            for (int j = 0; j < n; j++)
            {
                if (i != j) pr += 1.0 / (roots[i] - roots[j]);
            }
            complex<double> pp = Polyfun(roots[i], koeff, n) / PolyfunDerivative(roots[i], koeff, n);
            w[i] = pp / (1.0 - pp * pr);
        }

        // проверка максимального модуля чисел смещений
        err = -1.0;
        double err_max = 0.0;

        for (int i = 0; i < n; i++)
        {
            err_max = abs(w[i]);
            if (err_max > err) err = err_max;
            roots[i] -= w[i];
        }
        iter++;
#ifdef _DEBUG
        cout << "Кол-во итераций: " << iter << ". Значение погрешности вычислений: " << err << endl;
#endif

    }
    delete[] w;
}
