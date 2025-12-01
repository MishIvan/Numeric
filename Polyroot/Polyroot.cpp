#include <iostream>
#include <complex>
using namespace std;
#define EPS 1.0e-12


// Полином a[n]*x ^ n + a[n - 1] * x ^ (n - 1) + ... + a[0] 
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1

complex<double> Polyfun(const complex<double> &z, const double* koeff, const int n)
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

// Производная полинома n-ой степени

complex<double> PolyfunDerivative(const complex<double> &z, const double* koeff, const int n)
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
// Полином a[n]*x^n+a[n-1]*x^(n-1)+...+a[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n
// prev_iter - массив начальных значений корней полинома

void PolyrootsDC(const double* koeff, const int n, complex<double>* roots, complex<double> *prev_iter)
{
    for (int i = 0; i < n; i++)
        roots[i] = prev_iter[i];

    double err = 1.0;
    int iter = 0;
    while (err >= EPS && iter < 10000)
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

}

// Поиск всех корней полинома степени n с вещественными коэффициентами 
// методом Аберта-Эрлиха (https://en.wikipedia.org/wiki/Aberth_method)
// Полином a[n]*x^n+a[n-1]*x^(n-1)+...+a[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n
// w - массив начальных значений корней полинома, далее используется для чисел смещения

void PolyrootsAE(const double* koeff, const int n, complex<double>* roots, complex<double>* w)
{
    for (int i = 0; i < n; i++)
        roots[i] = w[i];

    double err = 1.0;
    int iter = 0;
    while (err >= EPS && iter < 10000)
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

}

int main()
{
    setlocale(LC_ALL, "");
    double a[8]{ 201,1,-16,-12,31,121,9,8 }; // массив вещественных коэффициентов полинома 
    int n = 7; // степень полинома
    complex<double> *roots = new complex<double>[n]; // массив корней полинома
    complex<double>* biter = new complex<double>[n]; // начальные значения корней полинома

    biter[0] = complex<double>(15, 0);
    biter[1] = complex<double>(0.5, 0.9);
    biter[2] = conj(biter[1]);
    biter[3] = complex<double>(0.7, 1.1);
    biter[4] = conj(biter[3]);
    biter[5] = complex<double>(0.35, 0.8);
    biter[6] = conj(biter[5]);


    PolyrootsAE(a, n, roots, biter); // поиск корней

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
    delete [] biter;
    
}


