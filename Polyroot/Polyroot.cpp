#include <iostream>
#include <complex>
using namespace std;
#define EPS 1.0e-12


// Полином x ^ n + a[n - 1] / a[n] * x ^ (n - 1) + ... + a[0] / a[n]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1

complex<double> Polyfun(const complex<double> z, const double* koeff, const int n)
{
    complex<double> poly(0, 0);
    poly += koeff[0] + koeff[1] * z + koeff[2] * z * z;
    double p = 0.0;
    for (int i = 3; i < n; i++)
    {
        p = (double)i;
        poly += koeff[i] * pow(z, p);

    }
    p = (double)n;
    poly /= koeff[n];
    return poly + pow(z, (double)n);
}

// Поиск всех корней полинома степени n с вещественными коэффициентами 
// методом Дюрана-Кернера (https://en.wikipedia.org/wiki/Durand–Kerner_method)
// Полином a[n]*x^n+a[n-1]*x^(n-1)+...+a[0]
// n - степень полинома
// koeff - массив вещественных коэффициентов полинома размерностью n+1
// roots - массив с корнями полинома размерностью n
// prev_iter - массив начальных значений корней полинома

void Polyroots(const double* koeff, const int n, complex<double>* roots, complex<double> *prev_iter)
{
    // Полином x^n+a[n-1]/a[n]*x^(n-1)+...+a[0]/a[n]
    auto Polyfun = [koeff, n](const complex<double> &z) -> complex<double>
    {
        complex<double> poly(0, 0);
        poly += koeff[0] + koeff[1] * z + koeff[2] * z * z;
        double p = 0.0;
        for (int i = 3; i < n; i++)
        {
            p = (double)i;
            poly += koeff[i] * pow(z, p);

        }
        p = (double)n;
        poly /= koeff[n];
        return poly + pow(z, (double)n);
    };

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
            roots[i] = prev_iter[i] - Polyfun(prev_iter[i]) / pr;            
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


    Polyroots(a, n, roots, biter); // поиск корней

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


