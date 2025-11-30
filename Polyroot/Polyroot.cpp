// Polyroot.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <complex>
using namespace std;

// Полином x^n+a[n-1]/a[n]*x^(n-1)+...+a[0]/a[n]
complex<double> Polyfun(complex<double> z, const double* koeff, const int n)
{
    complex<double> poly(0, 0);
    poly += koeff[0] / koeff[n - 1] + koeff[1] / koeff[n - 1] * z + koeff[2] / koeff[n - 1] * z * z;
    for (int i = 3; i < n - 1; i++)
    {
        double p = (double)i;
        poly += koeff[i] / koeff[n - 1] * pow(z, p);

    }
    return poly + pow(z, (double)n);
}
// Поиск всех корней полинома степени n с вещественными коэффициентами методом Дюрана-Кернера
// коэффициент при степени n равен единице
// Полином a[n]*x^n+a[n-1]*x^(n-1)+...+a[0]
// roots - массив с корнями полинома
void Polyroots(const double* koeff, const int n, complex<double>* roots)
{
    
}

int main()
{
    std::cout << "Hello World!\n";
}


