#include "PolyrootsTemplate.hpp"
#include <chrono>

int main()
{
    setlocale(LC_ALL, "");
    //double a[8]{ 201,1,-16,-12,0.031,121,9,8 }; // массив вещественных коэффициентов полинома 
    int n = 20; // степень полинома
    // подбор коэффициентов полинома генератором случайных чисел
    double* a = new double[n];
    srand(10);
    for (int i = 0; i <= n; i++)
        a[i] = randomDouble<double>();
    complex<double> *roots = new complex<double>[n]; // массив корней полинома

    for (int i = 0; i < 2; i++)
    {
         auto start = std::chrono::steady_clock::now();
        // поиск корней
        if (i == 0)
        {
            cout << "--- Нахождение корней методом Дюрана-Кернера" << endl;
            PolyrootsDC(a, n, roots); 
        }
        else
        {
            cout << " -- Нахождение корней методом Аберта-Эрлиха" << endl;
            PolyrootsAE(a, n, roots);
        }
        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        std::chrono::duration<double> secs = duration;
        cout << "Время решения " << secs.count() << " сек" << endl;
        
        cout << "*** Корни полинома и ошибка вычисления ***" << endl;
        for (int i = 0; i < n; i++)
        {
            cout << roots[i] << '\t' << abs(Polyfun(roots[i], a, n))  << endl;
        }

    }

    delete[] roots;
}


