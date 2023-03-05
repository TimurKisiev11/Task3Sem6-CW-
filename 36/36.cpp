// 36.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <cmath>
#include <windows.h>
#include <iomanip>
#include <cstdlib>

using namespace std;

//Ядро
double H(double x, double y)
{
    double a = -5.0 * y * cos(x);
    double b = x * y - 2.0;
    return(a / b);
}
//Функция правой части 
double f(double x)
{
    double a = exp(5 * x) - 1.0;
    double b = 2.0 * x;
    return(2.0 + a / b);
}
//Евклидова норма вектора 
double EVnorm_vec(double* vec, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
    {
        norm += fabs(pow(vec[i], 2));
    }
    norm = sqrt(norm);
    return norm;
}
//выводит матриицу A порядка n
void show_mat(double** A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << fixed << setprecision(6) << " \t" << A[i][j] << " \t";
        }
        std::cout << endl;
    }
}
//выводит вектор V
void show_vect(double* V, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << fixed << setprecision(6) << " \t" << V[i] << " \t";
        std::cout << endl;
    }
}
//строит нулевой фектор порядка n
double* get_null_vect(int n)
{
    double* V = new double[2 * n];
    for (int i = 0; i <= n; i++)
    {
        V[i] = 0.0;
    }
    return V;
}
//строит нулевую матрицу порядка n
double** get_null_mat(int n)
{
    double** A = new double* [2 * n];
    for (int i = 0; i <= n; i++)
    {
        A[i] = new double[2 * n];
        for (int j = 0; j <= n; j++)
            A[i][j] = 0.0;
    }
    return A;
}
//перемножение матриц
double** mul_mat(double** A, double** B, int n)
{
    double** C = new double* [n];
    for (int i = 0; i < n; i++)
    {
        C[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
            //cout << C[i][j] << " ";
        }
        //cout << endl;
    }
    // show(C,n);
    return C;
}
//для решения системы с матрицей "a" вектором правых частей "y" и порядком "n"
double* gauss2(double** a, double* y, int n)
{
    double* x, max;
    int k, index;
    const double eps = 0.00001;  // точность
    x = new double[n];
    k = 0;
    for (k = 0; k < n; k++)
    {
        // Поиск строки с максимальным a[i][k]
        max = fabs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (fabs(a[i][k]) > max)
            {
                max = fabs(a[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            return 0;
        }
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;
        // Нормализация уравнений
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (fabs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i] = y[i] - y[k];
        }
    }
    // обратная подстановка
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }
    return x;
}
//Многочлены Лежандра для поиска их корней
double Legendre(int n, double x)
{
    if (n == 0)
        return 1;
    else if (n == 1)
        return x;
    else
        return (Legendre(n - 1, x) * x * (2 * n - 1) / (n)) - (Legendre(n - 2, x) * (n - 1) / (n));
}
//Локализация корней многочленов Лежандра
double* loc_root(double a, double b, double N, int n)
{
    double* mas = new double[100];
    mas[0] = n;
    double h = (b - a) / N;
    double tmp = 0;
    int j = 2;
    //cout << "области единственности корня:" << endl;
    while (a + h <= b)
    {
        if ((Legendre(n, a) * Legendre(n, a + h)) <= 0)
        {

            //cout<<a<<" "<<a+h<< endl;
            mas[j] = a;
            mas[j + 1] = a + h;
            // cout << fixed << setprecision(25)<< "[" << mas[j] << " , " << mas[j + 1] << "]" << endl;
            j += 2;
        }
        a = a + h;
    }

    if (mas[0] > 0)
    {
        //cout << fixed << setprecision(0) << "количество корней ортоганального многочлена на отрезке (a,b) = " << mas[0] << endl;
        return mas;
    }

    else
    {
        cout << "нет корней на данном отрезке" << endl;
    }
    return 0;
}
//Коэффициенты для формулы Гаусса
double* koeff(double* z, int n)
{
    double* A = new double[100];
    for (int i = 0; i < n; i++)
    {
        A[i] = (2 * (1 - pow(z[i], 2)) / (pow(n, 2) * Legendre(n - 1, z[i]) * Legendre(n - 1, z[i])));
    }
    return A;
}
//Метод секущих 
double Secant(double a, double b, double EPS, int n)
{
    double fapprox, sapprox, x;
    int counter = 0;

    fapprox = a;
    sapprox = b;

    x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
    counter++;

    while (abs(x - sapprox) >= EPS)
    {
        fapprox = sapprox;
        sapprox = x;
        x = sapprox - (Legendre(n, sapprox) / (Legendre(n, sapprox) - Legendre(n, fapprox)) * (sapprox - fapprox));
        counter++;
    }

    return x;
}
//Вернет матрицу для каждого отрезка [a,b]; Первая строка - узлы КФ Гаусса, вторая - коэффициенты.
double** Gauss(double a, double b, int n)
{
    double** result = get_null_mat(2 * n);
    double* z = get_null_vect(2 * n);        //массив для узлов для формулы Гаусса
    int val = n;                             // количество узлов для формулы Гаусса
    double* MASS = loc_root(-1, 1, 100, val);//локализуем узлы многочлена Лежандра
    double counter = MASS[0];
    int cnt = 0;
    for (int j = 2; j <= 2 * counter + 1; j += 2)
    {
        z[cnt] = Secant(MASS[j], MASS[j + 1], 0.0000000001, val);//находим узлы (корни многочлена Лежандра)
        cnt++;
    }
    if (val % 2 > 0)
    {
        int m = val / 2;
        z[m] = 0.0;
    }
    double* A = koeff(z, val); //находим коэффициенты формулы гаусса
   // cout << "Узлы и коэффициенты для формулы Гаусса: " << endl;
   // cout << "      Узел      " << " <-> " << "    Коэффициент    " << endl;
    for (int k = 0; k < cnt; k++)
    {
        z[k] = (b - a) / 2.0 * z[k] + (b + a) / 2.0; //узлы для формулы гаусса
        A[k] = (b - a) / 2.0 * A[k];
        // cout << fixed << setprecision(2) << z[k] << "            <->            " << A[k] << endl;
    }
    result[0] = z;
    result[1] = A;
    // cout << "Соберем в матрицу узлы и коэффициенты" << endl;
    // show_mat(result, n);
    return result;
}
//Зная коэффициенты и узлы Гаусса и искомый вектор коэффициентов рассчитывает решение в данной точке x
double func_value_mech_quad(double** gauss, double* z, int n, double x)
{
    double value = 0.0;
    for (int k = 0; k < n; k++)
    {
        value += gauss[1][k] * H(x, gauss[0][k]) * z[k];
        //  cout << gauss[1][k] << "   " << H(x, gauss[0][k]) << "   " << z[k] <<"   "<<f(x)<< endl;
    }
    return (value + f(x));
}
//Метод механических квадратур, вернет значение решения в точке x
double Mech_quad(double a, double b, int n, double x)
{
    double** D = get_null_mat(2 * n);
    double* g = get_null_vect(2 * n);
    double** gauss = Gauss(a, b, n);
    double delta = 0.0;
    for (int j = 0; j < n; j++)
    {
        for (int k = 0; k < n; k++)
        {
            if (k == j)
                delta = 1.0;
            else
                delta = 0.0;
            D[j][k] = delta - gauss[1][k] * H(gauss[0][j], gauss[0][k]);
            //    cout << j и<< "   " << k << "   " << delta << "   " << gauss[0][j]<<"   "<< gauss[0][k]<<"   "<< H(gauss[0][j], gauss[0][k]) <<"   "<<f(gauss[0][j])<< endl;
        }
        g[j] = f(gauss[0][j]);
    }

    //    cout << "матрица D:" << endl;
    //    show_mat(D, n );
    //    cout << "Вектор g" << endl;
    //    show_vect(g, n);
    double* z = gauss2(D, g, n);
    //    cout << "Искомые коэффициенты:" << endl;
    //    show_vect(z, n);
    return func_value_mech_quad(gauss, z, n, x);
}
//Функции a_i из разложения Ядра, i - индекс функции (0,1,...), вернет значение i-ой функции в точке x
double a_i(int index_of_func, double x)
{
    switch (index_of_func)
    {
    case 0:
        return (2.5);
        break;
    case 1:
        return (1.25 * x);
        break;
    case 2:
        return (0.625 * pow(x, 2));
        break;
    case 3:
        return (0.3125 * pow(x, 3));
        break;
    case 4:
        return ((5 / 96) * pow(x, 4));
        break;
    case 5:
        return ((5 / 192) * pow(x, 5));
    }
}
//Функции b_i из разложения Ядра, i - индекс функции (0,1,...), вернет значение i-ой функции в точке x
double b_i(int index_of_func, double y)
{
    switch (index_of_func)
    {
    case 0:
        return (y);
        break;
    case 1:
        return (pow(y, 2));
        break;
    case 2:
        return (pow(y, 3) - 2 * y);
        break;
    case 3:
        return  (pow(y, 4) - 2 * pow(y, 2));
        break;
    case 4:
        return (3 * pow(y, 5) - 6 * pow(y, 3) + 2 * y);
        break;
    case 5:
        return (3 * pow(y, 6) - 6 * pow(y, 4) + 2 * pow(y, 2));
    }
}
//Зная искомый вектор коэффициентов рассчитывает решение в данной точке x
double func_value_ker_method(double* c, int n, double x)
{
    double value = 0.0;
    for (int i = 0; i < n; i++)
    {
        value += c[i] * a_i(i, x);
    }
    return (value + f(x));
}
//Метод замены ядра на вырожденное, вернет значение решения в точке x
double Ker_method(double a, double b, int n, double x)
{
    double** A = get_null_mat(2 * n);
    double** gamma = get_null_mat(2 * n);
    double* B = get_null_vect(2 * n);
    int val = 7;
    double** gauss = Gauss(a, b, val);
    double delta = 0.0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j == i)
                delta = 1.0;
            else
                delta = 0.0;
            for (int k = 0; k < val; k++)
            {
                gamma[i][j] += (gauss[1][k] * (b_i(i, gauss[0][k]) * a_i(j, gauss[0][k])));
                //               cout << i << "   " << j << "   " << delta << "   " << b_i(i, gauss[0][k]) << "   " << a_i(j, gauss[0][k]) << endl;
            }
            A[i][j] = delta - gamma[i][j];
            //               cout << i << " " << j << "  " << A[i][j] << endl;

        }
        for (int k = 0; k < val; k++)
            B[i] += (gauss[1][k] * (b_i(i, gauss[0][k]) * f(gauss[0][k])));
    }

    //  cout << "матрица A:" << endl;
    //  show_mat(A, n);
    //  cout << "Вектор B" << endl;
    //  show_vect(B, n);
    double* c = gauss2(A, B, n);
    //  cout << "Искомые коэффициенты:" << endl;
    //  show_vect(c, n);
    return func_value_ker_method(c, n, x);
}
int main()
{
    SetConsoleOutputCP(1251);
    int n = 0;
    cout << "Введите порядок задачи - >";
    cin >> n;
    int m = 9;
    cout << "       Задание 3" << endl;
    cout << "       Вариант 11" << endl;
    cout << "Уравнение вида: u(x) - 0.7\int^1_0 (y*cos(x))/(x*y - 2) * u(y)dy = 2 + (e^(5x) - 1)/(2*x)" << endl;
    double a = 0.00000000000000001;
    double b = 1.0;
    double mid = (a + b) / 2.0;
    double** mech_qaud_vect = get_null_mat(5 * n);
    double** ker_method_vect = get_null_mat(5 * n);
    double** differ_vect = get_null_mat(5 * n);
    cout << "|---------------|---------------|----------------|----------------|\n";
    cout << "|      x        |        a      |     (a+b)/2    |        b       | \n";
    cout << "|---------------|---------------|----------------|----------------|\n";
    for (int n = 3; n < m; n++)
    {
        mech_qaud_vect[n][0] = Mech_quad(a, b, n, a );
        mech_qaud_vect[n][1] = Mech_quad(a, b, n, mid);
        mech_qaud_vect[n][2] = Mech_quad(a, b, n, b);
        cout << fixed << setprecision(6) << "|     u^" << n << "(x)    |  " << mech_qaud_vect[n][0] << "  |  " << mech_qaud_vect[n][1] << "   |   " << mech_qaud_vect[n][2] << "   | \n";
        cout << "|---------------|---------------|----------------|----------------|\n";
    }
    cout << "|---------------|---------------|----------------|----------------|\n";
    cout << "|                   Метод замены ядра на вырожденное              | \n";
    cout << "|---------------|---------------|----------------|----------------|\n";
    for (int n = 3; n <= 6; n++)
    {
        ker_method_vect[n][0] = Ker_method(a, b, n, a );
        ker_method_vect[n][1] = Ker_method(a, b, n, mid);
        ker_method_vect[n][2] = Ker_method(a, b, n, b);
        cout << "|     u^" << n << "(x)    |  " << ker_method_vect[n][0] << "  |  " << ker_method_vect[n][1] << "   |   " << ker_method_vect[n][2] << "   | \n";
        cout << "|---------------|---------------|----------------|----------------|\n";
    }
    for (int i = 0; i < 3; i++)
        for (int n = 3; n <= 6; n++)
        {
            differ_vect[n][i] = fabs(mech_qaud_vect[n][i] - ker_method_vect[n][i]);
        }
    cout << fixed << setprecision(8) << "Норма вектора невязки двух методов для n = 3: " << EVnorm_vec(differ_vect[3], 3) << endl;
    cout << fixed << setprecision(8) << "Норма вектора невязки двух методов для n = 4: " << EVnorm_vec(differ_vect[4], 4) << endl;
    cout << fixed << setprecision(8) << "Норма вектора невязки двух методов для n = 5: " << EVnorm_vec(differ_vect[5], 5) << endl;
    cout << fixed << setprecision(8) << "Норма вектора невязки двух методов для n = 6: " << EVnorm_vec(differ_vect[6], 6) << endl;
    return 0;
}


