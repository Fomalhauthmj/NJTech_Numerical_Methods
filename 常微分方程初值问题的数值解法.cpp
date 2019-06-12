/*
    改进欧拉法
*/
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std;
double expected_1(double x)
{
    return -pow(2, 0.5) * pow(((2 * pow(x, 3) + 3 * x * x) / 6.0 + 7.0 / 6.0), 0.5);
}
double expected_2(double x)
{
    return x * x - 4 * exp(-x) - 2 * x + 5;
}
void Better_Euler(double x0, double y0, double h, int N, double f(double, double))
{
    int n = 1;
    double yp, yc, y_next;
    double x_next;
    do
    {

        x_next = x0 + h;
        yp = y0 + h * f(x0, y0);
        yc = y0 + h * f(x_next, yp);
        y_next = (yp + yc) / 2;
        cout<<setprecision(3)<<x_next<<"\t";
        cout<<setprecision(8)<<y_next<<"\t";
        cout<<setprecision(8)<<expected_1(x_next)<<"\t"<<endl;
        n++;
        x0 = x_next;
        y0 = y_next;
    } while (n <= N);
}
/*
    经典龙格-库塔法
        四阶精度
*/
void Runge_Kutta_4(double x0, double y0, double h, int N, double f(double, double))
{
    int n = 1;
    double y_next, x_next;
    double k1, k2, k3, k4;
    do
    {
        x_next = x0 + h;
        k1 = f(x0, y0);
        k2 = f(x0 + h / 2.0, y0 + h / 2.0 * k1);
        k3 = f(x0 + h / 2.0, y0 + h / 2.0 * k2);
        k4 = f(x_next, y0 + h * k3);
        y_next = y0 + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
        cout<<setprecision(3)<<x_next<<"\t";
        cout<<setprecision(10)<<y_next<<"\t";
        cout<<setprecision(10)<<expected_2(x_next)<<"\t"<<endl;
        n++;
        x0 = x_next;
        y0 = y_next;
    } while (n <= N);
}
void Task_7_1()
{
    auto f = [](double x, double t) {
        return (x * x + x) / t;
    };
    double h = 0.05;
    freopen("output.txt","w",stdout);
    Better_Euler(1, -2, h, 2 / h, f);
}
//ans =- 2 ^ (1 / 2) * ((x ^ 2 * (2 * x + 3)) / 6 + 7 / 6) ^ (1 / 2)
void Task_7_2()
{
    auto f = [](double x, double y) {
        return -y + x * x + 3;
    };
    double h = 0.05;
    freopen("output.txt","w",stdout);
    Runge_Kutta_4(0, 1, h, 3 / h, f);
}
//ans =x ^ 2 - 4 * exp(-x) - 2 * x + 5
int main()
{
    //Task_7_1();
    Task_7_2();
    system("pause");
}