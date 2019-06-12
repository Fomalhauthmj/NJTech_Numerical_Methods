#include <iomanip>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;
#define N 1000
/*
!:
初始近似值的搜索
    逐步搜索法：初步确定根的位置
    区间二分法：只能用于求单根，不能很好的处理重根和负根
*/
void Search_by_Step(double step, double left, double right, double func(double))
{
    cout << "逐步搜索法" << endl;
    time_t begin = clock();
    double pre = func(left);
    double temp;
    for (double i = left; i <= right + step; i += step)
    {
        temp = func(i);
        if (temp * pre < 0)
        {
            cout << "[" << i - step << "," << i << "]" << endl;
        }
        else if (fabs(temp - 0) < 1e-12)
        {
            cout << "x=" << i << endl;
        }
        pre = temp;
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
}
void Binary_Search(double e, double left, double right, double func(double))
{
    /*
    e：精度控制 
    left:左边界 
    right:右边界 
    func:函数
    */
    time_t begin = clock();
    double mid;
    double temp;
    int cnt = 0;
    while (right - left >= e)
    {
        mid = (left + right) / 2;
        temp = func(mid);
        cnt++;
        if (temp * func(left) > 0)
        {
            left = mid;
        }
        else
            right = mid;
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "x=" << (left + right) / 2 << " 二分次数:" << cnt << endl;
}
/*
!:
迭代法
    迭代法
    加权法加速
    Aitken method 埃特金加速法
    Steffensen's method 史蒂芬森迭代法

*/
void Iterative_Method(double x0, double e, double func(double x))
{
    /*
    x0：初值
    e:精度控制
    func：迭代函数
    */
    time_t begin = clock();
    int cnt = 1;
    double x1 = func(x0);
    while (fabs(x1 - x0) >= e)
    {
        if (cnt >= N)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            //cout << "x" << cnt << " = " << x1 << endl;
            cnt++;
            x0 = x1;
            x1 = func(x0);
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "x=" << x1 << " 迭代次数:" << cnt << endl;
}
void Iterative_Method(double x0, double e, double func(double x, int kind), int k)
{
    /*
    x0：初值
    e:精度控制
    func：迭代函数
    */
    time_t begin = clock();
    int cnt = 1;
    double x1 = func(x0, k);
    while (fabs(x1 - x0) >= e)
    {
        if (cnt >= N)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            cout << "x" << cnt << " = " << x1 << endl;
            cnt++;
            x0 = x1;
            x1 = func(x0, k);
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms";
    cout << " x=" << x1 << " 迭代次数:" << cnt << endl;
}
void Weight(double x0, double e, double func(double x, double c), double func_1(double))
{
    cout << "加权法" << endl;
    time_t begin = clock();
    int cnt = 1;
    double c = func_1(x0);
    double x1 = func(x0, c);
    while (fabs(x1 - x0) >= e)
    {
        if (cnt >= N)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            cnt++;
            x0 = x1;
            x1 = func(x0, c);
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "x=" << x1 << " 迭代次数:" << cnt << endl;
}
void Steffensen(double x0, double e, double func(double x))
{
    /*
    x0:初值
    e:精度控制
    func:迭代函数
    */
    time_t begin = clock();
    int cnt = 1;
    double x1 = x0 - pow(func(x0) - x0, 2) / (func(func(x0)) - 2 * func(x0) + x0);
    while (fabs(x1 - x0) >= e)
    {
        if (cnt >= N)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            cout << "x" << cnt << " = " << x1 << endl;
            cnt++;
            x0 = x1;
            x1 = x0 - pow(func(x0) - x0, 2) / (func(func(x0)) - 2 * func(x0) + x0);
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "x=" << x1 << " 迭代次数:" << cnt << endl;
}
/*
!:
牛顿迭代法：非线性方程线性化
    牛顿迭代法的修正
        简化牛顿迭代法 ：c
        牛顿下山法
        重根修正
        
*/
void Newton(double x0, double e, double func(double), double func_1(double))
{
    /*
    x0:初值
    e:精度控制
    func:函数
    func_1：func的一阶导数
    */
    time_t begin = clock();
    int cnt = 1;
    double x1 = x0;
    if (fabs(func_1(x0) - 0) < e)
    {
        cout << "奇异标记" << endl;
    }
    else
    {
        x1 = x0 - func(x0) / func_1(x0);
        while (fabs(x1 - x0) >= e)
        {
            if (cnt >= N)
            {
                cout << "迭代失败" << endl;
                return;
            }
            else
            {
                cout << "x" << cnt << " = " << x1 << endl;
                //cout << cnt << "\t" << x1 << endl;
                cnt++;
                x0 = x1;
                x1 = x0 - func(x0) / func_1(x0);
            }
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "x=" << x1 << " 迭代次数:" << cnt << endl;
}
void m_Newton(double x0, double e, double func(double), double func_1(double), int m)
{
    /*
    x0:初值
    e:精度控制
    func:函数
    func_1：func的一阶导数
    m:重根数
    */
    time_t begin = clock();
    int cnt = 1;
    double x1 = x0;
    if (fabs(func_1(x0) - 0) < e)
    {
        cout << "奇异标记" << endl;
    }
    else
    {
        x1 = x0 - m * func(x0) / func_1(x0);
        while (fabs(x1 - x0) >= e)
        {
            if (cnt >= N)
            {
                cout << "迭代失败" << endl;
                return;
            }
            else
            {
                cout << "x" << cnt << " = " << x1 << endl;
                //cout << cnt << "\t" << x1 << endl;
                cnt++;
                x0 = x1;
                x1 = x0 - m * func(x0) / func_1(x0);
            }
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "x=" << x1 << " 迭代次数:" << cnt << endl;
}
void Newton_sqrt(double c, double x0, double e)
{
    cout << "基于牛顿迭代法的计算平方根" << endl;
    time_t begin = clock();
    int cnt = 1;
    double x1 = (x0 + c / x0) / 2;
    while (fabs(x1 - x0) >= e)
    {
        if (cnt >= N)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            cnt++;
            x0 = x1;
            x1 = (x0 + c / x0) / 2;
        }
    }
    time_t end = clock();
    cout << "计算时间：" << end - begin << "ms" << endl;
    cout << "sqrt(x)=" << x1 << " 迭代次数:" << cnt << endl;
}
void Down_Hill(double x0, double ek, double e, double func(double), double func_1(double))
{
    /*
    x0:初值
    ek：下山因子控制
    e:精度控制
    func:函数
    func_1：func的一阶导数
    */
    double k = 1.0; //下山因子
    double x1 = x0 - k * func(x0) / func_1(x0);
    while (fabs(func(x1)) >= fabs(func(x0)))
    {
        if (k <= ek)
        {
            cout << "重选x0" << endl;
            return;
        }
        else
        {
            k /= 2;
            x1 = x0 - k * func(x0) / func_1(x0);
        }
    }
    cout << "下山成功,转入牛顿迭代法" << endl;
    cout << "下山因子:" << k << endl;
    //加快收敛速度 转入牛顿迭代法 保留x1为x0
    Newton(x1, e, func, func_1);
}
/*
!:
弦截法:避免计算导数
    单点弦法 导数项换成差商
    双点弦法
*/
double f(double x)
{
    return x * x + 2 * x * exp(x) + exp(2 * x);
}
//f的一阶导数
double f_1(double x)
{
    return 2 * x + exp(x) * 2 * (x + 1) + 2 * exp(2 * x);
}
double f_2(double x)
{
    return 4 - 4 * exp(2 * x);
}
double g(double x) //迭代形式
{
    return pow(x, 3) - 1;
}
double g_list(double x, int kind)
//迭代计算格式表
{
    switch (kind)
    {
    case 1:
        return (3 * x + 1) / (x * x);
    case 2:
        return (x * x * x - 1) / 3;
    case 3:
        return pow(3 * x + 1, 1.0 / 3);
    case 4:
        return 1 / (x * x - 3);
    case 5:
        return sqrt(3 + 1.0 / x);
    case 6:
        return x - (1.0 / 3) * (x * x * x - 3 * x - 1) / (x * x - 1);
    }
}
double g_1(double x)
{
    return 3 * pow(x, 2);
}
double g_weight(double x, double c)
{
    return (g(x) - c * x) / (1 - c);
}
double k(double x)
{
    //return x - log(x) - 2;
    return exp(x - 2);
}
double fH(double H)
{
    double B = 20;
    double n = 0.03;
    double S = 0.0002;
    double Q = 5;
    return (1 / n) * sqrt(S) * pow(B * H, 5.0 / 3) / pow(B + 2 * H, 2.0 / 3) - Q;
}
double fH_1(double H)
{
    double B = 20;
    double n = 0.03;
    double S = 0.0002;
    double Q = 5;
    return (1 / n) * sqrt(S) * ((5.0 / 3 * pow(B, 5.0 / 3) * pow(H, 2.0 / 3) * pow(B + 2 * H, -2.0 / 3)) - 4.0 / 3 * pow(B * H, 5.0 / 3) * pow(B + 2 * H, -5.0 / 3));
}
int main()
{

    /*
    for (int i = 1; i <= 6; i++)
    {
        cout << "迭代计算格式" << i << endl;
        Iterative_Method(-2.0, 1e-4, g_list, i);
        //注意收敛性对初值选取的影响
    }
    */
    /*
    double e = 1e-5;
    double ek = 0.01;
    double x0 = 0.0;
    Newton(x0, e, f, f_1);
    m_Newton(x0, e, f, f_1, 2);
    //二重根
    Down_Hill(x0, ek, e, f, f_1);
    */
    /*
    Iterative_Method(3.5,1e-5,k);
    Iterative_Method(0.5,1e-5,k);
    Steffensen(3.5, 1e-5, k);
    Steffensen(0.5, 1e-5, k);
    */

    Binary_Search(1e-5, 0, 10, fH);
    Newton(5, 1e-5, fH, fH_1);
    Down_Hill(5, 0.01, 1e-5, fH, fH_1);
    system("pause");
}