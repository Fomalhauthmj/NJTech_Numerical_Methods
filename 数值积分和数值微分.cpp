/*
    复化求积法
        复化梯形
        复化辛普森
        复化柯斯特
            精度提高，计算量增加 √
    变步长求积：解决计算量的增加
        变步长梯形求积
    龙贝格算法

*/
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
using namespace std;
/*
    积分区间左端点：a
    积分区间右端点: b
    积分函数：     f(x)
    步长:          h
    精度：         e
*/
//复化梯形
double T(double a, double b, double f(double), double h)
{
    double ans = 0;
    double k = a + h;
    while (k < b)
    {
        ans += f(k);
        k += h;
    }
    ans *= 2;
    ans += f(a) + f(b);
    return ans * h / 2.0;
}
//复化辛普森
double S(double a, double b, double f(double), double h)
{
    double ans = 0;
    double k = a + h;
    while (k < b + h / 2.0)
    {
        ans += 4 * f(k - h / 2.0) + 2 * f(k);
        k += h;
    }
    ans += f(a) - f(b);
    return ans * h / 6.0;
}
//复化柯斯特
double C(double a, double b, double f(double), double h)
{
    double ans = 0;
    double k = a;
    while (k < b)
    {
        ans += 32 * f(k + h / 4.0) + 12 * f(k + h / 2.0) + 32 * f(k + 3 * h / 4.0);
        k += h;
    }
    k = a + h;
    while (k < b)
    {
        ans += 14 * f(k);
        k += h;
    }
    ans += 7 * (f(a) + f(b));
    return ans * h / 90.0;
}
//变步长梯形法
double betterT(double a, double b, double e, double f(double))
{
    double h = b - a;
    double t1 = T(a, b, f, h);
loop:
    //cout << "cur t1:" << setprecision(10) << t1 << endl;
    double S = 0;
    double x = a + h / 2.0;
    do
    {
        S += f(x);
        x += h;
    } while (x < b);
    double t2 = (t1 + h * S) / 2.0;
    while (fabs(t2 - t1) >= e)
    {
        h = h / 2.0;
        t1 = t2;
        goto loop;
    }
    return t2;
}
double GetNextT(double a, double b, double f(double), double preT, int k)
{
    /*
        2^k:等分数
    */
    if (k == 0)
    {
        return (f(a) + f(b)) * (b - a) / 2.0;
    }
    else
    {
        double S = 0;
        double h = (b - a) / (1 << (k - 1));
        double x = a + h / 2;
        do
        {
            S += f(x);
            x += h;
        } while (x < b);
        return (preT + h * S) / 2.0;
    }
}
double Romberg(double a, double b, double e, double f(double))
{
    vector<double> s[4]; //计算序列
    int k;               //区间等分数
    for (k = 0; k < 4; k++)
    {
        if (k != 0)
            s[0].push_back(GetNextT(a, b, f, s[0].back(), k));
        else
            s[0].push_back(GetNextT(a, b, f, 0, k));
    }
    /*
        s[0]：梯形序列
        s[1]：辛普森序列
        s[2]：柯特斯序列
        s[3]：龙贝格序列
    */
    auto TtoS = [](double t1, double t2) {
        return (4.0 * t2 - t1) / 3.0;
    };
    auto StoC = [](double s1, double s2) {
        return (16.0 * s2 - s1) / 15.0;
    };
    auto CtoR = [](double c1, double c2) {
        return (64.0 * c2 - c1) / 63.0;
    };
    s[1].push_back(TtoS(s[0][0], s[0][1]));
    s[1].push_back(TtoS(s[0][1], s[0][2]));
    s[1].push_back(TtoS(s[0][2], s[0][3]));
    s[2].push_back(StoC(s[1][0], s[1][1]));
    s[2].push_back(StoC(s[1][1], s[1][2]));
    s[3].push_back(CtoR(s[2][0], s[2][1]));
    cout << s[0][0] << endl;
    cout << s[0][1] << " " << s[1][0] << endl;
    cout << s[0][2] << " " << s[1][1] << " " << s[2][0] << endl;
    cout << s[0][3] << " " << s[1][2] << " " << s[2][1] << " " << s[3][0] << endl;
    //得到第一个R值
    //需要计算：4次T值->3次S值->2次C值->1次R值
    //得到之后的R值
    //需要计算：1次T值->1次S值->1次C值->1次R值
    do
    {
        s[0].push_back(GetNextT(a, b, f, s[0].back(), k++));
        s[1].push_back(TtoS(s[0].at(s[0].size() - 2), s[0].back()));
        s[2].push_back(StoC(s[1].at(s[1].size() - 2), s[1].back()));
        s[3].push_back(CtoR(s[2].at(s[2].size() - 2), s[2].back()));
    } while (fabs(s[3].back() - s[3].at(s[3].size() - 2)) >= e);
    int pos = 4;
    while (pos < s[0].size())
    {
        cout << s[0][pos] << " " << s[1][pos - 1] << " " << s[2][pos - 2] << " " << s[3][pos - 3] << endl;
        pos++;
    }
    return s[3].back();
}
void Task_6_1()
{
    double a = 0;
    double b = 1;
    double e = 1e-10;
    auto f = [](double x) {
        return cos(x * x);
    };
    int bits = 10;
    double h = 0.001;
    cout << "T:" << setprecision(bits) << T(a, b, f, h) << endl;
    cout << "S:" << setprecision(bits) << S(a, b, f, h) << endl;
    cout << "C:" << setprecision(bits) << C(a, b, f, h) << endl;
    cout << "better T:" << setprecision(bits) << betterT(a, b, e, f) << endl;
}
void Task_6_2()
{
    double a = 1;
    double b = 3;
    double e = 1e-8;
    auto f = [](double x) {
        return 1.0 / x;
    };
    int bits = 8;
    cout << "Romberg:" << setprecision(bits) << Romberg(a, b, e, f) << endl;
}
int main()
{
    //Task_6_1();
    Task_6_2();
    system("pause");
}