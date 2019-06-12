#include <cstring>
#include <iostream>
#include <math.h>
using namespace std;
/*
    插值法 Interpolation
        Lagrange插值
        Newton插值
        Aitken插值:逐次线性插值
        样条插值
*/
void Lagrange(double *x, double *y, int N, double X)
{
    double t;
    double ans = 0;
    for (int i = 0; i < N; i++)
    {
        t = 1;
        for (int j = 0; j < N; j++)
        {
            if (j != i)
            {
                t *= (X - x[j]) / (x[i] - x[j]);
            }
        }
        ans += t * y[i];
    }
    cout << X << " 1= " << ans << endl;
}
void Aitken(double *x, double *y, int N, double X)
{
    for (int k = 1; k < N; k++)
    {
        for (int i = k; i < N; i++)
        {
            y[i] = (X - x[k - 1]) / (x[i] - x[k - 1]) * y[i] + (X - x[i]) / (x[k - 1] - x[i]) * y[k - 1];
        }
    }
    cout << X << " = " << y[N - 1] << endl;
}
#define N 5
//n-1 for spline_2_nature
double *Chase(double A[][N], double *B)
{
    double *l = new double[N];
    double *u = new double[N];
    l[0] = A[0][0];
    for (int i = 0; i < N - 1; i++)
    {
        u[i] = A[i][i + 1] / l[i];
        l[i + 1] = A[i + 1][i + 1] - A[i + 1][i] * u[i];
    }
    double *y = new double[N];
    y[0] = B[0] / l[0];
    for (int i = 1; i < N; i++)
    {
        y[i] = (B[i] - A[i][i - 1] * y[i - 1]) / l[i];
    }
    double *ans = new double[N];
    ans[N - 1] = y[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        ans[i] = y[i] - u[i] * ans[i + 1];
    }
    /*
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << ans[i] << endl;
    }
    */
    return ans;
}
double Spline_1(double *x, double *y, int n, double y0_1, double yn_1, double X)
{
    double h[n];
    for (int i = 0; i < n; i++)
    {
        h[i] = x[i + 1] - x[i];
        //cout << "Cur h" << i << ":" << h[i] << endl;
    }
    double a[n];
    double b[n];
    double d[n + 1];
    for (int i = 1; i <= n - 1; i++)
    {
        a[i] = h[i - 1] / (h[i - 1] + h[i]);
        b[i] = 1.0 - a[i];
        //cout << "Cur a" << i << ":" << a[i] << endl;
        //cout << "Cur b" << i << ":" << b[i] << endl;
    }
    for (int i = 1; i <= n - 1; i++)
    {
        d[i] = 6.0 / (h[i - 1] + h[i]) * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        //cout << "Cur d" << i << ":" << d[i] << endl;
    }
    d[0] = 6.0 / h[0] * ((y[1] - y[0]) / h[0] - y0_1);
    d[n] = 6.0 / h[n - 1] * (yn_1 - (y[n] - y[n - 1]) / h[n - 1]);
    double A[N][N]; //三弯矩
    memset(A, 0, sizeof(A));
    A[0][0] = 2;
    A[0][1] = 1;
    A[n][n] = 2;
    A[n][n - 1] = 1;
    for (int i = 1; i < n; i++)
    {
        A[i][i - 1] = a[i];
        A[i][i] = 2.0;
        A[i][i + 1] = b[i];
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    /*
    for (int i = 0; i < N; i++)
    {
        cout << d[i] << " ";
    }
    cout << endl;
    */
    double *M = Chase(A, d);
    //矩解完成
    //根据X确定i
    int i = 0;
    while (x[i] < X)
    {
        i++;
    }
    i--;
    //返回值由4项组成
    double temp1 = pow(x[i + 1] - x[i], 3) / (6.0 * h[i]) * M[i];
    double temp2 = pow(X - x[i], 3) / (6.0 * h[i]) * M[i + 1];
    double temp3 = (x[i + 1] - X) / h[i] * (y[i] - h[i] * h[i] / 6.0 * M[i]);
    double temp4 = (X - x[i]) / h[i] * (y[i + 1] - h[i] * h[i] / 6.0 * M[i + 1]);
    return temp1 + temp2 + temp3 + temp4;
}
double Spline_2_nature(double *x, double *y, int n, double X)
{
    double h[n];
    for (int i = 0; i < n; i++)
    {
        h[i] = x[i + 1] - x[i];
        //cout << "Cur h" << i << ":" << h[i] << endl;
    }
    double a[n];
    double b[n];
    double d[n + 1];
    for (int i = 1; i <= n - 1; i++)
    {
        a[i] = h[i - 1] / (h[i - 1] + h[i]);
        b[i] = 1.0 - a[i];
        //cout << "Cur a" << i << ":" << a[i] << endl;
        //cout << "Cur b" << i << ":" << b[i] << endl;
    }
    for (int i = 1; i <= n - 1; i++)
    {
        d[i] = 6.0 / (h[i - 1] + h[i]) * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        //cout << "Cur d" << i << ":" << d[i] << endl;
    }
    double A[N][N]; //三弯矩
    memset(A, 0, sizeof(A));
    for (int i = 0; i < N; i++)
    {
        A[i][i] = 2.0;
    }
    for (int i = 0; i < N - 1; i++)
    {
        A[i][i + 1] = b[i + 1];
    }
    for (int i = 1; i < N; i++)
    {
        A[i - 1][i] = a[i + 1];
    }
    double D[N];
    for (int i = 0; i < N; i++)
    {
        D[i] = d[i + 1];
    }
    double *m = Chase(A, D);
    double M[n + 1];
    M[0] = M[n] = 0.0;
    for (int i = 1; i <= n - 1; i++)
        M[i] = m[i - 1];
    //矩解完成
    //根据X确定i
    int i = 0;
    while (x[i] < X)
    {
        i++;
    }
    i--;
    //返回值由4项组成
    double temp1 = pow(x[i + 1] - X, 3) / (6.0 * h[i]) * M[i];
    double temp2 = pow(X - x[i], 3) / (6.0 * h[i]) * M[i + 1];
    double temp3 = (x[i + 1] - X) / h[i] * (y[i] - h[i] * h[i] / 6.0 * M[i]);
    double temp4 = (X - x[i]) / h[i] * (y[i + 1] - h[i] * h[i] / 6.0 * M[i + 1]);
    return temp1 + temp2 + temp3 + temp4;
}
void Task_4_1()
{
    double X[] = {2, 13, 27};
    for (auto ele : X)
    {
        double x[] = {0, 5, 10, 15, 20, 25, 30};
        double y[] = {1.57079, 1.56780, 1.55888, 1.54415, 1.52379, 1.49811, 1.46746};
        Aitken(x, y, 7, ele);
        //Lagrange(x,y,7,ele);
    }
}
void Task_4_2()
{
    double t[] = {0, 1, 2, 3, 4, 5, 6};
    double x[] = {1, 2, 3, 2, 1.2, 2, 2.7};
    double y[] = {1, 0, 1, 2.5, 3.4, 4, 3.2};
    for (double val = 0.0; val <= 6.0; val += 0.01)
    {
        cout << "[" << Spline_2_nature(t, x, 6, val) << "," << Spline_2_nature(t, y, 6, val) << "]," << endl;
        //cout << val << "\t" << Spline_2_nature(t, x, 6, val) << "\t" << Spline_2_nature(t, y, 6, val) << endl;
    }
}
int main()
{
    ios::sync_with_stdio(false);
    cin.tie(0);
    freopen("output.txt", "w", stdout);
    //Task_4_1();
    Task_4_2();
    system("pause");
}