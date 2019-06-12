#include <iostream>
#include <math.h>
#include <vector>
using namespace std;
/*
    Least squares
*/
#define N 2
//矩阵阶数
//函数空间中函数组数
double func(double x, int k)
{
    switch (k)
    {
    case 0:
        return 1;
    case 1:
        return x * x;
    }
}
//计算离散内积
double product(int k, int j, double *x, int m)
{
    double ans = 0;
    for (int i = 0; i < m; i++)
    {
        ans += func(x[i], k) * func(x[i], j);
    }
    return ans;
}
double product(int k, double *x, double *y, int m)
{
    double ans = 0;
    for (int i = 0; i < m; i++)
    {
        ans += func(x[i], k) * y[i];
    }
    return ans;
}
double product_weight(int k, int j, double *x, double *w, int m)
{
    double ans = 0;
    for (int i = 0; i < m; i++)
    {
        ans += w[i] * func(x[i], k) * func(x[i], j);
    }
    return ans;
}
double product_weight(int k, double *x, double *y, double *w, int m)
{
    double ans = 0;
    for (int i = 0; i < m; i++)
    {
        ans += w[i] * func(x[i], k) * y[i];
    }
    return ans;
}
//列选主元 多组数据成功测试
bool Maximal_Column(double A[][N], double *B, int k)
{
    //选主元
    double d = A[k][k];
    int line = k;
    for (int i = k + 1; i < N; i++)
    {
        if (fabs(A[i][k]) > fabs(d))
        {
            d = A[i][k]; //i行k列中选取
            line = i;
        }
    }
    if (fabs(d) < 1e-12) //奇异标志 列选主元失败
        return false;
    if (line == k)
        return true;
    else
    {
        //交换行
        for (int j = k; j < N; j++)
            swap(A[k][j], A[line][j]);
        swap(B[line], B[k]);
        return true;
    }
}
//AX=B 通过高斯消元法
double *Gauss(double A[][N], double *B)
{
    for (int k = 0; k < N - 1; k++) //这里只要N-1次消元
    {
        if (Maximal_Column(A, B, k)) //列选主元成功
        {
            for (int i = k + 1; i < N; i++)
            {
                A[i][k] = A[i][k] / A[k][k]; //计算系数
            }
            for (int i = k + 1; i < N; i++)
            {
                for (int j = k + 1; j < N; j++)
                {
                    A[i][j] = A[i][j] - A[i][k] * A[k][j]; //更新A
                }
            }
            for (int i = k + 1; i < N; i++)
            {
                B[i] = B[i] - A[i][k] * B[k]; //更新B
            }
        }
    }
    //回代求解
    B[N - 1] = B[N - 1] / A[N - 1][N - 1];
    //从最后一行逐步回代
    for (int i = N - 2; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            B[i] -= A[i][j] * B[j];
            //减去已解出元
        }
        B[i] /= A[i][i];
        //求得该元
    }
    return B;
}
double Least_squares(double *x, double *y, int m)
{
    double A[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = product(i, j, x, m);
        }
    }
    double B[N];
    for (int i = 0; i < N; i++)
    {
        B[i] = product(i, x, y, m);
    }
    /*
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            cout<<A[i][j]<<"\t";
        }
        cout<<endl;
    }
    for(int i=0;i<N;i++)
    {
        cout<<B[i]<<"\t";
    }
    cout<<endl;
    */
    double *ans = Gauss(A, B);
    for (int i = 0; i < N; i++)
    {
        cout << "a" << i << " = " << ans[i] << endl;
    }
}
double Least_squares_weight(double *x, double *y, double *w, int m)
{
    double A[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i][j] = product_weight(i, j, x, w, m);
        }
    }
    double B[N];
    for (int i = 0; i < N; i++)
    {
        B[i] = product_weight(i, x, y, w, m);
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
    for (int i = 0; i < N; i++)
    {
        cout << B[i] << "\t";
    }
    cout << endl;

    double *ans = Gauss(A, B);
    for (int i = 0; i < N; i++)
    {
        cout << "a" << i << " = " << ans[i] << endl;
    }
}
/*
    超定方程组的最小二乘解
    矩阵形式的求解
        矩阵的转置
        矩阵乘法
    步骤：AX=b
        A^TAx=A^Tb
    多变量、多项式的拟合在矩阵形式下只要得到Aα=y 即回到了超定方程组的求解问题上 殊途同归。
*/
struct Matrix
{
    vector<double> v;
    int r;
    int c;
    Matrix(int r, int c)
    {
        this->r = r;
        this->c = c;
        for (int i = 0; i < r * c; i++)
        {
            v.push_back(0);
        }
    }
    double &operator()(int i, int j)
    {
        return v.at(i * c + j);
    }
    void Disp()
    {
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                cout << v.at(i * c + j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    void Init(vector<double> &C)
    {
        int cnt = 0;
        for (int i = 0; i < r; i++)
        {
            for (int j = 0; j < c; j++)
            {
                v.at(i * c + j) = C[cnt++];
            }
        }
    }
};
Matrix Matrix_mul(Matrix m1, Matrix m2)
{
    Matrix ans(m1.r, m2.c);
    //矩阵乘法
    for (int i = 0; i < m1.r; i++)
    {
        for (int j = 0; j < m2.c; j++)
        {
            ans(i, j) = 0;
            for (int k = 0; k < m1.c; k++)
            {
                ans(i, j) += m1(i, k) * m2(k, j);
            }
        }
    }
    return ans;
}
vector<double> Matrix_mul(Matrix m1, vector<double> b)
{
    //返回一个向量
    vector<double> ans(m1.r);
    for (int i = 0; i < m1.r; i++)
    {
        ans[i] = 0;
        for (int j = 0; j < m1.c; j++)
        {
            ans[i] += m1(i, j) * b[j];
        }
    }
    return ans;
}
Matrix Matrix_trans(Matrix m1)
{
    Matrix ans(m1.c, m1.r);
    //矩阵的转置
    for (int i = 0; i < m1.c; i++)
    {
        for (int j = 0; j < m1.r; j++)
        {
            ans(i, j) = m1(j, i);
        }
    }
    return ans;
}
Matrix Inverse(Matrix m1)
{
}
bool Maximal_Column(Matrix &A, vector<double> &b, int k)
{
    //选主元
    double d = A(k, k);
    int line = k;
    for (int i = k + 1; i < A.r; i++)
    {
        if (fabs(A(i, k)) > fabs(d))
        {
            d = A(i, k); //i行k列中选取
            line = i;
        }
    }
    if (fabs(d) < 1e-12) //奇异标志 列选主元失败
        return false;
    if (line == k)
        return true;
    else
    {
        //交换行
        for (int j = k; j < A.r; j++)
            swap(A(k, j), A(line, j));
        swap(b[line], b[k]);
        return true;
    }
}
//AX=B 通过高斯消元法
vector<double> Gauss(Matrix A, vector<double> b)
{
    for (int k = 0; k < A.r - 1; k++) //这里只要N-1次消元
    {
        if (Maximal_Column(A, b, k)) //列选主元成功
        {
            for (int i = k + 1; i < A.r; i++)
            {
                A(i, k) = A(i, k) / A(k, k); //计算系数
            }
            for (int i = k + 1; i < A.r; i++)
            {
                for (int j = k + 1; j < A.r; j++)
                {
                    A(i, j) = A(i, j) - A(i, k) * A(k, j);
                }
            }
            for (int i = k + 1; i < A.r; i++)
            {
                b[i] = b[i] - A(i, k) * b[k]; //更新B
            }
        }
    }
    //回代求解
    //从最后一行逐步回代
    for (int i = A.r - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < A.r; j++)
        {
            b[i] -= A(i, j) * b[j];
            //减去已解出元
        }
        b[i] /= A(i, i);
        //求得该元
    }
    return b;
}
void Solve_byGauss(Matrix &A, vector<double> &x, vector<double> &b)
{
    Matrix left = Matrix_mul(Matrix_trans(A), A);
    //left.Disp();
    auto right = Matrix_mul(Matrix_trans(A), b);
    x = Gauss(left, right);
}
vector<double> Solve_byInverse(Matrix &A, vector<double> &b)
{
}
void Test()
{
    /*
    double x1[] = {19, 25, 31, 38, 44};
    double y1[] = {19.0, 32.3, 49.0, 73.0, 97.8};
    Least_squares(x1, y1, 5);
    double x2[] = {2, 4, 6, 8};
    double y2[] = {2, 11, 28, 40};
    double w2[] = {14, 27, 12, 1};
    Least_squares_weight(x2, y2, w2, 4);
    */
    /*
    Matrix matrix(5, 3);
    vector<double> init = {1, 2, 3, 1, 4, 5, 1, 5, 7, 1, 8, 9, 1, 9, 12};
    matrix.Init(init);
    vector<double> b = {48, 50, 51, 55, 56};
    vector<double> x = {0, 0, 0};
    Solve_byGauss(matrix, x, b);
    for (int i = 0; i < 3; i++)
    {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    */
    /*
    Matrix m(5, 3);
    vector<double> init = {0.1, 0.2, 1, 0.2, 0.3, 1, 0.4, 0.5, 1, 0.6, 0.7, 1, 0.9, 0.8, 1};
    m.Init(init);
    vector<double> z = {0.58, 0.63, 0.73, 0.83, 0.92};
    vector<double> x = {0, 0, 0};
    Solve_byGauss(m, x, z);
    */
    Matrix m(2, 2);
    vector<double> init = {5, 5327, 5327, 7277699};
    m.Init(init);
    vector<double> y = {271.4, 369321.5};
    vector<double> x = {0, 0};
    Solve_byGauss(m, x, y);
}
void Task_5_1_1()
{
    Matrix m(9, 3);
    vector<double> init = {
        1.0, 0.5, 2.0,
        1.0, 1.0, 4.0,
        1.0, 1.0, 5.0,
        1.0, 2.0, 2.0,
        1.0, 2.5, 4.0,
        1.0, 2.0, 5.0,
        1.0, 3.0, 2.0,
        1.0, 3.5, 4.0,
        1.0, 4.0, 5.0};
    m.Init(init);
    vector<double> x = {0, 0, 0};
    vector<double> b = {-0.19, -0.32, -1.00, 3.71, 4.49, 2.48, 6.31, 7.71, 8.51};
    Solve_byGauss(m, x, b);
    for (auto ele : x)
    {
        cout << ele << " ";
    }
    cout << endl;
    long double ans = 0;
    for (int i = 0; i < 9; i++)
    {
        long double temp = x[0] * m(i, 0) + x[1] * m(i, 1) + x[2] * m(i, 2);
        ans += powl(fabs(temp - b[i]), 2);
    }
    cout << ans << endl;
}
void Task_5_1_2()
{
    Matrix m(9, 4);
    vector<double> init = {
        1.0, 0.5, 2.0, 1.0,
        1.0, 1.0, 4.0, 4.0,
        1.0, 1.0, 5.0, 5.0,
        1.0, 2.0, 2.0, 4.0,
        1.0, 2.5, 4.0, 10.0,
        1.0, 2.0, 5.0, 10.0,
        1.0, 3.0, 2.0, 6.0,
        1.0, 3.5, 4.0, 14.0,
        1.0, 4.0, 5.0, 20.0};
    m.Init(init);
    vector<double> x = {0, 0, 0, 0};
    vector<double> b = {-0.19, -0.32, -1.00, 3.71, 4.49, 2.48, 6.31, 7.71, 8.51};
    Solve_byGauss(m, x, b);
    for (auto ele : x)
    {
        cout << ele << " ";
    }
    cout << endl;
    long double ans = 0;
    for (int i = 0; i < 9; i++)
    {
        long double temp = x[0] * m(i, 0) + x[1] * m(i, 1) + x[2] * m(i, 2) + x[3] * m(i, 3);
        ans += powl(fabs(temp - b[i]), 2);
    }
    cout << ans << endl;
}
void Task_5_2()
{
    Matrix m(6, 4);
    vector<double> init = {
        3, 5, 1, 4,
        3, -8, -2, 11,
        -1, 0, 0, -7,
        -10, -3, -2, -5,
        1, -5, -2, 5,
        8, -6, 2, -2};
    m.Init(init);
    vector<double> x = {0, 0, 0, 0};
    vector<double> b = {-7, 2, 8, 4, 10, 3};
    Solve_byGauss(m, x, b);
    for (auto ele : x)
    {
        cout << ele << " ";
    }
    cout << endl;
}
int main()
{
    //Test();
    //Task_5_1_1();
    //Task_5_1_2();
    Task_5_2();
    system("pause");
}