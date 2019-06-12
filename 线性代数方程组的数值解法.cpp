#include <cstring>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;
#define N 4
//全局宏定义
/*
    高斯消去法
        顺序高斯消去法：主元不为0
        列主元高斯消去法：全选主元/列选主元
        高斯-若尔当消去法  N次消元
            求非奇异矩阵的逆
*/
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
bool Gauss(double A[][N], double *B)
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
        else
            return false;
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
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << B[i] << endl;
    }
    return true;
}
//高斯-若尔当
bool G_J(double A[][N], double *B)
{
    for (int k = 0; k < N; k++) //注意要n次消元
    {
        if (Maximal_Column(A, B, k)) //列选主元成功
        {
            for (int i = k + 1; i < N; i++)
            {
                A[k][i] = A[k][i] / A[k][k]; //归一化
            }
            B[k] /= A[k][k];    //常数项归一化
            A[k][k] /= A[k][k]; //归一
            for (int i = 0; i < N; i++)
            {
                if (i != k)
                {
                    for (int j = k + 1; j < N; j++)
                    {
                        A[i][j] = A[i][j] - A[i][k] * A[k][j]; //更新A
                    }
                }
            }
            for (int i = 0; i < N; i++)
            {
                if (i != k)
                    B[i] = B[i] - A[i][k] * B[k]; //各常数项更新
            }
        }
        else
            return false;
    }
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << B[i] << endl;
    }
    return true;
}
void Disp(double A[][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}
//矩阵求逆
//等价于求系数阵相同的n个方程组
//! AB=E
bool Inverse(double A[][N])
{
    double B[N][N];
    double ans[N][N];
    memset(B, 0, sizeof(B));
    for (int i = 0; i < N; i++)
        B[i][i] = 1;
    //构造单位矩阵
    for (int i = 0; i < N; i++)
    {
        double C[N][N]; //复制临时阵
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                C[j][k] = A[j][k];
            }
        }
        if (G_J(C, B[i]))
        {
            for (int j = 0; j < N; j++)
            {
                ans[j][i] = B[i][j]; //更新至B
            }
        }
        else
            return false;
    }
    cout << "逆矩阵求解完成" << endl;
    Disp(ans);
    return true;
}
/*
    矩阵三角分解法
        矩阵的直接三角分解：LU分解
*/
//A=LU
void LU(double A[][N])
{
    for (int r = 0; r < N; r++)
    {
        for (int i = r; i < N; i++)
        {
            for (int k = 0; k <= r - 1; k++)
            {
                A[r][i] -= A[r][k] * A[k][i];
            }
        } //行计算
        for (int i = r + 1; i < N; i++)
        {
            for (int k = 0; k <= r - 1; k++)
            {
                A[i][r] -= A[i][k] * A[k][r];
            }
            A[i][r] /= A[r][r];
        } //列计算
    }
    cout << "LU分解完成" << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << 1 << endl;
    }
    cout << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j < i)
                cout << "0\t";
            else
                cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}
//利用杜利特尔分解解线性方程组Ax=b
void Doolittle(double A[][N], double *B)
{
    LU(A); //矩阵直接三角分解
    //求解Ly=b
    for (int i = 1; i < N; i++)
    {
        for (int j = 0; j <= i - 1; j++)
        {
            B[i] -= B[j] * A[i][j];
        }
    }
    //求解Ux=y
    B[N - 1] /= A[N - 1][N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            B[i] -= A[i][j] * B[j];
        }
        B[i] /= A[i][i];
    }
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << B[i] << endl;
    }
}
//AX=I
void Doolittle_Inverse(double A[][N])
{
    LU(A); //矩阵直接三角分解
    double B[N][N];
    memset(B, 0, sizeof(B));
    for (int i = 0; i < N; i++)
    {
        B[i][i] = 1;
    }
    //构造单位矩阵
    //求解Ly=b
    for (int k = 0; k < N; k++) //单位阵每一列
    {
        for (int i = 1; i < N; i++)
        {
            for (int j = 0; j <= i - 1; j++)
            {
                B[i][k] -= B[j][k] * A[i][j];
            }
        }
        //求解Ux=yj,j=0.N I的列
        B[N - 1][k] /= A[N - 1][N - 1];
        for (int i = N - 2; i >= 0; i--)
        {
            for (int j = i + 1; j < N; j++)
            {
                B[i][k] -= A[i][j] * B[j][k];
            }
            B[i][k] /= A[i][i];
        }
    }
    cout << "逆矩阵求解完成" << endl;
    Disp(B);
}
//追赶法
void Chase(double A[][N], double *B)
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
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << ans[i] << endl;
    }
}
//平方根法
void Cholesky(double A[][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k <= i - 1; k++)
        {
            A[i][i] -= A[i][k] * A[i][k];
        }
        A[i][i] = sqrt(A[i][i]);
        for (int j = i + 1; j < N; j++)
        {
            for (int k = 0; k <= i - 1; k++)
            {
                A[j][i] -= A[j][k] * A[i][k];
            }
            A[j][i] /= A[i][i];
        }
    }
    /*
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    */
}
//Ax=b
void Cholesky_Method(double A[][N], double *B)
{
    Cholesky(A);
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k <= i - 1; k++)
        {
            B[i] -= A[i][k] * B[k];
        }
        B[i] /= A[i][i];
    }
    for (int i = N - 1; i >= 0; i--)
    {
        for (int k = i + 1; k < N; k++)
        {
            B[i] -= A[k][i] * B[k];
        }
        B[i] /= A[i][i];
    }
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << B[i] << endl;
    }
}
void LDLt(double A[][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i - 1; j++)
        {
            for (int k = 0; k <= j - 1; k++)
            {
                A[i][j] -= A[i][k] * A[k][k] * A[j][k];
            }
            A[i][j] /= A[j][j];
        }
        for (int k = 0; k <= i - 1; k++)
        {
            A[i][i] -= A[k][k] * A[i][k] * A[i][k];
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (j != i)
                cout << A[i][j] << "\t";
            else
                cout << 1 << endl;
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                cout << A[i][j] << "\t";
            else
                cout << "0\t";
        }
        cout << endl;
    }
}
void LDLt_with_TLt(double A[][N])
{
    double t[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i - 1; j++)
        {
            t[i][j] = A[i][j];
            for (int k = 0; k <= j - 1; k++)
            {
                t[i][j] -= t[i][k] * A[j][k];
            }
            A[i][j] = t[i][j] / A[j][j];
        }
        for (int k = 0; k <= i - 1; k++)
        {
            A[i][i] -= t[i][k] * A[i][k];
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (j != i)
                cout << A[i][j] << " ";
            else
                cout << 1 << " ";
        }
        cout << endl;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                cout << A[i][j] << " ";
            else
                cout << " ";
        }
        cout << endl;
    }
}
void Better_Cholesky_Method(double A[][N], double *B)
{
    LDLt_with_TLt(A);
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k <= i - 1; k++)
        {
            B[i] -= A[i][k] * B[k];
        }
    }
    for (int i = 0; i < N; i++)
    {
        cout << "y" << i + 1 << " = " << B[i] << endl;
    }
    for (int i = N - 1; i >= 0; i--)
    {
        B[i] /= A[i][i];
        for (int k = i + 1; k < N; k++)
        {
            B[i] -= A[k][i] * B[k];
        }
    }
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << B[i] << endl;
    }
}
/*
    迭代法
        Jacobi迭代
        Gauss-Seidel迭代
        SOR
*/
//Ax=b
#define MAX 10000
double vector_norm_1(double *A)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += fabs(A[i]);
    }
    return sum;
}
bool Judge(double *x, double *y, double e)
{
    double max_temp = 0;
    double temp;
    for (int i = 0; i < N; i++)
    {
        temp = fabs(x[i] - y[i]);
        if (temp > max_temp)
        {
            max_temp = temp;
        }
    }
    if (max_temp < e)
        return true;
    else
        return false;
}
double *func_Jacobi(double A[][N], double *B, double *x)
{
    double *y = new double[N];
    for (int i = 0; i < N; i++)
    {
        y[i] = B[i];
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                y[i] -= A[i][j] * x[j];
            }
        }
        y[i] /= A[i][i];
    }
    return y;
}
void func_Gauss_Seidel(double A[][N], double *B, double *x)
{
    for (int i = 0; i < N; i++)
    {
        double temp = B[i];
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                temp -= A[i][j] * x[j];
            }
        }
        x[i] = temp / A[i][i];
    }
}
void func_SOR(double A[][N], double *B, double *x, double w)
{
    for (int i = 0; i < N; i++)
    {
        double temp = B[i];
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                temp -= A[i][j] * x[j];
            }
        }
        x[i] = (1.0 - w) * x[i] + w * temp / A[i][i];
    }
}
void Jacobi(double A[][N], double *B, double *x, double e)
{
    int k = 1;
    double *y = func_Jacobi(A, B, x);
    while (!Judge(x, y, e))
    {
        for (int i = 0; i < N; i++)
        {
            cout << "x" << i + 1 << " = " << y[i] << endl;
        }

        if (k >= MAX)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            k++;
            x = y;
            y = func_Jacobi(A, B, x);
        }
    }
    cout << "迭代次数：" << k << endl;
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << y[i] << endl;
    }
}
void Copy(double *a, double *b)
{
    for (int i = 0; i < N; i++)
    {
        a[i] = b[i];
    }
}
void Gauss_Seidel(double A[][N], double *B, double *x, double e)
{
    int k = 1;
    double *y = new double[N];
    Copy(y, x);
    func_Gauss_Seidel(A, B, x);
    while (!Judge(x, y, e))
    {
        for (int i = 0; i < N; i++)
        {
            cout << "x" << i + 1 << " = " << x[i] << endl;
        }
        if (k >= MAX)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            k++;
            Copy(y, x);
            func_Gauss_Seidel(A, B, x);
        }
    }
    cout << "迭代次数：" << k << endl;
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}
void SOR(double A[][N], double *B, double *x, double e, double w)
{
    int k = 1;
    double *y = new double[N];
    Copy(y, x);
    func_SOR(A, B, x, w);
    while (!Judge(x, y, e))
    {
        if (k >= MAX)
        {
            cout << "迭代失败" << endl;
            return;
        }
        else
        {
            k++;
            Copy(y, x);
            func_SOR(A, B, x, w);
        }
    }
    cout << "迭代次数：" << k << " 松弛因子：" << w << endl;
    for (int i = 0; i < N; i++)
    {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}
void Task_2_1()
{
    double A1[3][3] = {
        {2, 4, -6},
        {1, 5, 3},
        {1, 3, 2}};
    double B1[3] = {-4, 10, 5};
    double A2[3][3] = {
        {0.01, 2, -0.5},
        {-1, -0.5, 2},
        {5, -4, 0.5}};
    double B2[3] = {-5, 5, 9};
    double A3[3][3] = {
        {5, 2, 1},
        {2, 8, -3},
        {1, -3, -6}};
    double B3[3] = {8, 21, 1};
    double A4[4][4] = {
        {7, 2, 1, -2},
        {9, 15, 3, -2},
        {-2, -2, 11, 5},
        {1, 3, 2, 13}};
    double B4[4] = {4, 7, -1, 0};
    double x1[3] = {0, 0, 0};
    //SOR(A1, B1, x1, 1e-5,0.5);
    double x2[3] = {0, 0, 0};
    //SOR(A2, B2, x2, 1e-5,0.5);
    //Jacobi(A2, B2, x2, 1e-5);
    //Gauss_Seidel(A2, B2, x2, 1e-5);
    double x3[3] = {0, 0, 0};
    //SOR(A3, B3, x3, 1e-5,0.5);
    double x[4] = {0, 0, 0, 0};
    //SOR(A4, B4, x, 1e-5,0.5);
    //迭代的收敛性
}
void Task_2_2()
{
    double A[6][6] = {
        {4, -1, 0, -1, 0, 0},
        {-1, 4, -1, 0, -1, 0},
        {0, -1, 4, -1, 0, -1},
        {-1, 0, -1, 4, -1, 0},
        {0, -1, 0, -1, 4, -1},
        {0, 0, -1, 0, -1, 4}};
    double B[6] = {0, 5, 0, 6, -2, 6};
    /*
    收敛性，并求出使 || x(k + 1) – x(k) ||1 ≤0.0001 的近似解及相应的迭代次数。考虑用
    雅可比迭代法;赛德尔迭代法;超松驰迭代法(w取1.2,1.3,1.9,0.9)
    向量的1范数为各分量绝对值之和 <=1e-4
    */
    double x[6] = {0, 0, 0, 0, 0, 0};
    //SOR(A, B, x, 1e-4,0.9);
}
void matrix_mul(double A[][N], double B[][N])
{
    double ans[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ans[i][j] = 0;
            for (int k = 0; k < N; k++)
            {
                ans[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    Disp(ans);
}
void matrix_mul(double A[][N], double B[][N], double C[][N])
{
    double ans1[N][N];
    double ans2[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ans1[i][j] = 0;
            for (int k = 0; k < N; k++)
            {
                ans1[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            ans2[i][j] = 0;
            for (int k = 0; k < N; k++)
            {
                ans2[i][j] += ans1[i][k] * C[k][j];
            }
        }
    }
    Disp(ans2);
}
void Task_3_1()
{
    double A[9][9] = {
        {4, -1, 0, -1, 0, 0, 0, 0, 0},
        {-1, 4, -1, 0, -1, 0, 0, 0, 0},
        {0, -1, 4, 0, 0, -1, 0, 0, 0},
        {-1, 0, 0, 4, -1, 0, -1, 0, 0},
        {0, -1, 0, -1, 4, -1, 0, -1, 0},
        {0, 0, -1, 0, -1, 4, -1, 0, -1},
        {0, 0, 0, -1, 0, -1, 4, -1, 0},
        {0, 0, 0, 0, -1, 0, -1, 4, -1},
        {0, 0, 0, 0, 0, -1, 0, -1, 4}};
    //LDLt(A);
    //验证
    double L[9][9];
    memset(L, 0, sizeof(L));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j < i)
            {
                L[i][j] = A[i][j];
            }
            else if (j == i)
            {
                L[i][i] = 1;
            }
            else
            {
                L[i][j] = 0;
            }
        }
    }
    double D[9][9];
    memset(D, 0, sizeof(D));
    for (int i = 0; i < N; i++)
        D[i][i] = A[i][i];
    double Lt[9][9];
    memset(Lt, 0, sizeof(Lt));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (j > i)
                Lt[i][j] = A[j][i];
            else if (j == i)
                Lt[i][j] = 1;
            else
            {
                Lt[i][j] = 0;
            }
        }
    }
    //Disp(Lt);
    //matrix_mul(L, D, Lt);
}
void Task_3_2()
{
    double A[4][4] = {
        {6, 2, 1, -1},
        {2, 4, 1, 0},
        {1, 1, 4, -1},
        {-1, 0, -1, 3}};
    LU(A);
}
void Task_3_3()
{
    double A[4][4] = {
        {6, 2, 1, -1},
        {2, 4, 1, 0},
        {1, 1, 4, -1},
        {-1, 0, -1, 3}};
    Doolittle_Inverse(A);
}
int main()
{
    //freopen("output.txt", "w", stdout);
    //Task_3_3();
    /*
    double A[N][N] = {{4, 3, 0}, {3, 4, -1}, {0, -1, 4}};
    double B[N] = {24, 30, -24};
    double x[N] = {0, 0, 0};
    SOR(A, B, x, 1e-5, 1.22);
    */
    /*
    Jacobi(A, B, x, 1e-5);
    double C[N][N] = {{10, -2, -2}, {-2, 10, -1}, {-1, -2, 3}};
    double D[N] = {1, 0.5, 1};
    double y[N] = {0, 0, 0};
    Gauss_Seidel(C, D, y, 1e-5);
    */
    /*
    double A[N][N] = {
        {2, 1, 0},
        {0.5, 2, 0.5},
        {0, 1, 2}};
    double B[N] = {6, -12, 6};
    Chase(A, B);
    */
    /*
    double A[3][3] = {
        {16, 4, 8},
        {4, 5, -4},
        {8, -4, 22}};
    double B[3] = {-9, 3.25, -3};
    //Cholesky_Method(A, B);
    Better_Cholesky_Method(A, B);
    //LDLt_with_TLt(A);
    */
    
    double a[N][N] = {
        {5, 0, 10, 0},
        {0, 10, 0, 34},
        {10, 0, 34, 0},
        {0, 34, 0, 130}};
    double b[N] = {2.9, 4.2, 7, 14.4};
    Gauss(a, b);
    //Doolittle(a, b);
    /*
    double a[N][N] = {
        {1, 1, -1},
        {1, 2, -2},
        {-2, 1, 1}};
    Inverse(a);
    double b[N][N] = {
        {1, 1, -1},
        {1, 2, -2},
        {-2, 1, 1}};
    Doolittle_Inverse(b);
    */
    //Gauss(a, b);
    //G_J(a, b);
    system("pause");
}