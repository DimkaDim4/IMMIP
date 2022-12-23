#include "LOS_LU.h"

void LOS_LU::init(size_t * s_ig, size_t * s_jg, double * s_di, double * s_gu, double * s_gl, size_t s_n)
{
    ig = s_ig;
    jg = s_jg;
    gu = s_gu;
    gl = s_gl;
    di = s_di;
    n = s_n;

    precond();
}

void LOS_LU::precond()
{
    double sum_l, sum_u, sum_d; //Промежуточные переменные, для вычисления сумм


    size_t copy_end = ig[n];

    if(!Ll) Ll = new double [copy_end];
    if(!Uu) Uu = new double [copy_end];
    if(!Ld) Ld = new double [n];
    //Копируем старые в новые
    for(size_t i = 0; i < copy_end; i++)
    {
        Ll[i] = gl[i];
        Uu[i] = gu[i];
    }

    for(size_t i = 0; i < n; i++)
        Ld[i] = di[i];

    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        sum_d = 0;

        size_t i_s = ig[k1], i_e = ig[k];

        for(size_t m = i_s; m < i_e; m++)
        {

            sum_l = 0;
            sum_u = 0;
            size_t j_s = ig[jg[m]], j_e = ig[jg[m]+1];
            for(size_t i = i_s; i < m; i++)
            {
                for(size_t j = j_s ; j < j_e; j++)
                {
                    if(jg[i] == jg[j])
                    {
                        sum_l += Ll[i]*Uu[j];
                        sum_u += Ll[j]*Uu[i];
                        j_s++;
                    }
                }
            }
            Ll[m] = Ll[m] - sum_l;
            Uu[m] = (Uu[m] - sum_u) / Ld[jg[m]];
            sum_d += Ll[m]*Uu[m];
        }
        Ld[k1] = Ld[k1] - sum_d;
    }
}

double LOS_LU::dot_prod(double *a, double *b)
{
    double dp = 0;
    for(size_t i = 0; i < n; i++)
        dp += a[i]*b[i];
    return dp;
}

void LOS_LU::mull_A(double *f, double *&x)
{
    for(size_t i = 0; i < n; i++)
    {
        double v_el = f[i];
        x[i] = di[i]*v_el;
        for(size_t k = ig[i], k1 = ig[i+1]; k < k1; k++)
        {
            size_t j = jg[k];
            x[i] += gl[k]*f[j];
            x[j] += gu[k]*v_el;
        }
    }
}

void LOS_LU::solve_L(double *f, double *&x)
{
    for(size_t k = 1, k1 = 0; k <= n; k++, k1++)
    {
        double sum = 0;

        for(size_t i = ig[k1]; i < ig[k]; i++)
            sum += Ll[i]*x[jg[i]];

        x[k1] = (f[k1] - sum)/Ld[k1];
    }
}

void LOS_LU::solve_U(double *f, double *&x)
{

    double* f1 = new double [n];
    for(size_t i = 0; i < n; i++)
        f1[i] = f[i];

    for(size_t k = n, k1 = n-1; k > 0; k--, k1--)
    {

        x[k1] = f1[k1]/Ld[k1];
        double v_el = x[k1];

        for(size_t i = ig[k1]; i < ig[k]; i++)
            f1[jg[i]] -= Uu[i]*v_el;
    }

    delete[] f1;
}

void LOS_LU::solve(double * solution, double * s_rp, double epsilon)
{
    rp = s_rp;

    //Параметры решателя
    size_t max_iter = 2000;
    double eps = epsilon;
    double end_cycle = false;

    //Норма правой части, для выхода
    double rp_norm = sqrt(dot_prod(rp, rp));

    //Начинаем решение
    double* x0 = solution;

    double* r = new double [n]; //Вектор невязки
    double* z = new double [n];
    double* p = new double [n];
    double* s = new double [n]; //Вспомогательный вектор
    double* t = new double [n]; //Вспомогательный вектор

    //r0 = L^(-1) * (f - Ax0)
    mull_A(x0, s);
    for(size_t i = 0; i < n; i++)
        s[i] = rp[i] - s[i];
    solve_L(s, r);

    //z0 = U^(-1)r0
    solve_U(r, z);

    //p0 = L^(-1)Az0
    mull_A(z, s);
    solve_L(s, p);

    size_t iter;
    double discr;

    for(iter = 0; iter < max_iter && !end_cycle; iter++)
    {
        discr = sqrt(dot_prod(r, r)); // Абсолютная невязка
        if( discr / rp_norm > eps)  //Проверка условия выхода
        {
            if(iter%10 == 0)
            {
                printf("LOS_LU Residual:\t%5lu\t%.3e\n", (unsigned long)iter, discr / rp_norm);
                fflush(stdout);
            }

            double dot1 = dot_prod(p, p); //(p[k-1], p[k-1])
            double alpha = dot_prod(p ,r) / dot1; //a = (p[k-1], r[k-1]) / (p[k-1], p[k-1])

            for(size_t i = 0; i < n; i++)
            {
                x0[i] = x0[i] + alpha*z[i]; //x[k] = x[k-1] + a*z[k-1]
                r[i] = r[i] - alpha*p[i]; //r[k] = r[k-1] - a*p[k-1]
            }
            //betta = -(p[k-1], L^(-1)*A*U^(-1)r[k]) / (p[k-1], p[k-1])

            solve_U(r, s); // s = U^(-1)r[k]
            mull_A(s, t);
            solve_L(t, t);
            double betta = - dot_prod(p, t) / dot1;

            for(size_t i = 0; i < n; i++)
            {
                z[i] = s[i] + betta * z[i]; // z[k] = U^(-1)r[k] + b*z[k-1]
                p[i] = t[i] + betta * p[i]; // p[k] = L^(-1)*A*U^(-1)r[k] + b*p[k-1]
            }

            if(iter % n == 0)  //Обновление метода
            {
                //r0 = L^(-1) * (f - Ax0)
                mull_A(x0, s);
                for(size_t i = 0; i < n; i++)
                    s[i] = rp[i] - s[i];
                solve_L(s, r);

                //z0 = U^(-1)r0
                solve_U(r, z);

                //p0 = L^(-1)Az0
                mull_A(z, s);
                solve_L(s, p);
            }
        }
        else
        {
            end_cycle = true;
        }
    }

    printf("LOS_LU Residual:\t%5lu\t%.3e\n", (unsigned long)iter, discr / rp_norm);
    if(iter >= max_iter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    fflush(stdout);

    //И отчишаем память
    delete[] p;
    delete[] r;
    delete[] z;
    delete[] s;
    delete[] t;
}

LOS_LU::LOS_LU()
{
    Uu = Ll = Ld = NULL;
}

LOS_LU::~LOS_LU()
{
    delete [] Uu;
    delete [] Ll;
    delete [] Ld;
}
