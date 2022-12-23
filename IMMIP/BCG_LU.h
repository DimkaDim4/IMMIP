#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>

class BCG_LU
{
private:
    double *Mdi,*Mggl,*Mggu;
    double *Di,*Ggl,*Ggu,*F,*Res;
    size_t *Ig,*Jg;
    size_t N,N_ggu;
    size_t maxiter,iter;
    double eps;
    double nev_r,nev_r_;

public:
    //инициализация
    void init(size_t * s_ig, size_t * s_jg, double * s_di, double * s_gu, double * s_gl, size_t s_n);
    void solve(double * solution, double * s_rp, double epsilon); //Получение решения

    size_t flag;
    BCG_LU();
    ~BCG_LU();
    void bsg_LU();

    void LU();
    void slau_U(double *x, double *f);
    void slau_Ut(double *x, double *f);
    void slau_L(double *x, double *f);
    void slau_Lt(double *x, double *f);

    void MultVM(double *x,double *y);
    void MultVM_T(double *x,double *y);
    void copy(double *x,double *y,size_t n);
    double norma(double *v);
    double scalmult(double *x,double *y);
};

