#ifndef LOS_H_INCLUDED
#define LOS_H_INCLUDED

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>
#include <cstdlib>

//Решатель ЛОС, для СЛАУ с несимметричной матрицей в разреженном формате c LU - предобуславливанием

class LOS_LU
{
public:
    void init(size_t * s_ig, size_t * s_jg, double * s_di, double * s_gu, double * s_gl, size_t s_n); //инициализация
    void solve(double * solution, double * s_rp, double epsilon); //Получение решения

    LOS_LU();
    ~LOS_LU();
private:
    size_t n; //Размерность СЛАУ

    //Основные массивы
    size_t *ig, *jg;
    double *gu, *gl, *di;
    double *rp;

    //Массивы для преобславливателя
    double *Uu, *Ll, *Ld;

    void precond(); //Вычисление матриц L и U
    double dot_prod(double *a, double *b); //скалярное произведение
    void mull_A(double *f, double *&x); // x = Af
    void solve_L(double *f, double *&x); //Lx = f, прямой ход
    void solve_U(double *f, double *&x); //Ux = f, обратный ход
};

#endif // LOS_H_INCLUDED
