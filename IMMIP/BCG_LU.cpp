#include "BCG_LU.h"

//передаем указатели на нужные нам области памяти, где хранится матрица и правая часть
void BCG_LU::init(size_t * s_ig, size_t * s_jg, double * s_di, double * s_gu, double * s_gl, size_t s_n)
{
    N = s_n;
    maxiter = 10000;
    eps = 1e-20;

    Ig = s_ig;
    N_ggu = s_ig[s_n];
    Jg = s_jg;

    Ggl = s_gl;
    Ggu = s_gu;
    Di = s_di;

    flag=0;
    LU();
}

void BCG_LU::solve(double * solution, double * s_rp, double epsilon)
{
    F  = s_rp;
    eps = epsilon;
    Res  = solution;
    bsg_LU();
}

//Реализация метода
void BCG_LU::bsg_LU()
{
    double *r,*r_,*p,*p_,sc1,sc2;
    double *vec_1,*vec_2,*vec_3;
    double alpha,betta,old_nev_r=1.e+20;
    size_t i;
    r=new double [N];
    p=new double [N];
    r_=new double [N];
    p_=new double [N];
    vec_1=new double [N];
    vec_2=new double [N];
    vec_3=new double [N];

    iter=0;
    nev_r=old_nev_r/10;
    nev_r_=0;

    MultVM(Res,vec_1);
    for(i=0; i<N; i++) vec_1[i]=F[i]-vec_1[i];
    slau_L(r,vec_1);

    copy(r_,r,N);
    copy(p,r,N);
    copy(p_,r_,N);

    //Цикл алгоритма
    while(flag==0 && iter<maxiter)
    {
        sc1=scalmult(r,r_);

        slau_U(vec_1,p);

        MultVM(vec_1,vec_2); //input->output

        slau_L(vec_3,vec_2);

        sc2=scalmult(p_,vec_3);

        alpha=sc1/sc2;

        for(i=0; i<N; i++)	Res[i]=Res[i]+alpha*vec_1[i];

        for(i=0; i<N; i++)	r[i]=r[i]-alpha*vec_3[i];

        slau_Lt(vec_1,p_);

        MultVM_T(vec_1,vec_2);//input->output

        slau_Ut(vec_3,vec_2);

        for(i=0; i<N; i++)	r_[i]=r_[i]-alpha*vec_3[i];

        sc2=scalmult(r,r_);

        betta=sc2/sc1;

        old_nev_r=nev_r;
        nev_r=norma(r);

        if( nev_r < eps) flag=1;

        if(flag==0)
        {
            for(i=0; i<N; i++)	p[i]=r[i]+betta*p[i];

            for(i=0; i<N; i++)	p_[i]=r_[i]+betta*p_[i];

            iter++;
        }
/*
        if(iter%10 == 0)
        {
            printf("BCG_LU Residual:\t%5lu\t%.3e\n", (unsigned long)iter, nev_r);
            fflush(stdout);
        }
*/
    }
/*
    printf ("\nИтераций: %ld	Невязка:%.6e",iter,nev_r);
    if(flag==0) printf ("\nВыход по колебанию невязки\n");
    if(flag==1) printf ("\nРешение достигнуто\n");

    if(flag==2) printf ("\nРешение не достигнуто!\n- Измените вектор R_\n");
*/
    printf("BCG_LU Residual:\t%5lu\t%.3e\n", (unsigned long)iter, nev_r);
    if(iter >= maxiter)
        printf("Soulution can`t found, iteration limit exceeded!\n");
    fflush(stdout);

    delete[] r;
    delete[] p;
    delete[] r_;
    delete[] p_;
    delete[] vec_1;
    delete[] vec_2;
    delete[] vec_3;
}

void BCG_LU::copy(double *x,double *y,size_t n)
{
    for(size_t i=0; i<n; i++) x[i]=y[i];
}

double BCG_LU::scalmult(double *x,double *y)
{
    double s=0;
    for(size_t i=0; i<N; i++) s+=x[i]*y[i];
    return(s);
}

double BCG_LU::norma(double *v)
{
    return(sqrt(scalmult(v,v)));
}

//LU-разложение
void BCG_LU::LU()
{
    size_t a,b,i,j;
    if(!Mggl) Mggl=new double[N_ggu];
    if(!Mggu) Mggu=new double[N_ggu];
    if(!Mdi) Mdi=new double[N];
    copy(Mggl,Ggl,N_ggu);
    copy(Mggu,Ggu,N_ggu);
    copy(Mdi,Di,N);

    Mdi[0]=Di[0];
    for(i=1; i<N; i++)
    {
        Mdi[i]=Di[i];

        for(j=Ig[i]; j<Ig[i+1]; j++)

        {
            Mggl[j]=Ggl[j];
            Mggu[j]=Ggu[j];

            for(a=Ig[i]; a<j; a++)
            {
                for(b=Ig[Jg[j]]; b<=Ig[Jg[j]+1]-1; b++)
                {
                    if(Jg[a]==Jg[b])
                    {
                        Mggl[j]-=Mggl[a]*Mggu[b];
                        Mggu[j]-=Mggu[a]*Mggl[b];
                    }
                }
            }
            Mggu[j]/=Mdi[Jg[j]];
            Mdi[i]-=Mggu[j]*Mggl[j];
        }
    }
}

//решение СЛАУ с верхней треугольной матрицей
void BCG_LU::slau_U(double *x, double *g)
{
    size_t i,j,k;
    double *f;
    f=new double[N];
    copy(f,g,N);
    for(k=N; k>0; k--)
    {
        i = k - 1;
        x[i]=f[i];
        for(j=Ig[i]; j<Ig[i+1]; j++)
            f[Jg[j]]-=x[i]*Mggu[j];
    }
    delete[] f;
}

//решение СЛАУ с верхней треугольной транспонированной матрицей
void BCG_LU::slau_Ut(double *x, double *f)
{
    size_t i,j;
    for(i=0; i<N; i++)
    {
        x[i]=f[i];
        for(j=Ig[i]; j<Ig[i+1]; j++)
            x[i]-=x[Jg[j]]*Mggu[j];
    }
}

//решение СЛАУ с нижней треугольной матрицей
void BCG_LU::slau_L(double *x, double *f)
{
    size_t i,j;
    for(i=0; i<N; i++)
    {
        x[i]=f[i];
        for(j=Ig[i]; j<Ig[i+1]; j++)
            x[i]-=x[Jg[j]]*Mggl[j];
        x[i]/=Mdi[i];
    }
}

//решение СЛАУ с нижней треугольной транспонированной матрицей
void BCG_LU::slau_Lt(double *x, double *g)
{
    size_t i,j,k;
    double *f;
    f=new double[N];
    copy(f,g,N);
    for(k = N; k > 0; k--)
    {
        i = k - 1;
        x[i]=f[i]/Mdi[i];
        for(j=Ig[i]; j<Ig[i+1]; j++)
            f[Jg[j]]-=x[i]*Mggl[j];
    }
    delete[] f;
}

//умножение матрицы на вектор
void BCG_LU::MultVM(double *x,double *y)
{
    size_t i,j;
    for(i=0; i<N; i++)
        y[i]=x[i]*Di[i];
    for(i=0; i<N; i++)
        for(j=Ig[i]; j<Ig[i+1]; j++)

        {
            y[i]+=x[Jg[j]]*Ggl[j];
            y[Jg[j]]+=x[i]*Ggu[j];
        }

}

//умножение транспонированной матрицы на вектор
void BCG_LU::MultVM_T(double *x,double *y)
{
    size_t i,j;
    for(i=0; i<N; i++)
        y[i]=x[i]*Di[i];
    for(i=0; i<N; i++)
        for(j=Ig[i]; j<Ig[i+1]; j++)
        {
            y[i]+=x[Jg[j]]*Ggu[j];
            y[Jg[j]]+=x[i]*Ggl[j];
        }

}

BCG_LU::BCG_LU()
{
    Mdi = Mggl = Mggu = NULL;
    Di = Ggl = Ggu = F = Res = NULL;
    Ig = Jg = NULL;
}

BCG_LU::~BCG_LU()
{
    delete [] Mdi;
    delete [] Mggl;
    delete [] Mggu;
}
