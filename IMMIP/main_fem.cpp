#include "fem.h"

size_t TIMES_NUM = 201;
double TIMES_STEP = 1.0;//0.1;

double func_rp(const point & p, const quadrilateral * r, double t)
{
/*
    // -1/r d(r*lambda * du/dr)/dr -d(lambda * du/dz)/dz
    // 2x-y
    // 2x^2-y*x
    // 4x-y
    // 4-y/x
    return (-6.0 + p.y / p.x) * r->ph->lambda;
*/
    (void)p;
    (void)r;
    (void)t;
    return 0.0;
}

size_t get_type_of_bound(size_t gmsh_num, double t)
{
    (void)gmsh_num;
    (void)t;
    if(gmsh_num == 21) return 3;
    if(gmsh_num == 22) return 3;
    if(gmsh_num == 23) return 2;
    if(gmsh_num == 24) return 2;
    assert(0);
    return 1;
}

double get_bound_value(const point & p, const edge & e, double t)
{
/*
    return p.x * p.x + p.y * p.y - p.x * p.y;
*/
    (void)p;
    (void)e;
    (void)t;
    assert(0);
    return 0.0;
}

double get_theta(const point & p, const edge & e, double t)
{
/*
    if(e.gmsh_phys_num == 23) return 0.0;
    double l = e.fes[0]->ph->lambda;
    double du = 2.0 * p.y - p.x;
    //double du = 2.0 * p.x - p.y;
    return l * du;
*/
    (void)p;
    (void)e;
    (void)t;
    if(e.gmsh_phys_num == 23) return 0.0;

    double r1 = 0.0;
    double r2 = 0.08;
    double P = 1000;
    double S = M_PI * (r2 * r2 - r1 * r1);

    if(p.x > r2) return 0.0;
    if(p.x < r1) return 0.0;
    if(e.gmsh_phys_num == 24) return P / S;
    assert(0);
    return 0.0;
}

double get_u_beta(const point & p, const edge & e, double t)
{
/*
    double u = get_bound_value(p, e);
    double b = get_beta(e);
    double l = e.fes[0]->ph->lambda;
    double du = 2.0 * p.y - p.x;
    return u + du * l / b;
*/
    (void)p;
    (void)e;
    (void)t;
    return 20.0;
}

double get_beta(const edge & e, double t)
{
    (void)e;
    (void)t;
    return 50.0;
}

point get_V(const point & p, const quadrilateral * r, double t)
{
    return point(0,0);
    /*
      18
    19  17
      20
    */
    (void)p;
    (void)r;
    (void)t;
    double R = 0.08;
    double H = 0.09;
    double Vmax = 5e-2;
    double rz_translate = p.x / R;
    Vmax /= rz_translate;
    if(r->ph->gmsh_phys_num == 17)
    {
        //double d = (p.x - R / 2.0) / (R / 2.0);
        //return point(0, Vmax * d);
        //cout << Vmax * (p.x - R / 2.0) * (p.x - 3.0 * R / 2.0) << endl;
        return point(0, - Vmax * (p.x - R / 2.0) * (p.x - 3.0 * R / 2.0));
    }
    if(r->ph->gmsh_phys_num == 18)
    {
        //double d = (p.y - H / 2.0) / (H / 2.0);
        //return point(- Vmax * d, 0);
        return point(Vmax * (p.y - H / 2.0) * (p.y - 3.0 * H / 2.0), 0);
    }
    if(r->ph->gmsh_phys_num == 19)
    {
        //double d = (R / 2.0 - p.x) / (R / 2.0);
        //return point(0, - Vmax * d);
        //cout << Vmax * (p.x - R / 2.0) * (p.x - 0.0) << endl;
        return point(0, Vmax * (p.x - R / 2.0) * (p.x - 0.0));
    }
    if(r->ph->gmsh_phys_num == 20)
    {
        //double d = (H / 2.0 - p.y) / (H / 2.0);
        //return point(Vmax * d, 0);
        return point(- Vmax * (p.y - H / 2.0) * (p.y - 0.0), 0);
    }
    assert(0);
    return point();
}
