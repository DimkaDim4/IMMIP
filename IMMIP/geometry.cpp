#include "fem.h"

// ============================================================================

point::point()
{
    x = y = 0.0;
    num = 0;
}

point::point(double x, double y)
{
    this->x = x;
    this->y = y;
    num = 0;
}

point::point(double x, double y, size_t num)
{
    this->x = x;
    this->y = y;
    this->num = num;
}

point::point(const point & p)
{
    x = p.x;
    y = p.y;
    num = p.num;
}

double & point::operator [] (size_t i)
{
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    default:
        stringstream str;
        str << "Error: Unknown index " << i << " at point " << * this << endl;
        throw str.str();
    }
}

double point::operator [] (size_t i) const
{
    switch(i)
    {
    case 0:
        return x;
    case 1:
        return y;
    default:
        stringstream str;
        str << "Error: Unknown index " << i << " at point " << * this << endl;
        throw str.str();
    }
}

bool point::operator < (const point & t) const
{
    return num < t.num;
}

bool point::operator == (const point & t) const
{
    return num == t.num;
}

point & point::operator = (const point & other)
{
    if(this != & other)
    {
        this->x = other.x;
        this->y = other.y;
        this->num = other.num;
    }
    return * this;
}

ostream & operator << (ostream & os, point a)
{
    os << "{ " << a.x << ", " << a.y << " }";
    return os;
}

// ============================================================================

point quadrilateral::to_local(const point & p) const
{
    // Кирпич, страница 305

    double w = beta6 * (p.x - nodes[0]->x) - beta5 * (p.y - nodes[0]->y);

    double ksi = 0.0;
    double eta = 0.0;

    double eps = 1e-10;

    if(fabs(alpha1) < eps && fabs(alpha2) < eps) {
        ksi = (beta3 * (p.x - nodes[0]->x) - beta1 * (p.y - nodes[0]->y))/(beta2 * beta3 - beta1 * beta4);
        eta = (beta2 * (p.y - nodes[0]->y) - beta4 * (p.x - nodes[0]->x))/(beta2 * beta3 - beta1 * beta4);
    }
    else {
        if(fabs(alpha1) < eps) {
            ksi = (alpha2 * (p.x - nodes[0]->x) + beta1 * w)/(alpha2 * beta2 - beta5 * w);
            eta = -w/alpha2;
        }
        else {
            if(fabs(alpha2) < eps) {
                ksi = w / alpha1;
                eta = (alpha1 * (p.y - nodes[0]->y) - beta4 * w)/(alpha1 * beta3 + beta6 * w);
            }
            else {
                double a = beta5 * alpha2;
                double b = alpha2 * beta2 + alpha1 * beta1 + beta5 * w;
                double c = alpha1 * (nodes[0]->x - p.x) + beta2 * w;
                double D = b * b - 4.0 * a * c;
                double eta1 = (-b - sqrt(D)) / (2.0 * a);
                double eta2 = (-b + sqrt(D)) / (2.0 * a);
                double ksi1 = alpha2 / alpha1 * eta1 + w / alpha1;
                double ksi2 = alpha2 / alpha1 * eta2 + w / alpha1;
                if(eta1 + eps >= 0.0 && eta1 - eps <= 1.0 && ksi1 + eps >= 0.0 && ksi1 - eps <= 1.0) {
                    ksi = ksi1;
                    eta = eta1;
                }
                else {
                    if(eta2 + eps >= 0.0 && eta2 - eps <= 1.0 && ksi2 + eps >= 0.0 && ksi2 - eps <= 1.0) {
                        ksi = ksi2;
                        eta = eta2;
                    }
                    else
                    {
                        // Так как порой бывают случаи, когда значения бф лежат
                        // вне обычного интервала [0, 1], то тупое сравнение тут
                        // работать не будет. Введем типа функцию штрафа, которая
                        // будет квадратом расстояния отклонения бф от ожидаемого
                        // значения, затем выберем функцию с минимальным штрафом.

                        double penalty1 = 0.0;
                        if(eta1 + eps < 0.0)
                            penalty1 += eta1 * eta1;
                        if(ksi1 + eps < 0.0)
                            penalty1 += ksi1 * ksi1;
                        if(eta1 - eps > 1.0)
                            penalty1 += (eta1 - 1.0) * (eta1 - 1.0);
                        if(ksi1 - eps > 1.0)
                            penalty1 += (ksi1 - 1.0) * (ksi1 - 1.0);

                        double penalty2 = 0.0;
                        if(eta2 + eps < 0.0)
                            penalty2 += eta2 * eta2;
                        if(ksi2 + eps < 0.0)
                            penalty2 += ksi2 * ksi2;
                        if(eta2 - eps > 1.0)
                            penalty2 += (eta2 - 1.0) * (eta2 - 1.0);
                        if(ksi2 - eps > 1.0)
                            penalty2 += (ksi2 - 1.0) * (ksi2 - 1.0);

                        if(penalty1 < penalty2)
                        {
                            cerr << "Warning: Target point is outside of element, penalty = " << penalty1 << endl;
                            ksi = ksi1;
                            eta = eta1;
                        }
                        else
                        {
                            cerr << "Warning: Target point is outside of element, penalty = " << penalty2 << endl;
                            ksi = ksi2;
                            eta = eta2;
                        }
                    }
                }
            }
        }
    }

    return point(ksi, eta);
}

point quadrilateral::to_global(const point & p) const
{
    // Кирпич, страница 298, формулы 5.109-5.110
    double x = (1 - p.x) * (1 - p.y) * nodes[0]->x + p.x * (1 - p.y) * nodes[1]->x + (1 - p.x) * p.y * nodes[2]->x + p.x * p.y * nodes[3]->x;
    double y = (1 - p.x) * (1 - p.y) * nodes[0]->y + p.x * (1 - p.y) * nodes[1]->y + (1 - p.x) * p.y * nodes[2]->y + p.x * p.y * nodes[3]->y;
    return point(x, y);
}

double quadrilateral::det_J_local(const point & p) const
{
    // Кирпич, страница 300, формула после 5.116
    double jacobian = alpha0 + alpha1 * p.x + alpha2 * p.y;
    return jacobian;
}

void quadrilateral::init()
{
    // Кирпич, страница 301, формула 5.117
    alpha0 = (nodes[1]->x - nodes[0]->x) * (nodes[2]->y - nodes[0]->y) - (nodes[1]->y - nodes[0]->y) * (nodes[2]->x - nodes[0]->x);
    alpha1 = (nodes[1]->x - nodes[0]->x) * (nodes[3]->y - nodes[2]->y) - (nodes[1]->y - nodes[0]->y) * (nodes[3]->x - nodes[2]->x);
    alpha2 = (nodes[2]->y - nodes[0]->y) * (nodes[3]->x - nodes[1]->x) - (nodes[2]->x - nodes[0]->x) * (nodes[3]->y - nodes[1]->y);
    // Кирпич, страница 302, формула 5.118
    beta1 = nodes[2]->x - nodes[0]->x;
    beta2 = nodes[1]->x - nodes[0]->x;
    beta3 = nodes[2]->y - nodes[0]->y;
    beta4 = nodes[1]->y - nodes[0]->y;
    beta5 = nodes[0]->x - nodes[1]->x - nodes[2]->x + nodes[3]->x;
    beta6 = nodes[0]->y - nodes[1]->y - nodes[2]->y + nodes[3]->y;

    // Площадь через площадь двух треугольников
    double S1 = 0.5 * fabs((nodes[0]->x - nodes[2]->x) * (nodes[1]->y - nodes[2]->y) - (nodes[1]->x - nodes[2]->x) * (nodes[0]->y - nodes[2]->y));
    double S2 = 0.5 * fabs((nodes[1]->x - nodes[3]->x) * (nodes[2]->y - nodes[3]->y) - (nodes[2]->x - nodes[3]->x) * (nodes[1]->y - nodes[3]->y));
    S = S1 + S2;

    // https://ru.wikipedia.org/wiki/Список_квадратурных_формул
    static double g_a = sqrt((114.0 - 3.0 * sqrt(583.0)) / 287.0);
    static double g_b = sqrt((114.0 + 3.0 * sqrt(583.0)) / 287.0);
    static double g_c = sqrt(6.0 / 7.0);
    static double g_wa = 307.0 / 810.0 + 923.0 / (270.0 * sqrt(583.0));
    static double g_wb = 307.0 / 810.0 - 923.0 / (270.0 * sqrt(583.0));
    static double g_wc = 98.0 / 405.0;

    // Веса Гаусса
    for(size_t i = 0; i < 4; i++)
    {
        gauss_weights[i] = g_wc;
        gauss_weights[i+4] = g_wa;
        gauss_weights[i+8] = g_wb;
    }

    // Точки Гаусса в локальной с.к.
    double gauss_points_local[2][12] =
    {
        {-g_c, g_c, 0.0, 0.0, -g_a, g_a, -g_a, g_a, -g_b, g_b, -g_b, g_b},
        {0.0, 0.0, -g_c, g_c, -g_a, -g_a, g_a, g_a, -g_b, -g_b, g_b, g_b}
    };

    // Перевод с Гауссова мастер-элемента на базисный
    for(int i = 0; i < 12; i++)
    {
        // Базисный мастер-элемент у нас [0,1]x[0,1], Гауссов - [-1,1]x[-1,1]
        gauss_points[i] = point((gauss_points_local[0][i] + 1.0) / 2.0, (gauss_points_local[1][i] + 1.0) / 2.0);
    }
}

// Получние номеров одномерных БФ из номера двумерной
void quadrilateral::two2one(size_t two, size_t & one1, size_t & one2) const
{
    one1 = two % 3;
    one2 = two / 3;
}

// Одномерные лагранжевы базисные функции
double quadrilateral::bfunc_1d(size_t func_n, double ksi) const
{
    switch(func_n)
    {
    case 0:
        return 2.0 * (ksi - 0.5) * (ksi - 1.0);
    case 1:
        return - 4.0 * ksi * (ksi - 1.0);
    case 2:
        return 2.0 * ksi * (ksi - 0.5);
    };
    assert(func_n < 3);
    return 0.0;
}

// Производные одномерных лагранжевых базисных функций
double quadrilateral::dbfunc_1d(size_t func_n, double ksi) const
{
    switch(func_n)
    {
    case 0:
        return 4.0 * (ksi - 0.75);
    case 1:
        return 4.0 - 8.0 * ksi;
    case 2:
        return 4.0 * (ksi - 0.25);
    };
    assert(func_n < 3);
    return 0.0;
}

// Двумерные лагранжевы базисные функции
double quadrilateral::bfunc_2d_local(size_t num, const point & p) const
{
    size_t nx = 0, ny = 0;
    two2one(num, nx, ny);
    return bfunc_1d(nx, p.x) * bfunc_1d(ny, p.y);
}

double quadrilateral::bfunc_2d(size_t num, const point & p) const
{
    return bfunc_2d_local(num, to_local(p));
}

// Градиент двумерных лагранжевых базисных функций
point quadrilateral::grad_bfunc_2d_local(size_t num, const point & p) const
{
    size_t nx = 0, ny = 0;
    two2one(num, nx, ny);
    return point(dbfunc_1d(nx, p.x) * bfunc_1d(ny, p.y),
                 bfunc_1d(nx, p.x) * dbfunc_1d(ny, p.y));
}

point quadrilateral::grad_bfunc_2d(size_t num, const point & p) const
{
    return grad_bfunc_2d_local(num, to_local(p));
}

// Получение элемента матрицы жесткости
double quadrilateral::get_local_G(size_t i, size_t j) const
{
    double sign_alpha0 = 0.0;
    if(alpha0 > 0.0) sign_alpha0 = 1.0;
    if(alpha0 < 0.0) sign_alpha0 = -1.0;

    // Кирпич, страницы 302-303, формула 5.119
    double result = 0.0;
    for(int g = 0; g < 12; g++)
    {
        point grad_i = grad_bfunc_2d_local(i, gauss_points[g]);
        point grad_j = grad_bfunc_2d_local(j, gauss_points[g]);

        double J = det_J_local(gauss_points[g]);
        double ksi = gauss_points[g].x;
        double eta = gauss_points[g].y;
        double dfii_dksi = grad_i.x;
        double dfii_deta = grad_i.y;
        double dfij_dksi = grad_j.x;
        double dfij_deta = grad_j.y;

        double func = (dfii_dksi * (beta6 * ksi + beta3) - dfii_deta * (beta6 * eta + beta4)) *
                      (dfij_dksi * (beta6 * ksi + beta3) - dfij_deta * (beta6 * eta + beta4)) +
                      (dfii_deta * (beta5 * eta + beta2) - dfii_dksi * (beta5 * ksi + beta1)) *
                      (dfij_deta * (beta5 * eta + beta2) - dfij_dksi * (beta5 * ksi + beta1));
        func /= J;
        result += gauss_weights[g] * func * to_global(gauss_points[g]).x; // x это типа r
    }
    return sign_alpha0 * result / 4.0;
}

// Получение элемента матрицы массы
double quadrilateral::get_local_M(size_t i, size_t j) const
{
    double sign_alpha0 = 0.0;
    if(alpha0 > 0.0) sign_alpha0 = 1.0;
    if(alpha0 < 0.0) sign_alpha0 = -1.0;

    double result = 0.0;
    for(int g = 0; g < 12; g++)
    {
        result += gauss_weights[g] * det_J_local(gauss_points[g]) *
                  bfunc_2d_local(i, gauss_points[g]) * bfunc_2d_local(j, gauss_points[g])
                  * to_global(gauss_points[g]).x; // x это типа r
    }
    return sign_alpha0 * result / 4.0;
}

// Получение элемента матрицы массы
double quadrilateral::get_local_C(size_t i, size_t j, double t) const
{
    double sign_alpha0 = 0.0;
    if(alpha0 > 0.0) sign_alpha0 = 1.0;
    if(alpha0 < 0.0) sign_alpha0 = -1.0;

    // По аналогии с G
    double result = 0.0;
    for(int g = 0; g < 12; g++)
    {
        point point_glob = to_global(gauss_points[g]);
        point V = get_V(point_glob, this, t);
        point grad_j = grad_bfunc_2d_local(j, gauss_points[g]);

        double ksi = gauss_points[g].x;
        double eta = gauss_points[g].y;
        double dfij_dksi = grad_j.x;
        double dfij_deta = grad_j.y;

        double func = (V.x * (dfij_dksi * (beta6 * ksi + beta3) - dfij_deta * (beta6 * eta + beta4)) +
                       V.y * (dfij_deta * (beta5 * eta + beta2) - dfij_dksi * (beta5 * ksi + beta1))) *
                      bfunc_2d_local(i, gauss_points[g]);
        result += gauss_weights[g] * func * point_glob.x; // x это типа r
    }
    return sign_alpha0 * result / 4.0;
}

// Получение элемента правой части
double quadrilateral::get_local_rp(size_t i, double t) const
{
    double sign_alpha0 = 0.0;
    if(alpha0 > 0.0) sign_alpha0 = 1.0;
    if(alpha0 < 0.0) sign_alpha0 = -1.0;

    double result = 0.0;
    for(int g = 0; g < 12; g++)
    {
        result += gauss_weights[g] * det_J_local(gauss_points[g]) *
                  bfunc_2d_local(i, gauss_points[g]) * func_rp(to_global(gauss_points[g]), this, t)
                  * to_global(gauss_points[g]).x; // x это типа r
    }
    return sign_alpha0 * result / 4.0;
}

// Определение, внутри четырехугольника точка или нет
bool quadrilateral::inside(const point & p) const
{
    double eps = 1e-10;
    double S1 = 0.5 * fabs((nodes[0]->x - p.x) * (nodes[1]->y - p.y) - (nodes[1]->x - p.x) * (nodes[0]->y - p.y));
    double S2 = 0.5 * fabs((nodes[0]->x - p.x) * (nodes[2]->y - p.y) - (nodes[2]->x - p.x) * (nodes[0]->y - p.y));
    double S3 = 0.5 * fabs((nodes[1]->x - p.x) * (nodes[3]->y - p.y) - (nodes[3]->x - p.x) * (nodes[1]->y - p.y));
    double S4 = 0.5 * fabs((nodes[3]->x - p.x) * (nodes[2]->y - p.y) - (nodes[2]->x - p.x) * (nodes[3]->y - p.y));
    if(fabs(S - S1 - S2 - S3 - S4) < eps)
        return true;
    return false;
}

// ============================================================================

bool edge::operator < (const edge & t) const
{
    if(nodes[0]->num < t.nodes[0]->num) return true;
    if(nodes[0]->num > t.nodes[0]->num) return false;
    if(nodes[1]->num < t.nodes[1]->num) return true;
    return false;
}

bool edge::operator == (const edge &t) const
{
    if(nodes[0]->num == t.nodes[0]->num &&
       nodes[1]->num == t.nodes[1]->num) return true;
    return false;
}

edge::edge()
{
    nodes[0] = nodes[1] = NULL;
    gmsh_phys_num = 0;
    fes[0] = fes[1] = NULL;
}

edge::edge(const point * begin, const point * end)
{
    nodes[0] = begin;
    nodes[1] = end;
    gmsh_phys_num = 0;
    fes[0] = fes[1] = NULL;
}

point edge::get_freedom_position(size_t num) const
{
    if(num == 0) return point(*(nodes[0]));
    if(num == 2) return point(*(nodes[1]));
    return point((nodes[0]->x + nodes[1]->x) / 2.0,
            (nodes[0]->y + nodes[1]->y) / 2.0);
}

double edge::bfunc_1d(size_t func_n, double x) const
{
    switch(func_n)
    {
    case 0:
        return 2.0 * (x - 0.5) * (x - 1.0);
    case 1:
        return - 4.0 * x * (x - 1.0);
    case 2:
        return 2.0 * x * (x - 0.5);
    };
    assert(func_n < 3);
    return 0.0;
}

double edge::get_matrix_M(size_t i, size_t j)
{
    // https://en.wikipedia.org/wiki/Gaussian_quadrature
    static double gauss_points[5] =
    {
        0.0,
        1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)),
        - 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)),
        1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)),
        - 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0))
    };
    static double gauss_weights[5] =
    {
        128.0 / 225.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 + 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0,
        (322.0 - 13.0 * sqrt(70.0)) / 900.0
    };
    double dx = nodes[1]->x - nodes[0]->x;
    double dy = nodes[1]->y - nodes[0]->y;
    double h = sqrt(dx * dx + dy * dy);
    double jacobian = h / 2.0;

    bool need_r = false;
    const double ORIENT_EPS = 1e-10;
    assert(!(fabs(dx / h) > ORIENT_EPS && fabs(dy / h) > ORIENT_EPS));
    if(fabs(dx / h) > ORIENT_EPS) need_r = true;
    double r_min = nodes[0]->x;
    if(r_min > nodes[1]->x) r_min = nodes[1]->x;

    double result = 0.0;
    for(int g = 0; g < 5; g++)
    {
        double gauss_point_basis = (gauss_points[g] + 1.0) / 2.0;
        double add = gauss_weights[g] *
                bfunc_1d(i, gauss_point_basis) * bfunc_1d(j, gauss_point_basis);
        if(need_r) add *= r_min + gauss_point_basis * h;
        else       add *= r_min;
        result += add;
    }
    return result * jacobian;
}
