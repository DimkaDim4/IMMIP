#include "fem.h"

// ============================================================================

SLAE::SLAE()
{
    gu = gl = di = rp = x = NULL;
    ig = jg = NULL;
    n = 0;
}

SLAE::~SLAE()
{
    if(gu) delete [] gu;
    if(gl) delete [] gl;
    if(di) delete [] di;
    if(rp) delete [] rp;
    if(x)  delete [] x;
    if(ig) delete [] ig;
    if(jg) delete [] jg;
}

void SLAE::solve(double eps)
{
    cout << "Solving SLAE ..." << endl;
    solver.init(ig, jg, di, gu, gl, n);
    solver.solve(x, rp, eps);
}

void SLAE::alloc_all(size_t n_size, size_t gg_size)
{
    n = n_size;
    ig = new size_t [n + 1];
    jg = new size_t [gg_size];
    di = new double [n];
    gl = new double [gg_size];
    gu = new double [gg_size];
    rp = new double [n];
    x  = new double [n];

    memset(di, 0, sizeof(double) * n);
    memset(rp, 0, sizeof(double) * n);
    memset(x,  0, sizeof(double) * n);
    memset(gl, 0, sizeof(double) * gg_size);
    memset(gu, 0, sizeof(double) * gg_size);
}

void SLAE::add(size_t i, size_t j, double elem)
{
    if(i == j)
    {
        di[i] += elem;
        return;
    }

    size_t i_p = i, j_p = j;
    if(j_p > i_p)
        swap(i_p, j_p);
    size_t ind = 0;
    bool flag = false;
    for(size_t k = ig[i_p]; k < ig[i_p + 1] && !flag; k++)
    {
        if(jg[k] == j_p)
        {
            ind = k;
            flag = true;
        }
    }

    if(j < i)
        gl[ind] += elem;
    else
        gu[ind] += elem;
}

// ============================================================================

FEM::FEM()
{
    nodes = NULL;
    qls = NULL;
    phs = NULL;
    bounds = NULL;
    timed = NULL;

    nodes_num = 0;
    qls_num = 0;
    phs_num = 0;
    bounds_num = 0;
    timed_num = 0;
}

FEM::~FEM()
{
    if(nodes) delete [] nodes;
    if(qls) delete [] qls;
    if(phs) delete [] phs;
    if(bounds) delete [] bounds;
    if(timed)
    {
        for(size_t i = 0; i < timed_num; i++)
            delete [] timed[i];
        delete [] timed;
    }
}

void FEM::add_edge_freedom(const point * begin, const point * end, size_t ql_num)
{
    edge ed(begin, end);
    edge ei(end, begin);
    vector<edge>::iterator id = find(edges_freedom.begin(), edges_freedom.end(), ed);
    vector<edge>::iterator ii = find(edges_freedom.begin(), edges_freedom.end(), ei);
    if(id == edges_freedom.end() && ii == edges_freedom.end())
    {
        edges_freedom.push_back(ed);
        edges_freedom[edges_freedom.size()-1].fes[0] = qls + ql_num;
    }
    else
    {
        vector<edge>::iterator it;
        edge et;
        if(ii == edges_freedom.end())
        {
            it = id;
            et = ed;
        }
        else
        {
            it = ii;
            et = ei;
        }
        if(!it->fes[0])
            it->fes[0] = qls + ql_num;
        else if(!it->fes[1])
            it->fes[1] = qls + ql_num;
        else
            assert(!it->fes[0] || !it->fes[1]);
    }
}

void FEM::input()
{
    // У gmsh обход узлов обычно по часовой стрелке, но бывает и против часовой стрелки
    // А нам надо снизу-вверх-слева-направо, поэтому сделаем массивы конвертации
    size_t gmsh_convert_clock[4] = {2, 3, 1, 0};
    size_t gmsh_convert_conterclock[4] = {0, 1, 3, 2};

    ifstream ifs;
    cout << "Reading data..." << endl;

    ifs.open("data/phys.txt", ios::in);
    assert(ifs.good());
    if(!ifs.good())
    {
        cerr << "[ERROR] Physical areas file reading error!" << endl;
        exit(1);
    }

    ifs >> phs_num;
    cout << " > detected " << phs_num << " physical areas" << endl;
    phs = new phys_area [phs_num];
    for(size_t i = 0; i < phs_num; i++)
    {
        ifs >> phs[i].gmsh_phys_num;
        ifs >> phs[i].lambda;
        ifs >> phs[i].rho;
        ifs >> phs[i].cp;
        phs[i].num = i;
    }
    ifs.close();

    string line;
    double garbage;
    ifstream gmsh_file;
    gmsh_file.open("data/mesh.msh", ios::in);
    assert(gmsh_file.good());
    if(!gmsh_file.good())
    {
        cerr << "[ERROR] Mesh file reading error!" << endl;
        exit(1);
    }
    do
        getline(gmsh_file, line);
    while(line.find("$Nodes") == string::npos && gmsh_file.good());

    gmsh_file >> nodes_num;
    cout << " > detected " << nodes_num << " nodes" << endl;
    nodes = new point[nodes_num];
    for(size_t i = 0; i < nodes_num; i++)
    {
        gmsh_file >> garbage >> nodes[i].x >> nodes[i].y >> garbage;
        nodes[i].num = i;
    }

    do
        getline(gmsh_file, line);
    while(line.find("$Elements") == string::npos && gmsh_file.good());

    size_t num_elem;
    size_t type_elem;
    vector<edge> tmp_bounds;
    vector<quadrilateral> tmp_rects;
    gmsh_file >> num_elem;
    for(size_t i = 0; i < num_elem; i++)
    {
        gmsh_file >> garbage >> type_elem;

        if(type_elem == 1)
        {
            edge tmp_edge;
            gmsh_file >> garbage >> tmp_edge.gmsh_phys_num >> garbage;
            for(int j = 0; j < 2; j++)
            {
                size_t tmp_node;
                gmsh_file >> tmp_node;
                tmp_edge.nodes[j] = nodes + tmp_node - 1;
            }
            if(tmp_edge.nodes[0]->num > tmp_edge.nodes[1]->num)
                swap(tmp_edge.nodes[0], tmp_edge.nodes[1]);
            tmp_bounds.push_back(tmp_edge);
        }
        else if(type_elem == 3)
        {
            quadrilateral tmp_rect;

            size_t gmsh_phys_num;
            gmsh_file >> garbage >> gmsh_phys_num >> garbage;
            tmp_rect.ph = NULL;
            for(size_t j = 0; j < phs_num; j++)
                if(phs[j].gmsh_phys_num == gmsh_phys_num)
                    tmp_rect.ph = phs + j;

            point * tmp_nodes[4];
            for(int j = 0; j < 4; j++)
            {
                size_t tmp_node;
                gmsh_file >> tmp_node;
                tmp_nodes[j] = nodes + tmp_node - 1;
            }

            // Посчитаем z компоненту нормали, чтобы понять какой обход используется
            // по часовой стрелке, или же против оной
            double nz = (tmp_nodes[1]->x - tmp_nodes[0]->x) * (tmp_nodes[2]->y - tmp_nodes[1]->y) -
                        (tmp_nodes[2]->x - tmp_nodes[1]->x) * (tmp_nodes[1]->y - tmp_nodes[0]->y);
            // И применим соответствующую конвертацию
            if(nz < 0)
                for(int j = 0; j < 4; j++)
                    tmp_rect.nodes[gmsh_convert_clock[j]] = tmp_nodes[j];
            else
                for(int j = 0; j < 4; j++)
                    tmp_rect.nodes[gmsh_convert_conterclock[j]] = tmp_nodes[j];

            tmp_rects.push_back(tmp_rect);
        }
    }
    gmsh_file.close();

    qls_num = tmp_rects.size();
    cout << " > detected " << qls_num << " rectangles" << endl;
    qls = new quadrilateral [qls_num];
    for(size_t i = 0; i < qls_num; i++)
        qls[i] = tmp_rects[i];
    tmp_rects.clear();

    bounds_num = tmp_bounds.size();
    cout << " > detected " << bounds_num << " bounds" << endl;
    bounds = new edge [bounds_num];
    for(size_t i = 0; i < bounds_num; i++)
        bounds[i] = tmp_bounds[i];
    tmp_bounds.clear();

    // Вот тут будем делать степени свободы
    for(size_t i = 0; i < qls_num; i++)
    {
        add_edge_freedom(qls[i].nodes[0], qls[i].nodes[1], i);
        add_edge_freedom(qls[i].nodes[0], qls[i].nodes[2], i);
        add_edge_freedom(qls[i].nodes[1], qls[i].nodes[3], i);
        add_edge_freedom(qls[i].nodes[2], qls[i].nodes[3], i);
    }
    sort(edges_freedom.begin(), edges_freedom.end());
    vector<edge>::iterator it;
    it = unique(edges_freedom.begin(), edges_freedom.end());
    edges_freedom.resize(distance(edges_freedom.begin(), it));

    // Занумеруем ребра
    for(size_t i = 0; i < edges_freedom.size(); i++)
        edges_freedom[i].num = i;

    // У узлов пусть будут номера узлов
    for(size_t i = 0; i < edges_freedom.size(); i++)
    {
        edges_freedom[i].degree_of_freedom[0] = edges_freedom[i].nodes[0]->num;
        edges_freedom[i].degree_of_freedom[2] = edges_freedom[i].nodes[1]->num;
    }
    // Переменная, содержащая свободную степень свободы
    size_t curr_freedom = nodes_num;
    // Теперь заполним степени свободы на ребрах
    for(size_t i = 0; i < edges_freedom.size(); i++)
    {
        edges_freedom[i].degree_of_freedom[1] = curr_freedom;
        curr_freedom++;
    }
    // Дальше заполним центры и узлы, их сразу
    for(size_t i = 0; i < qls_num; i++)
    {
        qls[i].degree_of_freedom[4] = curr_freedom;
        curr_freedom++;
        qls[i].degree_of_freedom[0] = qls[i].nodes[0]->num;
        qls[i].degree_of_freedom[2] = qls[i].nodes[1]->num;
        qls[i].degree_of_freedom[6] = qls[i].nodes[2]->num;
        qls[i].degree_of_freedom[8] = qls[i].nodes[3]->num;
    }
    // Теперь и то, что на ребрах
    for (size_t i = 0; i < qls_num; i++)
    {
        vector<edge>::iterator it;
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[0], qls[i].nodes[1]));
        if(it == edges_freedom.end())
            it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[1], qls[i].nodes[0]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[1] = it->degree_of_freedom[1];
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[0], qls[i].nodes[2]));
        if(it == edges_freedom.end())
            it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[2], qls[i].nodes[0]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[3] = it->degree_of_freedom[1];
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[1], qls[i].nodes[3]));
        if(it == edges_freedom.end())
            it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[3], qls[i].nodes[1]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[5] = it->degree_of_freedom[1];
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[2], qls[i].nodes[3]));
        if(it == edges_freedom.end())
            it = find(edges_freedom.begin(), edges_freedom.end(), edge(qls[i].nodes[3], qls[i].nodes[2]));
        assert(it != edges_freedom.end());
        qls[i].degree_of_freedom[7] = it->degree_of_freedom[1];
    }

    // Заполним также инфу о степенях свободы в ребрах для краевых
    for(size_t i = 0; i < bounds_num; i++)
    {
        vector<edge>::iterator it;
        it = find(edges_freedom.begin(), edges_freedom.end(), bounds[i]);
        if(it == edges_freedom.end()) // Здравствуй жопа новый год
        {
            swap(bounds[i].nodes[0], bounds[i].nodes[1]);
            it = find(edges_freedom.begin(), edges_freedom.end(), bounds[i]);
        }
        assert(it != edges_freedom.end());
        for(size_t j = 0; j < 3; j++)
            bounds[i].degree_of_freedom[j] = it->degree_of_freedom[j];
        // За одно и инфу о КЭ
        for(size_t j = 0; j < 2; j++)
            bounds[i].fes[j] = it->fes[j];
    }

    degree_of_freedom_num = curr_freedom;
    cout << " > detected " << degree_of_freedom_num << " degree of freedom" << endl;

    // Перенумерация степеней свободы
    map<size_t, size_t> convert;
    size_t curr_deg = 0;
    for(size_t i = 0; i < qls_num; i++)
    {
        for(size_t j = 0; j < 9; j++)
        {
            if(convert.find(qls[i].degree_of_freedom[j]) == convert.end())
            {
                convert[qls[i].degree_of_freedom[j]] = curr_deg;
                curr_deg++;
            }
        }
    }

    // Применение перенумерованных значений
    // В четырехугольниках
    for(size_t i = 0; i < qls_num; i++)
        for(size_t j = 0; j < 9; j++)
            qls[i].degree_of_freedom[j] = convert[qls[i].degree_of_freedom[j]];
    // В ребрах
    for(size_t i = 0; i < edges_freedom.size(); i++)
        for(size_t j = 0; j < 3; j++)
            edges_freedom[i].degree_of_freedom[j] = convert[edges_freedom[i].degree_of_freedom[j]];
    // В краевых
    for(size_t i = 0; i < bounds_num; i++)
        for(size_t j = 0; j < 3; j++)
            bounds[i].degree_of_freedom[j] = convert[bounds[i].degree_of_freedom[j]];

    for(size_t i = 0; i < qls_num; i++)
        qls[i].init();
}

void FEM::make_portrait()
{
    // Формирование портрета
    cout << "Generating portrait ..." << endl;
    size_t gg_size = 0;
    // Создаем массив множеств для хранения связей
    set<size_t> * portrait = new set<size_t> [degree_of_freedom_num];

    // Связь есть, если степени свободы принадлежат одному КЭ
    // Поэтому обходим конечные элементы и добавляем в множество общие степени свободы
    for(size_t k = 0; k < qls_num; k++)
    {
        for(size_t i = 0; i < 9; i++)
        {
            size_t a = qls[k].degree_of_freedom[i];
            for(size_t j = 0; j < i; j++)
            {
                size_t b = qls[k].degree_of_freedom[j];
                if(b > a) portrait[b].insert(a);
                else      portrait[a].insert(b);
            }
        }
    }

    // Считем размер матрицы
    for(size_t i = 0; i < degree_of_freedom_num; i++)
        gg_size += portrait[i].size();
    slae.alloc_all(degree_of_freedom_num, gg_size);
    // Место для хранения решения по времени
    timed = new double * [timed_num];
    for(size_t i = 0; i < timed_num; i++)
        timed[i] = new double [degree_of_freedom_num];
    // Пусть начальное время будет u_beta
    for(size_t i = 0; i < degree_of_freedom_num; i++)
        timed[0][i] = get_u_beta(point(), edge(), 0);

    cout << " > slae.n_size = " << degree_of_freedom_num << endl;
    cout << " > slae.gg_size = " << gg_size << endl;

    // Заполнение портрета
    slae.ig[0] = 0;
    slae.ig[1] = 0;
    size_t tmp = 0;
    for(size_t i = 0; i < slae.n; i++)
    {
        for(set<size_t>::iterator j = portrait[i].begin(); j != portrait[i].end(); j++)
        {
            slae.jg[tmp] = *j;
            tmp++;
        }
        slae.ig[i + 1] = slae.ig[i] + portrait[i].size();
        portrait[i].clear();
    }

    delete [] portrait;
}

void FEM::assembling_global(size_t t)
{
    cout << "Assembling global matrix... [#" << t << "]" << endl;

    memset(slae.di, 0, sizeof(double) * slae.n);
    memset(slae.rp, 0, sizeof(double) * slae.n);
    memset(slae.gl, 0, sizeof(double) * slae.ig[slae.n]);
    memset(slae.gu, 0, sizeof(double) * slae.ig[slae.n]);

    double time = (double)t * delta_t;

    for(size_t k = 0; k < qls_num; k++)
    {
        double q_old[9];
        for(size_t i = 0; i < 9; i++)
            q_old[i] = timed[t - 1][qls[k].degree_of_freedom[i]];
        double rhocp = qls[k].ph->rho * qls[k].ph->cp;

        for(size_t i = 0; i < 9; i++)
        {
            double q_old_add = 0.0;
            for(size_t j = 0; j < 9; j++)
            {
                slae.add(qls[k].degree_of_freedom[i], qls[k].degree_of_freedom[j],
                         qls[k].get_local_G(i, j) * qls[k].ph->lambda);
                double Mij = qls[k].get_local_M(i, j) * rhocp / delta_t;
                slae.add(qls[k].degree_of_freedom[i], qls[k].degree_of_freedom[j], Mij);
                q_old_add += Mij * q_old[j];

                slae.add(qls[k].degree_of_freedom[i], qls[k].degree_of_freedom[j], rhocp * qls[k].get_local_C(i, j, time));
            }

            slae.rp[qls[k].degree_of_freedom[i]] += qls[k].get_local_rp(i, time) + q_old_add;
        }
    }
}

void FEM::applying_bounds(size_t t)
{
    cout << "Applying bounds... [#" << t << "]" << endl;

    double time = (double)t * delta_t;

    for(size_t bi = 0; bi < bounds_num; bi++)
    {
        if(get_type_of_bound(bounds[bi].gmsh_phys_num, time) == 2)
        {
            double theta[3] =
            {
                get_theta(bounds[bi].get_freedom_position(0), bounds[bi], time),
                get_theta(bounds[bi].get_freedom_position(1), bounds[bi], time),
                get_theta(bounds[bi].get_freedom_position(2), bounds[bi], time)
            };

            for(size_t i = 0; i < 3; i++)
            {
                double rp_add = 0.0;
                for(size_t j = 0; j < 3; j++)
                    rp_add += bounds[bi].get_matrix_M(i, j) * theta[j];
                slae.rp[bounds[bi].degree_of_freedom[i]] += rp_add;
            }
        }
        else if(get_type_of_bound(bounds[bi].gmsh_phys_num, time) == 3)
        {
            double beta = get_beta(bounds[bi], time);
            double u_beta[3] =
            {
                get_u_beta(bounds[bi].get_freedom_position(0), bounds[bi], time),
                get_u_beta(bounds[bi].get_freedom_position(1), bounds[bi], time),
                get_u_beta(bounds[bi].get_freedom_position(2), bounds[bi], time)
            };

            for(size_t i = 0; i < 3; i++)
            {
                size_t dof_i = bounds[bi].degree_of_freedom[i];
                double rp_add = 0.0;
                for(size_t j = 0; j < 3; j++)
                {
                    size_t dof_j = bounds[bi].degree_of_freedom[j];
                    double Mij = bounds[bi].get_matrix_M(i, j);
                    slae.add(dof_i, dof_j, beta * Mij);
                    rp_add += Mij * beta * u_beta[j];
                }
                slae.rp[bounds[bi].degree_of_freedom[i]] += rp_add;
            }
        }
    }

    for(size_t bi = 0; bi < bounds_num; bi++)
    {
        if(get_type_of_bound(bounds[bi].gmsh_phys_num, time) == 1)
        {
            for(size_t j = 0; j < 3; j++)
            {
                size_t freedom = bounds[bi].degree_of_freedom[j];
                double val = get_bound_value(bounds[bi].get_freedom_position(j), bounds[bi], time);
                // В диагональ пишем 1
                slae.di[freedom] = 1.0;
                // В правую часть пишем значение краевого
                slae.rp[freedom] = val;
                // Начальное приближение сразу знаем
                slae.x[freedom] = val;
                // Обнуляем строку
                for(size_t k = slae.ig[freedom]; k < slae.ig[freedom + 1]; k++)
                    slae.gl[k] = 0.0;
                for(size_t k = 0; k < slae.ig[slae.n]; k++)
                    if(slae.jg[k] == freedom)
                        slae.gu[k] = 0.0;
            }
        }
    }
}

const quadrilateral * FEM::get_fe(const point & p) const
{
    for(size_t i = 0; i < qls_num; i++)
    {
        if(qls[i].inside(p))
            return qls + i;
    }
    cerr << "Warning: Target point " << p << " is outside of area" << endl;
    return NULL;
}

double FEM::get_solution(const point & p, size_t time_curr) const
{
    return get_solution(p, get_fe(p), time_curr);
}

double FEM::get_solution(const point & p, const quadrilateral * fe, size_t time_curr) const
{
    if(fe)
    {
        double solution = 0.0;
        for(size_t i = 0; i < 9; i++)
            solution += fe->bfunc_2d(i, p) * timed[time_curr][fe->degree_of_freedom[i]];
        return solution;
    }
    return 0.0;
}

point FEM::get_grad(const point & p, size_t time_curr) const
{
    return get_grad(p, get_fe(p), time_curr);
}

point FEM::get_grad(const point & p, const quadrilateral * fe, size_t time_curr) const
{
    if(fe)
    {
        point solution(0.0, 0.0);
        for(size_t i = 0; i < 9; i++)
        {
            double q = timed[time_curr][fe->degree_of_freedom[i]];
            point gbf = fe->grad_bfunc_2d(i, p);
            solution.x += q * gbf.x;
            solution.y += q * gbf.y;
        }
        return solution;
    }
    return point(0.0, 0.0);
}

void FEM::set_timed()
{
    extern size_t TIMES_NUM;
    extern double TIMES_STEP;
    delta_t = TIMES_STEP;
    timed_num = TIMES_NUM;
}

void FEM::solve_timed()
{
    for(size_t i = 1; i < timed_num; i++)
    {
        assembling_global(i);
        applying_bounds(i);
        slae.solve(1e-20);
        memcpy(timed[i], slae.x, slae.n * sizeof(double));
    }
}

// ============================================================================

vector<edge>::const_iterator FEM::get_edge_by_nodes(const point & begin, const point & end) const
{
    vector<edge>::const_iterator it;
    it = find(edges_freedom.begin(), edges_freedom.end(), edge(&begin, &end));
    if(it == edges_freedom.end())
        it = find(edges_freedom.begin(), edges_freedom.end(), edge(&end, &begin));
    return it;
}

vector<edge>::const_iterator FEM::get_edge_by_fe(const quadrilateral & fe, size_t edge_num) const
{
    if(edge_num == 0) return get_edge_by_nodes(* fe.nodes[0], * fe.nodes[1]);
    if(edge_num == 1) return get_edge_by_nodes(* fe.nodes[0], * fe.nodes[2]);
    if(edge_num == 2) return get_edge_by_nodes(* fe.nodes[1], * fe.nodes[3]);
    if(edge_num == 3) return get_edge_by_nodes(* fe.nodes[2], * fe.nodes[3]);
    assert(edge_num < 4);
    return edges_freedom.end();
}

const point * FEM::get_nodes_by_edge(const edge & e, size_t node_num) const
{
    assert(node_num < 2);
    return e.nodes[node_num];
}

const quadrilateral * FEM::get_fe_by_edge(const edge & e, size_t fe_num) const
{
    assert(fe_num < 2);
    return e.fes[fe_num];
}

