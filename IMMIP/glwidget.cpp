

#include "glwidget.h"

// Конструктор
glwidget::glwidget()
{
    // На нулевом временном слое
    time_curr = 0;

    // Решаем МКЭ задачу
    fem.set_timed();
    fem.input();
    fem.make_portrait();
    fem.solve_timed();

    std::ofstream out;
    out.open("D:\\nodes.txt");
    out << fem.nodes_num << endl;
    for (size_t i = 0; i < fem.nodes_num; i++)
    {
        out << fem.nodes[i].x << "\t" << fem.nodes[i].y << endl;
    }
    out.close();

    out.open("D:\\elements.txt");
    out << fem.qls_num * 2 << endl;
    for (size_t i = 0; i < fem.qls_num; i++)
    {
        out << fem.qls[i].nodes[0]->num << "\t" << fem.qls[i].nodes[1]->num << "\t" << fem.qls[i].nodes[2]->num << endl;
        out << fem.qls[i].nodes[1]->num << "\t" << fem.qls[i].nodes[2]->num << "\t" << fem.qls[i].nodes[3]->num << endl;
    }
    out.close();
}

// Изменение количества сегментов, на которые разбивается каждый КЭ
void glwidget::set_div_num(size_t num)
{
    std::ofstream out;
    string s1;
    string s2;
    string s3;
    string s;

    s1 = "D:\\solution";
    s2 = to_string(time_curr);
    s3 = ".txt";
    s = s1 + s2 + s3;
    out.open(s);
    out << fem.nodes_num << endl;
    for (size_t i = 0; i < fem.nodes_num; i++)
    {
        out << fem.get_solution(fem.nodes[i], time_curr) << endl;
    }
    out.close();
}

// Задание временного слоя
void glwidget::set_time_num(size_t num)
{
    time_curr = num;
    set_div_num(div_num);
}

// Подгонка осей под реальность и вычисление шагов координатной сетки
void glwidget::adjustAxis(double & min, double & max, size_t & numTicks)
{
    static const double axis_epsilon = 1.0 / 10000.0;
    if(max - min < axis_epsilon)
    {
        min -= 2.0 * axis_epsilon;
        max += 2.0 * axis_epsilon;
    }

    static const size_t MinTicks = 10;
    double grossStep = (max - min) / MinTicks;
    double step = pow(10, floor(log10(grossStep)));

    if (5 * step < grossStep)
        step *= 5;
    else if (2 * step < grossStep)
        step *= 2;

    numTicks = (size_t)(ceil(max / step) - floor(min / step));
    min = floor(min / step) * step;
    max = ceil(max / step) * step;
}

// Инициализация сцены
void glwidget::initializeGL()
{
    
}

// Действие при изменении размеров виджета
void glwidget::resizeGL(int nWidth, int nHeight)
{

}

// Отрисовка сцены
void glwidget::paintGL()
{

}
