#ifndef GLWIDGET_H
#define GLWIDGET_H

#include "fem.h"
#include <fstream>

// Класс треугольник
class triangle
{
public:
    point nodes[3];
    double solution[3];
    triangle() {}
    triangle(const point & node1, const point & node2, const point & node3)
    {
        nodes[0] = node1;
        nodes[1] = node2;
        nodes[2] = node3;
    }
};

// Класс OpenGL виджет
class glwidget
{
public:
    // Инициализация сцены
    void initializeGL();
    // Действие при изменении размеров виджета
    void resizeGL(int nWidth, int nHeight);
    // Отрисовка сцены
    void paintGL();
    // Конструктор
    glwidget();

    // Изменение количества сегментов, на которые разбивается каждый КЭ
    void set_div_num(size_t num);
    // Задание временного слоя
    void set_time_num(size_t num);
private:
    // Класс МКЭ
    FEM fem;

    // Минимальные и максимальные значения геометрии + размер
    double min_x, max_x, size_x;
    double min_y, max_y, size_y;

    // Количество шагов координатной сетки
    size_t num_ticks_x, num_ticks_y;
    // Подгонка осей под реальность и вычисление шагов координатной сетки
    void adjustAxis(double & min, double & max, size_t & numTicks);

    // Вспомогательные шаги по цвету для закраски
    double step_u_big, step_u_small;

    // Треугольники, которые будем рисовать
    vector<triangle> triangles;

    // Текущий временной слой
    size_t time_curr;

    size_t div_num;
    size_t isol_num;
};

#endif // GLWIDGET_H
