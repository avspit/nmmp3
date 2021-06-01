import finiteElementMethod as fem
import constants.variables as const
import utils.logger as logger
import utils.util as util
from matplotlib.widgets import Slider
import pylab


def update_grafs():
    # Будем использовать A и C, установленные с помощью слайдеров
    global slider_A
    global slider_C
    global graph_axes

    # Используем атрибут val, чтобы получить значение слайдеров
    A = slider_A.val
    C = slider_C.val

    graph_axes.clear()

    # Производим вычисления х и у
    for n in const.N_ARR:
        h = util.init_h(n)
        logger.log(text='Начинаем вычисление для h:', value=str(h), force=True)
        x = util.init_x(h, n)
        y = fem.solve(A, C, n, x)
        graph_axes.plot(x, y, label='h='+str(h))
        logger.log(text='Вычисление завершено для h:', value=str(h), force=True)
    graph_axes.legend()
    pylab.draw()
    pylab.show()


def onChangeValue(value):
    # Обработчик события изменения значений слайдеров
    update_grafs()


if __name__ == '__main__':
    # Создадим окно с графиком
    fig, graph_axes = pylab.subplots()

    # Оставим снизу от графика место для виджетов
    fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.4)
    # Создаем слайдер для задания параметра A
    axes_slider_A = pylab.axes([0.05, 0.25, 0.85, 0.04])
    slider_A = Slider(axes_slider_A,
                          label='A',
                          valmin=0,
                          valmax=5,
                          valinit=1,
                          valfmt='%1.2f')
    # Подпишемся на событие при изменении значения слайдера.
    slider_A.on_changed(onChangeValue)

    # Создаем слайдер для задания параметра С
    axes_slider_C = pylab.axes([0.05, 0.17, 0.85, 0.04])
    slider_C = Slider(axes_slider_C,
                       label='C',
                       valmin=-10,
                       valmax=10,
                       valinit=0,
                       valfmt='%1.2f')
    # Подпишемся на событие при изменении значения слайдера.
    slider_C.on_changed(onChangeValue)

    update_grafs()



