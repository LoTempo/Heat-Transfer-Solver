# Heat Transfer Solver
Простое приложение на Python для решения одномерной задачи теплопереноса в пористой среде в соответствии с двухтемпературной моделью при течении несжимаемой жидкости.

Математическая постановка решаемой задачи представлена в файле problem.pdf. Задача решается методом конечных разностей по явной схеме аппроксимации.

В приложении реализовано построение графика изменения температуры твердотельного каркаса и жидкости в пределах рассматриваемой области. При помощи слайдера возможна визуализация изменения температурного состояния пористой среды с течением времени.

Интерфейс приложение выполнен при помощи модуля Tkinter.
