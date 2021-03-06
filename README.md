# Pathplanning-2018
Проектная работа на ФКН 2018. Выполнил Козлов Василий, студент 2-го курса группы БПМИ173.
#### Аннотация
В данной работе мною будут реализованы алгоритмы поиска кратчайшего пути: A-star, Дейстра, Theta*, Jump Point Search. Данный список может быть дополнен в ходе работы над проектом.
#### Цели:
* Улучшить навыки в сфере написания и разработки алгоритмов
* Овладеть основами ООП
* Попрактиковаться в изучении материала на английском языке
#### Как запустить проект:
Мною будет рассмотрен запуск проекта на Mac OS. Необходимым условием является наличие Qt-creator, в противном случае, его можно скачать по ссылке https://www.qt.io/download.
В качестве входным данных моя програма принимает .xml файл(формат данного файла описан ниже).
* скачать файлы(все кроме .md) и поместить их в одну папку.
* запустить командную строку
* написать команду: 
    `cd путь к вашей папке`
* написать команду: 
    `qmake ASearch.pro`
* написать команду: 
    `make`
* написать команду: 
    `./ASearch путь до вашего файла .xml`   
    
    
В результате проделаных действий в папке с проектом будет создан .xml файл, в котором будет описан искомый путь, а при его отсутсвие, будет сказано об этом. 
#### Формат входного и выходного файлов:
Здесь описан формат входного и выходного файлов. Входной файл должен быть формата .xml. Приведём тэги и опишем, как составить входной файл. Каждый блок закрывается `<\”тэг открывающий этот блок”>`.(например, `<map>` и `<\map>`)
*	`<map>` обозначает начало блока, в котором описано начальное пространство.
*	`<width>` находиться внутри блока, обозначенного <map>, отвечает за ширину (кол-во столбцов) поля
*	`<height>` находиться внутри блока, обозначенного <map>, отвечает за высоту (кол-во строк) поля
*	`<startx>` находиться внутри блока, обозначенного <map>, обозначает блок, задающий x – координату начальной вершины (положительное целое число).
*	`<starty>` находиться внутри блока, обозначенного <map>, обозначает блок, задающий y – координату начальной вершины (положительное целое число).
*	`<finishx>` находиться внутри блока, обозначенного <map>, обозначает блок, задающий x – координату конечной вершины (положительное целое число).
*	`<finishy>` находиться внутри блока, обозначенного <map>, обозначает блок, задающий y – координату конечной вершины (положительное целое число).
*	`<grid>` находиться внутри блока, обозначенного <map>, обозначает начало блока, задающего описания ячеек пространство.
*	`<row>` находиться внутри блока, обозначенного <grid>, обозначает начало блока, в котором описан ряд поля, т.е. задана последовательность из 1 и 0, разделенных пробелами. Символ 1 в i – ом ряду j – ом столбце обозначает, что эта ячейка поля не проходима, символ 0 – проходима.
*	`<algorithm>` обозначает начало блока, задающего параметры, с которыми запускается алгоритм.
*	`<searchtype>` находится в блоке, заданном <algorithm>, определяет используемый алгоритм (“astar”, “dijkstra”, “jp_search”, “theta”, “bfs”).
*	`<metrictype>` находится в блоке, заданном <algorithm>, определяет используемую эвристику (“diagonal”, “manhattan”, “euclidean”, “chebyshev”).
*	`<breakingties>` находится в блоке, заданном <algorithm>, определяет выбор вершины, при равенстве F + G (“g-max”, “f-max”).
*	`<hweight>` находится в блоке, заданном <algorithm>, определяет вес эвристики (целое число).
*	`<allowdiagonal>` находится в блоке, заданном <algorithm>, определяет разрешено ли ходить по диагоналям (“true”, “false”).
*	`<cutcorners>` находится в блоке, заданном <algorithm>, определяет разрешено ли срезать углы (“true”, “false”).
*	`<allowsqueeze>` находится в блоке, заданном <algorithm>, определяет разрешено ли просачиваться (“true”, “false”).
*	`<logpath>` обозначает блок, в котором описан путь для выходного файла.
Выходной файл содержит информацию о пространстве и алгоритме, описанную в входном файле.
*	`<log>` о том, как отработал алгоритм.
###### Перечисленные далее тэги лежат в `<log>`.
*	`<path>` визуализировано поле (аналогично <grid>). Клетки, по которым пролегает оптимальный путь помечены как "*".
*	`<lplevel>` по прядку перечисленные вершины, по которым пролегает оптимальный путь (<node x = “x - координата”, y = “y - координата”, number = “номер вершины”,>).
*	`<hplevel>` по прядку перечисленные вершины, в которых меняется направление (аналогично <lplevel>).
*	`<summary>` записано количество шагов, созданных вершин и длина пути

#### Дедлайны:
1. 17.12.2018. 
* выложить в LMS проект технического задания
* реализация Дейксты и А-стар
* сдать в учебный офис титульный лист ТЗ
2. 15.02.2018. 
* выложить в LMS проект пояснительной записки
* сдать в учебный офис титульный лист пояснительной записки
* предоставление доступа к github учебному офису
* реализация Theta-star
* реализация JPS
3. 25.03.2018. 
* выложить в LMS и сдать в УО пояснительную записку
* выложить в LMS и сдать в УО ТЗ
* выложить в LMS и сдать в УО отзыв руководителя практики
* выложить в LMS и сдать в УО отчёт о проверке пояснительной записки на плагиат
