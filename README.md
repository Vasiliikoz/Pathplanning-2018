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
Входной файл должен быть формата .xml. Тэг `<map>` обозначает начало блока, в котором описано начальное простраство.
Тэги `<width>` и `<height>` обозначают блоки, в которых описаны размеры поля. Тэги `<startx>` и `<starty>`, `<finishx>` и `<finishy>` обозначают координаты начальной и конечной точек соответсвенно. Тэг `<grid>` означает начало блока, в котором описана карта. Если в i-ом ряду(`<row>`) на j-ом месте стоит 0, то данная клетка проходима, и если стоит 1 - то не проходима. Тэг `<algorithm>` обозначает начало блока, задающего параметры, с которыми запускается алгоритм, в данном блоке должны находиться такие тэги: `<searchtype>`(определяет алгоритм поиска), `<metrictype>`(определяет используемую эвристику), `<breakingties>`(определяет выбор вершин при равентсве F + G выражений), `<hweight>`(вез эврестической функции), `<allowdiagonal>`(разрешено ли ходить по диагоналям), `<cutcorners>`(разрешено ли срезание углов), `<allowsqueeze>`(разрешено ли просачивание). Тэг `<logpath>` обозначает блок, в котором описан путь для выходного файла.
#### Формат выходного и выходного файлов:
Выходной файл содержит информацию о пространстве и алгоритме, описанную в входном файле. Также в блоке, соответсвующему тэгу `<log>`, он содержит информацию о том, как отработал алгоритм. В блоке тэга `<path>` визуализировано поле(аналогично `<grid>`), клетки, по которым пролегает оптимальный путь помечены как "*". В блоке тэга `<lplevel>` по прядку перечисленные вершины, по которым пролегат оптимальный путь. В блоке тэга `<summary>` записано количество шагов, созданных вершин и длина пути.
#### Дедлайны:
17.12.2018   

* выложить в LMS проект технического задания
* реализация Дейксты и А-стар
* сдать в учебный офис титульный лист ТЗ
15.02.2018.  


* выложить в LMS проект пояснительной записки
* сдать в учебный офис титульный лист пояснительной записки
* предоставление доступа к github учебному офису
* реализация Theta-star
* реализация JPS
25.03.2018   


* выложить в LMS и сдать в УО пояснительную записку
* выложить в LMS и сдать в УО ТЗ
* выложить в LMS и сдать в УО отзыв руководителя практики
* выложить в LMS и сдать в УО отчёт о проверке пояснительной записки на плагиат
