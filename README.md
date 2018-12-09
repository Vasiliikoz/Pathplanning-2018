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
* написать команду: \\
    `cd путь к вашей папке`
* написать команду: \\
    `qmake ASearch.pro`
* написать команду: \\
    `make`
* написать команду: \\
    `./ASearch путь до вашего файла .xml`
В результате проделаных действий в папке с проектом будет создан .xml файл, в котором будет описан искомый путь, а при его отсутсвие, будет сказано об этом. 
#### Формат входного и выходного файлов:
