# ABC-lab4
# Оптимизация доступа к памяти
1. На языке С/С++/C# реализовать функцию DGEMM BLAS последовательное умножение двух
квадратных матриц с элементами типа double. Обеспечить возможность задавать
размерности матриц в качестве аргумента командной строки при запуске программы.
Инициализировать начальные значения матриц случайными числами.
2. Провести серию испытаний и построить график зависимости времени выполнения
программы от объёма входных данных. Например, для квадратных матриц с числом
строк/столбцов 1000, 2000, 3000, ... 10000.
3. Оценить предельные размеры матриц, которые можно перемножить на вашем
вычислительном устройстве.
4. Реализовать дополнительную функцию DGEMM_opt_1, в которой выполняется
оптимизация доступа к памяти, за счет построчного перебора элементов обеих матриц.
5. * Реализовать дополнительную функцию DGEMM_opt_2, в которой выполняется
оптимизация доступа к памяти, за счет блочного перебора элементов матриц. Обеспечить
возможность задавать блока, в качестве аргумента функции.
6. ** Реализовать дополнительную функцию DGEMM_opt_3, в которой выполняется
оптимизация доступа к памяти, за счет векторизации кода.
7. Оценить ускорение умножения для матриц фиксированного размера, например,
1000х1000, 2000х2000, 5000х5000, 10000х10000.
* Для блочного умножения матриц определить размер блока, при котором достигается
максимальное ускорение.
8. С помощью профилировщика для исходной программы и каждого способа оптимизации
доступа к памяти оценить количество промахов при работе к КЭШ памятью (cache-misses).
9. Подготовить отчет отражающий суть, этапы и результаты проделанной работы.
