#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <locale.h>

using namespace std;

//2.Проиллюстрировать возможность инициализации локальных переменных, объявленных как threadprivate,
//начальным значением одноименной глобальной переменной из главного потока, при помощи опции copyin
int n;
#pragma omp threadprivate(n)
void task_2()
{
	//setlocale(LC_ALL, ".ACP");
	n = 1;
#pragma omp parallel num_threads(4) copyin(n)
	{
		printf("Значение n: %d\n", n);
	}
}

//глобальная переменная вносится в качестве локальной для потоков



//4. На примере задачи нахождения суммы векторов достаточно большой размерности, проиллюстрировать распределение итераций цикла 
//между потоками в параллельной области при помощи директивы omp for.
#define NN 35
void task_4()
{
	printf("\n");
	//setlocale(LC_ALL, ".ACP");
	int A[NN], B[NN], C[NN], i, n;
	/* Заполним исходные массивы */
	for (i = 0; i < NN; i++) { A[i] = i; B[i] = 2 * i; C[i] = 0; }
#pragma omp parallel num_threads(4) shared(A, B, C) private(i, n)
	{
		/* Получим номер текущей нити */
		n = omp_get_thread_num();
#pragma omp for
		for (i = 0; i < NN; i++)
		{
			C[i] = A[i] + B[i];
			printf("Поток %d сложил элементы с номером %d\n",
				n, i);
		}
	}
}
//свободный поток выполняет итерации последовательно по их освобождению


//6. Разработать программу, иллюстрирующую возможность реализации конечного (неитеративного) параллелизма при помощи директив omp sections и omp section
void task_6()
{
	printf("\n");
	//setlocale(LC_ALL, ".ACP");
	int n;
#pragma omp parallel num_threads(3) private(n)
	{
		n = omp_get_thread_num();
#pragma omp sections
		{
#pragma omp section
			{
				printf("Первая секция, поток %d\n", n);
			}
#pragma omp section
			{
				printf("Вторая секция, поток %d\n", n);
			}
#pragma omp section
			{
				printf("Третья секция, поток %d\n", n);
			}
#pragma omp section
			{
				printf("Четвертая секция, поток %d\n", n);
			}
#pragma omp section
			{
				printf("Пятая секция, поток %d\n", n);
			}
		}
		printf("Параллельная область, поток %d\n", n);
	}

}
//секции распределяются между потоками, пока есть секция для выполнения, свободный поток ее выполняет, пока все не выполнятся секции, потоки не закроются

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Lab 3\n";
	cout << "Task 2\n";
	task_2();
	cout << "Task 4\n";
	task_4();
	cout << "Task 6\n";
	task_6();
}

