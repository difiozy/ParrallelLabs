﻿#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <locale.h>

using namespace std;

//2.  Проиллюстрировать возможность директивы single, осуществив из параллельной области печать 
//сообщения лишь одним потоком(не обязательно главным)
void task_2()
{
	//setlocale(LC_ALL, ".ACP");

#pragma omp parallel num_threads(2)
	{
		printf("Печать сообщения 1\n");
#pragma omp single
		{
			printf("Один поток\n");
		}
		printf("Печать сообщения 2\n");
	}

}
//Мы задаем 3 потока, но single выполняется только один раз в первом потоке, остальные игнорирует. 
//Код в директиве выполняется 1 раз вне зависимости от числа потоков


//4. Проиллюстрировать возможность синхронизации в параллельной области значений локальных для потоков переменных
// (т.е.явного задания нужного значения) посредством опции copyprivate директивы single.
void task_4()
{
	//setlocale(LC_ALL, ".ACP");
	int n;
	printf("\n");
#pragma omp parallel num_threads(3) private(n)
	{
		n = omp_get_thread_num(); //Каждый поток печатает свой номер
		printf("Значение n (начало): %d\n", n);
#pragma omp single copyprivate(n)
		{
			n = 100; //На выходе из секции single всем 
		} //локальным переменным будет присвоено значение 100
		printf("Значение n (конец): %d\n", n);
	}
}


//6. Считается, что в большинстве реализаций опция private действует следующим образом: 1) значения локальных переменных, 
//создаваемых в параллельной области, не инициализируются значением
//одноименной переменной, созданной в последовательной области до
//входа в параллельную; 2) значение одноименной переменной, созданной
//в последовательной области до входа в параллельную, не зависит от значений локальных переменных и остается неизменным после выхода из
//параллельной области.Разработать соответствующую программу
int n;
#pragma omp threadprivate(n)
void task_6()
{
	//setlocale(LC_ALL, ".ACP");
	n = 1;
	printf("\nn в последовательной области (начало): %d\n", n);
#pragma omp parallel num_threads(4) copyin(n)
	{
		printf("Значение n на нити (на входе): %d\n", n);
		/* Присвоим переменной n номер текущей нити */
		n = omp_get_thread_num();
		printf("Значение n на нити (на выходе): %d\n", n);
	}
	printf("n в последовательной области (конец): %d\n", n);
}

//Как видно, при межпроцедурной оптимизации некоторые компиляторы: 1) сопоставляют локальной переменной основного потока 
//в параллельной области одноименную переменную, созданную в последовательной области до входа в параллельную;
//2) инициализируют в параллельной области значения локальных для потоков переменных значением
//одноименной переменной, созданной в последовательной области до входа в параллельную.


int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Lab 2\n";
	cout << "Task 2\n";
	task_2();
	cout << "Task 4\n";
	task_4();
	cout << "Task 6\n";
	task_6();
}
