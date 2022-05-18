#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <locale.h>

using namespace std;

//2. . При помощи опции ordered и директивы omp ordered
//внутри параллельного цикла гарантировать выполнения операций в том
//порядке, как они выполнялись бы в последовательном цикле.
#define NNN 12
void task_2()
{
	int i, n;
#pragma omp parallel private (i, n)
	{
		n = omp_get_thread_num(); //номер выполняемого потока
#pragma omp for ordered
		for (i = 0; i < NNN; i++)
		{
#pragma omp ordered
			{
				printf("Поток %d, итерация %d\n", n, i);
			}
		}
	}
}



//4. Частным случаем критических секций является директива omp atomic, непосредственно воздействующая на следующий за 
//ней оператор присваивания.Используя данную директиву, написать программу, которая подсчитывает число потоков в параллельной области.

void task_4()
{
	printf("\n");
	int count = 0;
#pragma omp parallel
	{
#pragma omp atomic //на каждой итерации вызывает новый поток
		count++;
	}
	printf("Число потоков: %d\n", count);
}



//6.  
// При помощи блокировки (обеспечивается стандартными функциями OpenMP omp_init_lock, omp_set_lock, omp_unset_lock), 
//добиться того, чтобы сообщения, печатаемые в разных потоках в параллельной секции, не перемежались друг с другом.
//Для имитации вычислений использовать функцию sleep(Sleep в ОС Windows).
// То же самое, но с попыткой установить блокировку, и,
//если секция уже заблокирована, с печатью ожидающими потоками с интервалом 
//в 1 секунду сообщения о том, что секция блокирована. Использовать функцию omp_test_lock

#ifdef _WIN32
#include <windows.h>
void sleep(int Secnds)
{
	Sleep(1000 * Secnds);
}
#endif

omp_lock_t lock;
void task_6()
{
	printf("\n");
	int n;
	omp_init_lock(&lock);
#pragma omp parallel num_threads(4) private (n)
	{
		n = omp_get_thread_num();
		while (!omp_test_lock(&lock))
		{
			printf("Секция закрыта, поток %d\n", n);
			sleep(2);
		}
		printf("Начало закрытой секции, поток %d\n", n);
		sleep(2 + 3 * ((n + 1) % 2));
		printf("Конец закрытой секции, поток %d\n", n);
		omp_unset_lock(&lock);
	}
	omp_destroy_lock(&lock);
}
//если секция занята потоком, то другой поток ее не трогает, а выводит сообщение о ее занятости


int main()
{
	setlocale(LC_ALL, "Russian");
	cout << "Lab 4\n";
	cout << "Task 2\n";
	task_2();
	cout << "Task 4\n";
	task_4();
	cout << "Task 6\n";
	task_6();
}

