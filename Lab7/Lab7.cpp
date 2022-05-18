#include "mpi.h"
#include <iostream>
void lab7_tast_2
{
    double a = 0, b = 0, c = 0;
    int const BufSize = sizeof(double) + MPI_BSEND_OVERHEAD;
    MPI_Init(NULL, NULL);
    int Size;
    MPI_Comm_size(MPI_COMM_WORLD, &Size);
    int Rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    int Posl = (Rank + 1) % Size; //Номер последующего процесса
    int Pred = Rank ? Rank - 1 : Size - 1; //Номер предшествующего процесса
    a = Rank + 0.7;
    //void *Addr=NULL;
    void* Buf = malloc(BufSize);
    MPI_Buffer_attach(Buf, BufSize);
    MPI_Bsend(&a, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD);
    MPI_Status St;
    MPI_Recv(&b, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &St);
    int Tmp;
    MPI_Buffer_detach(Buf, &Tmp);
    MPI_Buffer_attach(Buf, BufSize);
    MPI_Bsend(&a, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD);
    MPI_Recv(&c, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &St);
    MPI_Buffer_detach(Buf, &Tmp);
    if (!Buf) free(Buf);
    cout << "Process " << Rank << " a=" << a << " b=" << b << " c=" << c << endl;
    MPI_Finalize();
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Lab 7\n";
    cout << "Task 2\n";
    task_2();
    cout << "Task 4\n";
    task_4();
    cout << "Task 6\n";
    task_6();
}