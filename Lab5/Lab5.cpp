#include <iostream>
#include <cmath>
#include <time.h>

#define N 1000
using namespace std;

int main()
{
    int i, j, k;
    double h = 1.0 / N;
    double A = 0, B = 0;
    double Tms;
    Tms = clock();
    printf("cycle 0\n");
    //---------------------------2
#pragma omp parallel for reduction(+: A)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += exp(-h * (i + 0, 5)) / (1 + pow(sin(h * (i + 0, 5) + h * (j + 0, 5) + h * (k + 0, 5)), 2));
        }
    }
    printf("cycle 1\n");
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (static)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += exp(-h * (i + 0, 5)) / (1 + pow(sin(h * (i + 0, 5) + h * (j + 0, 5) + h * (k + 0, 5)), 2));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time static = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (dynamic)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += exp(-h * (i + 0, 5)) / (1 + pow(sin(h * (i + 0, 5) + h * (j + 0, 5) + h * (k + 0, 5)), 2));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time dynamic = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

    // ---------------------- 4
#pragma omp parallel for reduction(+: A)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += 1 / (h * (k + 0, 5) + exp(-h * (k + 0, 5)) * sin(pow(h * (i + 0, 5), 2) - pow(h * (j + 0, 5), 2)));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (static)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += 1 / (h * (k + 0, 5) + exp(-h * (k + 0, 5)) * sin(pow(h * (i + 0, 5), 2) - pow(h * (j + 0, 5), 2)));

        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time static = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (dynamic)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += 1 / (h * (k + 0, 5) + exp(-h * (k + 0, 5)) * sin(pow(h * (i + 0, 5), 2) - pow(h * (j + 0, 5), 2)));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time dynamic = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

    // ------------------------- 6
#pragma omp parallel for reduction(+: A)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += (cos(h * (i + 0, 5))) / (1 + exp(-(pow(h * (i + 0, 5), 2) + (pow(h * (j + 0, 5), 2) - pow(h * (k + 0, 5), 2)))));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (static)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += (cos(h * (i + 0, 5))) / (1 + exp(-(pow(h * (i + 0, 5), 2) + (pow(h * (j + 0, 5), 2) - pow(h * (k + 0, 5), 2)))));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time static = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: A) schedule (dynamic)
    for (i = 0; i < N - 1; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < N - 1; k++)
                A += (cos(h * (i + 0, 5))) / (1 + exp(-(pow(h * (i + 0, 5), 2) + (pow(h * (j + 0, 5), 2) - pow(h * (k + 0, 5), 2)))));
        }
    }
    A *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time dynamic = " << Tms << " sec" << endl;
    cout << "A = " << A << endl;
    A = 0;

#pragma omp parallel for reduction(+: B)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
                B += exp(h * i) / (1 + pow(sin(h * i + h * j + h * k), 2));
        }
    }
    B *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time = " << Tms << " sec" << endl;
    cout << "B = " << B << endl;
    B = 0;

#pragma omp parallel for reduction(+: B) schedule (static)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
                B += exp(h * i) / (1 + pow(sin(h * i + h * j + h * k), 2));
        }
    }
    B *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time static = " << Tms << " sec" << endl;
    cout << "B = " << B << endl;
    B = 0;

#pragma omp parallel for reduction(+: B) schedule (dynamic)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N - 1; j++) {
            for (k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
                B += exp(h * i) / (1 + pow(sin(h * i + h * j + h * k), 2));
        }
    }
    B *= pow(h, 3);
    Tms = (clock() - Tms) / CLOCKS_PER_SEC;
    cout << "Time dynamic = " << Tms << " sec" << endl;
    cout << "B = " << B << endl;
    B = 0;
    printf("cycle last");
}