/*PThread. Разработать, используя средства многопотокового программирования, параллельную программу решения одномерной нестационарной краевой задачи методом конечных разностей (МКР) с использованием явной вычислительной схемы.
Дан цилиндрический стержень длиной L и площадью поперечного сечения S. Цилиндричекая поверхность стержня теплоизолирована. На торцевых поверхностях стержня слева и справа могут иметь место граничные условия первого, второго или третьего родов. Распределение поля температур по длине стержня описывается уравнением теплопроводности

dT/dt = aT*d2T/dx2+gT

, где

aT=lambda/(CT*p) - коэффициент температуропроводности;

lambda - коэффициент теплопроводности среды;

CT - удельная теплоемкость единицы массы;

p - плотность среды;

gT=GT/(CT*p) - приведенная скорость взаимного превращения тепловой энергии в другие виды энергии, в нашем случае gT=0.

В явной вычислительной схеме МКР для аппроксимации производной температуры по времени в узле, принадлежащем i-ому временному и j-ому пространственному слоям, используется "разница вперед":

dT/dt|ij = (Ti+1j-Ti j)/ht
, где ht - шаг дискретизации по оси времени.

Для аппроксимации второй производной температуры по пространственной координате x используется "центральная разница":

d2T/dx2|ij = (Tij +1-2*Tij+Tij-1)/hx 2
, где hx - шаг дискретизации по пространственной координате.

Тогда алгебраизированное уравнение теплопроводности для узла, принадлежащего i-ому временному и j-ому пространственному слоям, принимает следующий вид (gT=0):

(Ti+1j-Tij)/ht = aT*(Tij+1-2*Tij+Ti j-1)/h2x

Такой вид уравнения позволяет в явном виде выразить единственную неизвестную:

Ti+1j = aT*(Ti j+1-2*Tij+Ti j-1)*ht/h2x+Tij

Что, в свою очередь, дает возможность просто организовать вычислительный процесс в виде "цикл в цикле" (без деталей, связанных с граничными условиями различного рода):

for (i=0; i<n; i++)
for (j=0; j<m; j++)
Ti+1j = ...
Напомним, что значения T0j известны из начальных условий. Здесь n=Tкон/ht, m - количество пространственных слоёв расчетной сетки.
Требуется разработать параллельную программу, в которой каждый поток управления ответственнен за расчеты для "полосы" расчетной сетки шириной в m/N пространственных узлов, где N - число потоков управления. Программа должна демонстрировать ускорение по сравнению с последовательным вариантом. Предусмотреть визуализацию результатов посредством утилиты gnuplot.*/

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#define L
#define lambda
#define C
#define ro

typedef struct data {
    pthread_t ptid;
    int first;
    int last;
} ThreadRecord;

int done = 0;
pthread_barrier_t barr1, barr2;
double **T, **TT;
double ht = 0.1, hx = 1.0, at = 1.0;
int Ny, Nx;
int n;
FILE *f;

void *
mysolver(void *arg_p)
{
    int i, j;
    ThreadRecord *arg = (ThreadRecord *) arg_p;
    while(!done) {
        pthread_barrier_wait(&barr1);

        for (i = 0; i < Ny - 1; i++)
            for (j = arg->first; j < arg->last; j++)
                if (j != 0 && j != Nx - 1)
                    T[i + 1][j] = (at * (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) * ht / (hx * hx)) + T[i][j];

        pthread_barrier_wait(&barr2);
    }
}
int main(int argc, char *argv[])
{
    int i, j, k, m, t, h, p, q;
    int nt = atoi(argv[1]);
    double Tl = 0.0, Tr = 100.0, T0 = 0.0;
    double tmp = 0.0;
    double delimiter = 0.0;
    ThreadRecord *threads;
    pthread_attr_t pattr;
    struct timeval tim; 

//    at = lamda / (C * ro);
//    hx = L / (Nx - 1);
//    ht = 0.00000000001;
    if ( Nx % nt ) {
        fprintf(stderr, "Размерность должна быть кратна количеству потоков\n");
        exit (2);
    };
    Nx = atoi(argv[2]);
    Ny = atoi(argv[3]);
    T = (double **) malloc( sizeof(double*) * Ny);
    TT = (double **) malloc( sizeof(double*)* Ny);

    for(i = 0; i < Ny; i++)
        T[i] = (double *) malloc(sizeof(double) * Nx);
    for(i = 0; i < Ny; i++)
        TT[i]=(double *) malloc(sizeof(double) * Nx);

    for(i = 0; i < Ny; i++)
        for(j = 0; j < Nx; j++){
            T[i][j]=T0;
            TT[i][j]=T0;
        }
    for(i = 0, j = 0; i < Ny; i++ ){
        T[i][j] = Tl;
        TT[i][j] = Tl;
    }
    for(i = 0, j = Nx - 1; i < Ny; i++ ) {
        T[i][j] = Tr;
        TT[i][j] = Tr;
    }
    pthread_attr_init (&pattr);
    pthread_attr_setscope (&pattr, PTHREAD_SCOPE_SYSTEM);
    pthread_attr_setdetachstate (&pattr,PTHREAD_CREATE_JOINABLE);

    threads = (ThreadRecord *) malloc(nt * sizeof(ThreadRecord));

    pthread_barrier_init(&barr1, NULL, nt + 1);
    pthread_barrier_init(&barr2, NULL, nt + 1);

    j = Nx/nt;
    for (i = 0; i < nt; i++) {
        threads[i].first = tmp;
        threads[i].last = j + tmp;
        tmp += j;
        if (pthread_create(&(threads[i].ptid), NULL, mysolver, (void *) &threads[i]))
            perror("Error pthread create\n");
    }
	gettimeofday(&tim, NULL);  
        double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
    pthread_barrier_wait(&barr1);
    for (k = 0; k < nt; k++) {
        pthread_barrier_wait(&barr2);

        for (t = 0; t < Ny; t++)
            for (h = 0; h < Nx; h++)
                TT[t][h] = T[t][h];

        pthread_barrier_wait(&barr1);
    }
	char c[255];
    for (i = 0; i < Ny; i++) {
	sprintf(c, "step%d", i);
	f = fopen(c, "w+");
        for (j = 0; j < Nx; j++)
            fprintf(f, "%d %lf\n", j, T[i][j]);
    }
	fclose(f);
	gettimeofday(&tim, NULL);  
        double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	printf("Time= %lf\n",t2-t1);

	fclose(f);
    done = 1;
	free(T);
	free(TT);
    return 0;
}
