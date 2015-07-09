/*MPI.  Разработать средствами MPI параллельную программу решения двухмерной нестационарной краевой задачи методом конечных разностей с использованием неявной вычислительной схемы. Объект моделирования - прямоугольная пластина постоянной толщины. Подробности постановки подобной задачи даны ниже. Возможны граничные условия первого и второго рода в различных узлах расчетной сетки. Временной интервал моделирования и количество узлов расчетной сетки - параметры программы. Программа должна демонстрировать ускорение по сравнению с последовательным вариантом. Предусмотреть визуализацию результатов посредством утилиты gnuplot.
 Примечание. Использование неявной вычислительной схемы связано с решением системы ЛАУ на каждом временном слое. Следовательно, распараллеливание решения краевой задачи неявным методом сводится к распараллеливанию решения системы ЛАУ методом Гаусса. Распараллеливание же метода Гаусса проще всего реализуется в ситуации, когда матрица коэффициентов системы имеет блочно-диагональный с окаймлением вид. Матрица же коэффициентов автоматически получит такой вид, если нумерацию узлов расчетной сетки вести по такой простой схеме: сначала "внутренние" узлы всех фрагментов, на которые разбивается стержень, а в последнюю очередь - "соединительные" узлы. Кстати, сказанное выше справедливо и для метода прогонки (поскольку этот метод - частный случай метода Гаусса).*/
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>
#define _REENTRANT

#define N_PER_PROC 100	// Кол-во строк матрицы на процесс
#define PERVOGO_RODA_TOP 30.0
#define PERVOGO_RODA_BOTTOM 10.0
#define PERVOGO_RODA_LEFT 20.0
#define COUNT_TRHEADS 3 //Заглушка




int main(int argc, char **argv) {
    
    int myrank, total;
    if (argc != 5){
        write(1,"Error\n",7);
        exit(-1);
    }
    double *A, *B, *T;
    double *a, *b, *t;  
    int *smesh;
    int *smeshForB;
    int countPoints;
    int *lastBlock;
    int *sizeBlock;
    int *sizeBlockForB;
    int intBuf[1];
    int i, j;
    int startPosition;
    int endPosition;
    int m;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &total);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    lastBlock=(int*)malloc(sizeof(int)*total);
    smesh=(int*)malloc(sizeof(int)*total);
    smeshForB=(int*)malloc(sizeof(int)*total);
    sizeBlock=(int*)malloc(sizeof(int)*total);
    sizeBlockForB=(int*)malloc(sizeof(int)*total);
    //printf ("Total=%d, rank=%d\n", total, myrank);

    struct timeval tim;
    double t1, t2;
    if(!myrank){
        int countX, countY;
        double deltaX, deltaY;
        
        countX=atoi(argv[1]);
        deltaX=atof(argv[2]);

        countY=atoi(argv[3]);
        deltaY=atof(argv[4]);

     
        countPoints=countX*countY;
       
        A=(double*)malloc(sizeof(double)*countPoints*countPoints);
        B=(double*)malloc(sizeof(double)*countPoints);
        for(i=0; i < countPoints; ++i){
            B[i]=0.0;
        }
        T=(double*)malloc(sizeof(double)*countPoints);
        for(i=0; i < countPoints; ++i){
            T[i]=0.0;
        }
        int tempCountPoints = countPoints - (total-1)*countY;
        int tempCountInOneBlock = tempCountPoints/total;
        
        j=0;
        while(1){                                                      
            smesh[j]=j*tempCountInOneBlock*countPoints;
            smeshForB[j]=j*tempCountInOneBlock;
            if(j+1 == total || tempCountPoints <= tempCountInOneBlock){
                lastBlock[j] = j*tempCountInOneBlock + tempCountPoints - 1; 
                sizeBlock[j] = (lastBlock[j]+1)*countPoints - smesh[j];
                sizeBlockForB[j] = sizeBlock[j]/countPoints;
                break;
            }
            lastBlock[j] = j*tempCountInOneBlock + tempCountInOneBlock - 1; 
            sizeBlock[j] = (lastBlock[j]+1)*countPoints - smesh[j];
            sizeBlockForB[j] = lastBlock[j]-smeshForB[j]+1;
            tempCountPoints -= tempCountInOneBlock;
            j++;
        }

        for (i = 0; i < countPoints; ++i) {
            if(i%countY == 0){
                B[i]=PERVOGO_RODA_TOP;
                A[i*countPoints+i]=1.0;
                continue;
            }
            if(i%countY == (countY-1)){
                B[i]=PERVOGO_RODA_BOTTOM;
                A[i*countPoints+i]=1.0;
                continue;
            } 
            if (i < countY){
                B[i]=PERVOGO_RODA_LEFT;
                 A[i*countPoints+i]=1.0;
                continue;
            } 
            B[i]=0.0;
            if( i > (countPoints-total*countY) && i < (countPoints-(total-1)*countY)){
                A[i*countPoints+i]=(-1.0)*((1.0/deltaX)+1.0);
                A[i*countPoints+i-countY]=(1.0)/(deltaX);
                continue;
            }
            int currentBlock;
            for( currentBlock = 0; currentBlock <= total; currentBlock++){
                if (currentBlock == total)
                    break;
                if(i < lastBlock[currentBlock]){
                    break;
                }
            }
            A[i*countPoints+i-1]=((1.0)/(deltaX*deltaX));   //Заполнение по Y
            A[i*countPoints+i]=((-2.0)/(deltaX*deltaX))+((-2.0)/(deltaY*deltaY));
            A[i*countPoints+i+1]=((1.0)/(deltaX*deltaX));
            if( i > lastBlock[total-1] ){            //Игнорируем окаймление
                continue;
            }
            if(currentBlock == 0){
                A[i*countPoints+i-countY]=((1.0)/(deltaY*deltaY));
                /* Заполнения окаймления */
                if(i > (lastBlock[currentBlock]-countY)){
                    A[i*countPoints+lastBlock[total-1]+1+i%countY]=((1.0)/(deltaY*deltaY));
                    A[(lastBlock[total-1]+1+i%countY)*countPoints+i]=((1.0)/(deltaY*deltaY));
                } else { 
                    A[i*countPoints+i+countY]=((1.0)/(deltaY*deltaY));
                }
                continue;
            }
            if (currentBlock == (total-1)){
                A[i*countPoints+i+countY]=((1.0)/(deltaY*deltaY));
                if ( i > lastBlock[currentBlock-1] && i < (lastBlock[currentBlock-1]+countY)){
                    A[i*countPoints+lastBlock[total-1]+1+countY*(total-2)+i%countY]=((1.0)/(deltaY*deltaY));
                    A[(lastBlock[total-1]+1+countY*(total-2)+i%countY)*countPoints+i]=((1.0)/(deltaY*deltaY));
                } else {
                    A[i*countPoints+i-countY]=((1.0)/(deltaY*deltaY));
                }
                continue;
            }
            if ( i > lastBlock[currentBlock-1] && i < (lastBlock[currentBlock-1]+countY)){      //Первый столбец блока
                A[i*countPoints+lastBlock[total-1]+1+countY*(currentBlock-1)+i%countY]=((1.0)/(deltaY*deltaY));
                A[(lastBlock[total-1]+1+countY*(currentBlock-1)+i%countY)*countPoints+i]=((1.0)/(deltaY*deltaY));
            } else {
                A[i*countPoints+i-countY]=((1.0)/(deltaY*deltaY));
            }
            if(i > (lastBlock[currentBlock]-countY)){                                           //Последний столбец блока
                A[i*countPoints+lastBlock[total-1]+1+countY*(currentBlock)+i%countY]=((1.0)/(deltaY*deltaY));
                A[(lastBlock[total-1]+1+countY*(currentBlock)+i%countY)*countPoints+i]=((1.0)/(deltaY*deltaY));
            } else{
                A[i*countPoints+i+countY]=((1.0)/(deltaY*deltaY));
            }
        }
    }

    if (!myrank) {  
        intBuf[0] = countPoints;
    }




    MPI_Bcast((void *)intBuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    countPoints=intBuf[0];
    MPI_Bcast((void *)lastBlock, total, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)smesh, total, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)sizeBlock, total, MPI_INT, 0, MPI_COMM_WORLD);
    

    
    m = sizeBlock[myrank]/countPoints;
    //printf("%d\n",countPoints*m);

    a = (double*)malloc(sizeof(double)*countPoints*m);
    b=(double *) malloc (sizeof(double)*countPoints);

    MPI_Scatterv((void *)A,sizeBlock,smesh,MPI_DOUBLE,
          (void *)a, countPoints*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv((void *)B,sizeBlockForB,smeshForB, MPI_DOUBLE,
        (void*) b, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int startBlockIndex = lastBlock[myrank]-m+1;
    if(!myrank){
        gettimeofday(&tim, NULL);  
        t1=tim.tv_sec+(tim.tv_usec/1000000.0);
    }
    for(j = 0; j < m; ++j){
        for(i = startBlockIndex; i <= lastBlock[myrank]; ++i){
            //printf("%1.0lf|",a[j*countPoints+i]);
            if(j != (i-startBlockIndex)) {
                double temp = a[(i-startBlockIndex)*countPoints+j+startBlockIndex]/a[j*countPoints+j+startBlockIndex];
                int k;
                for(k = startBlockIndex; k <= lastBlock[myrank]; ++k){
                    a[(i-startBlockIndex)*countPoints+k]-=a[(j*countPoints+k)]*temp;
                }
                for(k = (lastBlock[total-1]+1); k < countPoints; ++k){
                        //printf("i=%d, j=%d, k=%d, start=%d\n", i, j,k, startBlockIndex );
                        a[(i-startBlockIndex)*countPoints+k]-=a[(j*countPoints+k)]*temp;
                }
                b[(i-startBlockIndex)]-=b[(j)]*temp;
            }
        }
        //printf("\n");
    }
    MPI_Gatherv((void *)a, m*countPoints, MPI_DOUBLE,(void *)A,
        sizeBlock,smesh,MPI_DOUBLE,0, MPI_COMM_WORLD);
     MPI_Gatherv((void *)b, m, MPI_DOUBLE,(void *)B,
        sizeBlockForB,smeshForB,MPI_DOUBLE,0, MPI_COMM_WORLD);
    if(!myrank){
        for(i = (lastBlock[total-1]+1); i < countPoints; ++i){
            for(j = 0; j < i; ++j){
                if (A[i*countPoints+j] == 0.0)
                    continue;
                double temp = A[i*countPoints+j]/A[j*countPoints+j];
                int k;
                for(k = j; k < countPoints; ++k){
                    A[(i*countPoints+k)]-=A[(j*countPoints+k)]*temp;
                }
                B[i]-=B[j]*temp;
            }
        }
       
        for( i = (countPoints-1); i >= 0; i--){
            double temp=0.0;
            for( j = (countPoints-1); j >= 0; j--){
                if(i != j){
                    temp+=T[j]*A[i*countPoints+j];
                }
                T[i]=(B[i]-temp)/(A[i*countPoints+i]);
            }

	}
	gettimeofday(&tim, NULL);  
       t2=tim.tv_sec+(tim.tv_usec/1000000.0);
	printf("Time= %lf\n",t2-t1);
       // for (i = 0; i < countPoints; ++i){
           // printf("T%d=%lf\n", i, T[i]);
       // }
    }
    MPI_Finalize();
    exit(0);
}
