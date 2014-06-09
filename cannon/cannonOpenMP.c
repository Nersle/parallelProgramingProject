/*
*
*名称:	cannonMPI.c
*功能:	OpenMP实现的cannon算法
*作者：	LH
*时间：2014-05-29
*
*/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

 
int main (int argc, char *argv[]) 
{
	int matrixScale = atoi(argv[1]); 
    int	tid, nthreads, i, j, k, chunk;
    double	matrixA[matrixScale][matrixScale], matrixB[matrixScale][matrixScale],matrixC[matrixScale][matrixScale];          
    double start, end;
	chunk = 10;                  

	#pragma omp parallel shared(matrixA,matrixB,matrixC,nthreads,chunk) private(tid,i,j,k)
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Starting matrix multiple example with %d threads\n",nthreads);
            printf("Initializing matrices...\n");
        }
		
		#pragma omp for schedule (static, chunk) 
        for (i=0; i<matrixScale; i++)
            for (j=0; j<matrixScale; j++)
                matrixA[i][j]= i+j;
		#pragma omp for schedule (static, chunk)
        for (i=0; i<matrixScale; i++)
            for (j=0; j<matrixScale; j++)
                matrixB[i][j]= i*j;
		#pragma omp for schedule (static, chunk)
        for (i=0; i<matrixScale; i++)
            for (j=0; j<matrixScale; j++)
                matrixC[i][j]= 0;

		start = omp_get_wtime();
        printf("Thread %d starting matrix multiply...\n",tid);
		#pragma omp for schedule (static, chunk)
        for (i=0; i<matrixScale; i++)    
        {
            printf("Thread=%d did row=%d\n",tid,i);
			#pragma  omp for
            for(j=0; j<matrixScale; j++) 	
				#pragma  omp for      
                for (k=0; k<matrixScale; k++)
                    matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
        }
    }  
		end =omp_get_wtime();	
	/*
    printf("Result Matrix:\n");
    for (i=0; i<matrixScale; i++)
    {
        for (j=0; j<matrixScale; j++) 
            printf("%6.2f   ", matrixC[i][j]);
        printf("\n"); 
    }
	*/
	printf("Total time of running cannonOpenMP.exe is %lf\n\n", end-start ); 

}
