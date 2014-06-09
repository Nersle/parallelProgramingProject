/*
*
*名称:	cannonSerial.c
*功能:	cannon算法串行实现
*作者：	LH
*时间：2014-05-29
*要求：测试数据文件第一个数n为方程组的阶数，
*	   后续n*n个数为方程组Ax=b的A矩阵，接着n
*	   个数为B向量
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define IDEBUG 0	//调试信息输出控制

double  aa[1002][1002],bb[1002],xx[1002];
int  matrix_print_off (int nr, int nc, double A[1002][1002]);
int  vector_print_off (int nr, double x[1002]);
void gauss(double a[1002][1002], double b[1002], double x[1002], int n);


int main (void)
{
    float  aij,bi;
    int    i,j,n;
    FILE   *finput;

    finput = fopen("gauss.dat","r");
    if (finput == NULL)
    {
        printf("Data file gauss.dat not found\n");
        return(-1);
    }
    fscanf(finput, "%d",&n);
    printf("\nDimension of matrix = %d\n",n);
    for (i=1; i<=n; i++)
    {
        for (j=1; j<=n; j++)
        {
            fscanf(finput,"%f ",&aij);
            aa[i][j] = (double) aij;
        }
    }
    for (i=1; i<=n; i++)
    {
        fscanf(finput,"%f ",&bi);
        bb[i] = (double) bi;
    }
    fclose(finput);
	if(IDEBUG == 1)
	{
	   printf("\nMatrix A\n");
	   matrix_print_off (n,n,aa);
	   printf("\nVector b\n");
	   vector_print_off (n,bb);
	}

	clock_t start,end;//定义CPU时钟
	double time=0;
	start=clock();

    gauss(aa,bb,xx,n);

	end=clock();
	time=(double)(end-start)/CLOCKS_PER_SEC;

    printf("\nThis time gaussOpenMP.exe runtime=%lf\n\n",time);
	if(IDEBUG == 1)
	{
		vector_print_off (n,xx);
	}
    return(0);
}
/*
*功能：高斯消元解线性方程组主函数
*
*/
void gauss(double a[1002][1002], double b[1002], double x[1002], int n)
{
    int   i,j,k,m,rowx;
    double xfac,temp,temp1,amax;
    rowx = 0;
	#pragma omp parallel
    for (k=1; k<=n-1; ++k)
    {
        amax = (double) fabs(a[k][k]) ;
        m = k;
        for (i=k+1; i<=n; i++)
        {
            xfac = (double) fabs(a[i][k]);
            if(xfac > amax)
            {
                amax = xfac;
                m=i;
            }
        }
        if(m != k)
        {
            rowx = rowx+1;
            temp1 = b[k];
            b[k]  = b[m];
            b[m]  = temp1;
            for(j=k; j<=n; j++)
            {
                temp = a[k][j];
                a[k][j] = a[m][j];
                a[m][j] = temp;
            }
        }
    #pragma omp parallel for
    for (k=1; k<=n-1; ++k)
	{
        for (i=k+1; i<=n; ++i)
        {
            xfac = a[i][k]/a[k][k];

            for (j=k+1; j<=n; ++j)
            {
                a[i][j] = a[i][j]-xfac*a[k][j];
            }
            b[i] = b[i]-xfac*b[k];
        }
	}
        if(IDEBUG == 1)
        {
            printf("\n A after decomposition step %d\n\n",k);
            matrix_print_off (n, n, a);
        }
    }
    for (j=1; j<=n; ++j)
    {
        k=n-j+1;
        x[k] = b[k];
        for(i=k+1; i<=n; ++i)
        {
            x[k] = x[k]-a[k][i]*x[i];
        }
        x[k] = x[k]/a[k][k];
    }
    if(IDEBUG == 1)
		printf("\nNumber of row exchanges = %d\n",rowx);
}

/*
*功能：输出矩阵
*
*/
int matrix_print_off (int nr, int nc, double A[1002][1002])
{
    int i,j;
    if ( nr <= 0 ) return (-1);
    if ( nc <= 0 ) return (-2);
    for (i = 1; i <= nr; i++)
    {

        for (j = 1; j <= nc; j++)
        {
            printf ("%9.4f  ", A[i][j]);
        }
        printf("\n");
    }
    return (0);
}


/*
*功能：输出向量
*
*/
int vector_print_off (int nr, double x[])
{
    int i;
    if ( nr <= 0 ) return (-1);
    for (i = 1; i <= nr; i++)
    {
        printf ("%9.4f  \n", x[i]);
    }
    printf("\n");
    return (0);
}
