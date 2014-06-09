/*
*
*名称:	quickSortOpenMP.c
*功能:	quickSort算法OpenMP实现
*作者：	LH
*时间：2014-05-29
*
*/
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <malloc.h>

/*
*功能：输出数组
*
*/
void print(double *x,int n)
{
    int i=0;
    for(;i<n;i++)
    {
        printf("%d ",x[i]);
    }
    printf("\n");
}

/*
*功能：根据枢纽元素划分数组
*
*/
int partions(double l[],int low,int high)
{
    //scanf("%d",&x);
    int middle=low;
    double flag=l[middle];
    while((low<high))
    {
        while(low<high&&l[high]>flag)
            high--;
        while(low<high&&l[low]<=flag)
            low++;
        if(low<high)
        {
            double temp=l[low];
            l[low]=l[high];
            l[high]=temp;
        }
    }
    double tp2=l[high];
    l[high]=l[middle];
    l[middle]=tp2;
    return low;
}

/*
*功能：递归快速排序
*
*/
void qsort(double l[],int low,int high)
{
    int prvotloc;
    if(low<high)
    {
        prvotloc=partions(l,low,high);    //将第一次排序的结果作为枢轴
        #pragma omp parallel 
	    //printf("%d\\%d run!\n",omp_get_thread_num(),omp_get_num_threads());
          {
          #pragma omp sections nowait
            {
                #pragma omp section//section
                 qsort(l,low,prvotloc-1); //递归调用排序 由low 到prvotloc-1
                #pragma omp section // section
                 qsort(l,prvotloc+1,high); //递归调用排序 由 prvotloc+1到 high

            }

          }

    }
}

void quicksort(double l[],int n)
{
    qsort(l,0,n-1); //第一个作为枢轴 ，从第一个排到第n个
}

int main()
{

    FILE   *finput;
    finput = fopen("sortData.dat","r");
    if (finput == NULL)
    {
        printf("Data file sortData.dat not found\n");
        return(-1);
    }
    int n;
    fscanf(finput, "%d",&n);
    printf("scale=%dx\n",n);
    double *ax;
    int j;
    ax=malloc(n*sizeof(double));
    double input;
    for (j=0; j<n; j++)
    {
        fscanf(finput,"%d",&input);
        ax[j] = input;
    }

    int b,c; 

	clock_t start,end; 
	double time=0;
	start=omp_get_wtime();
    printf("\n");

    quicksort(ax,n);

    end=omp_get_wtime();
	time=(double)(end-start) ;

	printf("\nThis time quickSortOMP.exe runtime=%lf\n\n",time);

}
