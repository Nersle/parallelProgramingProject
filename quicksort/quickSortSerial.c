/*
*
*名称:	quickSortSerial.c
*功能:	quickSort算法串行实现
*作者：	LH
*时间：2014-05-29
*
*/
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>

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
void qsortMe(double l[],int low,int high)
{
    int prvotloc;
    if(low<high)
    {
        prvotloc=partions(l,low,high);    //将第一次排序的结果作为枢轴
        qsortMe(l,low,prvotloc-1); //递归调用排序 由low 到prvotloc-1
        qsortMe(l,prvotloc+1,high); //递归调用排序 由 prvotloc+1到 high
    }
}

void quicksort(double l[],int n)
{
    qsortMe(l,0,n-1); //第一个作为枢轴 ，从第一个排到第n个
}

void main()
{
    FILE   *finput;
    finput = fopen("sortData.dat","r");
    if (finput == NULL)
    {
        printf("Data file gauss.dat not found\n");
        return(-1);
    }
    int n;
    fscanf(finput, "%d",&n);
    printf("sacle=%dx\n",n);
    double *ax;
    int j;
    ax=malloc(n*sizeof(double));
    double input;
    for (j=0; j<n; j++)
    {
        fscanf(finput,"%d",&input);
        ax[j] = input;
    }

	clock_t start,end; 
	double time=0;
	start=clock();
    printf("\n");

    quicksort(ax,n);

    end=clock();
	time=(double)(end-start)/CLOCKS_PER_SEC;

	printf("\nThis time quickSortSerial.exe runtime=%lf\n\n",time);
}
