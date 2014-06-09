/*
*
*名称:	cannonSerial.c
*功能:	cannon算法串行实现
*作者：	LH
*时间：2014-05-29
*
*/
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

float **matrixA, **matrixB, **matrixC;             
int matrixScale ,myRank;
double starttime,endtime,alltime;

int main(int argc, char *argv[])
{ 
    int i; 
    matrixScale = atoi(argv[1]);  
    clock_t start,end ;//定义CPU时钟
    double time=0;
    start=clock();
    printf("\n");
    alltime=0;

    matrixA = (float **)malloc( matrixScale * sizeof(float*) );
    matrixB = (float **)malloc( matrixScale * sizeof(float*) );
    matrixC = (float **)malloc( matrixScale * sizeof(float*) );

    for(i=0; i<matrixScale; i++)
    {
        matrixA[i] = (float *)malloc( matrixScale * sizeof(float) );
        matrixB[i] = (float *)malloc( matrixScale * sizeof(float) );
        matrixC[i] = (float *)malloc( matrixScale * sizeof(float) );
    }

    randomAB(); 

    start=clock();	
    mainShift();    
    end=clock();
    time=(double)(end-start)/CLOCKS_PER_SEC; 
    printf("\nThis time cannonSerial.exe runtime=%lf\n\n",time);             
    return 0;
}


/*
*功能：从文件中读取测试数据
*
*/
void randomABsave()
{
    int i,j;   
    FILE   *finput;
    finput = fopen("cannonData.dat","r");
    if (finput == NULL)
    {
        printf("Data file cannonData.dat not found\n");
        return(-1);
    }	
    int input;
    /*随机生成matrixA,matrixB,并初始化matrixC*/
    for(i=0; i<matrixScale ; i++)
        for(j=0; j<matrixScale ; j++)
        {
            fscanf(finput,"%d",&input); 
            matrixA[i][j] =input;
            fscanf(finput,"%d",&input); 
            matrixB[i][j] = input;
            matrixC[i][j] = 0.0;
        } 
}

/*
*功能：临时生成测试数据
*
*/
void randomAB()
{
    int i,j;
    srand((unsigned int)time(NULL));   
    /*随机生成matrixA,matrixB,并初始化matrixC*/
    for(i=0; i<matrixScale ; i++)
        for(j=0; j<matrixScale ; j++)
        {
            matrixA[i][j] = rand()%10;
            matrixB[i][j] = rand()%10;
            matrixC[i][j] = 0.0;
        }
}

/*
*功能：分块矩阵左移和上移，并计算分块c
*
*/
void mainShift()
{
    int i,j,k; 
    for(i=0; i<matrixScale; i++) 
        for(j=0; j<matrixScale; j++) 
            for(k=0; k<matrixScale; k++)
                matrixC[i][j] += matrixA[i][k]*matrixB[k][j];

}

void print(float **m,char *str)
{
    int i,j;
    printf("%s",str);
    /*打印矩阵m*/
    for(i=0;i<matrixScale;i++)
    {
        for(j=0;j<matrixScale;j++)
            printf("%15.0f    ",m[i][j]);
        printf("\n");
    }
    printf("\n");
}


