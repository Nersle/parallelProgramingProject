/*
*
*名称:	cannonMPI.c
*功能:	MPI实现的cannon算法
*作者：	LH
*时间：2014-05-29
*
*/
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

float **A, **B, **C;
float *a, *b, *c, *tempA, *tempB;
int n, locN, locN2,process, sqrtNum;
int myRank, myRow, myCol;
MPI_Status status;
double starttime,endtime;

int main(int argc, char *argv[])
{
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &process);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   starttime=MPI_Wtime();
   sqrtNum = sqrt(process);
   if (sqrtNum*sqrtNum != process)
   {
      if (myRank == 0)
	  printf("Number of processors is not a quadratic number!\n");
      MPI_Finalize();
      exit(1);
   }
   int i;
   n = atoi(argv[1]);
   locN  = n / sqrtNum;
   locN2 = locN * locN;
   myCol =  myRank % sqrtNum ;
   myRow = (myRank-myCol) / sqrtNum ;

   a = (float *)malloc( locN2 * sizeof(float) );
   b = (float *)malloc( locN2 * sizeof(float) );
   c = (float *)malloc( locN2 * sizeof(float) );

   for(i=0; i<locN2 ; i++)
     c[i] = 0.0;

   tempA = (float *)malloc( locN2 * sizeof(float) );
   tempB = (float *)malloc( locN2 * sizeof(float) );

   if (myRank == 0)
   {
      A = (float **)malloc( n * sizeof(float*) );
      B = (float **)malloc( n * sizeof(float*) );
      C = (float **)malloc( n * sizeof(float*) );

      for(i=0; i<n; i++)
      {
         A[i] = (float *)malloc( n * sizeof(float) );
         B[i] = (float *)malloc( n * sizeof(float) );
         C[i] = (float *)malloc( n * sizeof(float) );
      }
      randomAB();
	  starttime=MPI_Wtime();
      scatterAB();
   } else
   {
       MPI_Recv(a, locN2, MPI_FLOAT, 0 , 1, MPI_COMM_WORLD, &status);
       MPI_Recv(b, locN2, MPI_FLOAT, 0 , 2, MPI_COMM_WORLD, &status);
   }
   initAlignment();
   mainShift();
   if(myRank == 0)
   {
     togetherResult();
	 /*
	 *输出运算结果，数据量较大时注释该语句
     print(A,"random matrix A : \n");
	 print(B,"random matrix B : \n");
	 print(C,"Matrix C = A * B : \n");
	 */
   } else
   {
      MPI_Send(c,locN2,MPI_FLOAT,0,1,MPI_COMM_WORLD);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   endtime=MPI_Wtime();
   printf("Total time of %d is %lf\n",myRank,endtime-starttime);

   MPI_Finalize();
   return 0;
}

/*
*功能：处理器逻辑阵列坐标至rank号的转换 
*
*/
int getIndex(int row, int col, int sqrtNum)
{
   return ((row+sqrtNum)%sqrtNum)*sqrtNum + (col+sqrtNum)%sqrtNum;
}

/*
*功能：从文件中读取测试数据
*
*/
void randomABSave()
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
	/*随机生成A,B,并初始化C*/
    for(i=0; i<n ; i++)
      for(j=0; j<n ; j++)
	  {
		fscanf(finput,"%d",&input);
	    A[i][j] =input;
		fscanf(finput,"%d",&input);
        B[i][j] = input;
        C[i][j] = 0.0;
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

	/*随机生成A,B,并初始化C*/
    for(i=0; i<n ; i++)
      for(j=0; j<n ; j++)
	  {
	    A[i][j] = rand()%10;
        B[i][j] = rand()%10;
        C[i][j] = 0.0;
	  }
}


/*
*功能：:rank为0的处理器向其他处理器发送A、B矩阵的相关块 
*
*/
void scatterAB()
{
   int i,j,k,l;
   int p_imin,p_imax,p_jmin,p_jmax;
   
   for(k=0; k<p; k++)
   {
	  p_jmin = (k % sqrtNum    ) * locN;
  	  p_jmax = (k % sqrtNum + 1) * locN-1;
	  p_imin = (k - (k % sqrtNum))/sqrtNum * locN;
	  p_imax = ((k - (k % sqrtNum))/sqrtNum +1) *locN -1;
      l = 0;

      for(i=p_imin; i<=p_imax; i++)
      {
      	  for(j=p_jmin; j<=p_jmax; j++)
      	  {
              tempA[l] = A[i][j];
			  tempB[l] = B[i][j];
			  l++;
          }
      }

      if(k==0)
      {
         memcpy(a, tempA, locN2 * sizeof(float));
		 memcpy(b, tempB, locN2 * sizeof(float));
      } else
      {
          MPI_Send(tempA, locN2, MPI_FLOAT, k, 1, MPI_COMM_WORLD);
		  MPI_Send(tempB, locN2, MPI_FLOAT, k, 2, MPI_COMM_WORLD);
      }
   }
}

/*
*功能：接收初始测试数据
*
*/
void initAlignment()
{

   MPI_Sendrecv(a, locN2, MPI_FLOAT, getIndex(myRow,myCol-myRow,sqrtNum), 1,
            tempA, locN2, MPI_FLOAT, getIndex(myRow,myCol+myRow,sqrtNum), 1, MPI_COMM_WORLD, &status);
   memcpy(a, tempA, locN2 * sizeof(float) );

   MPI_Sendrecv(b, locN2, MPI_FLOAT, getIndex(myRow-myCol,myCol,sqrtNum), 1,
            tempB, locN2, MPI_FLOAT, getIndex(myRow+myCol,myCol,sqrtNum), 1, MPI_COMM_WORLD, &status);
   memcpy(b, tempB, locN2 * sizeof(float) );
}

/*
*功能：分块矩阵左移和上移，并计算分块c
*
*/
void mainShift()
{
   int i,j,k,l;

   for(l=0; l<sqrtNum; l++)
   {

     for(i=0; i<locN; i++)
       for(j=0; j<locN; j++)
         for(k=0; k<locN; k++)
           c[i*locN+j] += a[i*locN+k]*b[k*locN+j];

      MPI_Send(a , locN2, MPI_FLOAT, getIndex(myRow, myCol-1, sqrtNum), 1, MPI_COMM_WORLD);
      MPI_Recv(a , locN2, MPI_FLOAT, getIndex(myRow, myCol+1, sqrtNum), 1, MPI_COMM_WORLD, &status);
      MPI_Send(b , locN2, MPI_FLOAT, getIndex(myRow-1, myCol, sqrtNum), 1, MPI_COMM_WORLD);
      MPI_Recv(b , locN2, MPI_FLOAT, getIndex(myRow+1, myCol, sqrtNum), 1, MPI_COMM_WORLD, &status);
   }
}

/*
*功能：主处理器接收各个处理器的运算结果并整合
*
*/
void togetherResult()
{
   int i,j,i2,j2,k;
   int p_imin,p_imax,p_jmin,p_jmax;

   for (i=0;i<locN;i++)
	 for(j=0;j<locN;j++)
	   C[i][j]=c[i*locN+j];

   for (k=1;k<p;k++)
   {
       MPI_Recv(c, locN2, MPI_FLOAT, k, 1, MPI_COMM_WORLD, &status);

       p_jmin = (k % sqrtNum    ) *locN;
       p_jmax = (k % sqrtNum + 1) *locN-1;
       p_imin =  (k - (k % sqrtNum))/sqrtNum     *locN;
       p_imax = ((k - (k % sqrtNum))/sqrtNum +1) *locN -1;

       i2=0;

       for(i=p_imin; i<=p_imax; i++)
       {
           j2=0;
           for(j=p_jmin; j<=p_jmax; j++)
           {
               C[i][j]=c[i2*locN+j2];
               j2++;
           }
           i2++;
       }
   }
}


/*
*功能：打印矩阵
*
*/
void print(float **m,char *str)
{
   int i,j;
   printf("%s",str);
   for(i=0;i<n;i++)
   {
       for(j=0;j<n;j++)
           printf("%15.0f    ",m[i][j]);
       printf("\n");
   }
   printf("\n");
}


