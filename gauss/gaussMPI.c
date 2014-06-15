/*
*
*名称:	cannonMPI.c
*功能:	cannon算法MPI实现
*作者：	LH
*时间：2014-05-29
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

 

int M,N,m;
float *A,*B;
double starttime,time1,time2;
int my_rank,p,l;
MPI_Status status;


int main(int argc, char **argv)
{
    int i,j,t,k,my_rank,group_size;
    int i1,i2,v,w,tem;
    float temp,lmax;
    float *sum,*f;
    float *a,*b,*x;
    int *shift;
    FILE *fdA,*fdB;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&group_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    p=group_size;

    if (my_rank==0)
    {
        starttime=MPI_Wtime(); 
        fdA=fopen("Gauss.data","r");
        fscanf(fdA,"%d %d", &M, &N);
        if (M != N-1)
        {
            printf("the input is wrong\n");
            exit(1);
        }

        A=(float *)malloc(sizeof(float)*M*M);
        B=(float *)malloc(sizeof(float)*M);

        for(i = 0; i < M; i++)
        {
            for(j = 0; j < M; j++)
            {
                fscanf(fdA,"%f", A+i*M+j);
            }
            fscanf(fdA,"%f", B+i);
        }
        fclose(fdA);
    }

    MPI_Bcast(&M,1,MPI_INT,0,MPI_COMM_WORLD);    
    m=M/p;
    if (M%p!=0) m++;

    f=(float*)malloc(sizeof(float)*(M+1));        
    a=(float*)malloc(sizeof(float)*m*M);          
    b=(float*)malloc(sizeof(float)*m);            
    sum=(float*)malloc(sizeof(float)*m);
    x=(float*)malloc(sizeof(float)*M);
    shift=(int*)malloc(sizeof(int)*M);

    if (a==NULL||b==NULL||f==NULL||sum==NULL||x==NULL||shift==NULL)
	{
	    printf("allocate error\n");
		exit(1);
	} 

    for(i=0;i<M;i++)
        shift[i]=i;


    if (my_rank==0)
    {
        for(i=0;i<m;i++)
            for(j=0;j<M;j++)
                a[i*M+j]=A[i*p*M+j];

        for(i=0;i<m;i++)
            b[i]=B[i*p];
    }

    if (my_rank==0)
    {
        for(i=0;i<M;i++)
            if ((i%p)!=0)
        {
            i1=i%p;
            i2=i/p+1;

            MPI_Send(&A[i*M],M,MPI_FLOAT,i1,i2,MPI_COMM_WORLD);
            MPI_Send(&B[i],1,MPI_FLOAT,i1,i2,MPI_COMM_WORLD);
        }
    }                                             /*  my_rank==0 */
    else                                          /*  my_rank !=0 */
    {
        for(i=0;i<m;i++)
        {
            MPI_Recv(&a[i*M],M,MPI_FLOAT,0,i+1,MPI_COMM_WORLD,&status);
            MPI_Recv(&b[i],1,MPI_FLOAT,0,i+1,MPI_COMM_WORLD,&status);
        }
    }

    time1=MPI_Wtime();                           

    for(i=0;i<m;i++)                            
        for(j=0;j<p;j++)
    {
        if (my_rank==j)                         
        {
            v=i*p+j;                           
            lmax=a[i*M+v];
            l=v;

            for(k=v+1;k<M;k++)                   
                if (fabs(a[i*M+k])>lmax)
            {
                lmax=a[i*M+k];
                l=k;
            }

            if (l!=v)                           
            {
                for(t=0;t<m;t++)
                {
                    temp=a[t*M+v];
                    a[t*M+v]=a[t*M+l];
                    a[t*M+l]=temp;
                }

                tem=shift[v];
                shift[v]=shift[l];
                shift[l]=tem;
            }

            for(k=v+1;k<M;k++)                 
                a[i*M+k]=a[i*M+k]/a[i*M+v];

            b[i]=b[i]/a[i*M+v];
            a[i*M+v]=1;

            for(k=v+1;k<M;k++)
                f[k]=a[i*M+k];
            f[M]=b[i];
            MPI_Bcast(&f[0],M+1,MPI_FLOAT,my_rank,MPI_COMM_WORLD); 
			MPI_Bcast(&l,1,MPI_INT,my_rank,MPI_COMM_WORLD);
        }
        else
        {
            v=i*p+j;
            MPI_Bcast(&f[0],M+1,MPI_FLOAT,j,MPI_COMM_WORLD);
            MPI_Bcast(&l,1,MPI_INT,j,MPI_COMM_WORLD);

            if (l!=v)
            {
                for(t=0;t<m;t++)
                {
                    temp=a[t*M+v];
                    a[t*M+v]=a[t*M+l];
                    a[t*M+l]=temp;
                }

                tem=shift[v];
                shift[v]=shift[l];
                shift[l]=tem;
            }
        }

        if (my_rank<=j)
            for(k=i+1;k<m;k++)
        {
            for(w=v+1;w<M;w++)
                a[k*M+w]=a[k*M+w]-f[w]*a[k*M+v];
            b[k]=b[k]-f[M]*a[k*M+v];
        }

        if (my_rank>j)
            for(k=i;k<m;k++)
        {
            for(w=v+1;w<M;w++)
                a[k*M+w]=a[k*M+w]-f[w]*a[k*M+v];
            b[k]=b[k]-f[M]*a[k*M+v];
        }
    }                                             

    for(i=0;i<m;i++)
        sum[i]=0.0;
    for(i=m-1;i>=0;i--)                            
        for(j=p-1;j>=0;j--)
            if (my_rank==j)
            {
                x[i*p+j]=(b[i]-sum[i])/a[i*M+i*p+j];

                MPI_Bcast(&x[i*p+j],1,MPI_FLOAT,my_rank,MPI_COMM_WORLD);

                for(k=0;k<i;k++)
                    sum[k]=sum[k]+a[k*M+i*p+j]*x[i*p+j];
            }
            else
            {
        MPI_Bcast(&x[i*p+j],1,MPI_FLOAT,j,MPI_COMM_WORLD);

        if (my_rank>j)
            for(k=0;k<i;k++)
                sum[k]=sum[k]+a[k*M+i*p+j]*x[i*p+j];

        if (my_rank<j)
            for(k=0;k<=i;k++)
                sum[k]=sum[k]+a[k*M+i*p+j]*x[i*p+j];
    }

    if (my_rank!=0)
        for(i=0;i<m;i++)
            MPI_Send(&x[i*p+my_rank],1,MPI_FLOAT,0,i,MPI_COMM_WORLD);
    else
        for(i=1;i<p;i++)
            for(j=0;j<m;j++)
                MPI_Recv(&x[j*p+i],1,MPI_FLOAT,i,j,MPI_COMM_WORLD,&status);

    if (my_rank==0)
    {
        printf("Input of file \"dataIn.txt\"\n");
        printf("%d\t%d\n", M, N);
        for(i=0;i<M;i++)
        {
            for(j=0;j<M;j++) printf("%f\t",A[i*M+j]);
            printf("%f\n",B[i]);
        }
        printf("\nOutput of solution\n");
        for(k=0;k<M;k++)
        {
            for(i=0;i<M;i++)
            {
                if (shift[i]==k) printf("x[%d]=%f\n",k,x[i]);
            }
        }
    }

    time2=MPI_Wtime();

    if (my_rank==0)
    { 
		printf("\nThis time gaussMPI.exe runtime=%lf\n\n",time2-time1); 
    }

    MPI_Finalize(); 
	free(a);
    free(b);
    free(x);
    free(f);
    return(0);
}