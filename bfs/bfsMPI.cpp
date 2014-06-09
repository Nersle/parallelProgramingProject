/*
*
*名称:	bfsMPI.c
*功能:	宽度优先搜索算法MPI实现
*作者：	LH
*时间： 2014-05-29
*
*/
#include <iostream>
#include <Queue>
#include <stdlib.h>

using namespace std;

#define N 400
#define DEBUG 0

int adjacency[N][N]={0};
/*
*功能：随机生成邻接矩阵，25%的概率两个顶点之间不存在一条边
*
*/
int myRand(int row, int col)
{
    if(row==col)
        return 0;
    int temp=rand()%4;
    if(temp<3)
        return 1;
    else
        return 0;
}
/*
*功能：随机生成邻接矩阵
*
*/
void initMap()
{
    int i,j,rands;
	#pragma omp parallel for
    for(i=0; i<N; i++)
		#pragma omp parallel for
        for(j=i; j<N; j++)
        {
            rands=myRand(i,j);
            adjacency[i][j]=rands;
            adjacency[j][i]=rands;
        }
	//调试状态下输出邻接矩阵
    if(DEBUG==1)
    {
        for(i=0; i<N; i++)
        {
            cout<<i<<" ";
            for(j=0; j<N; j++)
            {
                cout<<adjacency[i][j]<<" ";
            }
            cout<<endl;
        }
    }


}
/*
*功能：宽度优先搜索
*
*/
void boradSearch(int start)
{
    queue<int> nodeQueue;;
    int flag[N]={0};
    nodeQueue.push(start);
    flag[start]=1;
    int head,i   ;
    while(!nodeQueue.empty())
    {
        head=nodeQueue.front();
        nodeQueue.pop();
		#pragma omp parallel for
        for(i=0; i<N; i++)
        {
            if(adjacency[head][i]==1&&flag[i]==0)
            {
                flag[i]=1;
                nodeQueue.push(i);
            }
        }
        cout<<head<<" -> ";
    }
}

int main(int argc, char* argv)
{
    initMap();
    boradSearch(0);//从顶点0开始遍历随机图
    return 0;
}
