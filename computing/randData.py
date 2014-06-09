#-*- coding:utf-8 -*-
#
#  本程序主要用于生成测试数据集
#
#

import random

def getGaussData():
    f=open("gauss.dat","w")
    x=input("__Gauss Num__")
    num=int(x)
    f.write("Gauss solution!\n")
    f.write(x+"\n")
    for  i in range (num):
        for j in range (num):
            p=int(random.random()*50)
            f.write(str(p)+"   ")
        f.write("\n")

    for i in range (num):
        p=int(random.random()*100)
        f.write(str(p)+"\n")
    f.close()
def getSortData():
    f=open("sortData.dat","w")
    x=input("_Sort Num_")
    num=int(x)
    f.write(str(x)+"\n")
    for i in range (num):
        p=int(random.random()*100000)
        f.write(str(p)+"   ")
    f.close()


if __name__=="__main__":
    #getGaussData()
    getSortData()
