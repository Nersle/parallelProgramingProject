#-*- coding:utf-8 -*-
#
#  本程序主要用于计算加速比和CPU效率
#
#


def calc(a,b,c):
  return float(a)/(b*c)

#x=[12,78,132, 823, 1852 ,12897 ,32323 ]
#x=[3,44,390,5270,96380]
#x=[12,68,289,1032]
x=[13,28,138,289,1489]
print x
#y=[35,122,186,938,1753,12596,27001]
#y=[6,54,352,4321,73429]
#y=[19,95,258,815]
y=[18,37,125,254,1162]

#z=[7,52,402,5417,107382]
#z=[27,97,172,890,1880,13128,35992]
#z=[16,87,332,1329]
z=[16,32,146,316,1728]
lens=len(x)
i=0
print "x1:"
while i<lens:
  print float(x[i])/y[i]," "
  i=i+1
print "x2:"
i=0
while i<lens:
  print float(x[i])/z[i]," "
  i=i+1
 

#chuan=[12897,12897,12897]
#chuan=[289,289,289]
chuan=[1489,1489,1489]
#mpi=[13632,14443,12596]
#mpi=[318,268,258]
mpi=[1338,1290,1162]
#openmp=[13871,13766,13128]
#openmp=[318,351,332]
openmp=[1691,1756,1728]
lens=len(chuan)
print "MPI"
i=0
while i<lens:
  print calc(chuan[i],2**i,mpi[i])
  i=i+1
print "OpenMP"
i=0
while i<lens:
  print calc(chuan[i],2**i,openmp[i])
  i=i+1
