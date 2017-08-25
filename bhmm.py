# -*- coding: utf-8 -*-

"""
  @author: Necva Bolucu (@necvabolucu)
  @author: Salih Tuc (@salihtuc)
"""

from collections import defaultdict
from math import exp
class BHMM(object):
    
    def __init__(self):
        self.frequencyTrans = defaultdict(float)
        self.frequencyEmit = defaultdict(float)
        self.frequencyTransExp = defaultdict(float)
        self.frequencyEmitExp = defaultdict(float)
        
        self.minTrans=float(-8.390173385616183e-07)
        self.minEmit=float(-2.765386569039885e-10)
    def __change_count(self,matrix, x, y, i,a):
#        print(x,y)
        matrix["%s|%s" % (x, y)] += (i +a)
        matrix["%s" % x] +=(i +a)
#        print( matrix["%s|%s" % (x, y)], matrix["%s" % x])
    
    def __change_countExp(self,matrix, x, y, a,i):
    #        print(x,y)
            matrix["%s|%s" % (x, y)] += exp(i-a+0.0001)
            matrix["%s" % x] +=exp( i -a+0.0001 )
                  
    def __get_value(self,matrix, x,y):
#        print(x, y)
        return matrix["%s|%s" % (x, y)]/matrix["%s" % x]
    
    def __get_valueExp(self,matrix, x,y):
#        print(x, y)
        return exp(matrix["%s|%s" % (x, y)]-matrix["%s" % x])
        
    def run(self):
        count=1
        file=open("result.txt")
        file1=open("out.txt","w")
        
        for line in file.readlines():
            
            splitted=line.strip("\n").strip().split()
            splittedValue=splitted[0].split("-",1)

            if count<170:
#            print(splittedValue)
                self.__change_count(self.frequencyTrans, splittedValue[0], splittedValue[1], float(splitted[1]),(-1*self.minTrans+0.001))
#                self.__change_countExp(self.frequencyTransExp, splittedValue[0], splittedValue[1],self.minTrans,float(splitted[1]))
            else:
                self.__change_count(self.frequencyEmit, splittedValue[0], splittedValue[1],float(splitted[1]),(-1*self.minEmit+0.001))
#                self.__change_countExp(self.frequencyEmitExp, splittedValue[0], splittedValue[1],self.minEmit,float(splitted[1]))
            count=count+1

        print(min([y for x,y in self.frequencyTrans.items() if "|" in x]),max([y for x,y in self.frequencyTrans.items() if "|" in x]))
        print(min([y for x,y in self.frequencyEmit.items() if "|" in x]),max([y for x,y in self.frequencyEmit.items() if "|" in x]))

#        print(len(self.frequency.keys()))
        for key in self.frequencyTrans.keys():
            if "|" in key:
                splittedValue=key.split("|")
#                print(key)
#                print(key,self.__get_value(self.frequency,splittedValue[0],splittedValue[1]))
                file1.write(key+" "+str(self.__get_value(self.frequencyTrans,splittedValue[0],splittedValue[1]))+"\n")
#                file2.write(key+" "+str(self.__get_value(self.frequencyTransExp,splittedValue[0],splittedValue[1]))+"\n")
                
                
        for key in self.frequencyEmit.keys():
            if "|" in key:
                splittedValue=key.split("|")
        #                print(key)
        #                print(key,self.__get_value(self.frequency,splittedValue[0],splittedValue[1]))
                file1.write(key+" "+str(self.__get_value(self.frequencyEmit,splittedValue[0],splittedValue[1]))+"\n")
#                file2.write(key+" "+str(self.__get_value(self.frequencyEmitExp,splittedValue[0],splittedValue[1]))+"\n")
                
#        print(min(self.frequencyTrans.values()),min(self.frequencyEmit.values()))

        
    
    
