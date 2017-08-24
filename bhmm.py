# -*- coding: utf-8 -*-

"""
  @author: Necva Bolucu (@necvabolucu)
  @author: Salih Tuc (@salihtuc)
"""

from collections import defaultdict
from math import exp
class BHMM(object):
    
    def __init__(self):
        self.frequency = defaultdict(float)
        self.frequencyExp = defaultdict(float)
    def __change_count(self,matrix, x, y, i):
#        print(x,y)
        matrix["%s|%s" % (x, y)] += i
        matrix["%s" % x] += i
#        print( matrix["%s|%s" % (x, y)], matrix["%s" % x])
    
    def __change_countExp(self,matrix, x, y, i):
    #        print(x,y)
            matrix["%s|%s" % (x, y)] = +exp(i)
            matrix["%s" % x] +=exp( i  )
                  
    def __get_value(self,matrix, x,y):
#        print(x, y)
        return matrix["%s|%s" % (x, y)]/matrix["%s" % (x)]
    
    def __get_valueExp(self,matrix, x,y):
#        print(x, y)
        return exp(matrix["%s|%s" % (x, y)]-matrix["%s" % (x)])
        
    def run(self):
        file=open("result.txt")
        file1=open("out.txt","w")
        file2=open("outExp.txt","w")
        for line in file.readlines():
            splitted=line.strip("\n").strip().split()
            splittedValue=splitted[0].split("-",1)
#            print(splittedValue)
            self.__change_count(self.frequency, splittedValue[0], splittedValue[1],float(splitted[1]))
            self.__change_countExp(self.frequencyExp, splittedValue[0], splittedValue[1],float(splitted[1]))
            
        print(len(self.frequency.keys()))
        for key in self.frequency.keys():
            if "|" in key:
                splittedValue=key.split("|")
#                print(key)
#                print(key,self.__get_value(self.frequency,splittedValue[0],splittedValue[1]))
                file1.write(key+" "+str(self.__get_value(self.frequency,splittedValue[0],splittedValue[1]))+"\n")
                file2.write(key+" "+str(self.__get_value(self.frequencyExp,splittedValue[0],splittedValue[1]))+"\n")
        
    
    
