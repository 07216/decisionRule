# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 16:44:41 2016

@author: Zhan
"""

from gurobipy import *
import numpy as np

class Input:
    def __init__(self):
        self.m = Model("rALP")
        
    def readIn(self,filename):        
        with open(filename) as file:
            next(file)
            self.n = int(next(file))
            
            for i in range(0,3):
                next(file)
            self.f = int(next(file))
            self.flight = []
            for i in range(0,self.f):
                self.flight.append([int(x) for x in next(file).split()])
                
            for i in range(0,3):
                next(file)
            self.i = int(next(file))
            self.iti = []
            for i in range(0,self.i):
                t = next(file).split()
                self.iti.append((int(t[0]),int(t[1]),int(t[2]),float(t[3])))
                
            for i in range(0,3):
                next(file)
            self.prob = []
            self.bdic = {}
            for i in range(0,self.n):
                t = next(file).split()
                r = []
                r.append(int(t[0]))
                j = 2
                while(j<len(t)):
                    r.append((int(t[j]),int(t[j+1]),int(t[j+2]),float(t[j+4])))
                    #bdic -- index (time,from,to,class) and return 
                    self.bdic[(i,int(t[j]),int(t[j+1]),int(t[j+2]))] = float(t[j+4])
                    j = j+6
                self.prob.append(r)
        
    def product(self):
        self.pdic = {}
        self.prdic = {}
        self.pval = []
        for i in range(0,self.i):
            t = self.iti[i]
            #pdic -- index (from,to,class) and return id of the product
            #its inverse prdic
            #self.pval[pdic] store the value
            self.pdic[(t[0],t[1],t[2])] = i
            self.prdic[i] = (t[0],t[1],t[2])
            self.pval.append(t[3])
    
    def leg(self):
        self.ldic = {}
        self.lrdic = {}
        self.leg = []
        for i in range(0,self.f):
            t = self.flight[i]
            #ldic -- index (from,to) and return id of the leg
            #its inverse lrdic
            #self.flight[i][2] store numbers of legs
            self.ldic[(t[0],t[1])] = i
            self.lrdic[i] = (t[0],t[1])

    def lam(self):
        self.lmda = []
        for t in range(0,self.n):
            lmda = {}
            
    
    def constrRefforIJ(self,a,b):
        b=b-1
        if self.dicAforI.has_key(a):
            self.dicAforI[a].append(b)
        else:
            self.dicAforI[a]=[b]
        if self.dicAforJ.has_key(b):
            self.dicAforJ[b].append(a)
        else:
            self.dicAforJ[b] = [a]
    
    def construct(self):
        #self.a = np.zeros((20,2000),dtype=np.int8)
        self.val = []
        self.id = 0
        self.listA = []
        self.ref = []
        #self.lmda = np.zeros((1000,200),dtype=np.float)
        self.dicAforI = {}
        self.dicAforJ = {}
        for j in range(0,self.i):
            t = self.prdic[j]
            if self.ldic.has_key((t[0],t[1])):
                self.id = self.id + 1
                self.ref.append(j)
                self.val.append(self.pval[j])
                #for p in range(0,self.n):
                    #self.lmda[p,self.id] = self.bdic[(p,t[0],t[1],t[2])]
                #self.a[self.ldic[(t[0],t[1])],self.id] = 1
                self.listA.append((self.ldic[(t[0],t[1])],self.id))
                self.constrRefforIJ(self.ldic[(t[0],t[1])],self.id)
            else:
                for k in range(0,self.f/2+1):
                    if k != t[0] and k != t[1] and self.ldic.has_key((t[0],k)) and self.ldic.has_key((k,t[1])):
                        self.id = self.id+1
                        self.ref.append(j)
                        #self.a[self.ldic[(t[0],k)],self.id] = 1
                        #self.a[self.ldic[(k,t[1])],self.id] = 1
                        self.val.append(self.pval[j])
                        #for p in range(0,self.n):
                            #self.lmda[p,self.id] = self.bdic[(p,t[0],t[1],t[2])]
                        self.listA.append((self.ldic[(t[0],k)],self.id))
                        self.listA.append((self.ldic[(k,t[1])],self.id))                
                        self.constrRefforIJ(self.ldic[(t[0],k)],self.id)                
                        self.constrRefforIJ(self.ldic[(k,t[1])],self.id)        

    def echoVal(self):
        print "id:%d"%self.id
        print "f:%d"%self.f
        print "i:%d"%self.i
        print self.prob
        print self.pval
        print self.flight
    
    def addVAR(self):
        self.q = {}
        self.y = {}
        self.z = {}
        for t in range(0,self.n):
            for i in range(0,self.id):
                self.q[t,i] = self.m.addVar(vtype=GRB.CONTINUOUS,name = 'q_%s_%s' % (t,i))
            for i in range(0,self.f):
                for k in range(0,self.flight[i][2]):
                    self.y[t,i,k] = self.m.addVar(vtype=GRB.CONTINUOUS,name = 'y_%s_%s_%s' % (t,i,k))
            for j in range(0,self.id):
                for i in range(0,self.f):
                    for k in range(0,self.flight[i][2]+1):
                        #print self.n*self.id*self.f*self.flight[i][2]
                        self.z[t,j,i,k] = self.m.addVar(vtype=GRB.CONTINUOUS,name = 'z_%s_%s_%s' % (t,j,i,k))
        self.m.update()
        
    def addOpt(self):
        return 