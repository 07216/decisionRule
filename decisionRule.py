# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:23:11 2016

@author: Zhan
"""

from gurobipy import *
import numpy as np
import scipy as sp

class decisionRule:
    def __init__(self):
        self.m = Model("Decision Rule Approch")
        
        
        #Sparse P
        #self.P = np.zeros((20000,20000),dtype=np.float)
        #self.prob = np.zeros((self.t*self.T*self.j+1,1),dtype=np.float)
        #Sparse W
        #self.W = np.zeros((2*(self.t*self.j+1),self.t*self.j+1),dtype=np.float)

    def inputDemand(self,recorder):
        self.i = recorder.i
        self.j = recorder.j
        self.t = recorder.t
        self.A = recorder.A
        self.xi = recorder.xi
        self.v = recorder.v
        self.c = recorder.c
        self.h = recorder.h
        self.refJ = recorder.refJ
        self.t = recorder.t
        self.limt  = recorder.limt
        self.d = recorder.d
        self.r = recorder.r
        self.w = recorder.w
        self.a = recorder.a
        self.b = recorder.b
        self.r = recorder.r
        self.row = sp.shape(self.w)[0]
        self.col = sp.shape(self.w)[1]
        self.lessThanXi = recorder.mesh
           
    def echoInput(self):
        print "i:%d"%self.i
        print "j:%d"%self.j
        print "t:%d"%self.t
        print "legs:"
        print self.c
        print "value:"
        print self.v
        print "h:"
        print self.h
        print "A:"
        print self.A
        print "xi:"
        print self.xi
    
    def echoVal(self):
        for x in self.h:
            print x
    
    def addVar(self):
        #add X
        self.x = np.empty((self.t,self.j,self.a,self.b),dtype=object)
        #first limt time
        for t in range(0,self.t):
            for j in range(0,self.j):
                for a in range(0,self.a):
                    for b in range(self.b):
                        self.x[t,j,a,b] = self.m.addVar(lb=0, ub=self.lessThanXi[t,j][0][a])        
        self.l = np.empty((self.i,self.row),dtype=object)
        for i in range(self.i):
            for j in range(self.row):
                self.l[i,j] = self.m.addVar(lb=-GRB.INFINITY,ub=0)
        
        self.m.update()
        
        #x size for first limt t, (self.j , t*self.j+1)
    
    def addOpt(self):
        obj  = LinExpr()
        for t in range(0,self.t):
            for j in range(0,self.j):
                for a in range(self.a):
                    for b in range(self.b):
                        obj += self.v[j,0] * self.x[t,j,a,b] * self.xi[t,j,a,b]
                    
        self.m.setObjective(obj,GRB.MAXIMIZE)
    
    def addConstr(self):
        tmp = sp.dot(self.l,self.h)
        for item in tmp:
            self.m.addConstr(item, GRB.LESS_EQUAL, 0)
        tmp = np.empty((self.i,self.col),dtype=object)
        for i in range(self.i):
            for j in range(self.col):
                tmp[i,j] = LinExpr()
        (row,col)  = self.w.nonzero()
        for i in range(len(row)):
            for j in range(self.i):
                tmp[j,col[i]] += self.l[j,row[i]] * self.w[row[i],col[i]]
        z = np.empty((self.i,self.col),dtype=object)
        for i in range(self.i):
            for j in range(self.col):
                z[i,j] = LinExpr()
        for t in range(self.t):
            for j in range(self.j):
                for i in self.refJ[j]:
                    for a in range(self.a):
                        for b in range(self.b):
                            z[i,self.t * self.j + (t*self.j + j) * self.r + a *self.b + b] += self.x[t,j,a,b]
                            '''
                        start = self.t * self.j + (t*self.j + j) * self.r + a *self.b
                        end = start + self.b
                        z[i,start:end] += self.x[t,j,a]
                                '''
        for i in range(self.i):
            z[i,self.col-1] += -2*self.c[i]
        for t in range(self.t):
            for j in range(self.j):
                for i in self.refJ[j]:
                    for a in range(self.a):
                        for b in range(self.b):
                            z[i,self.col-1] += self.x[t,j,a,b] * self.xi[t,j,a,b]
        for i in range(tmp.shape[0]):
            for j in range(tmp.shape[1]):
                self.m.addConstr(tmp[i,j],GRB.EQUAL,z[i,j])
                
    def solve(self):
        self.m.optimize()
    
    def echoOpt(self):
        for j in range(self.j):
            for a in range(self.a):
                print j,a,self.x[0,j,a,self.b-1].x
                
    
    def expectedLeftDemand(self):
        c = self.c.copy()
        for t in range(self.t):
            for j in range(self.j):
                for i in range(self.refJ[j]):
                    for a in range(self.a):
                        for b in range(self.b):
                            c[i] -= self.x[t,j,a,b] * self.xi[t,j,a,b]
        print c
    
    def writeMPS(self):
        self.m.write("out.lp")
    
#    def echoSolution(self):
