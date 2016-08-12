# -*- coding: utf-8 -*-
"""
Created on Sun Aug 07 13:12:01 2016

@author: Zhan
"""

from gurobipy import *
import numpy as np

class rALP:
    def __init__(self,relrALP):
        self.m = Model("Decision Rule Approch")
        self.i = relrALP.f
        self.j = relrALP.i
        self.t = relrALP.n
        
        
        self.t = 20#customized periods
        
        
        #Sparse P
        #self.P = np.zeros((20000,20000),dtype=np.float)
        self.A = np.zeros((self.i,self.j),dtype=np.float)
        self.xi = np.zeros((self.t*self.j+1,1),dtype=np.float)
        self.c = np.zeros((self.i,1),dtype=np.float)
        self.h = np.zeros((2*(self.t*self.j+1),1),dtype=np.float)
        #Sparse W
        #self.W = np.zeros((2*(self.t*self.j+1),self.t*self.j+1),dtype=np.float)
        self.v = np.zeros((self.j,1),dtype=np.float)
        self.rALP = relrALP
        self.refI = relrALP.dicAforI
        self.refJ = relrALP.dicAforJ
        
    def construct(self):
        #construct A
        for item in self.rALP.listA:
            #Hint Item[1] is from 1 to ...
            #print item
            self.A[item[0],item[1]-1] = 1
        #construct c
        for item in range(0,self.i):
            self.c[item] = self.rALP.flight[item][2]
            #self.c[item] = 10
        #construct v
        for j in range(0,self.j):
            self.v[j] = self.rALP.pval[j]
        #construct Exi
        self.xi[0] = 1
        for t in range(0,self.t):
            for j in range(0,self.j):
                index = self.rALP.prdic[j]
                self.xi[1+t*self.j+j] = self.rALP.bdic[(t,index[0],index[1],index[2])]
                
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
    
    def addVar(self):
        #add y
        self.y = {}
        for t in range(0,self.t):
            for i in range(0,self.i):
                for k in range(0,self.c[i]):
                    self.y[t,i,k] = self.m.addVar(name = 'y %d %d %d' %(t,i,k) )
        self.z = {}
        for t in range(0,self.t):
            for j in range(0,self.j):
                for i in range(0,self.i):
                    for k in range(0,int(self.c[i])):
                        self.z[t,j,i,k] = self.m.addVar(name = 'z %d %d %d %d' %(t,j,i,k) )
                    self.z[t,j,i,int(self.c[i,0])] = self.m.addVar(ub=0,lb=0,name = 'z %d %d %d %d' %(t,j,i,self.c[i,0]) )
        self.q = {}
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.q[t,j] = self.m.addVar(name = 'q %d %d' %(t,j) )
        self.m.update()
    
    def addOpt(self):
        obj  = LinExpr()
        for t in range(0,self.t):
            for j in range(0,self.j):
                obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.q[t,j]
        self.m.setObjective(obj,GRB.MAXIMIZE)
    
    def addConstr(self):
        #y
        for i in range(0,self.i):
            for k in range(0,self.c[i]):
                self.m.addConstr(self.y[0,i,k],GRB.EQUAL,1,name = 'Constr y %d %d %d' % (0,i,k) )
        for t in range(1,self.t):
            for i in range(0,self.i):
                for k in range(0,self.c[i]):
                    rhs = LinExpr()
                    for j in self.refI[i]:
                        rhs += self.xi[1+(t-1)*self.j+j,0] * (self.z[t-1,j,i,k] - self.z[t-1,j,i,k+1])
                    self.m.addConstr(self.y[t-1,i,k]-self.y[t,i,k],GRB.EQUAL,rhs,name = 'Constr y %d %d %d' %(t,i,k))
        #q
        for t in range(0,self.t):
            for j in range(0,self.j):
                for i in self.refJ[j]:
                    self.m.addConstr(self.q[t,j],GRB.EQUAL,self.z[t,j,i,1],name = 'Constr q %d %d %d' %(t,j,i))
        #z
        for t in range(0,self.t):
            for j in range(0,self.j):
                for i in self.refJ[j]:
                    for k in range(0,self.c[i]):
                        self.m.addConstr(self.z[t,j,i,k+1],GRB.LESS_EQUAL,self.z[t,j,i,k],name = 'Constr z %d %d %d %d'%(t,j,i,k+1))
                        self.m.addConstr(self.z[t,j,i,k],GRB.LESS_EQUAL,self.y[t,i,k],name = 'Constr z<=y %d %d %d %d'%(t,j,i,k))
                        
    def solve(self):
        self.m.optimize()
    
    def echoOpt(self):
        print self.m.getObjective().getValue()
    
    def writeMPS(self):
        self.m.write("out.mps")
    
#    def echoSolution(self):
