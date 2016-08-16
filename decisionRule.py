# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:23:11 2016

@author: Zhan
"""

from gurobipy import *
import numpy as np

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
        self.x = {}
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j+1):
                    self.x[t,j,p] = self.m.addVar(name = 'X %d %d %d' % (t,j,p))
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                for p in range(0,self.limt*self.j+1):
                    self.x[t,j,p] = self.m.addVar(name = 'X %d %d %d' % (t,j,p))
        #add Lambda
        self.l = {}
        for p in range(0,self.i):
            for i in range(0,2*(self.t*self.j+1)):
                self.l[p,i] = self.m.addVar(lb=0, name = 'Lambda %d %d' % (p,i))
                '''
        #add Gamma
        self.g = {}
        for t in range(0,self.t):
            for p in range(0,self.j):
                for i in range(0,2*(self.t*self.j+1)):
                    self.g[t,p,i] = self.m.addVar(ub=0, name = 'Gamma %d %d %d' % (t,p,i))
        '''
        #add Omega
        self.o = {}
        for t in range(0,self.t):
            for p in range(0,self.j):
                for i in range(0,2*(self.t*self.j+1)):
                    self.o[t,p,i] = self.m.addVar(lb=0, name = 'Omega %d %d %d' %(t,p,i))

        self.m.update()
        
        #x size for first limt t, (self.j , t*self.j+1)
    
    def addOpt(self):
        obj  = LinExpr()
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j+1):
                    if self.v[j,0] * self.xi[p,0] ==0:
                        continue
                    obj += self.v[j,0] * self.x[t,j,p] * self.xi[p,0]
                    #obj += self.v[j,0] * self.x[j,p]
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                #first element of xi = 1
                obj += self.v[j,0] * self.x[t,j,0]
                for p in range(1,self.limt*self.j+1):
                    if self.v[j,0] * self.xi[(t-self.limt)*self.j+p,0] ==0:
                        continue
                    obj += self.v[j,0] * self.x[t,j,p] * self.xi[(t-self.limt)*self.j+p,0]
                    
        self.m.setObjective(obj,GRB.MAXIMIZE)
    
    def addConstr(self):
        #Lambda 
        for p in range(0,self.i):
            l = LinExpr()
            for i in range(0,2*(self.t*self.j+1)):
                if self.h[i,0] != 0:
                    l += self.h[i,0] * self.l[p,i]
            self.m.addConstr(l,GRB.LESS_EQUAL,0,'h^T Lambda <=0')
        #Gamma
            '''
        for t in range(0,self.t):
            for p in range(0,self.j):
                g = LinExpr()
                for i in range(0,2*(self.t*self.j+1)):
                    if self.h[i,0] != 0:
                        g += self.h[i,0] * self.g[t,p,i]
                self.m.addConstr(g,GRB.GREATER_EQUAL,0,'h^T Gamma >=0')
                '''
        #Omega
        for t in range(0,self.t):
            for p in range(0,self.j):
                o = LinExpr()
                for i in range(0,2*(self.t*self.j+1)):
                    if self.h[i,0] != 0:
                        o += self.h[i,0] * self.o[t,p,i]
                self.m.addConstr(o,GRB.LESS_EQUAL,0,'h^T Omega <= 0')
        #Lambda for w^T Lambda = Z1
        for p in range(0,self.i):
            lhs = {}
            rhs = {}
            for i in range(0,self.t*self.j+1):
                lhs[i] = LinExpr()
                rhs[i] = LinExpr()
            for i in range(0,self.t*self.j+1):
                lhs[i] += self.l[p,2*i]-self.l[p,2*i+1]
            #-CQ
            #Q 1*(t*j+1) (1,0,0,...)
            #C i*1 (c1,c2,c3,...,ci)^T
            rhs[0] += -self.c[p,0]
            #from 0 to limt
            for ii in range(0,self.j):
                if self.A[p,ii] ==0 :
                    continue
                for t in range(0,self.limt):
                    for j in range(0,t*self.j+1):
                        rhs[j] += self.A[p,ii] * self.x[t,ii,j]
            #from limt to final time
            for ii in range(0,self.j):
                if self.A[p,ii] ==0:
                    continue
                for t in range(self.limt,self.t):
                    rhs[0] += self.A[p,ii] * self.x[t,ii,0]
                    #rhs[0] += self.A[p,ii] * self.x[t,ii,0]
                    for j in range((t-self.limt)*self.j+1,t*self.j+1):
                        rhs[j] += self.A[p,ii] * self.x[t,ii,j-(t-self.limt)*self.j]
                        #rhs[j] += self.A[p,ii] * self.x[t,ii,j-(t-self.limt)*self.j]
            for j in range(0,self.t*self.j+1):
                self.m.addConstr(rhs[j],GRB.EQUAL,lhs[j],'Z1 %d %d' %(p,j))
        #Gamma for w^T Gamma = Z2
        #t before limt
                '''
        for t in range(0,self.limt):
            for p in range(0,self.j):
                lhs = {}
                rhs = {}
                for i in range(0,self.t*self.j+1):
                    lhs[i] = LinExpr()
                    rhs[i] = LinExpr()
                for i in range(0,self.t*self.j+1):
                    lhs[i] += self.g[t,p,2*i]-self.g[t,p,2*i+1]
                for j in range(0,t*self.j+1):
                    rhs[j] += self.x[t,p,j]
                for j in range(0,self.t*self.j+1):
                    self.m.addConstr(lhs[j],GRB.EQUAL,rhs[j],'Z2 time %d, %d %d' %(t,p,j))
        #after first limt time
        for t in range(self.limt,self.t):
            for p in range(0,self.j):
                lhs = {}
                rhs = {}
                for i in range(0,self.t*self.j+1):
                    lhs[i] = LinExpr()
                    rhs[i] = LinExpr()
                for i in range(0,self.t*self.j+1):
                    lhs[i] += self.g[t,p,2*i]-self.g[t,p,2*i+1]
                rhs[0] += self.x[t,p,0]
                for j in range((t-self.limt)*self.j+1,t*self.j+1):
                    rhs[j] += self.x[t,p,j-(t-self.limt)*self.j]
                for j in range(0,self.t*self.j+1):
                    self.m.addConstr(lhs[j],GRB.EQUAL,rhs[j],'Z2 time %d, %d %d' %(t,p,j))
                    '''
        #Omega for w^T Omega = Z3
        #before limt
        for t in range(0,self.limt):
            for p in range(0,self.j):
                lhs = {}
                rhs = {}
                for i in range(0,self.t*self.j+1):
                    lhs[i] = LinExpr()
                    rhs[i] = LinExpr()
                for i in range(0,self.t*self.j+1):
                    lhs[i] += self.o[t,p,2*i] - self.o[t,p,2*i+1]
                for j in range(0,t*self.j+1):
                    rhs[j] += self.x[t,p,j]
                #for j in range(t*self.j+1,(t+1)*self.j+1):
                #    rhs[j] += -1
                rhs[t*self.j+p+1] += -1
                for j in range(0,self.t*self.j+1):
                    self.m.addConstr(lhs[j],GRB.EQUAL,rhs[j],'Z3 time %d,%d %d' %(t,p,j))
        #after limt
        for t in range(self.limt,self.t):
            for p in range(0,self.j):
                lhs = {}
                rhs = {}
                for i in range(0,self.t*self.j+1):
                    lhs[i] = LinExpr()
                    rhs[i] = LinExpr()
                for i in range(0,self.t*self.j+1):
                    lhs[i] += self.o[t,p,2*i] - self.o[t,p,2*i+1]
                rhs[0] += self.x[t,p,0]
                for j in range(1,self.limt*self.j+1):
                    rhs[j+(t-self.limt)*self.j] += self.x[t,p,j]
                #for j in range(t*self.j+1,(t+1)*self.j+1):
                    #rhs[j] += -1
                rhs[t*self.j+p+1] += -1
                for j in range(0,self.t*self.j+1):
                    self.m.addConstr(lhs[j],GRB.EQUAL,rhs[j],'Z3 time %d,%d %d' %(t,p,j))
                    
    def solve(self):
        self.m.optimize()
        return 0
    
    def echoOpt(self):
        
        obj  = 0
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j+1):
                    if self.xi[1+t*self.j+j,0] * self.v[j,0] * self.xi[p,0] ==0:
                        continue
                    obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.x[t,j,p].X * self.xi[p,0] 
                    #obj += self.v[j,0] * self.x[j,p]
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                #first element of xi = 1
                obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.x[t,j,0].X
                for p in range(1,self.limt*self.j+1):
                    if self.xi[1+t*self.j+j,0] * self.v[j,0] * self.xi[(t-self.limt)*self.j+p,0] ==0:
                        continue
                    obj += self.xi[1+t*self.j+j,0] * self.v[j,0] * self.x[t,j,p].X * self.xi[(t-self.limt)*self.j+p,0]
        print obj
    
    def writeMPS(self):
        self.m.write("out.mps")
    
#    def echoSolution(self):
