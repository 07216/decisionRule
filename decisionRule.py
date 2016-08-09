# -*- coding: utf-8 -*-
"""
Created on Thu Aug 04 15:23:11 2016

@author: Zhan
"""

from gurobipy import *
import numpy as np

class decisionRule:
    def __init__(self,relrALP):
        self.m = Model("Decision Rule Approch")
        self.i = relrALP.f
        self.j = relrALP.i
        self.t = relrALP.n
        self.t = 2
        self.T = 10
        self.limt = 2
        #Sparse P
        #self.P = np.zeros((20000,20000),dtype=np.float)
        self.A = np.zeros((self.i,self.j),dtype=np.float)
        self.xi = np.zeros((self.t*self.j+1,1),dtype=np.float)
        self.prob = np.zeros((self.t*self.T*self.j+1,1),dtype=np.float)
        self.c = np.zeros((self.i,1),dtype=np.float)
        self.h = np.zeros((2*(self.t*self.j+1),1),dtype=np.float)
        #Sparse W
        #self.W = np.zeros((2*(self.t*self.j+1),self.t*self.j+1),dtype=np.float)
        self.v = np.zeros((self.j,1),dtype=np.float)
        self.rALP = relrALP
        #look resources I required by J
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
            #self.c[item] = 2
        #construct h
        self.h[0] = 1
        self.h[1] = -1
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.h[2*(1+t*self.j+j)] = self.cal(t,j,1e-5)
                self.h[2*(1+t*self.j+j)+1] = 0
        #construct v
        for j in range(0,self.j):
            self.v[j] = self.rALP.pval[j]
        #construct Exi
        self.xi[0] = 1
        self.prob[0] = 1
        for t in range(0,self.t * self.T):
            for j in range(0,self.j):
                index = self.rALP.prdic[j]
                self.prob[1+t*self.j+j] = self.rALP.bdic[(t,index[0],index[1],index[2])]
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.xi[1+t*self.j+j] = 0
                for T in range(0,self.T):
                    self.xi[1+t*self.j+j] += self.prob[1+(t*self.T+T)*self.j+j]
        
        for t in range(0,self.t):
            for j in range(0,self.j):
                self.xi[1+t*self.j+j] *= 100
        print sum(self.xi)
                
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
        
    def cal(self,t,j,minp):
        s = 0
        for T in range(0,self.T):
            index = self.rALP.prdic[j]
            p = self.rALP.bdic[(t*self.T+T,index[0],index[1],index[2])]
            i = 0
            q = p
            tmp = 1
            while(q>minp):
                i = i+1
                tmp = tmp*(i+1)
                q = (p**(i+1))/tmp
            s += i
        return s
    
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
        #add Gamma
        self.g = {}
        '''
        for t in range(0,self.limt):
            for p in range(0,self.i):
                for i in range(0,2*(t*self.j+1)):
                    self.g[t,p,i] = self.m.addVar(lb=0, name = 'Gamma %d %d %d' % (t,p,i))
        '''
        for t in range(0,self.t):
            for p in range(0,self.i):
                for i in range(0,2*(self.t*self.j+1)):
                    self.g[t,p,i] = self.m.addVar(ub=0, name = 'Gamma %d %d %d' % (t,p,i))
            '''
        #add Beta
        for i in range(0,2*(self.t*self.j+1)):
            self.b[i] = self.m.addVar(lb=0, name = 'Beta %d' %i)
            '''
        self.m.update()
        
        #x size for first limt t, (self.j , t*self.j+1)
    
    def addOpt(self):
        obj  = LinExpr()
        #first limt time
        for t in range(0,self.limt):
            for j in range(0,self.j):
                for p in range(0,t*self.j+1):
                    if self.prob[t*self.T*self.j+j+1,0] * self.v[j,0] * self.xi[p,0] ==0:
                        continue
                    obj += self.prob[t*self.T*self.j+j+1,0] * self.v[j,0] * self.x[t,j,p] * self.xi[p,0]
                    #obj += self.v[j,0] * self.x[j,p]
        #after first limt time
        for t in range(self.limt,self.t):
            for j in range(0,self.j):
                #first element of xi = 1
                obj += self.prob[t*self.T*self.j+j+1,0] * self.v[j,0] * self.x[t,j,0]
                for p in range(1,self.limt*self.j+1):
                    if self.prob[t*self.T*self.j+j+1,0] * self.v[j,0] * self.xi[(t-self.limt)*self.j+p,0] ==0:
                        continue
                    obj += self.prob[t*self.T*self.j+j+1,0] * self.v[j,0] * self.x[t,j,p] * self.xi[(t-self.limt)*self.j+p,0]
                    
        self.m.setObjective(obj,GRB.MAXIMIZE)
    
    def addConstr(self):
        #Lambda 
        for p in range(0,self.i):
            l = LinExpr()
            for i in range(0,2*(self.t*self.j+1)):
                l += self.h[i,0] * self.l[p,i]
            self.m.addConstr(l,GRB.LESS_EQUAL,0,'h^T Lambda <=0')
        #Gamma
        for t in range(0,self.t):
            for p in range(0,self.i):
                g = LinExpr()
                for i in range(0,2*(self.t*self.j+1)):
                    g += self.h[i,0] * self.g[t,p,i]
                self.m.addConstr(g,GRB.GREATER_EQUAL,0,'h^T Gamma >=0')
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
            for t in range(0,self.limt):
                for ii in range(0,self.j):
                    if self.A[p,ii] ==0 :
                        continue
                    for j in range(0,t*self.j+1):
                        rhs[j] += self.prob[t*self.T*self.j+ii+1,0] * self.A[p,ii] * self.x[t,ii,j]
                        #rhs[j] += self.A[p,ii] * self.x[t,ii,j]
            #from limt to final time
            for t in range(self.limt,self.t):
                for ii in range(0,self.j):
                    if self.A[p,ii] ==0:
                        continue
                    rhs[0] += self.prob[t*self.T*self.j+ii+1,0] * self.A[p,ii] * self.x[t,ii,0]
                    #rhs[0] += self.A[p,ii] * self.x[t,ii,0]
                    for j in range((t-self.limt)*self.j+1,t*self.j+1):
                        rhs[j] += self.prob[t*self.T*self.j+ii+1,0] * self.A[p,ii] * self.x[t,ii,j-(t-self.limt)*self.j]
                        #rhs[j] += self.A[p,ii] * self.x[t,ii,j-(t-self.limt)*self.j]
            for j in range(0,self.t*self.j+1):
                self.m.addConstr(rhs[j],GRB.EQUAL,lhs[j],'Z1 %d %d' %(p,j))
        #Gamma for w^T Gamma = Z2
        for t in range(0,self.limt):
            for p in range(0,self.i):
                lhs = {}
                rhs = {}
                for i in range(0,self.t*self.j+1):
                    lhs[i] = LinExpr()
                    rhs[i] = LinExpr()
                for i in range(0,self.t*self.j+1):
                    lhs[i] += self.g[t,p,2*i]-self.g[t,p,2*i+1]
                for j in range(0,t*self.j+1):
                    rhs[j] += self.x[t,p,j]
                for j in range(t*self.j+1,self.t*self.j+1):
                    rhs[j] += 0
                for j in range(0,self.t*self.j+1):
                    self.m.addConstr(lhs[j],GRB.EQUAL,rhs[j],'Z2 time %d, %d %d' %(t,p,j))
        #after first limt time
        for t in range(self.limt,self.t):
            for p in range(0,self.i):
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
                
    def solve(self):
        self.m.optimize()
    
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
