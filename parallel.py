# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 19:13:59 2016

@author: Zhan
"""

from gurobipy import *
import numpy as np

def paraLambda(pt,pj,pd):
    return 1
    '''
    print "First",pt,pj,pd
    base = 2+(pt*ins.j+pj)*(ins.d+1)+pd
    lhs = {}
    rhs = {}
    for i in range(0,ins.i):
        lhs[i] = LinExpr()
        rhs[i] = LinExpr()
    if ins.seg[pt,pj][1] != ins.seg[pt,pj][0]:
        for i in range(0,ins.i):
            lhs[i] += ins.l[i,base] * -1.0 / (ins.seg[pt,pj][pd+1] - ins.seg[pt,pj][pd])
            lhs[i] += ins.l[i,base+1] * 1.0 / (ins.seg[pt,pj][pd+1] - ins.seg[pt,pj][pd])
    else:
        for i in range(0,ins.i):
            lhs[i] += ins.l[i,base] * -1.0 
            lhs[i] += ins.l[i,base+1] * 1.0                        
    #Ax xi
    col = 1+(pt*ins.j+pj)*ins.d+pd
    for t in range(pt+1,ins.limt):
        for j in range(0,ins.j):
            for i in ins.refJ[j]:
                rhs[i] += ins.x[t,j,col]
    for t in range(max(pt+1,ins.limt),min(ins.t,pt+1+ins.limt)):
        for j in range(0,ins.j):
            for i in ins.refJ[j]:
                rhs[i] += ins.x[t,j,1+((pt-(t-ins.limt))*ins.j+pj)*ins.d+pd]
    #Axx xi
    for i in ins.refJ[pj]:
        rhs[i] += ins.xx[pt,pj,pd]
    return (lhs,rhs)
    '''
    
def paraGamma(ins,t,pt,pj,pd):
    print "Second",t,pt,pj,pd
    base = 2+(pt*ins.j+pj)*(ins.d+1)+pd
    lhs = {}
    rhs = {}
    for i in range(0,ins.j):
        lhs[i] = LinExpr()
        rhs[i] = LinExpr()
    if ins.seg[pt,pj][1] != ins.seg[pt,pj][0]:
        for i in range(0,ins.j):
            lhs[i] += ins.g[t,i,base] * -1.0 / (ins.seg[pt,pj][pd+1] - ins.seg[pt,pj][pd])
            lhs[i] += ins.g[t,i,base+1] * 1.0 / (ins.seg[pt,pj][pd+1] - ins.seg[pt,pj][pd])
    else:
        for i in range(0,ins.j):
            lhs[i] += ins.g[t,i,base] * -1.0 
            lhs[i] += ins.g[t,i,base+1] * 1.0                        
    #Ax xi
    col = 1+(pt*ins.j+pj)*ins.d+pd
    if t < ins.limt:
        if pt < t:                              
            for j in range(0,ins.j):
                rhs[j] += ins.x[t,j,col]
    else:
        if pt < t and pt >= t - ins.limt:                            
            for j in range(0,ins.j):
                rhs[j] += ins.x[t,j,1+((pt-(t-ins.limt))*ins.j+pj)*ins.d+pd]                            
    #Axx xi
    if pt == t:
        rhs[pj] += ins.xx[pt,pj,pd]
    return (lhs,rhs)
        
def paraOmega(ins,t,pt,pj,pd):
    print "Third",t,pt,pj,pd
    base = 2+(pt*ins.j+pj)*(ins.d+1)+pd
    lhs = {}
    rhs = {}
    for i in range(0,ins.j):
        lhs[i] = LinExpr()
        rhs[i] = LinExpr()
    if ins.seg[pt,pj][1] != ins.seg[pt,pj][0]:
        for i in range(0,ins.j):
            lhs[i] += ins.o[t,i,base] * -1.0 / (ins.seg[pt,pj][pd+1] - ins.seg[pt,pj][pd])
            lhs[i] += ins.o[t,i,base+1] * 1.0 / (ins.seg[pt,pj][pd+1] - ins.seg[pt,pj][pd])
    else:
        for i in range(0,ins.j):
            lhs[i] += ins.o[t,i,base] * -1.0 
            lhs[i] += ins.o[t,i,base+1] * 1.0                        
    #Ax xi
    col = 1+(pt*ins.j+pj)*ins.d+pd
    if t < ins.limt:
        if pt < t:                              
            for j in range(0,ins.j):
                rhs[j] += ins.x[t,j,col]
    else:
        if pt < t and pt >= t - ins.limt:                            
            for j in range(0,ins.j):
                rhs[j] += ins.x[t,j,1+((pt-(t-ins.limt))*ins.j+pj)*ins.d+pd]                            
    #Axx xi
    if pt == t:
        rhs[pj] += ins.xx[pt,pj,pd] - 1
    return (lhs,rhs)
            