# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 16:50:33 2016

@author: Zhan
"""


import Input
import sys
import decisionRule
import rALP
import time
import simulation
import CustomizeDemand
import gene
import numpy as np
#from guppy import hpy

reductionALP=0
firstCaseInReSolve=1
secondCaseInReSolve=2

choose,t,limt,d,T = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5])
demander = CustomizeDemand.CustomizeDemand(choose,t,limt,d,T)

first = 0

if first == 0:#Decision Rule Approch 
    
    decisionSolver = decisionRule.decisionRule()
    #   decisionSolver.echoInput()
    #decisionSolver.echoVal()
    decisionSolver.inputDemand(demander)
    start = time.clock()
    print "Input End",start
    decisionSolver.addVar()
    print "AddVar End",time.clock()
    decisionSolver.addOpt()
    print "AddOpt End",time.clock()
    decisionSolver.addConstr()
    end = time.clock()
    print "AddCons End",end
    #decisionSolver.writeMPS()
    decisionSolver.solve()
    decisionSolver.expectedLeftDemand()
    #decisionSolver.echoOpt()    
    
    print t,limt,d,end-start
    
    simulator = simulation.simulation(decisionSolver,demander)
    '''
    realDemand = {}
    for i in range(0,3):
        realDemand[i] = demander.sim()
    
    opt = gene.gene(decisionSolver,simulator)
    opt.initX()
    opt.setSTD(realDemand)
    opt.evolve()
    '''
    simulator.initX()
    simulator.initXX()
    #simulator.echoXX()
    simulator.bookLimLeft()
    print simulator.bookLimRun(np.round,1000)
    print simulator.bookLimRun(np.ceil,1000)
    print simulator.bookLimRun(np.floor,1000)
    print simulator.bookLimRun(simulator.atLeastOne,1000)
    print simulator.run(1000)

elif first == 1:#reduction of Approximate Linear Programming
    
    reduction = rALP.rALP(demander.reader)
    reduction.construct()
    reduction.addVar()
    reduction.addOpt()
    reduction.addConstr()
    reduction.solve()

elif first == 2:#Gene Algorithm
    
    decisionSolver = decisionRule.decisionRule()
    decisionSolver.inputDemand(demander)    
    decisionSolver.x = {}
    simulator = simulation.simulation(decisionSolver,demander)

    realDemand = {}
    for i in range(0,100):
        realDemand[i] = demander.sim()

    opt = gene.gene(decisionSolver,simulator)
    opt.zeroX()
    opt.setSTD(realDemand)
    opt.evolve()
    #simulator.run(1000)
    print simulator.bookLimRun(100,opt.X)
