# -*- coding: utf-8 -*-
"""
Created on Wed Aug 03 16:50:33 2016

@author: Zhan
"""


import Input
import decisionRule
import rALP
import time
import simulation
import CustomizeDemand
#from guppy import hpy

reader = Input.Input()
reader.readIn('data/rm_200_4_1.0_4.0.txt')
reader.product()
reader.leg()
reader.construct()

demander = CustomizeDemand.CustomizeDemand(1)

first = 1

if first:

    
    decisionSolver = decisionRule.decisionRule()
    #decisionSolver.echoInput()
    #decisionSolver.echoVal()
    decisionSolver.inputDemand(demander)
    decisionSolver.addVar()
    decisionSolver.addOpt()
    decisionSolver.addConstr()
    decisionSolver.solve()
    #decisionSolver.echoOpt()
    
    simulator = simulation.simulation(decisionSolver,demander)
    simulator.run(1000)

else:
    
    reduction = rALP.rALP(reader)
    reduction.construct()
    reduction.addVar()
    reduction.addOpt()
    reduction.addConstr()
    reduction.solve()