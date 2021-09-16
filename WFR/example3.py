#import necessary modules
import wfr
import numpy as np
from matplotlib import pyplot as plt#
import time

#init WFR object
foo = wfr.WFR()

#define materials
n1=3.065
n2=1.00027
losstan1=15E-4
losstan2=0
foo.addMaterial("Al2O3",np.power(n1,2),losstan1,'disc')
foo.addMaterial("Air",np.power(n2,2),losstan2,'spacing')

#output list of defined materials
foo.listMaterial()

#define a stack element with/without a random normal distributed error (e.g for manufacturing tolerances)
useD1Error=True
useD2Error=True
d1spread=0.127E-3/3
d2spread=0.1E-3/3
incidentAngle=0
d1=1E-3
d2=1E-3
stackElement=foo.defineStackElement(material=["Al2O3","Air"],thickness=[d1,d2], spread=[d1spread,d2spread], useError=[useD1Error,useD2Error], absSpread=[True,True])

#Method: defineStack
#-------------------

#Description:
#define stack consisting of 5 times the above defined stack element, arbitrary combinations possible
stack1=foo.defineStack(elem=[stackElement],qty=[5])
stack2=foo.defineStack(elem=[stackElement],qty=[10])

#Method: listSpecs
#-----------------

#Description:
#list filter specs
foo.listSpecs()

#Method: compareStacks
#------------------

#Description:
#performs the calculation for multiple stacks and compares them for the points of interest defined in parameter "pois"
#plots the result

#optional parameters in square brackets, default values are shown
#Usage: obj.compareStacks(fRange=[fStart,fStop,fStep),[errorAnalysis=True],pois=[frequencies in units of GHz],[numCalc=100],[parallelize=False],[numCores=2],[chunksize=1])

foo.compareStacks(stacks=[stack1,stack2],fRange=[128E9,143E9,100E6],poi=[128.0,133.0,135.0,140.0,137.3,142.8])

