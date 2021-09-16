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

#Method: listSpecs
#-----------------

#Description:
#list filter specs
foo.listSpecs()

#Method: findStackElement
#------------------

#Description:
#performs a parameter search to find the best stack (element) which fulfills defined requirements

#optional parameters in square brackets, default values are shown
#Usage: obj.findStackElement(goals,sweep,stackElement=elem,[numCores=2],[parallelize=True],[numCalc=100],[chunksize=1])

#Parameter:goals
#It is an array containing dictionaries of the following form:
#{'f':frequency,'qty':array containing one or more of the following keywords 'T','devT','R','devR','A','devA', 
#'target': array containing one or more of the following keywords 'min'/'max',
#'threshold': array containg as much values as defined quantities}
#
#Each dictionary defines a requirement the stack has to fulfill, devT e.g. stands for the standard deviation of the transmission
#
#e.g. goals=[{'f':140.0E9,'qty':['T','devT'],'taget':['min','min'], 'threshold':[-40,6]}]
#defines that for 140GHz we want to minimize the transmision and the error of the transmission, so that it is smaller than -40 or 6 dB respectively

#Parameter:sweep
#defines the search space
#It is a dictionary of the following form
#{'layer': one or more indices of the layer of the given stack element,
#'range': array containing each an array with values to try e.g. thickness of the layer with index 0,
#'count': array containg quantities for given stack element e.g. [5,10]}

#Parameter:stackElement
#expects a stack element to be optimized

result=foo.findStackElement( goals=[{'f':140.0E9,'qty':['T','devT'],'taget':['min','min'], 'threshold':[-40,6]}],sweep={'layer':[0,1],'range':[d1,d2],'count':lset},stackElement=stackElement)

labels=[str(['{:1.2e}'.format(elem['stackElement'][0]['thickness'][0]),'{:1.2e}'.format(elem['stackElement'][0]['thickness'][1]),elem['quantity'][0]]) for elem in result['stackElement']]

#plot candidates found 
foo.compareCandidates(result, labels, qty='T')

