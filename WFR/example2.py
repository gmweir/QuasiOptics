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
foo.defineStack(elem=[stackElement],qty=[5])

#Method: listSpecs
#-----------------

#Description:
#list filter specs
foo.listSpecs()

#Method: errorStudy
#------------------

#Description:
#performs an error study: repeat calculation of the frequency response x times to obtain error bars for the e.g. transmission curves
#and to answer the question how much the filter characteristics are changed by manufacturing tolerances

#optional parameters in square brackets, default values are shown
#Usage: obj.errorStudy(fRange=[fStart,fStop,fStep),[incidentAngle=0],[Quiet=False],[extractPOIs=True],pois=[frequencies in units of GHz],[numCalc=100],[parallelize=False],[numCores=2],[chunksize=1])

#Parameter: parallelize, numCores
#This calculation is parallelized, define whether you would like to use parallelization by parameter "parallelize=True" 
#and number of cores by "numCores=4"

#Parameter: numCalc
#the parameter "numCalc" defines how often the calculation of the frequency response is repeated with a slightly different stack taking into account normal distributed errors for the distances.

#Parameter: extractPOIs, pois
#If you would like to extract and display directly the transmisson/reflection etc. for interesting frequencies, 
#define them in parameter "pois" (points of interest) and set "extractPOIs" to True

#Parameter: incidentAngle
#You can define the incident angle of the beam when it enters the stack in units of degrees
[studyResult,poiResult]=foo.errorStudy(fRange=[127E9,144E9,100E6],incidentAngle=incidentAngle,extractPOIs=True,pois=[128.00,130.00,140.0,133.0,137.3,142.8,135.00, 140.0, 142.0],numCalc=500, parallelize=False)
#[studyResult,poiResult]=foo.errorStudy(fRange=[127E9,144E9,100E6],incidentAngle=incidentAngle,extractPOIs=True,pois=[128.00,130.00,140.0,133.0,137.3,142.8,135.00],numCalc=500, parallelize=True,numCores=4,chunksize=1)

#plot the result with markers for interesting frequencies
foo.plot(result=studyResult, \
marker=[{'f':140.02,'color':'red','width':2,'style':"-."},\
{'f':137.3,'color':'green','width':2,'style':"-."},\
{'f':133.0,'color':'red','width':2,'style':"-."},\
{'f':142.8,'color':'green','width':2,'style':"-."},\
{'f':130.0,'color':'red','width':2,'style':"-."},\
{'f':135.0,'color':'red','width':2,'style':"-."},\
{'f':128.0,'color':'red','width':2,'style':"-."}],\
plot=[{'key':'T','label':'Transmission [dB]','color':'blue','width':2,'style':'-'}],sizepx=[640,480])
