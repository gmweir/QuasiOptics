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

#define stack consisting of 5 times the above defined stack element, arbitrary combinations possible
foo.defineStack(elem=[stackElement],qty=[5])

#build the stack
foo.buildStack(useError=False)

#list filter specs
foo.listSpecs()

#Method: calculateFrequencyResponse
#----------------------------------

#optional parameters in square brackets, default values shown
#Usage: 
#obj.calculateFrequencyResponse(fStart,fStop,fStep,[incidentAngle=0],[quiet=False])

#Description:
#perform a single calculation with console output
calcResult=foo.calculateFrequencyResponse(120E9,150E9,20E6)

#plot the result with markers for interesting frequencies
foo.plot(result=calcResult, \
marker=[{'f':140.02,'color':'red','width':2,'style':"-."},\
{'f':137.3,'color':'green','width':2,'style':"-."},\
{'f':133.0,'color':'red','width':2,'style':"-."},\
{'f':142.8,'color':'green','width':2,'style':"-."},\
{'f':130.0,'color':'red','width':2,'style':"-."},\
{'f':135.0,'color':'red','width':2,'style':"-."},\
{'f':128.0,'color':'red','width':2,'style':"-."}],\
plot=[{'key':'T','label':'Transmission [dB]','color':'blue','width':2,'style':'-'}],sizepx=[640,480])
