from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import sys
import copy
from progressbar import progressBar
import multiprocessing
import time
from cycler import cycler
import itertools
import re

class WFR(object):
	"""
	The class WFR encapsulates all methods to calculate and plot the frequency response for
	a multi disc window. It is a port of a MATLAB script 
	by H.Oosterbeek, based on script by Hartfuss, Geist, Nickel, Thumm
	
	Fabian Wilde, IPP-HGW, E3, 6/2017

	Version 0.9 - 13.6.2017
	-----------------------
	-Parallelization is still missing
	
	Version 1.0 - 15.6.2017
	-----------------------
	-Parallelization for error study achieved
	
	"""

	def __init__(self):
		self.layers = dict({'thickness':[],'material':[],'epsilonR':[],'lossTangent':[]})
		self.stack = dict({'stackElement':[],'quantity':[],'label':'','length':[]})
		self.stackElement = dict({'thickness':[],'material':[],'useError':[],'spread':[],'absSpread':[],'label':''})
		self.materials = dict()
		self.WFR = dict({'freqs':[],'T':[],'R':[],'A':[]})

	def clear(self):
		self.__init__()

	def s2t(self, s):
		"""
		converts the S-matrix to a T-matrix
		
		Parameters:
		-----------
		s : numpy.ndarray - 2x2 numpy array containing the S-matrix
		"""

		# if not isinstance(s, np.ndarray):
		# 	raise ValueError("expected numpy array for input argument 's'.")
		# else:
		# 	print(s.size)
		# 	if not s.size == (2, 2):
		# 		raise ValueError("expected numpy array with size (2,2).")

		t = np.zeros((2, 2), dtype=complex)
		t[0, 0] = s[0, 1] - ((s[0, 0] * s[1, 1]) / (s[1, 0]))
		t[0, 1] = s[0, 0] / s[1, 0]
		t[1, 0] = -s[1, 1] / s[1, 0]
		t[1, 1] = 1 / s[1, 0]

		return t

	def t2s(self, t):
		"""
		converts the T-matrix to a S-matrix
		
		Parameters:
		-----------
		t : numpy.ndarray - 2x2 numpy array containing the T-matrix
		"""

		# if not isinstance(t, np.ndarray):
		# 	raise ValueError("expected numpy array for input argument 's'.")
		# else:
	# 		if not t.size == (2, 2):
		# 		raise ValueError("expected numpy array with size (2,2).")

		s = np.zeros((2, 2), dtype=complex);
		s[0, 0] = t[0, 1] / t[1, 1];
		s[0, 1] = t[0, 0] - ((t[0, 1] * t[1, 0]) / (t[0, 1]));
		s[1, 0] = 1 / t[1, 1]
		s[1, 1] = -t[1, 0] / t[1, 1]

		return s

	def Scir(self, epsilonR, waveguideR, lossTangent, layerThickness, modeEigenvalue, frequency, freeWave=False):
		"""
		returns S-matrix for a layer with given thickness, relative permitivity, 
		loss tangent for a given mode in a waveguide with defined radius at a certain frequency
		
		Nickel, Thumm (Rho and S-matrix)
		Hartfuss, Geist (T-matrix concatenation)
		"""
		e0 = 8.849E-12
		u0 = 4 * np.pi * 1E-7
		c = 1 / np.sqrt(e0 * u0)

		# Complex permitivitty OK
		er = epsilonR * (1 - 1j * lossTangent)

		# Propagation constant gamma
		k_0 = 2 * np.pi * (frequency / c)  # phase constant in vacuum

		if freeWave:
			k_cmn = k_0
		# modification for free wave TEM
			rho_mn = (1 - np.sqrt(epsilonR)) / (1 + np.sqrt(epsilonR))
		else:
			# modification phase constant resulting from waveguide walls OK
			k_cmn = modeEigenvalue / waveguideR
			# Reflection coefficient rho for TE OK
			rho_mn = (np.sqrt(1 - np.power((k_cmn / k_0), 2)) - np.sqrt(epsilonR - np.power((k_cmn / k_0), 2))) / (np.sqrt(1 - np.power((k_cmn / k_0), 2)) + np.sqrt(epsilonR - np.power((k_cmn / k_0), 2)))
			
		#gm_mn now OK!
		gm_mn = np.sqrt((np.power(k_cmn, 2) - er*np.power(k_0, 2)))
		
		# Scattering matrix S_m - PROBLEM HERE
		ab0 = 1 - np.power(rho_mn, 2) * np.exp(-2 * gm_mn * layerThickness)
		#ab0 OK
		a1 = rho_mn * (1 - np.exp(-2 * gm_mn * layerThickness))
		b1 = (1 - np.power(rho_mn, 2)) * np.exp(-1 * gm_mn * layerThickness)
		b2 = (1 - np.power(rho_mn, 2)) * np.exp(-1 * gm_mn * layerThickness)
		a2 = rho_mn * (1 - np.exp(-2 * gm_mn * layerThickness))
		S_m = np.zeros((2, 2), dtype=complex)
		S_m[0, 0] = a1 / ab0;
		S_m[0, 1] = b2 / ab0;
		S_m[1, 0] = b1 / ab0;
		S_m[1, 1] = a2 / ab0;

		return S_m

	def calculateFrequencyResponse(self, startF, stopF, stepF, buildStack=False,  useError=False, quiet=False,  incidentAngle=0, modeEigenvalue=1.8412, waveguideR=20E-3, freeWave=False, output=""):			
		
		if buildStack:
			self.buildStack(useError=useError)

		if len(self.layers['thickness']) == 0:
			raise ValueError("no layers defined. Cannot calculate frequency response. Add layers first.")
			
		if not incidentAngle == 0:
			angleRadians = incidentAngle * (np.pi/180)
			opticalLengths=[]
			#assume first border is Air/Material
			for i in range(0,len(self.layers['thickness'])):
				if i == 0:
					alpha=np.arcsin((np.sqrt(self.layers['epsilonR'][i])/1.00027)*np.sin(angleRadians))
				else:
					alpha=np.arcsin((np.sqrt(self.layers['epsilonR'][i])/np.sqrt(self.layers['epsilonR'][i-1]))*np.sin(angleRadians))
				opticalLengths.append(self.layers['thickness'][i]/np.cos(alpha))
			self.layers['thickness']=opticalLengths
				
		if not quiet:
			print("Calculating frequency response...")
				
		epsilon0 = 8.854E-12
		mu0 = 4 * np.pi * 1E-7
		c = 1 / np.sqrt(mu0 * epsilon0)

		# eigenvalue for TE11 mode 1.8412

		# preallocate
		if startF == stopF:
			fPoints = 1
			S_11 = np.zeros(1)
			S_21 = np.zeros(1)
			freqs = np.zeros(1)
		else:
			fPoints = int((stopF - startF) / stepF);
			S_11 = np.zeros(fPoints, dtype=complex)
			S_21 = np.zeros(fPoints, dtype=complex)
			freqs = np.zeros(fPoints)

		# calculate S11 and S22 for each frequency
		m = len(self.layers['thickness'])
		Sm = np.zeros((2, 2, m), dtype=complex)
		Tm = np.zeros((2, 2, m), dtype=complex)
		f = startF
		fCount = 0
		if startF == stopF:
			freqs = [startF]
		else:
			freqs = np.linspace(startF,stopF,int((stopF-startF)/stepF)+1)
			freqs = freqs[:-1]
		
		# compute cascaded S-matrix for each frequency
		for f in freqs:
			for n in range(0, m):
				materialProps=self.getMaterial(label=self.layers['material'][n])
				Sm[:, :, n] = self.Scir(epsilonR=materialProps['epsilonR'], layerThickness=self.layers['thickness'][n], lossTangent=materialProps['lossTangent'], modeEigenvalue=modeEigenvalue, waveguideR=waveguideR, freeWave=freeWave, frequency=f)
				Tm[:, :, n] = self.s2t(Sm[:, :, n])

			Tm_cas = Tm[:, :, m-1]
					
			# cascade matrices
			for n in range(1, m):
				Tm_cas = np.dot(Tm[:, :, m - 1 - n],Tm_cas)

			# cascaded S-matrix for each specific frequency
			Sm_cas = self.t2s(Tm_cas)

			S_21[fCount] = Sm_cas[1, 0]
			S_11[fCount] = Sm_cas[0, 0]
			fCount += 1
			
		self.WFR['T'] = 10*np.log10(np.power(np.abs(S_21), 2))
		self.WFR['R'] = 10*np.log10(np.power(np.abs(S_11), 2))
		self.WFR['A'] = 10*np.log10(1 - np.power(np.abs(S_21), 2) - np.power(np.abs(S_11), 2))
		self.WFR['freqs'] = np.divide(freqs, 1E9)
		
		if len(output) > 0:
			np.save(output,self.WFR)

		return self.WFR
		
	def calculateFrequencyResponseMinimal(self, f=140E9, stack=[], useError=False, modeEigenvalue=1.8412, waveguideR=20E-3, freeWave=False):			
		
		if not len(stack)==0:
			self.stack = stack
			
		self.buildStack(useError=useError)

		if len(self.layers['thickness']) == 0:
			raise ValueError("no layers defined. Cannot calculate frequency response. Add layers first.")
				
		epsilon0 = 8.854E-12
		mu0 = 4 * np.pi * 1E-7
		c = 1 / np.sqrt(mu0 * epsilon0)

		# eigenvalue for TE11 mode 1.8412

		# preallocate
		fPoints = 1
		S_11 = np.zeros(1,dtype=complex)
		S_21 = np.zeros(1,dtype=complex)
		freqs = np.zeros(1)

		# calculate S11 and S22 for each frequency
		m = len(self.layers['thickness'])
		Sm = np.zeros((2, 2, m), dtype=complex)
		Tm = np.zeros((2, 2, m), dtype=complex)
		
		# compute cascaded S-matrix for each frequency
		for n in range(0, m):
			Sm[:, :, n] = self.Scir(epsilonR=self.layers['epsilonR'][n], layerThickness=self.layers['thickness'][n], lossTangent=self.layers['lossTangent'][n], modeEigenvalue=modeEigenvalue, waveguideR=waveguideR, freeWave=freeWave, frequency=f)
			Tm[:, :, n] = self.s2t(Sm[:, :, n])

		Tm_cas = Tm[:, :, m-1]
				
		# cascade matrices
		for n in range(1, m):
			Tm_cas = np.dot(Tm[:, :, m - 1 - n],Tm_cas)

		# cascaded S-matrix for each specific frequency
		Sm_cas = self.t2s(Tm_cas)

		S_21[0] = Sm_cas[1, 0]
		S_11[0] = Sm_cas[0, 0]

		T = 10*np.log10(np.power(np.abs(S_21), 2))
		R = 10*np.log10(np.power(np.abs(S_11), 2))
		A = 10*np.log10(1 - np.power(np.abs(S_21), 2) - np.power(np.abs(S_11), 2))
		freqs = np.divide(freqs, 1E9)

		return {'T':T,'R':R,'A':A,'freqs':freqs}
		
	def calculateFrequencyResponse_helper(self,args):
		return self.calculateFrequencyResponse(*args)
		
	def extractPOI(self, pois=[], result=[], quiet=True, output=""):
		if (len(result) == 0):
			result=self.WFR
		elif not type(result) is dict:
			raise ValueError("Invalid frequency response result. Expected dict for input argument 'result'.")
		
		#check input argument result
		#check if all expected keys are present in the dict
		requiredKeys = ['freqs','T','R','A']
		for elem in requiredKeys:
			if not elem in list(result.keys()):
				raise ValueError("One or more required keys are missing in the dict given for input argument 'result'.")
		
		if not ((type(pois) is np.array) | (type(pois) is list)):
			pois = [pois]
		if not type(pois) is np.array:
			pois = np.array(pois)
			
		poiFreqs = []
		poiT = []
		poiR = []
		poiA = []
		poiDevT = []
		poiDevR = []
		poiDevA = []
		
		if 'devT' in list(result.keys()):
			poiWithErrors = True
		else:
			poiWithErrors = False
			
		for poi in pois:
			cand=np.where(result['freqs'] == poi)
			if len(cand[0]) > 0:
				poiFreqs.append(poi)
				poiT.append(result['T'][cand[0]][0])
				poiR.append(result['R'][cand[0]][0])
				poiA.append(result['A'][cand[0]][0])
				if poiWithErrors:
					poiDevT.append(result['devT'][cand[0]][0])
					poiDevR.append(result['devR'][cand[0]][0])
					poiDevA.append(result['devA'][cand[0]][0])
				if not quiet:
					print("Frequency:"+str(poi)+" GHz")
					if poiWithErrors:
						print("Transmission:"+str(poiT[-1])+" dB +/- "+str(poiDevT[-1])+" dB")
					else:
						print("Transmission:"+str(poiT[-1])+" dB")
					print("")
					
		outputResult=dict({'poiFreqs':poiFreqs,'poiT':poiT,'poiR':poiR,'poiA':poiA,'poiDevT':poiDevT,'poiDevR':poiDevR,'poiDevA':poiDevA})
					
		if len(output)>0:
			np.save(output,outputResult)
		
		return outputResult
		
	def plot(self,marker=[],result=[],plot=[{'key':'T','label':'Transmission [dB]','color':'red','style':'.-'}],labelSize=18,tickSize=14, sizepx=[640,480], dpi=75, output=""):
		if len(result) == 0:
			result=self.WFR
		elif not type(result) is dict:
			raise ValueError("Invalid data type given for input argument 'result'. Expected dict.")
			
		requiredKeys = ['freqs','T','R','A']
		for elem in requiredKeys:
			if not elem in list(result.keys()):
				raise ValueError("One or more required keys are missing in the dict given for input argument 'result'.")
				
		fig=plt.figure(figsize=(sizepx[0]/dpi,sizepx[1]/dpi),dpi=dpi)
		
		for elem in plot:
			if not elem['key'] in list(result.keys()):
				raise ValueError("Invalid key given for result to plot.")
			plt.plot(result['freqs'],result[elem['key']],linewidth=elem['width'],linestyle=elem['style'],color=elem['color'],label=elem['label'])
		if 'devT' in list(result.keys()):
			plt.errorbar(result['freqs'], result['T'], yerr=result['devT'], fmt='.', color="blue", ecolor="red")
			lb=np.nanmin(result['T'])-np.nanmax(result['devT'])
		else:
			lb=np.nanmin(result['T'])
		for poi in marker:
			plt.plot(np.full(10,poi['f']),np.linspace(1.1*lb,0,10),color=poi['color'],linestyle=poi['style'],linewidth=poi['width'])
		plt.ylim([1.1*lb,0])
		plt.xticks(fontsize=tickSize)
		plt.yticks(fontsize=tickSize)
		plt.xlabel("Frequency [GHz]",fontsize=labelSize)
		if len(plot) > 1:
			plt.legend()
		else:
			plt.ylabel(plot[0]['label'],fontsize=labelSize)
		plt.grid()
		plt.show()
		
	def clearLayers(self):
		self.layers = dict({'thickness':[],'material':[],'epsilonR':[],'lossTangent':[]})
		self.clearWFR()
	
	def clearWFR(self):
		self.WFR = dict({'freqs':[],'T':[],'R':[],'A':[]}) 

	def addLayer(self, thickness, materialLabel, useError=False, absSpread=False, spread=0):
		"""
		adds a layer to the multi disc window
		"""
		if not materialLabel in self.materials.keys():
			print("Material undefined.")
		else:
			self.layers['epsilonR'].append(self.materials[materialLabel]['epsilonR'])
			self.layers['lossTangent'].append(self.materials[materialLabel]['lossTangent'])
			self.layers['material'].append(materialLabel)
			if useError:
				if not absSpread:
					scalef=(spread * thickness) / 3
				else:
					scalef=(spread)
				self.layers['thickness'].append(np.random.normal(loc=thickness, scale=scalef))
			else:
				self.layers['thickness'].append(thickness)
	
	def clearStackElement(self):
		self.stackElement = dict({'thickness':[],'material':[],'useError':[],'spread':[],'absSpread':[]})
			
	def addStackElement(self,num,useError=False,elem=[]):
		if len(elem) == 0:
			elem=self.stackElement
		for num in range(0,num):
			for i in range(0,len(elem['useError'])):
				if elem['useError'][i] & useError:
					self.addLayer(elem['thickness'][i],elem['material'][i],useError=True,absSpread=elem['absSpread'][i],spread=elem['spread'][i])
				else:
					self.addLayer(elem['thickness'][i],elem['material'][i])
		
	def addMaterial(self,label,epsilonR,lossTangent,groupLabel):
		"""
		adds a material
		"""
		
		materialLabels = list(self.materials.keys())
		if not label in materialLabels:
			materialProperties={'epsilonR':epsilonR,'lossTangent':lossTangent,'groupLabel':groupLabel}
			material={label : materialProperties}
			self.materials.update(material)
		else:
			print("Material '"+label+"' already exists. Overwriting properties.")
			self.materials[label]['epsilonR'] = epsilonR
			self.materials[label]['lossTangent'] = lossTangent
			self.materials[label]['groupLabel'] = groupLabel
			
	def getMaterial(self,label="",groupLabel=""):
		if not len(label)==0:
			if label in self.materials.keys():
				return self.materials[label]
			else:
				raise KeyError("Invalid material label '"+label+"'.")
		if not len(groupLabel) == 0:
			result=[]
			for key in self.materials.keys():
				if groupLabel == self.materials[key]['groupLabel']:
					result.append(dict({key:self.materials[key]}))
			return result
			
	def listMaterial(self):
		for label in self.materials.keys():
			print("Material:"+label)
			print("Group:"+self.materials[label]['groupLabel'])
			print("n:"+str(np.sqrt(self.materials[label]['epsilonR'])))
			print("EpsilonR:"+str(self.materials[label]['epsilonR']))
			print("LossTangent:"+str(self.materials[label]['lossTangent']))
			print("")
	
	def getLayers(self,material=""):
		if len(material)==0:
			return self.layers['thickness']
		else:
			if not type(material) is list:
				material = [material]
			result=[]
			for label in material:
				cand=np.where(np.array(self.layers['material'])==label)[0]
				result.extend(list(np.array(self.layers['thickness'])[cand]))
			return result
			
	def listSpecs(self,materialGroup='disc'):
		print("Filter length:"+str(np.sum(self.layers['thickness'])))
		sheetMaterial=self.getMaterial(groupLabel=materialGroup)
		sheetMaterial = [list(elem.keys())[0] for elem in sheetMaterial]
		print("No of layers:"+str(len(self.getLayers(material=sheetMaterial))))
		print("Layer thickness:"+str(self.getLayers(material=sheetMaterial)))
		print("Airgaps:"+str(self.getLayers(material="Air")))
		
	def defineStackElement(self,material=[],thickness=[], spread=[], useError=[], absSpread=[],label=""):
		if len(material) == 0:
			raise ValueError("No materials defined for stack element layers.")
		if len(thickness) == 0:
			raise ValueError("No thickness defined for stack element layers.")
		if not len(useError) == 0:
			if len(spread) == 0:
				raise ValueError("No spreads defined for thickness errors.")
			if len(absSpread) == 0:
				raise ValueError("Not defined whether absolute or relative spreads should be used.")
		else:
			withError = np.full(len(thickness),False)
			spread = np.full(len(thickness),0)
			absSpread = np.full(len(thickness),False)
		
		self.stackElement['thickness']=thickness
		self.stackElement['material']=material
		self.stackElement['useError']=useError
		self.stackElement['absSpread']=absSpread
		self.stackElement['spread']=spread
		self.stackElement['label']=label
		
		return self.stackElement.copy()
		
	def defineStack(self,elem=[],qty=[],label=""):
		if len(elem) == 0 | len(qty) == 0:
			raise ValueError("Stackelement list or quantity list is empty.")
		if not len(elem) == len(qty):
			raise ValueError("Stackelement list and quantity list have different size.")
	
		self.stack['stackElement'] = elem
		self.stack['quantity'] = qty
		self.stack['label'] = label
		self.stack['length'] = 0
		for i in range(0,len(qty)):
			self.stack['length'] += qty[i]*np.sum(elem[i]['thickness']) 
		
		return self.stack.copy()
		
	def buildStack(self,useError=False):
		if len(self.stack['stackElement']) == 0:
			raise ValueError("You have to define a stack first.")
			
		self.clearLayers()
		for elem in range(0,len(self.stack['stackElement'])):
			self.addStackElement(elem=self.stack['stackElement'][elem],num=self.stack['quantity'][elem],useError=useError)

	def errorStudyWithConstraints_helper(self,args):
		return self.errorStudyWithConstraintsMinimal(*args)
		
	def errorStudy(self,fRange=[120E9,150E9,50E6],stack=[],quiet=False,incidentAngle=0,extractPOIs=False,pois=[128.00,130.00,140.0,133.0,137.3,142.8,135.00],numCalc=10, parallelize=False, numCores=2, output="",chunksize=1):
		
		"""
		Runs the calculation for the filter characteristic $numCalc times with normal distributed random errors for layer
		thickness and then averages over transmission, reflection etc.
		In particular the standard deviation is calculated for every frequency in the transmission curve to 
		quantify the sensitivity of the filter design to manufacturing tolerances
		"""		
				
		if not len(stack) == 0:
			self.stack = stack
			
		sigmas=[]
		means=[]
		freq=[]
		
		if fRange[0] == fRange[1]:
			stepWidth=1
			steps=len(np.linspace(fRange[0],fRange[0],stepWidth))
		else:
			stepWidth=int((fRange[1]-fRange[0])/fRange[2])
			steps=len(np.linspace(fRange[0],fRange[1],stepWidth))
			
		raw=dict({'T':np.zeros((numCalc,steps)), 'R':np.zeros((numCalc,steps)), 'A':np.zeros((numCalc,steps))})
		
		avail_cores = multiprocessing.cpu_count()
		
		if not quiet:
			print(str(avail_cores)+" cores available.")
		if parallelize & (not quiet):
			if numCores > avail_cores:
				print("Warning: demanded "+str(numCores)+" cores but only "+str(avail_cores)+" cores available.")
				numCores = avail_cores
				print("Using "+str(numCores)+" cores...")
				results = Parallel(n_jobs=numCores)(delayed(calculateFrequencyResponse)(i) for i in inputs)
		else:
			numCores=1
								
		if not quiet:
			print("Using "+str(numCores)+" core(s).")
		
		if not quiet:
			print("Performing random error sweep...")		
			pb=progressBar(total=numCalc)		
			pb.show()	
			
		if parallelize:
				
			job_args = [(fRange[0],fRange[1],fRange[2],True, True, True, incidentAngle) for i in range(numCalc)] 
			with multiprocessing.Pool(processes=numCores) as pool:
				poolResult=pool.map_async(self.calculateFrequencyResponse_helper,job_args,chunksize)
				while (True):
					if (poolResult.ready()): break
					if not quiet:
						tasksLeft=pool._cache[list(pool._cache.keys())[0]]._number_left
						chunkSize=pool._cache[list(pool._cache.keys())[0]]._chunksize
						pb.update(numCalc - chunkSize*tasksLeft)
					time.sleep(1)
				pool.close()
				pool.join()
			if not quiet:
				pb.finish(timeElapsed=True)
			poolResult=poolResult.__dict__['_value']
				
			for itervar in range(0,numCalc):
				raw['T'][itervar,:] = poolResult[itervar]['T']
				raw['R'][itervar,:] = poolResult[itervar]['T']
				raw['A'][itervar,:] = poolResult[itervar]['T']
			result=poolResult[0]

		else:
			itervar=0
			while itervar < numCalc:
				result = self.calculateFrequencyResponse(fRange[0],fRange[1],fRange[2],useError=True, buildStack=True, quiet=True, incidentAngle=incidentAngle)
				raw['T'][itervar,:] = result['T']
				raw['R'][itervar,:] = result['R']
				raw['A'][itervar,:] = result['A']
				if not quiet:
					pb.update(itervar)
				itervar += 1
			if not quiet:
				pb.finish(timeElapsed=True)

		freq=result['freqs']	
		result=dict({'freqs':freq,'devT':np.zeros(steps),'devA':np.zeros(steps),'devR':np.zeros(steps),\
		'T':np.zeros(steps),'A':np.zeros(steps),'R':np.zeros(steps)})
			
		if not quiet:
			print("Calculating mean and deviation...")
			pb=progressBar(total=steps)		
			pb.show()
		for stepvar in range(0,steps):
			result['T'][stepvar]=np.mean(raw['T'][:,stepvar])
			result['devT'][stepvar]=np.std(raw['T'][:,stepvar])
			result['R'][stepvar]=np.mean(raw['R'][:,stepvar])
			result['devR'][stepvar]=np.std(raw['R'][:,stepvar])
			result['A'][stepvar]=np.mean(raw['A'][:,stepvar])
			result['devA'][stepvar]=np.std(raw['A'][:,stepvar])
			#if not quiet:
			#	pb.update(stepvar)
			stepvar += 1
			
		if not quiet:
			pb.finish()
			
		if extractPOIs:
			poiResult=self.extractPOI(pois=pois,result=result,quiet=quiet)
		else:
			poiResult=dict({})
		
		if len(output)>0:
			np.save(output,[result,poiResult])
		
		return [result,poiResult]
		
	def errorStudyMinimal(self,f,stack,numCalc=1000):
		
		"""
		minimal version for parallelization
		"""		
				
		if not len(stack) == 0:
			self.stack = stack
		
		if not type(f) is list:
			f = [f]
			
		sigmas=[]
		means=[]
		freq=[]
		raw=dict({'T':np.zeros((numCalc,len(f))), 'R':np.zeros((numCalc,len(f))), 'A':np.zeros((numCalc,len(f)))})
		
		for fIndex in range(0,len(f)):
			itervar=0
			while itervar < numCalc:
				result = self.calculateFrequencyResponseMinimal(f[fIndex],stack, useError=True)
				raw['T'][itervar,fIndex] = result['T']
				raw['R'][itervar,fIndex] = result['R']
				raw['A'][itervar,fIndex] = result['A']
				itervar+=1

		result=dict({'freqs':f,'devT':np.zeros(len(f)),'devA':np.zeros(len(f)),'devR':np.zeros(len(f)), 'T':np.zeros(len(f)),'A':np.zeros(len(f)),'R':np.zeros(len(f))})
		
		for fIndex in range(0,len(f)):
			result['T'][fIndex]=np.mean(raw['T'][:,fIndex])
			result['devT'][fIndex]=np.std(raw['T'][:,fIndex])
			result['R'][fIndex]=np.mean(raw['R'][:,fIndex])
			result['devR'][fIndex]=np.std(raw['R'][:,fIndex])
			result['A'][fIndex]=np.mean(raw['A'][:,fIndex])
			result['devA'][fIndex]=np.std(raw['A'][:,fIndex])
		
		return result
		
	def errorStudyWithConstraintsMinimal(self,f,stack,goals,numCalc=1000):
		
		"""
		minimal version for parallelization
		"""		
				
		if not len(stack) == 0:
			self.stack = stack
		
		if not type(f) is list:
			f = [f]
			
		sigmas=[]
		means=[]
		freq=[]
		raw=dict({'T':np.zeros((numCalc,len(f))), 'R':np.zeros((numCalc,len(f))), 'A':np.zeros((numCalc,len(f)))})
		
		for fIndex in range(0,len(f)):
			itervar=0
			while itervar < numCalc:
				result = self.calculateFrequencyResponseMinimal(f[fIndex],stack, useError=True)
				raw['T'][itervar,fIndex] = result['T']
				raw['R'][itervar,fIndex] = result['R']
				raw['A'][itervar,fIndex] = result['A']
				itervar+=1

		result=dict({'freqs':f,'devT':np.zeros(len(f)),'devA':np.zeros(len(f)),'devR':np.zeros(len(f)), 'T':np.zeros(len(f)),'A':np.zeros(len(f)),'R':np.zeros(len(f))})
		
		for fIndex in range(0,len(f)):
			result['T'][fIndex]=np.mean(raw['T'][:,fIndex])
			result['devT'][fIndex]=np.std(raw['T'][:,fIndex])
			result['R'][fIndex]=np.mean(raw['R'][:,fIndex])
			result['devR'][fIndex]=np.std(raw['R'][:,fIndex])
			result['A'][fIndex]=np.mean(raw['A'][:,fIndex])
			result['devA'][fIndex]=np.std(raw['A'][:,fIndex])
			
		#check goals
		goalResult=[]
		for goalIndex in range(0,len(goals)):
			elem=goals[goalIndex]
			for qtyIndex in range(0,len(elem['qty'])):
				if elem['target'][qtyIndex] == 'max': 
					if result[elem['qty'][qtyIndex]][goalIndex] > elem['threshold'][qtyIndex]:
						goalResult.append(True)
					else:
						goalResult.append(False)
				if elem['target'][qtyIndex] == 'min':
					if result[elem['qty'][qtyIndex]][goalIndex] < elem['threshold'][qtyIndex]:
						goalResult.append(True)
					else:
						goalResult.append(False)
		#check if all goals met = candidate found
		if np.all(goalResult):
			return [stack,result]
		else:
			return []
		
	#Description:
	#------------
	#Compares multiple defined stacks with respect to the points of interest
	#
	#Parameters:
	#----------
	#errorAnalysis: True/False (Default: True)
	#If True, an error analysis is performed for every given stack, so error bars are also obtained for the points of interest
		
	def compareStacks(self,stacks=[], fRange=[125E9,145E9,100E6], errorAnalysis=True, numCalc=100, poi=[142.8,137.3,140.0], parallelize=True, numCores=2, plot=True):
		stackResults = []
		print("Comparing stacks...")
		if not errorAnalysis:
			pb = progressBar(total=len(stacks))
			pb.show()
			itervar = 0
		for stack in stacks:
			self.clearLayers()
			self.clearWFR()
			self.defineStack(elem=stack['stackElement'],qty=stack['quantity'])
			if errorAnalysis:
				self.buildStack(useError=True)
				[studyResult,poiResult]=self.errorStudy(fRange=fRange, pois=poi, numCalc=numCalc, extractPOIs=True, parallelize=parallelize, numCores=numCores, quiet=False)
			else:
				result = self.calculateFrequencyResponse(fRange[0],fRange[1],fRange[2],useError=True, buildStack=True, quiet=True)
				poiResult=self.extractPOI(pois=poi,result=result,quiet=True)
				itervar += 1
				pb.update(itervar)
			stackResults.append(poiResult)
		if not errorAnalysis:
			pb.finish()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.set_prop_cycle(cycler('color', ['r','g','b','y','m','c']))
		matplotlib.rcParams['errorbar.capsize'] = 3
		for stackIndex in range(0,len(stackResults)):
			if errorAnalysis:	
				ax.errorbar(stackResults[stackIndex]['poiFreqs'],stackResults[stackIndex]['poiT'],yerr=stackResults[stackIndex]['poiDevT'],fmt="x",label=stacks[stackIndex]['label'])
			else:
				ax.scatter(stackResults[stackIndex]['poiFreqs'],stackResults[stackIndex]['poiT'],label=stacks[stackIndex]['label'])								
		ax.set_xlabel("Frequency [GHz]")
		ax.set_ylabel("Transmission [dB]")
		plt.legend()
		plt.grid()
		plt.show()
		return stackResults
		
	def compareCandidates(self, result, labels=[], qty=''):
		fig = plt.figure()
		freqs = result['result'][0]['freqs']
		numFreqs = len(freqs)
		ax = fig.add_subplot(111)
		matplotlib.rcParams['errorbar.capsize'] = 3
		for freqIndex in range(0,numFreqs):
			#x=[int(elem['label']) for elem in result['stackElement']]
			x=list(np.linspace(0,len(result['stackElement']),num=len(result['stackElement'])))
			y=[elem[qty][freqIndex] for elem in result['result']]
			yerr=[elem['dev'+qty][freqIndex] for elem in result['result']]
			#candidates=[labels[elem] for elem in x]
			ax.errorbar(y,x,xerr=yerr,fmt="o",label=str(freqs[freqIndex]/1E9)+" GHz")
		ax.set_ylabel("Configuration #")
		if qty == "T": xlbl="Transmission [dB]"
		if qty == "R": xlbl="Reflection [dB]"
		if qty == "A": xlbl="Absorption [dB]"
		ax.set_xlabel(xlbl)
		if not len(labels) == 0:
			plt.yticks(x, labels, rotation='horizontal')
		plt.legend()
		plt.grid(which="both")
		plt.show()
	
	#Description:
	#Time consuming method! Tries to optimize given stack with respect to given goals (constraints)
	#
	#Example for sweep definition:
	#-----------------------------
	#To define parameter sweeps, a list of strings is given
	#	
	#Each dict contains a command-like structure:
	#{'label':'elem1','layer':0,'thickness':[start,stop,step]}
	#to explicitly sweep absolute thickness of a single layer in the stack element
	#
	# -or-
	#
	#{'elem1','count',[start,stop,step]}
	#to use different number of stack elements in stack
	
	def findStackElement(self, numCores=4, parallelize=True, numCalc=100, stackElement=[], goals=[], sweep=[],chunksize=1):
		
		if len(goals) == 0:
			raise ValueError("You have to define goals for the stack optimization.")
		if (not type(goals) is list) | (not type(goals[0]) is dict):
			raise TypeError("Wrong data type given for input parameter goals. Expected list of dicts.")		
		if (len(stackElement) == 0) & (len(self.stackElement) == 0):
			raise ValueError("No stack element given for optimization.")
		elif len(stackElement) == 0:
			stackElement = self.stackElement
		if not type(sweep) is dict:
			raise ValueError("Expected dict for sweep definition.")

		print("Generating test configurations...")
		
		tmp=[]
		for r in sweep['range']:
			tmp.append(r)
		r=sweep['count']
		tmp.append(r)
		
		sweepSpace = list(itertools.product(*tmp))
		stackSpace = []

		c=0
		for config in sweepSpace:
			foo=np.array(list(config))
			elem=self.defineStackElement(material=stackElement['material'], thickness=foo[np.array(sweep['layer'])],  spread=[0.1E-3/3,0.1E-3/3], useError=[True,True], absSpread=[True,True])
			stack=self.defineStack(elem=[elem],qty=[int(foo[-1])], label=str(c))
			stackSpace.append(stack.copy())
			c+=1
			
		print(str(len(stackSpace))+" configurations generated.")
		print(str(len(stackSpace)*len(goals)*numCalc)+" calculations in total.")
		
		avail_cores = multiprocessing.cpu_count()
		print(str(avail_cores)+" cores available.")
		if parallelize:
			if numCores > avail_cores:
				print("Warning: demanded "+str(numCores)+" cores but only "+str(avail_cores)+" cores available.")
				numCores = avail_cores
		else:
			numCores=1
								
		print("Using "+str(numCores)+" core(s).")
		
		#calculate frequency response for all given configurations
		finalists=[]
		finalistsResult=[]
		configResults=[]
		freqs = [elem['f'] for elem in goals]

		pb=progressBar(total=len(stackSpace))
		i=0
		total=len(stackSpace)
		
		if parallelize:
			t1=time.time()
			job_args=[]
			#prepare job arguments
			for config in stackSpace:
				job_args.append((freqs,config,goals,numCalc)) 
			with multiprocessing.Pool(processes=numCores) as pool:
				poolResult=pool.map_async(self.errorStudyWithConstraints_helper,job_args,chunksize)
				while (True):
					if (poolResult.ready()): break
					tasksLeft=pool._cache[list(pool._cache.keys())[0]]._number_left
					chunkSize=pool._cache[list(pool._cache.keys())[0]]._chunksize
					pb.update(total - chunkSize*tasksLeft)
					time.sleep(1)
				pool.close()
				pool.join()
			pb.finish(timeElapsed=True)
				
			poolResult=poolResult.__dict__['_value']
			
			#check if all goals met = candidate found
			for elem in poolResult:
				if len(elem) == 2:
					finalists.append(elem[0])	
					finalistsResult.append(elem[1])
				
			pb.finish(timeElapsed=True)
		else:
			t1=time.time()
			for config in stackSpace:
				#calculate frequency response for POI and check goals
				goalResult = []
				configResult = []
				#frequency response is only calculated for points of interest since otherwise it would consume too much time
				studyResult=self.errorStudyWithConstraintsMinimal(freqs, config, goals, numCalc)
				
				if len(studyResult) == 2:
					finalists.append(studyResult[0])	
					finalistsResult.append(studyResult[1])
				
				i += 1
				pb.update(i)
				
		timeElapsed=time.time()-t1
		totalCalculations=len(stackSpace)*len(freqs)*numCalc
		
		print("")
		print(str(len(finalists))+" candidates found with given constraints.")
		print("Total time elapsed:"+str(timeElapsed)+" seconds.")
		print("Configurations tested:"+str(len(stackSpace)))
		print("Total number of calculations:"+str(totalCalculations))
		print(str(np.round(totalCalculations/timeElapsed,2))+" calculations per second")
				
		return {'stackElement':finalists,'result':finalistsResult}