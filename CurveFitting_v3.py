# Curve-fitting: 
# 1. use least-squares approach to find curves for each animal across temperatures
# 2. extend this approach to different events?
# 3. determine if there are trends in the normalized data (see if they are stat sig different from vertical) -> can try F statistic for these

import math, os, re, string, sys
import matplotlib as mpl
import pylab as pyl
import numpy as np
import scipy.stats as st

def SxySqrd(coordinates):
	n = len(coordinates)
	b,a = findSlopeAndIntercept(coordinates)
	runningSum = 0
	for (x,y) in coordinates:
		runningSum += (y-(b*x+a))**2
	return runningSum/(n-2)

def findSlopeAndIntercept(coordinates):
	n = len(coordinates)
	vectorX = []
	vectorY = []
	runningProduct = 0
	runningSquare = 0
	for (x,y) in coordinates:
		vectorX.append(x)
		vectorY.append(y)
		runningProduct += x*y
		runningSquare += x**2
	b = (n*runningProduct - sum(vectorX) * sum(vectorY))/(n*runningSquare-(sum(vectorX))**2)
	a = (sum(vectorY) * runningSquare - sum(vectorX) * runningProduct)/(n*runningSquare-(sum(vectorX))**2)
	return b, a
	
def findVariance(vectorA):
	n = len(vectorA)
	runningSum = 0
	for i in vectorA:
		runningSum += i**2
	return (runningSum-(sum(vectorA)**2/n))/(n-1)

## Enter in the data
speciesNames = [['Mj', 'Dmoj'],['An', 'Dana'],['Er','Dere'],['Me', 'Dmel'],['Me18','Dm18'],['Ps','Dpse'],['Pe','Dper'],['Si','Dsim'],['Se','Dsec'],['Wi','Dwil'],['Vi','Dvir'],['Ya','Dyak']]
speciesOrdering = ['Dsim', 'Dsec','Dmel','Dm18','Dere','Dyak','Dana','Dper','Dpse','Dwil','Dmoj','Dvir']

def changeName(i):
	for (shortName, fullName) in speciesNames:
		if i == shortName:
			return fullName
	newName = ''
	if '-merge' not in sys.argv:
		newName = i
	if i in melStrains:
		return 'Dmel'+newName
	elif i in yakStrains:
		return 'Dyak'+newName
	elif i in simStrains:
		return 'Dsim'+newName
def changeTemp(j):
	if len(j) == 4: return j
	elif len(j) == 2: return j+'.0'
	
	
RaleighStrains = ['R303','R329','R379','R380','R437','R705']
HeterochronicStrains = ['12414','12592','12622','16532']
melStrains = RaleighStrains + HeterochronicStrains
CYyakStrains = ['CY01A','CY04B','CY08A','CY15B4','CY17C','CY22B','CY23A','CY21B','CY27A9']
NYyakStrains = ['NY20','NY48','NY55','NY56','NY60','NY62','NY65','NY73','NY81','NY85']
yakStrains = CYyakStrains + NYyakStrains 
NSsimStrains = ['NS37','NS67','NS79','NS113','NS137']
MDsimstrains = ['MD15','MD63','MD73','MD105','MD106','MD199','MD221','MD233','MD251']
simStrains = NSsimStrains + MDsimstrains

otherSpecies = ['An','Er','Me','MeAcc','Mj','Pe','Ps','Se','Si','Vi','Wi','Ya']

speciesListDefault = melStrains + yakStrains + simStrains + otherSpecies	
temperatureListDefault = ['12.5','15','17.5','17','20','22.5','25','27.5','30','32.5','35']

SpList = []; Tlist = []

for species in speciesOrdering:
	for i in sys.argv:
		# if i in speciesListDefault:
		if changeName(i) == species:
			# print changeName(i)
			SpList.append(changeName(i))
for i in sys.argv:
	if i in temperatureListDefault:
		Tlist.append(changeTemp(i))
		
comboDict = {}
animalList = []

timeLapseInfoFile = open('TimeLapseInstructions_run_' + str(sys.argv[1]) + '.txt','r')
for i in timeLapseInfoFile:
	[date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temp, scaling, zStack, filler] = re.split("[\s]",i)
	if changeName(strain) in SpList:
		speciesFull = changeName(strain) # change the species name
		if changeTemp(temp) in Tlist:
			temperature = changeTemp(temp)
			if date+position not in animalList:
				animalList.append(date+position)
				if speciesFull not in comboDict:
					comboDict[speciesFull] = {}
				if temperature not in comboDict[speciesFull]:
					comboDict[speciesFull][temperature] = []
				comboDict[speciesFull][temperature].append((float(scaling)*(int(tfTime)-int(mryTime)))/60) # record time in hours

## Curve-fitting
# set up tValues (log transform the temperature)

fig = pyl.figure(figsize=(25.5,33))
counter = 0
for species in SpList:
	counter +=1
	tValues, yValues, xValues = [], [], []
	xOriginals = []
	for T in Tlist:
		for t in comboDict[species][T]:
			tValues.append((math.log((float(T)-15)),t))
			xValues.append(math.log((float(T)-15)))
			yValues.append(t)
			xOriginals.append(float(T))
	b,a = findSlopeAndIntercept(tValues)	
	n = len(tValues)
	Rsquared = 1 - ((n-2)*SxySqrd(tValues)/((n-1)*findVariance(yValues))) # Pearson Product-Moment Correlation Coefficient
	pyl.subplot(4,3,counter)
	pyl.plot(yValues, xOriginals, marker = 's', markersize = 5, c = 'cyan', linestyle = '')
	pyl.title(species, fontsize = 30)
	pyl.ylabel('Temperature', fontsize = 20)
	pyl.xlabel('time (hours post-cellularization)', fontsize = 20)
	pyl.xlim(0,60)
	pyl.ylim(17,28)
	# find confidence interval for future estimates (lower line and upper line)
	z = np.arange(17, 28.5,.1)
	w = []
	confidenceGap50pos, confidenceGap95pos, confidenceGap50neg, confidenceGap95neg = [], [], [], []
	talpha = st.t.isf(0.05/2,(n-2)) # isf is one-tailed, so divide by 2 for 2-tailed critical value
	tbeta = st.t.isf(0.5/2,(n-2))
	standardErrorOfTheEstimate = math.sqrt((findVariance(yValues)-findVariance(xValues)*b**2)*(n-1)/(n-2))
	for u in z:
		w.append(b*math.log(u-15)+a)
		confidenceGap95pos.append(b*math.log(u-15)+a + talpha*standardErrorOfTheEstimate*math.sqrt(1.+(1./n)+(math.log(u-15)-sum(xOriginals)/n)**2/(findVariance(xValues)*(n-1))))
		confidenceGap95neg.append(b*math.log(u-15)+a - talpha*standardErrorOfTheEstimate*math.sqrt(1.+(1./n)+(math.log(u-15)-sum(xOriginals)/n)**2/(findVariance(xValues)*(n-1))))
		confidenceGap50pos.append(b*math.log(u-15)+a + tbeta*standardErrorOfTheEstimate*math.sqrt(1.+(1./n)+(math.log(u-15)-sum(xOriginals)/n)**2/(findVariance(xValues)*(n-1))))
		confidenceGap50neg.append(b*math.log(u-15)+a - tbeta*standardErrorOfTheEstimate*math.sqrt(1.+(1./n)+(math.log(u-15)-sum(xOriginals)/n)**2/(findVariance(xValues)*(n-1))))
	pyl.plot(w, z, marker = '',linestyle = '-', lw = 1.5, c = 'green')
	# try shading between these lines
	# pyl.fill_between(z,confidenceGap95neg, confidenceGap95pos, facecolor = 'yellow')
#	pyl.fill_between(confidenceGap50pos, z, confidenceGap50neg, facecolor = 'orange', edgecolor = 'orange')
	pyl.plot(confidenceGap95pos, z, marker = '', lw = 0.5, linestyle = '--', c = 'orange')
	pyl.plot(confidenceGap95neg, z, marker = '', lw = 0.5, linestyle = '--', c = 'orange')
	pyl.plot(confidenceGap50pos, z, marker = '', lw = 0.5, linestyle = '-', c = 'orange')
	pyl.plot(confidenceGap50neg, z, marker = '', lw = 0.5, linestyle = '-', c = 'orange')

	print '$t_{',species,'} =', '%.2f' % a,'%.2f' % b,' \ln({T-15})$ & ','%.3f' %Rsquared,'&$t_{',species,'} \plusminus', '%.3f' % (talpha*standardErrorOfTheEstimate/math.sqrt((findVariance(xValues)*(n-1)))), '\sqrt{','%.2f' % (1.+(1./n)*(findVariance(xValues)*(n-1))),'+ (\ln(T-15)-','%.2f' % (sum(xOriginals)/n),')^{2}}$'

pyl.savefig('curveFitting/curveFitting_allSp_logTemp.pdf', transparent = True)
pyl.close()

## compare all vs. all (statistical backing to independence between climate groups)

# go through each species and generate the relevant statistics
for species1 in SpList:
	for species2 in SpList:
		tValues1, tValues2, tValuesAll = [], [], []
		for T in Tlist:
			for t in comboDict[species1][T]:
				tValues1.append((math.log((float(T)-15)),t))
				tValuesAll.append((math.log((float(T)-15)),t))
			for t in comboDict[species2][T]:
				tValues2.append((math.log((float(T)-15)),t))
				tValuesAll.append((math.log((float(T)-15)),t))
		n1 = len(tValues1)
		n2 = len(tValues2)
		var_pooled = ((n1-2)*SxySqrd(tValues1) + (n2-2)*SxySqrd(tValues2))/(n1+n2-4)
		var_common = SxySqrd(tValuesAll)
		var_improvement = ((n1 + n2 - 2)*var_common - (n1 + n2 -4)* var_pooled)/2
		dof_common = 2
		dof_improvement = n1 + n2 -4
		F = var_improvement/var_pooled
		F_critical = st.f.isf(0.05,dof_common, dof_improvement)
		if F > F_critical:
			# print species1, species2,'F=',F
			pass
		elif species1 == species2:
			pass
		else:
			print species1, species2, 'similar'
## Are there trends in the normalized data?
# read in data from the dotLineNumbers outputs

#yes = no

#for species in SpList:
#	normalizedNumbersFile = open('dotLinePlots/dotLineNumbers_'+ str(sys.argv[1]) +'_'+species+'.txt', 'r')
#	for i in normalizedNumbersFile:
#		[T, sp, devEvent, t, stdev] = re.split("[\s]",i)
		









