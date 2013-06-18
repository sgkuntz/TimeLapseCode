# example command: python devEventAdjDataProcessing_vLinePlot1.py trial15 Pe Me 27.5 30 white
# must be run from Dropbox folder
# analyzes only manually verified data, including mry, tf, and GUI-verified data

import os, sys, re, string, math, array
import scipy.stats as st
import pylab as pyl
import matplotlib as mpl
import numpy as np

### Part 0 - Can read in info from the TimeLapseInstructions file and the TimeLapseCorrections files
global backgroundSetting
if 'white' in sys.argv:
	backgroundSetting = 'white'
	markerShape = 's'
	lineColor = 'k'
	faceColor = 'w'
elif 'black' in sys.argv:
	backgroundSetting = 'black'
	markerShape = 'o'
	lineColor = 'w'
	faceColor = 'k'


speciesNames = [['Mj', 'Dmoj'],['An', 'Dana'],['Er','Dere'],['Me', 'Dmel'],['Me18','Dm18'],['Ps','Dpse'],['Pe','Dper'],['Si','Dsim'],['Se','Dsec'],['Wi','Dwil'],['Vi','Dvir'],['Ya','Dyak']]
tempNames = [['12.5','12.5'],['15','15.0'],['17.5','17.5'],['20','20.0'],['22.5','22.5'],['25','25.0'],['27.5','27.5'],['30','30.0'],['32.5','32.5'],['17','17.0']]


RaleighStrains = ['R303','R329','R379','R380','R437','R705']
HeterochronicStrains = ['12414','12592','12622','16532']
melStrains = RaleighStrains + HeterochronicStrains
CYyakStrains = ['CY01A','CY04B','CY08A','CY15B4','CY17C','CY22B','CY23A','CY21B','CY27A9']
NYyakStrains = ['NY20','NY48','NY55','NY56','NY60','NY62','NY65','NY73','NY81','NY85']
yakStrains = CYyakStrains + NYyakStrains 
NSsimStrains = ['NS37','NS67','NS79','NS113','NS137']
MDsimstrains = ['MD15','MD63','MD73','MD105','MD106','MD199','MD221','MD233','MD251']
simStrains = NSsimStrains + MDsimstrains

otherSpecies = ['An','Er','Me','MeAcc','Me18','Mj','Pe','Ps','Se','Si','Vi','Wi','Ya']

speciesListDefault = melStrains + yakStrains + simStrains + otherSpecies
temperatureListDefault = ['12.5','15','17.5','17','20','22.5','25','27.5','30','32.5','35']
possibleEventList = ['posterior_gap','pole_bud_appears','nuclei_at_periphery','pole_cells_form','yolk_contraction','cellularization_begins','membrane_reaches_yolk',
				'pole_cells_migrate','cephalic_furrow_formation','pole_cell_invagination','cephalic_furrow_reclines','amnioproctodeal_invagination',
				'transversal_fold_formation','anterior-midgut_primordium','stomodeal_plate_forms','stomodeum_invagination','clypeolabral_lobe_forms',
				'germband_maxima','clypeolabrum_rotates','posterior_gap_before_germband_shortens','germband_retraction_starts',
				'amnioserosa','germband_retracted','dorsal_divot','clypeolabrum_retracts','anal_plate_forms','midgut_unified','clypeolabrum_ventral_lobe_even',
				'gnathal_lobes_pinch','head_involution_done','heart-shaped_midgut','convoluted_gut','muscle_contractions','trachea_fill','hatch']

SpList = []; Tlist = []


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
	for (temp, tempInt) in tempNames:
		if j == temp:
			return tempInt

for i in sys.argv:
	if i in speciesListDefault:
		print changeName(i)
		SpList.append(changeName(i))
	elif i in temperatureListDefault:
		Tlist.append(changeTemp(i))

if SpList == []:
	for species in speciesListDefault:
		SpList.append(changeName(species))
if Tlist == []:
	for temp in temperatureListDefault:
		Tlist.append(changeTemp(temp))

comboDict = {} ; animalDict = {}; eventList = ['membrane_reaches_yolk','trachea_fill']

for (dirpath, dirnames, filenames) in os.walk("."):
	for filename in filenames:
		if 'TimeLapseCorrections_' in filename:
			if 'deleted' not in filename:
				try:
					[title, trial, species, temperature] = re.split("_",filename)
				except:
					species = 'fail'
				if changeName(species) in SpList:
					if temperature[:4] in Tlist:
						correctionsFile = open(dirpath+'/'+filename,'r')
						for i in correctionsFile:
							line = re.split("[\s]",i)
							[condition, event, datePosition, timepoint, dilation, filler] = re.split("[\s]",i)
							if condition <> 'conditions':
								if timepoint <> 'NaN':
									if event not in eventList:
										eventList.append(event)
									if '-merge' in sys.argv:
										if condition[4:] in yakStrains:
											condition = condition[:4] + 'Ya'
										elif condition[4:] in melStrains:
											condition = condition[:4] + 'Me'
										elif condition[4:] in simStrains:
											condition = condition[:4] + 'Si'
									if condition[:4] in Tlist:
										condition = condition[:4] + changeName(species)
									elif condition[:2] in temperatureListDefault:
										condition = changeTemp(condition[:2]) + changeName(species)
									# print condition
									if condition not in comboDict:
										comboDict[condition] = {}
									if datePosition not in comboDict[condition]:
										comboDict[condition][datePosition] = {}
									if event not in comboDict[condition][datePosition]:
										comboDict[condition][datePosition][event] = [changeName(species), int(dilation), timepoint]
									if event not in eventList:
										eventList.append(event)
										
timeLapseInfoFile = open('TimeLapseInstructions_run_' + str(sys.argv[1]) + '.txt','r')
for i in timeLapseInfoFile:
	[date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temp, scaling, zStack, filler] = re.split("[\s]",i)
	if changeName(strain) in SpList:
		speciesFull = changeName(strain) # change the species name
		if changeTemp(temp) in Tlist:
			temperature = changeTemp(temp)
			if temperature+speciesFull not in comboDict:
				comboDict[temperature+speciesFull] = {}
			if date+'_'+position not in comboDict[temperature+speciesFull]:
				comboDict[temperature+speciesFull][date+'_'+position] = {}
			if 'membrane_reaches_yolk' not in comboDict[temperature+speciesFull][date+'_'+position]:
				comboDict[temperature+speciesFull][date+'_'+position]['membrane_reaches_yolk'] = [speciesFull, int(scaling), mryTime]
			if 'trachea_fill' not in comboDict[temperature+speciesFull][date+'_'+position]:
				comboDict[temperature+speciesFull][date+'_'+position]['trachea_fill'] = [speciesFull, int(scaling), tfTime]

comboListPrime = comboDict.keys()
comboListPrime.sort()

speciesOrdering = ['Dsim', 'Dsec','Dmel','Dere','Dyak','Dana','Dper','Dpse','Dwil','Dmoj','Dvir']

comboList = []
for species in speciesOrdering:
	for combo in comboListPrime:
		if species == combo[4:8]:
			comboList.append(combo)

print comboList
			
meanDict = {}

numberOfTemperatures = len(Tlist) # number of Temperatures
numberOfSpecies = len(SpList) # number of Species

#outputFile = open('Rinput_event_'+str(sys.argv[1])+'.txt','w')	# want to label the header as (event, t, T, Sp, sd, se, ci)
#outputFile.write('condition' + '\t' + 'event' +'\t' + 'mean' +'\t' + 'meanRev' +'\t' + 'median' +'\t' + 'stDev' +'\t' + 'stDevRev' +'\t' + 'MAD' +'\t' + 'errMin' +'\t' + 'errMax' +'\t' + 'numberOfAnimals' +'\n')
#outputTable = open('RinputTable_event_'+str(sys.argv[1])+'.txt','w')

# make SpList unique

# SpList = list(set(SpList))
# SpList.sort()

# if 'dotPlot' in sys.argv:


### Part 0b - Dot Line Graph (event determines color, condition determines y-value, time determines x-value
fig = pyl.figure(figsize=(12, 8)) # [8, 14)) for poster] # 15,10))  8, 14)) [ 12, 8)) for figures] [15, 15)) for other figures]
ax = fig.add_subplot(111)
yCounter = 0

colorDict = {}
colorCounter = 0
for event in eventList:
	colorCounter += 1
	colorDict[event] = colorCounter/10.

conditionInterval = 150

for conditions in comboList: # cycle through each T+Sp combo
	# print conditions
	normalizedTiming_allEvents, normalizedTimingRev_allEvents = [], []
	numberOfAnimals_allEvents = []
	animalListTemp = comboDict[conditions].keys()
	animalListTemp.sort()
	yCounter += conditionInterval
	posCounter = 0
	normalizedTiming = []
	for animal in animalListTemp:
		normalizedTiming = []
		posCounter += 3
		eventListTemp = comboDict[conditions][animal].keys()
		eventListTemp.sort()
		for event in eventList:
			if event <> 'membrane_reaches_yolk':
				try:
					normalizedTiming = (comboDict[conditions][animal][event][1]*(float(comboDict[conditions][animal][event][2])-float(comboDict[conditions][animal]['membrane_reaches_yolk'][2])))
					yCoor = (int(posCounter + yCounter-3*len(animalListTemp)/2))
					pyl.plot(normalizedTiming, yCoor, marker = markerShape, linestyle = '', markersize = 7, c=mpl.cm.gist_ncar(colorDict[event]))# spectral(colorDict[event]))
				except:
					pass
indexVector = []
for x in range(len(comboList)+2):
	indexVector.append(x*conditionInterval)
if len(Tlist) > 1:
	labelList = ['']
	for combo in comboList:
		labelList.append((combo[:4]+ unichr(176)+'C'))
	labelList.append('')
	figTitle = SpList[0]
if len(SpList) > 1:
	labelList = ['']
	for combo in comboList:
		labelList.append(combo[4:])
	labelList.append('')
	figTitle = str(Tlist[0]) + unichr(176)+'C'
index = np.array(indexVector)
pyl.xticks(color = lineColor, fontsize = 20)
pyl.yticks(index, labelList, color = lineColor, fontsize = 25)
pyl.xlabel('time (min post-cellularization)', color = lineColor, fontsize = 25)
pyl.title(figTitle, color = lineColor, fontsize = 30)
ax.spines['bottom'].set_color(lineColor)
ax.spines['left'].set_color(lineColor)
ax.yaxis.label.set_color(lineColor)
ax.tick_params(axis='x', colors=lineColor)
ax.tick_params(axis='y', colors=lineColor)
pyl.savefig('dotLinePlots/dotLinePlot_' + str(sys.argv[1]) + '_' + figTitle + '_' + backgroundSetting +  '.pdf', transparent = True, facecolor = faceColor, linecolor = 'w') # to change from black to white background, make non-transparent
pyl.close()

# Normalized graph
textOutput = open('dotLinePlots/dotLineNumbers_' + str(sys.argv[1]) + '_' + figTitle + '.txt', 'w')

fig = pyl.figure(figsize=(14,6)) # [14, 6)) for poster]
ax = fig.add_subplot(111)
yCounter = 0
statsDict = {}
for conditions in comboList: # cycle through each T+Sp combo
	#if conditions not in meanDict:
	#	meanDict[conditions] = {}
	conditionDict = {}
	normalizedTiming_allEvents, normalizedTimingRev_allEvents = [], []
	numberOfAnimals_allEvents = []
	animalListTemp = comboDict[conditions].keys()
	animalListTemp.sort()
	yCounter += conditionInterval
	posCounter = 0
	normalizedTiming = []
	for animal in animalListTemp:
		normalizedTiming = []
		posCounter += 3
		eventListTemp = comboDict[conditions][animal].keys()
		eventListTemp.sort()
		for event in eventList:
			try:
				normalizedTiming = (comboDict[conditions][animal][event][1]*(float(comboDict[conditions][animal][event][2])-float(comboDict[conditions][animal]['membrane_reaches_yolk'][2])))
				yCoor = (int(posCounter + yCounter-3*len(animalListTemp)/2))
				scaledNormalizedTiming = (normalizedTiming/(comboDict[conditions][animal][event][1]*(float(comboDict[conditions][animal]['trachea_fill'][2])-float(comboDict[conditions][animal]['membrane_reaches_yolk'][2]))))
				pyl.plot(scaledNormalizedTiming, yCoor, marker = markerShape, linestyle = '', markersize = 7, c=mpl.cm.spectral(colorDict[event]))
				if event not in conditionDict:
					conditionDict[event] = [0,0, []]
				conditionDict[event][0] += scaledNormalizedTiming
				conditionDict[event][1] += 1
				conditionDict[event][2].append(scaledNormalizedTiming)
				# textOutput.write(str(conditions[:4]) + '\t' str(conditions[4:]) + '\t' + str(event) + '\t' + str(scaledNormalizedTiming) + '\n')
				if event not in statsDict: # write to a dictionary of vectors for statistical analysis (one key for each event)
					statsDict[event] = []
				# statsDict[event].append((scaledNormalizedTiming, float(conditions[:4]))) # vector needs two dimensions: time and temperature/species -> record conditions
				statsDict[event].append((float(conditions[:4]), scaledNormalizedTiming)) # vector needs two dimensions: time and temperature/species -> record conditions
			except:
				pass
	for event in eventList:
		try:
			meanValue = conditionDict[event][0]/conditionDict[event][1]
			varianceValue = []
			for scaledValue in conditionDict[event][2]:
				varianceValue.append((scaledValue-meanValue)**2)
			stDevValue = math.sqrt(sum(varianceValue)/(len(varianceValue)-1))
			textOutput.write(str(conditions[:4]) + '\t' + str(conditions[4:]) + '\t' + str(event) + '\t' + str(meanValue) + '\t' + str(stDevValue) + '\n') 
		except:
			pass
index = np.array(indexVector)
pyl.yticks(index, labelList, color = lineColor, fontsize = 20)
pyl.xticks(color = lineColor, fontsize = 20)
pyl.title(figTitle, color = lineColor, fontsize = 25)
pyl.xlim(-0.2,1.)
ax.spines['bottom'].set_color(lineColor)
ax.spines['left'].set_color(lineColor)
ax.yaxis.label.set_color(lineColor)
ax.tick_params(axis='x', colors=lineColor)
ax.tick_params(axis='y', colors=lineColor)
pyl.xlabel('proportion of development', color = lineColor, fontsize = 20)
pyl.savefig('dotLinePlots/dotLinePlot_' + str(sys.argv[1]) + '_' + figTitle + '_' + backgroundSetting +  '_Normalized.pdf', transparent = True, facecolor = faceColor, linecolor = 'w') # to change from black to white background, make non-transparent
pyl.close()

statsList = statsDict.keys()
statsList.sort()
exciseList = ['membrane_reaches_yolk', 'trachea_fill', 'germband_maxima', 'clypeolabrum_retracts', 'clypeolabrum_ventral_lobe_even']
for event in exciseList:
	try:
		statsList.pop(statsList.index(event))
	except:
		pass

if 'distribution' in sys.argv:
	# make a distribution of the occurrences, x axis is time, y axis is number of occurrences
	pyl.figure(figsize=(14,6))
	for event in statsList:
		xValues = []
		for (temperature, time) in statsDict[event]:
			xValues.append(time)
		n = len(xValues)
		pyl.hist(xValues, n, edgecolor = mpl.cm.spectral(colorDict[event]), color = mpl.cm.spectral(colorDict[event])) #, histtype = 'stepfilled')
	pyl.xlim(-.2, 1.)
	pyl.xlabel('proportion of development', color = lineColor, fontsize = 20)
	pyl.ylabel('number of occurrences', color = lineColor, fontsize = 20)
	pyl.savefig('curveFitting/eventDistribution.pdf',transparent = True)
	pyl.close()



# Statistical analysis
if 'stats' in sys.argv:

	def SyxSqrd(coordinates):
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
		vectorB = []
		for (i,j) in vectorA:
			runningSum += i**2
			vectorB.append(i)
		return (runningSum-(sum(vectorB)**2/n))/(n-1)

	numberOfComparisons = 0

	colorDict = {}
	colorCounter = 0
	for event in eventList:
		colorCounter += 1
		colorDict[event] = colorCounter/10.


	redundancyList = []
	for event1 in statsList:
		for event2 in statsList:
			if event1 <> event2:
				numberOfComparisons +=1
	for event1 in statsList:
		redundancyList.append(event1)
		for (timing, temperature) in statsDict[event1]:
			pyl.figure(0)
			pyl.plot(timing, temperature, marker = markerShape, linestyle = '', markersize = 3, c=mpl.cm.spectral(colorDict[event1]))
			pyl.figure(1)
			pyl.plot(temperature, timing, marker = markerShape, linestyle = '', markersize = 3, c=mpl.cm.spectral(colorDict[event1]))
		for event2 in statsList:
			if event1 <> event2 and event2 not in redundancyList:
				# find the linear regression
				slope1, intercept1 = findSlopeAndIntercept(statsDict[event1])
				slope2, intercept2 = findSlopeAndIntercept(statsDict[event2])
				n1 = len(statsDict[event1])
				n2 = len(statsDict[event2])
				# find if their slopes are significantly different (assume the intercepts will be)
				Syx1 = SyxSqrd(statsDict[event1])
				Syx2 = SyxSqrd(statsDict[event2])
				var_pooled = (Syx1*(n1-2)+Syx2*(n2-2))/(n1+n2-4)
				var1 = findVariance(statsDict[event1])
				var2 = findVariance(statsDict[event2])		
				tTest = (slope1-slope2)/(math.sqrt((var_pooled/((n1-1)*var1))+(var_pooled/((n2-1)*var2))))
				dof = n1+n2-4
				criticalValue = st.t.isf(0.001/2,dof)
				pValue = (numberOfComparisons)*(1.0 - (st.t.cdf(math.fabs(tTest),dof) - st.t.cdf(-math.fabs(tTest), dof)))
			
				print event1,'slope=',slope1,'int=',intercept1,event2,'slope=',slope2,'int=',intercept2,'t=',math.fabs(tTest),'dof=',dof,'tcritical0.001=',criticalValue,'P=', pValue

	for event in statsList:
		z = np.arange(17,35,2)
		slope, intercept = findSlopeAndIntercept(statsDict[event])
		pyl.figure(0)
		pyl.plot(z, z*slope+intercept, marker = '', linestyle = '--', lw = 1, c=mpl.cm.spectral(colorDict[event]))
		pyl.figure(1)
		pyl.plot(z*slope+intercept, z, marker = '', linestyle = '--', lw = 1, c=mpl.cm.spectral(colorDict[event]))
	
		n = len(statsDict[event])
		talpha = st.t.isf(0.05/2,(n-2))
		var = findVariance(statsDict[event])
		runningSum = 0
		stError = []
		for (x,y) in statsDict[event]:
			runningSum += x
		mean = runningSum/n
		upGap = []
		downGap = []
		for u in z:
			stError.append(talpha*math.sqrt(SyxSqrd(statsDict[event]))*math.sqrt((1/n)+((u-mean)**2)/((n-1)*var)))
			upGap.append((u*slope+intercept - talpha*math.sqrt(SyxSqrd(statsDict[event]))*math.sqrt((1/n)+((u-mean)**2)/((n-1)*var))))
			downGap.append((u*slope+intercept + talpha*math.sqrt(SyxSqrd(statsDict[event]))*math.sqrt((1/n)+((u-mean)**2)/((n-1)*var))))
		try:
			pyl.figure(0)
			pyl.fill_between(z, downGap, upGap, edgecolor = mpl.cm.spectral(colorDict[event]), facecolor = mpl.cm.spectral(colorDict[event]))
		except:
			print "error in figure 1"
		# try:
#		pyl.figure(1)
#		pyl.plot(downGap,z, marker = '', linestyle = '-', lw = 1, c = mpl.cm.spectral(colorDict[event]))
#		pyl.plot(upGap,z, marker = '', linestyle = '-', lw = 1, c = mpl.cm.spectral(colorDict[event]))
#		pyl.fill_between(downGap, upGap, z, edgecolor = 'g', facecolor = 'g') # pattern 2 # mpl.cm.spectral(colorDict[event]), facecolor = mpl.cm.spectral(colorDict[event]))
#		pyl.fill_between(z, upGap, downGap, edgecolor = 'g', facecolor = 'g') # nothing
#		pyl.fill_between(z, downGap, upGap, edgecolor = 'b', facecolor = 'b') # nothing
#		pyl.fill_between(upGap, downGap, z, edgecolor = 'r', facecolor = 'r') # pattern 1
#		pyl.fill_between(upGap, z, downGap, edgecolor = 'y', facecolor = 'y') # pattern 1
#		pyl.fill_between(downGap, z, upGap, edgecolor = 'cyan', facecolor = 'cyan') # pattern 2 # where=upGap>=z
		#except:
		#	print "error in figure 2"

	pyl.figure(0)
	pyl.plot(z,z*0., marker = '', linestyle = '--', lw = 1, c = 'k')
	pyl.figure(1)
	pyl.plot(z*0.,z, marker = '', linestyle = '--', lw = 1, c = 'k')

	# establish confidence intervals for the regression lines


	spNames = ''
	Tnames = ''
	for species in SpList:
		spNames += species+'_'
	for temperature in Tlist:
		Tnames +=temperature+'_'
	pyl.figure(0)
	pyl.ylim(-0.2,1.)
	pyl.figure(1)
	pyl.ylim(17,33)
	pyl.xlim(-0.2,1.)
	pyl.figure(0)
	pyl.savefig('curveFitting/curveFitting_'+spNames+Tnames+'eventTrends_inv.pdf', transparent = True)
	pyl.close()
	pyl.figure(1)
	pyl.savefig('curveFitting/curveFitting_'+spNames+Tnames+'eventTrends.pdf', transparent = True)
	pyl.close()
		
		
		