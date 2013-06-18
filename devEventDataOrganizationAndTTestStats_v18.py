# example command: python devEventDataOrganizationAndTTestStats_v12.py trial8 -s Pe Me -T 27.5 30
# must be run from Dropbox folder


import os, sys, re, string, math, array
import scipy.stats as st
import pylab as pyl
import matplotlib as mpl
import numpy as np

### Part 0 - First read in list of positions to be analyzed
fileNamesData = []
fileNamesFiles = os.listdir("./ImageProcessing/")
for filename in fileNamesFiles:
	if '_EventTimesTrial_avg_' in filename:
		line = re.split("_", filename)
		if len(sys.argv) > 2:
			if '-s' and '-t' in sys.argv: # find which species and temperatures are listed
				speciesAndTempList = []
				for i in sys.argv:
					if i in ['An','Er','Me','MeAcc','Me18','Mj','Pe','Ps','Se','Si','Vi','Wi','Ya','NY20','NY62','CY01A','CY21B','12414','12592','12622','16532','12.5','15','17.5','17','20','22.5','25','27.5','30','32.5']:
						speciesAndTempList.append(i)
				if line[6] in speciesAndTempList:
					if line[7] in speciesAndTempList:
						fileNamesData.append((filename, line[1], line[2], line[3], line[6], line[7], line[8][:1]))
			elif '-s' in sys.argv: # find which species are listed
				speciesList = []
				for i in sys.argv:
					if i in ['An','Er','Me','MeAcc','Me18','Mj','Pe','Ps','Se','Si','Vi','Wi','Ya','NY20','NY62','CY01A','CY21B','12414','12592','12622','16532']:
						speciesList.append(i)
				if line[6] in speciesList:
					fileNamesData.append((filename, line[1], line[2], line[3], line[6], line[7], line[8][:1]))
			elif '-T' in sys.argv: # find which temperatures are listed
				tempList = []
				for j in sys.argv:
					if j in ['12.5','15','17.5','17','20','22.5','25','27.5','30','32.5']:
						tempList.append(j)
				try:
					if line[7] in tempList:
						fileNamesData.append((filename, line[1], line[2], line[3], line[6], line[7], line[8][:1]))
				except: pass
		elif len(sys.argv) == 2:
			try:
				fileNamesData.append((filename, line[1], line[2], line[3], line[6], line[7], line[8][:1]))
			except: pass
		else:
			print 'error'
			sys.exit()
			
comboDict = {}; meanDict = {}; Tdict = {}; SpDict = {}; focusWeight = {}; totalData = []
speciesNames = [['An', 'Dana'],['Er','Dere'],['Me', 'Dmel'],['Me18','Dm18'],['Mj', 'Dmoj'],['Ps','Dpse'],['Pe','Dper'],['Si','Dsim'],['Se','Dsec'],['Wi','Dwil'],['Vi','Dvir'],['Ya','Dyak'],['CY01A','DyakCY01A'],['CY21B','DyakCY21B'],['NY20','DyakNY20'],['NY62','DyakNY62']]
CYyakStrains = ['CY01A','CY04B','CY08A','CY15B4','CY17C','CY22B','CY23A','CY21B','CY27A9']
NYyakStrains = ['NY20','NY48','NY55','NY56','NY60','NY62','NY65','NY73','NY81','NY85']
taiyakStrains = ['tai18E2']
yakStrains = CYyakStrains + NYyakStrains + taiyakStrains
NSsimStrains = ['NS37','NS67','NS79','NS113','NS137']
MDsimstrains = ['MD15','MD63','MD73','MD105','MD106','MD199','MD221','MD233','MD251']
simStrains = ['Si'] + NSsimStrains + MDsimstrains


# Now read in the data from all those positions
for (fileName, date, pos, z, speciesShort, T, dilation) in fileNamesData:
	for (shortName, longName) in speciesNames:
		if speciesShort == shortName:
			species = longName
			if 'Ya' in sys.argv:
				if speciesShort in yakStrains:
					species = 'Dyak'
			if 'Si' in sys.argv:
				if speciesShort in simStrains:
					species = 'Dsim'
	entryData = []
	entryData2 = []
	matlabFile = open('ImageProcessing/'+fileName,'r') # these are the matlab generated event time files
	for i in matlabFile:
		line = re.split("[\s]",i)
		for value in line:
			if value <> '': entryData.append(value)
	totalData.append((date + pos, z, species, T, dilation, entryData))
	if T + species not in comboDict: comboDict[T+species] = []
	if T + species not in meanDict: meanDict[T+species] = []
	if T not in Tdict: Tdict[T] = []	
 	if species not in SpDict: SpDict[species] = []
	if date + pos + z not in focusWeight: focusWeight[date+pos+z] = []

comboList = comboDict.keys()
Tlist = Tdict.keys()
SpList = SpDict.keys()
animalList = focusWeight.keys()
Tlist.sort()
SpList.sort()
comboList.sort()

numberOfTemperatures = len(Tlist) # number of Temperatures
numberOfSpecies = len(SpList) # number of Species

for (Pos, z, sp, T, dilation, data) in totalData:
	comboDict[T+sp].append((Pos, z, int(dilation), data))
	
for zVideo in animalList:
	for filename in fileNamesFiles:
		if zVideo[:8] in filename: # match the date
			if zVideo[8:-1] in filename: # match the position
				line = re.split("_", filename)
				if zVideo[-1:] == line[3]: # match the z-stack
					try:
						focusfile = open('FocusFiles/FocusScoresUnaligned_'+zVideo[:8]+'_'+zVideo[8:-1]+'_'+zVideo[-1:]+'.txt', 'r')
						for i in focusfile:
							focusScores = re.split("[\s]", i)
							focusWeight[zVideo].append(focusScores[1]) # 0 is normedVariance, 1 is autocorrelation, 2 is stdevcorrelation
					except:
						print 'focus file missing for', zVideo[:8], zVideo[8:-1], zVideo[-1:]

outputFile = open('Rinput_event_'+str(sys.argv[1])+'.txt','w')	# want to label the header as (event, t, T, Sp, sd, se, ci)
outputFile.write('condition' + '\t' + 'event' +'\t' + 'mean' +'\t' + 'meanRev' +'\t' + 'median' +'\t' + 'stDev' +'\t' + 'stDevRev' +'\t' + 'MAD' +'\t' + 'errMin' +'\t' + 'errMax' +'\t' + 'numberOfAnimals' +'\n')
outputTable = open('RinputTable_event_'+str(sys.argv[1])+'.txt','w')
# outputVariance = open('EventVariance_event_'+str(sys.argv[1])+'.txt','w')

### Part 1 - Generate statistical summaries (mean and standard deviation) of each group
for conditions in comboList: # cycle through each T+Sp combo
	print conditions
	tempData = []
	positionDict = {}
	positionList = []
	for (pos, z, dilation, data) in comboDict[conditions]: # can try first going through and figuring out which are unique. If not unique, then can cycle through and average in the components
		if pos not in positionDict: 
			positionDict[pos] = []
		positionDict[pos].append((dilation, data)) # first generate list of positions
	positionList = positionDict.keys()
	numberOfEvents = len(positionDict[pos][0][1]) # this length should be equivalent to 35
	for position in positionList:
		averageData = []
		focusAdjustmentRecord = 0
		for event in range(numberOfEvents): # go through each time event in the data
			focusAdjustmentTracker = 0
			sumPos = 0
			sumWeight = 0
			z= 0
			try:
				correctTime = []
				for (dilation, zPos) in positionDict[position]: # go through each z-stack
					# for zPos[event], which is a time, need to go through all the focus score of all the focus stacks
					weightRankPrep = []
					eventTime = int(float(zPos[event])) - 1 # event times are called in matlab, which is not 0-counting
					for i in range(len(positionDict[position])):
						weightRankPrep.append(float(focusWeight[position+str(i)][eventTime]))
					# now need to find the maximum within weightRankPrep
					isitfirst = weightRankPrep.index(max(weightRankPrep))
					if z == isitfirst:
						correctTime.append(eventTime)
					z += 1
				# if correctTime == []:
					# print 'error, no points in focus'
				for timeOfEvent in correctTime:
					sumPos += timeOfEvent
				averageData.append(sumPos/len(correctTime)) # record the selected position
			except:
				sumPos = 0
				for (dilation, zPos) in positionDict[position]: # go through each z-stack
					sumPos += float(zPos[event])
				averageData.append(sumPos/len(positionDict[position])) # record the average values for each time event for an animal
				focusAdjustmentTracker = 2
			if focusAdjustmentTracker == 0:
				pass
				# print 'focus adjusted for', position, event
			elif focusAdjustmentTracker == 2:
				focusAdjustmentRecord +=1
				# print 'focus not adjusted for', position, event
		tempData.append((dilation, averageData))
		# print averageData
		if focusAdjustmentRecord >= 1:
			print 'focus not adjusted for', position, 'in', focusAdjustmentRecord, 'events'
	numberOfAnimals = len(positionList)
#	median, MAD, meanRev, adjustmentFactor, stDevRev = [0.]*numberOfEvents, [0.]*numberOfEvents, [0.]*numberOfEvents, [0]*numberOfEvents, [0.]*numberOfEvents
#	tempSum, tempVar, mean, stDev, ymin, ymax = [0.0]*numberOfEvents, [0.0]*numberOfEvents, [0.0]*numberOfEvents, [0.0]*numberOfEvents, [0.0]*numberOfEvents, [0.0]*numberOfEvents
	mean_allEvents, median_allEvents, meanRev_allEvents, MAD_allEvents = [], [], [], []
	stDev_allEvents, stDevRev_allEvents, ymin_allEvents, ymax_allEvents = [],[], [],[]
	normalizedTiming_allEvents, normalizedTimingRev_allEvents = [], []
	adjustmentFactor = [0]*numberOfEvents
#	normalizedTiming = [0.0]*numberOfEvents
	for event in range(numberOfEvents): # go through each of the 35 events to determine the mean and median
		normalizedTiming = []
		for (dilation, position) in tempData: # go through each of the animals (positions)
			normalizedTiming.append(dilation*(float(position[event])-float(position[6])))
		mean = sum(normalizedTiming)/len(tempData) # generate a mean across positions for a single event
		median = np.ma.median(normalizedTiming)
		outputTable.write(str(conditions) + '\t' + str(event+1) + '\t' + str(normalizedTiming) + '\t' + str(mean) + '\n') 

		# determine the MAD
		tempMAD = []
		for incident in normalizedTiming:
			tempMAD.append(math.fabs(incident-median))
		if len(normalizedTiming) > 1: 
			MAD = np.ma.median(tempMAD)
		else:
			MAD = 1000.	
			
		# determine the revised mean	
		tempSumRev = 0
		normalizedTimingRev = []
		for incident in normalizedTiming:
			zScore = math.fabs((incident-mean)/(0.0001+1.4826*MAD))
			if zScore < 3.5:
				normalizedTimingRev.append(incident)
			elif math.fabs(mean-incident) > math.fabs(median-incident): # if the 'outlier' is closer to the median, it is not a real outlier
				normalizedTimingRev.append(incident)
			else:
				adjustmentFactor[event] += 1
		meanRev = sum(normalizedTimingRev)/len(normalizedTimingRev)

		# determine the standard deviation and revised standard deviation
		var = []
		varRev = []
		for incident in normalizedTiming:
			var.append((incident-mean)**2)
		for incident in normalizedTimingRev:
			zScore = math.fabs((incident-meanRev)/(0.0001+1.4826*MAD))
			if zScore < 3.5:
				varRev.append((incident-meanRev)**2) # changed this to revised time list
			elif math.fabs(mean-incident) > math.fabs(median-incident):
				varRev.append((incident-meanRev)**2)
		if len(normalizedTiming) > 1:
			stDev = math.sqrt(sum(var)/(len(normalizedTiming)-1))
			stDevRev = math.sqrt(sum(varRev)/(len(normalizedTimingRev)))
			ymin = mean - stDev
			ymax = mean + stDev
		else:
			stDev = 1000.
			stDevRev = 1000
		incidentOutput = normalizedTimingRev
		
		if max(normalizedTiming) - min(normalizedTiming) < 5:
			meanRev = mean
			stDevRev = stDev
			if max(normalizedTiming) - min(normalizedTiming) == 0:
				if len(normalizedTiming) > 2:
					stDevRev = 2.0
				
		# redo if failed					
		if event <> 6:
			MAD2 = MAD
			mean2 = mean			
			while stDevRev == 0.:
				# first remove the most problematic entry
				print event, normalizedTiming, mean2, MAD2, median
				if MAD2 < mean2:
					normalizedTiming.pop(normalizedTiming.index(max(normalizedTiming)))
				elif MAD2 > mean2:
					normalizedTiming.pop(normalizedTiming.index(min(normalizedTiming)))
				# now recalculate the variance (sum the tempVarComponents)
				if len(normalizedTiming) <= 1:
					print 'too few instances'
					stDevRev = 1000.
					break
				mean2 = sum(normalizedTiming)/len(normalizedTiming)
				median2 = np.ma.median(normalizedTiming)
				tempMAD2 = []
				var2 = []
				for incident in normalizedTiming:
					var2.append((incident-mean2)**2)
					tempMAD2.append(math.fabs(incident-median2))
				MAD2 = np.ma.median(tempMAD2)
				# remove all remaining outliers again
				normalizedTimingRev2 = []
				for incident in normalizedTiming:	
					zScore = math.fabs((incident-mean2)/(0.0001+1.4826*MAD))
					if zScore < 3.5:
						normalizedTimingRev2.append(incident)
					elif math.fabs(mean2-incident) > math.fabs(median2-incident):
						normalizedTimingRev2.append(incident)
				if len(normalizedTimingRev2) <= 1:
					print 'too few instances'
					stDevRev = 1000.
					break
				meanRev = sum(normalizedTimingRev2)/len(normalizedTimingRev2)
				varRev = []	
				for incident in normalizedTimingRev2:	
					zScore = math.fabs((incident-meanRev)/(0.0001+1.4826*MAD))
					if zScore < 3.5:
						varRev.append((incident-meanRev)**2)
					elif math.fabs(meanRev-incident) > math.fabs(median2-incident):
						varRev.append(incident)
				if len(normalizedTimingRev2) > 1:
					stDevRev = math.sqrt(sum(varRev)/(len(normalizedTimingRev2)))
				incidentOutput = normalizedTimingRev2
				if max(normalizedTiming) - min(normalizedTiming) == 0:
					if len(normalizedTiming) > 2:
						stDevRev = 2.0
				if stDevRev <> 0:
					print 'success'
		outputFile.write(str(conditions) + '\t' + str(event+1) + '\t' + str(mean) + '\t' + str(meanRev) + '\t' + str(median) + '\t' + str(stDev) + '\t' + str(stDevRev) + '\t' + str(MAD) + '\t' + str(ymin) + '\t' + str(ymax) + '\t' + str(numberOfAnimals) + '\t' + str(normalizedTiming) + '\t' + str(incidentOutput)+ '\n')
		mean_allEvents.append(mean)
		median_allEvents.append(median)
		meanRev_allEvents.append(meanRev)
		MAD_allEvents.append(MAD)
		stDev_allEvents.append(stDev)
		stDevRev_allEvents.append(stDevRev)
		ymin_allEvents.append(ymin)
		ymax_allEvents.append(ymax)
		normalizedTiming_allEvents.append(normalizedTiming)
		normalizedTimingRev_allEvents.append(incidentOutput)
		
		
	if adjustmentFactor[33] > 0:
		print conditions, adjustmentFactor[33]
#	print conditions, 'tempVar', tempVar
#	print conditions, 'stDev', stDev
#	print conditions, 'stDevRev', stDevRev  
	meanDict[conditions] = (mean_allEvents, stDev_allEvents, numberOfAnimals, normalizedTimingRev_allEvents, median_allEvents, MAD_allEvents, meanRev_allEvents, stDevRev_allEvents) # save the mean and stdev for each event for use later
	outputFile.write('\n')

### Part 2 - t-Test statistics
tTestOutput = open('tTestOutputs_event_'+str(sys.argv[1])+'.txt','w')
tTestOutput.write('event'+'\t'+'speciesA'+'\t'+'speciesB'+'\t'+'temperatureA'+'\t'+'temperatureB'+'\t'+'t-score'+'\t'+'DoF'+'\t'+'p-value' +'\t'+'cutoff'+'\n')
tTestOutput_s = open('tTestOutputs_event_'+str(sys.argv[1])+'_crossSpecies.txt','w')
tTestOutput_s.write('event'+'\t'+'speciesA'+'\t'+'speciesB'+'\t'+'temperatureA'+'\t'+'temperatureB'+'\t'+'t-score'+'\t'+'DoF'+'\t'+'p-value' +'\t'+'cutoff'+'\n')
tTestOutput_T = open('tTestOutputs_event_'+str(sys.argv[1])+'_crossTemperature.txt','w')
tTestOutput_T.write('event'+'\t'+'speciesA'+'\t'+'speciesB'+'\t'+'temperatureA'+'\t'+'temperatureB'+'\t'+'t-score'+'\t'+'DoF'+'\t'+'p-value' +'\t'+'cutoff'+'\n')

confidenceLevelCutoff = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05]

# first calculate the actual number of comparisons in preparation for the Bonferroni correction
numberOfComparisons = 0
comparisonDict = {}
for event in range(numberOfEvents):	
	for species in SpList:
		for TA in Tlist:
			for TB in Tlist:
				if TA <> TB:
					try:
						(meanA, stDevA, nA, normTimingA, medianA, madA, meanRA, stDevRA) = meanDict[TA+species]
						(meanB, stDevB, nB, normTimingB, medianB, madB, meanRB, stDevRB) = meanDict[TB+species]
						if species + TA + TB + str(event) not in comparisonDict: 
							if species + TB + TA + str(event) not in comparisonDict:
								comparisonDict[species+TA+TB+str(event)] = []
								numberOfComparisons += 1
					except:
						pass
	for T in Tlist:
		for SpA in SpList:
			for SpB in SpList:
				if SpA <> SpB:
					try:
						(meanA, stDevA, nA, normTimingA, medianA, madA, meanRA, stDevRA) = meanDict[T+SpA]
						(meanB, stDevB, nB, normTimingB, medianB, madB, meanRB, stDevRB) = meanDict[T+SpB]
						if SpA + SpB + T+str(event) not in comparisonDict: 
							if SpB + SpA + T+str(event) not in comparisonDict: 
								comparisonDict[SpA + SpB + T+str(event)] = []
								numberOfComparisons += 1
					except:
						pass

def findpValue(mean1, mean2, stDev1, stDev2, n1, n2, constant, eventStage):
	# if ((constant*stDev1[eventStage])**2/n1)+((constant*stDev2[eventStage])**2/n2) == 0:
	#	print 'flag', mean1[eventStage], mean2[eventStage], stDev1[eventStage], stDev2[eventStage], n1, n2, constant, eventStage, stDev1, stDev2
	tTest = (mean1[eventStage]-mean2[eventStage])/math.sqrt(((constant*stDev1[eventStage])**2/n1)+((constant*stDev2[eventStage])**2/n2))
	degreesOfFreedom = ((constant*stDev1[eventStage])**2/n1 + (constant*stDev2[eventStage])**2/n2)**2/((((constant*stDev1[eventStage])**2/n1)**2/(n1-1))+(((constant*stDev2[eventStage])**2/n2)**2/(n2-1)))
	pValue = (numberOfComparisons)*(1.0 - (st.t.cdf(math.fabs(tTest),degreesOfFreedom) - st.t.cdf(-math.fabs(tTest), degreesOfFreedom))) # 2-tailed p-value
	result = []
	if pValue < 0.05:
		for cutoff in confidenceLevelCutoff:
			if pValue <= cutoff:
				result = [tTest, degreesOfFreedom, pValue, cutoff]
				break
	elif pValue > 0.05:
		result = []
	return result

comparisonDict = {}
# Need to cycle through all events	
for event in range(numberOfEvents):	
	# First cycle through the species to compare all temperatures to each other 
	for species in SpList:
		for TA in Tlist:
			for TB in Tlist:
				if TA <> TB:
					try:
						(meanA, stDevA, nA, normTimingA, medianA, madA, meanRA, stDevRA) = meanDict[TA+species]
						(meanB, stDevB, nB, normTimingB, medianB, madB, meanRB, stDevRB) = meanDict[TB+species]
						if str(event) + species + TA + TB not in comparisonDict:
							if str(event) + species + TB + TA not in comparisonDict:
								# if event <> 6:
									# if stDevRA[event] == 0:
									#	print 'stDevRA', species, TA, nA, event, normTimingA[event]
									# if stDevA[event] == 0:
									#	print 'stDevA', species, TA, nA, event, normTimingA[event]
									# if madA[event] == 0:
									#	print 'madA', species, TA, nA, event, normTimingA[event]
									# if madB[event] ==0:
									#	print 'madB', species, TB, nB, event, normTimingB[event]
									# if stDevB[event] == 0:
									#	print 'stDevB', species, TB, nB, event, normTimingB[event]
									# if stDevRB[event] == 0:
									#	print 'stDevRB', species, TB, nB, event, normTimingB[event]
								comparisonDict[str(event)+species+TA+TB] = []
								meanPvalues = findpValue(meanA, meanB, stDevA, stDevB, nA, nB, 1, event)
								if meanPvalues <> []:
									tTestOutput.write(str(event+1)+'\t'+str(species) + '\t' + str(species) + '\t' + str(TA) + '\t' + str(TB) + '\t' + str(meanPvalues[0]) + '\t'  + str(meanPvalues[1]) + '\t' + str(meanPvalues[2])  + '\t' + str(meanPvalues[3]) + '\t' + 'mean' + '\n')
								meanRevPvalues = findpValue(meanRA, meanRB, stDevRA, stDevRB, nA, nB, 1, event)
								if meanRevPvalues <> []:
									tTestOutput.write(str(event+1)+'\t'+str(species) + '\t' + str(species) + '\t' + str(TA) + '\t' + str(TB) + '\t' + str(meanRevPvalues[0]) + '\t'  + str(meanRevPvalues[1]) + '\t' + str(meanRevPvalues[2])  + '\t' + str(meanRevPvalues[3]) + '\t' + 'meanRev' + '\n')
									tTestOutput_T.write(str(event+1)+'\t'+str(species) + '\t' + str(species) + '\t' + str(TA) + '\t' + str(TB) + '\t' + str(meanRevPvalues[0]) + '\t'  + str(meanRevPvalues[1]) + '\t' + str(meanRevPvalues[2])  + '\t' + str(meanRevPvalues[3]) + '\n')
								medianPvalues = findpValue(medianA, medianB, madA, madB, nA, nB, 1.4826, event)
								if medianPvalues <> []:
									tTestOutput.write(str(event+1)+'\t'+str(species) + '\t' + str(species) + '\t' + str(TA) + '\t' + str(TB) + '\t' + str(medianPvalues[0]) + '\t'  + str(medianPvalues[1]) + '\t' + str(medianPvalues[2])  + '\t' + str(medianPvalues[3]) + '\t' + 'median' + '\n')
					except:
						pass
						
	# Second cycle through the temperatures to compare all species to each other
	for T in Tlist:
		for SpA in SpList:
			for SpB in SpList:
				if SpA <> SpB:
					try:
						(meanA, stDevA, nA, normTimingA, medianA, madA, meanRA, stDevRA) = meanDict[T+SpA]
						(meanB, stDevB, nB, normTimingB, medianB, madB, meanRB, stDevRB) = meanDict[T+SpB]
						if str(event) + SpA + SpB + T not in comparisonDict:
							if str(event) + SpB + SpA + T not in comparisonDict:
								# if event <>6:
								#	if stDevRA[event] == 0:
								#		print 'stDevRA', T, SpA, nA, event, normTimingA[event]
								#	if stDevA[event] == 0:
								#		print 'stDevA', T, SpA, nA, event, normTimingA[event]
								#	if stDevB[event] == 0:
								#		print 'stDevB', T, SpB, nB, event, normTimingB[event]
									# if madA[event] == 0:
									#	print 'madA', species, TA, nA, event, normTimingA[event]
									# if madB[event] ==0:
									#	print 'madB', species, TB, nB, event, normTimingB[event]
								#	if stDevRB[event] == 0:
								#		print 'stDevRB', T, SpB, nB, event, normTimingB[event]
								comparisonDict[str(event)+SpA+SpB+T] = []
								meanPvalues = findpValue(meanA, meanB, stDevA, stDevB, nA, nB, 1, event)
								if meanPvalues <> []:
									tTestOutput.write(str(event+1)+'\t'+str(SpA) + '\t' + str(SpB) + '\t' + str(T) + '\t' + str(T) + '\t' + str(meanPvalues[0]) + '\t'  + str(meanPvalues[1]) + '\t' + str(meanPvalues[2])  + '\t' + str(meanPvalues[3]) + '\t' + 'mean' + '\n')
								meanRevPvalues = findpValue(meanRA, meanRB, stDevRA, stDevRB, nA, nB, 1, event)
								if meanRevPvalues <> []:
									tTestOutput.write(str(event+1)+'\t'+str(SpA) + '\t' + str(SpB) + '\t' + str(T) + '\t' + str(T) + '\t' + str(meanRevPvalues[0]) + '\t'  + str(meanRevPvalues[1]) + '\t' + str(meanRevPvalues[2])  + '\t' + str(meanRevPvalues[3]) + '\t' + 'meanRev' + '\n')
									tTestOutput_s.write(str(event+1)+'\t'+str(SpA) + '\t' + str(SpB) + '\t' + str(T) + '\t' + str(T) + '\t' + str(meanRevPvalues[0]) + '\t'  + str(meanRevPvalues[1]) + '\t' + str(meanRevPvalues[2])  + '\t' + str(meanRevPvalues[3]) + '\n')
								medianPvalues = findpValue(medianA, medianB, madA, madB, nA, nB, 1.4826, event)
								if medianPvalues <> []:
									tTestOutput.write(str(event+1)+'\t'+str(SpA) + '\t' + str(SpB) + '\t' + str(T) + '\t' + str(T) + '\t' + str(medianPvalues[0]) + '\t'  + str(medianPvalues[1]) + '\t' + str(medianPvalues[2])  + '\t' + str(medianPvalues[3]) + '\t' + 'median' + '\n')
					except:
						pass
						
					
### Part 3: Box and Whisker graphical analysis 
# for event in range(numberOfEvents):
for event in range(1):
	event = 33
	if len(Tlist) > 1:
		for species in SpList:
			plottingList = []
			labelList = []
			meanList = []
			for T in Tlist:
				try:
					(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
					meanList.append(meanR[event])
					plottingList.append(normTiming[event])
					labelList.append(T+unichr(176)+'C')
				except:
					pass
			fig = pyl.figure()
			bp = pyl.boxplot(plottingList, patch_artist = True)
			ax = fig.add_subplot(111)
			pyl.setp(bp['boxes'], color = 'w', linewidth = 3) # to switch from black background to white background, simply change the 'w' to 'k'
			pyl.setp(bp['whiskers'], color = 'w', linewidth = 3)
			pyl.setp(bp['caps'], color = 'w', linewidth = 3)
			pyl.setp(bp['medians'], color = 'b', linewidth = 3)
			pyl.setp(bp['fliers'], color = 'w', marker = 'x', markersize = 10, linewidth = 5)
			pyl.xticks(range(1,len(labelList)+1), labelList, color = 'w', fontsize = 15)
			pyl.suptitle(species, color = 'w', fontsize = 20)
			pyl.xlabel('Temperature', color = 'w', fontsize = 15)
			pyl.plot(range(1,len(labelList)+1), meanList, linestyle = '--', marker='D', color = 'r', markeredgecolor='r')
			ax.spines['bottom'].set_color('w')
			ax.spines['left'].set_color('w')
			ax.yaxis.label.set_color('w')
			ax.tick_params(axis='x', colors='w')
			ax.tick_params(axis='y', colors='w')
			pyl.ylabel('time (minutes post cellularization)', color = 'w', fontsize = 15)
			pyl.savefig('boxPlots/'+str(event+1) + '_' + str(species) + 'boxAndWhiskerPlot_' + str(sys.argv[1]) + '.pdf', transparent = True, facecolor = 'k', linecolor = 'w') # to change from black to white background, make non-transparent
			pyl.clf()
			pyl.cla()
			pyl.close(fig)
			pyl.close()
		
	if len(SpList) > 1:	
		for T in Tlist:
			plottingList = []
			labelList = []
			meanList = []
			for species in SpList:
				try:
					(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
					meanList.append(meanR[event])
					plottingList.append(normTiming[event])
					labelList.append(species)
				except:
					pass
			fig = pyl.figure()
			bp = pyl.boxplot(plottingList, patch_artist = True)
			ax = fig.add_subplot(111)
			pyl.setp(bp['boxes'], color = 'w', linewidth = 3)
			pyl.setp(bp['whiskers'], color = 'w', linewidth = 3)
			pyl.setp(bp['caps'], color = 'w', linewidth = 3)
			pyl.setp(bp['medians'], color = 'b', linewidth = 3)
			pyl.setp(bp['fliers'], color = 'w', marker = 'x', markersize = 10, linewidth = 5)
			pyl.xticks(range(1,len(labelList)+1), labelList, color = 'w', fontsize = 15)
			pyl.suptitle(T+unichr(176)+'C', color = 'w', fontsize = 20)
			pyl.xlabel('Species', color = 'w', fontsize = 15)
			pyl.plot(range(1,len(labelList)+1), meanList, linestyle = 'none', marker='D', color = 'r', markeredgecolor='r')
			ax.spines['bottom'].set_color('w')
			ax.spines['left'].set_color('w')
			ax.yaxis.label.set_color('w')
			ax.tick_params(axis='x', colors='w')
			ax.tick_params(axis='y', colors='w')
			pyl.ylabel('time (minutes post cellularization)', color = 'w', fontsize = 15)
			pyl.savefig('boxPlots/'+str(event+1) + '_' + str(T) + 'boxAndWhiskerPlot_' + str(sys.argv[1]) + '.pdf', transparent = True, facecolor = 'k', linecolor = 'w')
			pyl.clf()
			pyl.cla()
			pyl.close(fig)
			pyl.close()
			
### Part 4 - Stacked Bar Graphs

# chronologicalEventsList = [0,1,2,3,4,5,6,7,8,10,9,12,11,13,14,15,20,17,19,16,18,21,22,26,23,24,30,25,27,31,28,29,32,33,34]
# chronologicalEventsList = [0,1,2,3,4,5,6,7,8,10,9,12,11,13,14,15,17,20,19,16,18,21,22,23,26,24,30,25,27,31,28,29,32,33,34]
chronologicalEventsList = [1,33,30,21,11,9,6] # ,9,11,21,30,33]


# first graph across temperatures for each species
if len(Tlist) > 1:
	try:
		for species in SpList:
			labelList = []
			meanList = []
			stDevList = []
			eventDict = {}
			for event in chronologicalEventsList: #range(numberOfEvents):
				if event not in eventDict:
					eventDict[event] = []
			for event in chronologicalEventsList: # range(numberOfEvents):
				meanList = []
				stDevList = []
				labelList = []
				for T in Tlist:
					try:
						(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
						meanList.append(meanR[event])
						stDevList.append(stDevR[event])
						labelList.append(T+unichr(176)+'C')
					except:
						pass
				eventDict[event].append((meanList,stDevList))
			fig = pyl.figure()
			ax = fig.add_subplot(111)
			pyl.suptitle(species, color = 'w', fontsize = 20)
			pyl.xlabel('Temperature', color = 'w', fontsize = 15)
			index = np.arange(len(labelList))
			for e in chronologicalEventsList: #range(numberOfEvents):
				# want to cycle through colors so that it hits r, g, b and c, m, y : r-y-g-m-b-c
				# can use math.fabs(math.cos(event/pi)), math.fabs(math.sin(event/pi)), math.fabs(math.cos(event/pi))
				if e < 6:
					event = e
					pyl.bar(index+1, eventDict[event][0][0], 0.75, color = (math.fabs(math.cos(5*event*3.14/17)), math.fabs(math.sin(2*event*3.14/17)), math.fabs(math.sin(3.4*event*3.14/17))), yerr = eventDict[event][0][1])		
				elif e >= 6:
					# print event, eventDict[event][0][0], eventDict[event][0][1]
					event = e
					pyl.bar(index+1, eventDict[event][0][0], 0.75, color = (math.fabs(math.cos(5*event*3.14/17)), math.fabs(math.sin(2*event*3.14/17)), 0), yerr = eventDict[event][0][1])
			ax.spines['bottom'].set_color('w')
			ax.spines['left'].set_color('w')
			pyl.xticks(index+1+0.75/2.0, labelList, color = 'w', fontsize = 15)
			ax.yaxis.label.set_color('w')
			ax.tick_params(axis='x', colors='w')
			ax.tick_params(axis='y', colors='w')
			pyl.ylabel('time (minutes post cellularization)', color = 'w', fontsize = 15)
			pyl.savefig('barPlots/'+str(event+1) + '_' + str(species) + 'StackedBarGraphPlot_' + str(sys.argv[1]) + '.pdf', transparent = True, facecolor = 'k', linecolor = 'w') # to change from black to white background, make non-transparent
			pyl.close()
	except:
		print species

# now graph across species for each temperature	
if len(SpList) > 1:
	try:
		for T in Tlist:
			labelList = []
			meanList = []
			stDevList = []
			eventDict = {}
			for event in chronologicalEventsList: #range(numberOfEvents):
				if event not in eventDict:
					eventDict[event] = []
			for event in chronologicalEventsList: #range(numberOfEvents):
				meanList = []
				stDevList = []
				labelList = []
				for species in SpList:
					try:
						(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
						meanList.append(meanR[event])
						stDevList.append(stDevR[event])
						labelList.append(species)
					except:
						pass
				eventDict[event].append((meanList,stDevList))
			fig = pyl.figure()
			ax = fig.add_subplot(111)
			pyl.suptitle(T+unichr(176)+'C', color = 'w', fontsize = 20)
			pyl.xlabel('Species', color = 'w', fontsize = 15)
			index = np.arange(len(labelList))
			for e in chronologicalEventsList: #range(numberOfEvents):
				# want to cycle through colors so that it hits r, g, b and c, m, y : r-y-g-m-b-c
				# can use math.fabs(math.cos(event/pi)), math.fabs(math.sin(event/pi)), math.fabs(math.cos(event/pi))
				if e < 6:
					event = e
					pyl.bar(index+1, eventDict[event][0][0], 0.75, color = (math.fabs(math.cos(5*event*3.14/17)), math.fabs(math.sin(2*event*3.14/17)), math.fabs(math.sin(3.4*event*3.14/17))), yerr = eventDict[event][0][1])		
				elif e >= 6:
					event = e
					pyl.bar(index+1, eventDict[event][0][0], 0.75, color = (math.fabs(math.cos(5*event*3.14/17)), math.fabs(math.sin(2*event*3.14/17)), math.fabs(math.sin(3.4*event*3.14/17))), yerr = eventDict[event][0][1])
			ax.spines['bottom'].set_color('w')
			ax.spines['left'].set_color('w')
			pyl.xticks(index+1+0.75/2.0, labelList, rotation = 90, color = 'w', fontsize = 15)
			ax.yaxis.label.set_color('w')
			ax.tick_params(axis='x', colors='w')
			ax.tick_params(axis='y', colors='w')
			pyl.ylabel('time (minutes post cellularization)', color = 'w', fontsize = 15)
			pyl.savefig('barPlots/'+str(event+1) + '_' + str(T) + 'StackedBarGraphPlot_' + str(sys.argv[1]) + '.pdf', transparent = True, facecolor = 'k', linecolor = 'w') # to change from black to white background, make non-transparent
			pyl.close()			
	except:
		print T
			
### Part 5 - Heat Map graphical analysis
# convert all times to fractional or normalized time (t/(tracheaFill-mry))
# times will be non-dimensional and thus easier to compare
# need to set color boundaries

# goodEventsList = [4,9,21,22,29,33,34]
# goodEventsList = [0,1,2,3,4,5,8,9,10,11,21,22,28,29,30,31,32,33,34]
# goodEventsList = range(numberOfEvents)
goodEventsList = [1,9,11,21,30,33]

if len(Tlist) > 1:
	for species in SpList:
		eventList = []
		normEventList = []
		for event in goodEventsList:
			labelList = []
			propMeanList = []
			normPropMeanList = []
			runningSum = 0
			for T in Tlist:
				try:
					(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
					proportionalMean = meanR[event]/(meanR[33]-meanR[6]) # unclear if these should be the revised estimates in the denominator
					propMeanList.append(math.fabs(proportionalMean))
					labelList.append(T+unichr(176)+'C')
					runningSum += proportionalMean
				except:
					# print 'error in event', event, 'at temperature', T, 'when finding the average time'
					pass
			averageTime = runningSum/len(labelList)
			if averageTime == 0.0:
				averageTime = 0.01
			for T in Tlist:
				try:
					(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
					normalizedProportionalMean = (meanR[event]/(meanR[33]-meanR[6]))/averageTime # unclear if these should be the revised estimates in the denominator
					normPropMeanList.append(normalizedProportionalMean)
				except:
					# print 'error in event', event, 'at temperature', T, 'when finding the proportional time'
					pass
			eventList.append(propMeanList)
			normEventList.append(normPropMeanList)
		mappingList = zip(*normEventList)
		Z = np.array(mappingList)
		pyl.figure(figsize=(11, 7.5))
		pyl.pcolor(Z,cmap=pyl.get_cmap("Greens_r"))
		pyl.colorbar()
		pyl.xlabel('Developmental Event', color = 'w', fontsize = 15)
		indexx = np.arange(len(goodEventsList))
		pyl.xticks(indexx+0.5, map(lambda x:x+1, goodEventsList), color = 'w', fontsize=8)
		pyl.ylabel('Temperature', color = 'w', fontsize = 15)
		pyl.title(species, color = 'w', fontsize = 20)
		indexy = np.arange(len(labelList))
		pyl.yticks(indexy+0.5,labelList, color = 'w', fontsize=15)
		pyl.savefig('heatMaps/'+str(event+1) + '_' + str(species) + 'Heatmap_' + str(sys.argv[1]) + '.pdf', transparent = True, facecolor = 'k', linecolor = 'w')
		print 'checkpoint'
		pyl.close()

backgroudSetting = 'black'
if len(SpList) > 1:
	for T in Tlist:
		eventList = []
		normEventList = []
		for event in goodEventsList:
			labelList = []
			propMeanList = []
			normPropMeanList = []
			runningSum = 0
			for species in SpList:
				try:
					(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
					proportionalMean = meanR[event]/(meanR[33]-meanR[6]) # unclear if these should be the revised estimates in the denominator
					propMeanList.append(math.fabs(proportionalMean))
					labelList.append(species)
					runningSum += proportionalMean
				except:
					# print 'error in event', event, 'for species', species, 'when finding the average time'
					pass
			averageTime = runningSum/len(labelList)
			if averageTime == 0.0:
				averageTime = 0.01
			for species in SpList:
				try:
					(mean, stDev, num, normTiming, median, mad, meanR, stDevR) = meanDict[T+species]
					normalizedProportionalMean = (meanR[event]/(meanR[33]-meanR[6]))/averageTime # unclear if these should be the revised estimates in the denominator
					# normalizedProportionalMean = mean[event]/averageTime
					normPropMeanList.append(normalizedProportionalMean)
				except:
					# print 'error in event', event, 'for species', species, 'when finding the proportional time'
					pass
			eventList.append(propMeanList)
			normEventList.append(normPropMeanList)
		mappingList = zip(*normEventList)
		Z = np.array(mappingList)
		pyl.pcolor(Z,cmap=pyl.get_cmap("PRGn")) # "gist_rainbow_r"))
		pyl.colorbar()
		if backgroudSetting == 'black':
			pyl.xlabel('event', color = 'w', fontsize = 15)
			indexx = np.arange(len(goodEventsList))
			pyl.xticks(indexx+0.5, map(lambda x:x+1, goodEventsList), color = 'w', fontsize=8)
			pyl.ylabel('species', color = 'w', fontsize = 15)
			pyl.title(T, color = 'w', fontsize = 20)
			indexy = np.arange(len(labelList))
			pyl.yticks(indexy+0.5,labelList, color = 'w', fontsize=15)
			pyl.savefig('heatMaps/'+str(T) + 'Heatmap_' + str(sys.argv[1]) + '.pdf', transparent = True, facecolor = 'k', linecolor = 'w')
		elif backgroudSetting == 'white':
			pyl.xlabel('event', color = 'k', fontsize = 15)
			indexx = np.arange(len(goodEventsList))
			pyl.xticks(indexx+0.5, map(lambda x:x+1, goodEventsList), color = 'k', fontsize=8)
			pyl.ylabel('species', color = 'k', fontsize = 15)
			pyl.title(T+unichr(176)+'C', color = 'k', fontsize = 20)
			indexy = np.arange(len(labelList))
			pyl.yticks(indexy+0.5,labelList, color = 'k', fontsize=15)
			pyl.savefig('heatMaps/'+str(T) + 'Heatmap_' + str(sys.argv[1]) + '.pdf', transparent = True, linecolor = 'w') # , facecolor = 'k', linecolor = 'w')
		print 'checkpoint'
		pyl.close()
