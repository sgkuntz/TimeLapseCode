# sample command: Dropbox Hypsiglena $ python ImageProcessingWrapper_v6.py EmbryoSeries#_u.txt trial#
# must be run from Dropbox folder on the iMac

# this program will consolidate information from EmbryoSeries.ods
# EmbryoSeries is transposed from the original source file that Image ProcessingWrapper was written for
# necessary outputs will be the date, position, filename, orientation, filename structure, membrane reaches yolk, trachea fills, species, temperature, timeDilation
# output will be used for Matlab analysis

import os, sys, re, string, math, array, shutil

# read in list of dates, positions, and file names
fileNamesData = []
counter = 0
# for each file name, fullDate is the first folder, abbDate is this without '-', pos is the second folder with letters replaced and Pos gone, fileName is the file name sans folders
for (dirpath, dirnames, filenames) in os.walk("../../../Volumes/Agkistrodon/Imaging/"): # go to ../../Volumes/Agkistrodon/Imaging
# for (dirpath, dirnames, filenames) in os.walk("../Desktop/Imaging/"):
	if 'aligned' not in dirpath:
		for filename in filenames:
			if '_t000_' in filename: # search for all file names with */*/*_t000_* and */*/*_t0000_*
				# print dirpath
				fullDate = dirpath[37:45] 
				abbDate = fullDate[:2]+fullDate[3:5]+fullDate[6:]
				pos = dirpath[49:] # should this be 50?
				fileNamesData.append((fullDate, abbDate, pos, filename))
			elif '_t0000_' in filename:
				fullDate = dirpath[37:45]
				abbDate = fullDate[:2]+fullDate[3:5]+fullDate[6:]
				pos = dirpath[49:]
				fileNamesData.append((fullDate, abbDate, pos, filename))
			elif '_t000.' in filename:
				fullDate = dirpath[37:45]
				abbDate = fullDate[:2]+fullDate[3:5]+fullDate[6:]
				pos = dirpath[49:] # should this be 49?
				fileNamesData.append((fullDate, abbDate, pos, filename))
			elif '_t0000.' in filename:
				fullDate = dirpath[37:45]
				abbDate = fullDate[:2]+fullDate[3:5]+fullDate[6:]
				pos = dirpath[49:]
				fileNamesData.append((fullDate, abbDate, pos, filename))
	counter += 1
	print "path", counter,"done"

print "images screened"
# print fileNamesData

# next, read in the time-lapse spreadsheet [date, Position, translatedPosition, Species, Temperature, timeDilation, orientation, membrane reaches yolk, trachea fills]
embryoSeriesFile = open(str(sys.argv[1]),'r') # this should be a text file generated from EmbryoSeries.ods
positionData = []
for k in embryoSeriesFile:
	line = re.split("[\s]",k)
	try:
		positionData.append((line[0], line[1], line[2], line[5], line[6], line[7], line[8], line[16], line[44]))
	except:
		try:
			positionData.append((line[0], line[1], line[2], line[5], line[6], line[7], line[8], '0', '0'))
		except:
			error = 0

# can generate the already completed list from the dropbox list or the imaging list -> which is better to use? -> can check both
dirList_dropbox = os.listdir("../Dropbox/") # check the dropbox directory
dirList_imaging = os.listdir("../Desktop/Imaging/")
dirList = dirList_dropbox + dirList_imaging
alreadyCompletedList = []
alreadyCompletedDates = []

print "analysis files screened"

for filename in dirList:
	if 'EventTimesTrial' in filename:
		alreadyCompletedList.append(filename)
		line = re.split("_", filename)
		alreadyCompletedDates.append((line[1], line[2][3:])) # full date, Position
		
		
outputFile = open('TimeLapseInstructions_all_' + str(sys.argv[2]) + '.txt', 'w') # this should be the TimeLapseInstructions.txt file
outputFileNotTested = open("../Desktop/Imaging/"+'TimeLapseInstructions_notRun_' + str(sys.argv[2]) + '.txt','w')
outputFileTested = open("../Desktop/Imaging/"+'TimeLapseInstructions_run_' + str(sys.argv[2]) + '.txt','w')
outputFileNotAnnotated = open("../Desktop/Imaging/"+'TimeLapseInstructions_notAnnotated_' + str(sys.argv[2]) + '.txt','w')
outputData = []
outputFile.write("Date" + '\t' + "Position" + '\t' + "FileName" + '\t' + "Orientation" + '\t' + "Structure" + '\t' + "mry" + '\t' + "trachea" + '\t' + "Species" + '\t' + "temperature" + '\t' + "scaling"  + '\t' + "zPosition" + '\n')

for (fullDate, shortDate, translatedPosition, fileName) in fileNamesData: # info from the image files
	for (abbrDate, Pos, numericalPos, species, temperature, timeDilation, orientation, mry, trachea) in positionData: # info from EmbryoSeries file
		if len(abbrDate) == 5:
			abbrDate = '0'+abbrDate
		if shortDate == abbrDate:
			if translatedPosition == Pos:
				z = 0
				# need to redefine the fileName to get the filenameStructure
				# always of the form PosXXX_SXXX_CropXXX_tXXX*.jpg, where * can include an extra time digit (1), _zX (3), _chXX (5)
				if len(fileName) == 28 or len(fileName) == 29: # no z or channel info
					filenameStructure = 'Pos*_S*_Crop*_t*.jpg'
				elif len(fileName) == 31 or len(fileName) == 33 or len(fileName) == 36: # t is 3 digits, need to incorporate z information
					filenameStructure = 'Pos*_S*_Crop*_t*'+fileName[24:]
					if len(fileName) == 31 or len(fileName) == 36:
						z = fileName[26:27]
				elif len(fileName) == 32 or len(fileName) == 34 or len(fileName) == 37: # t is 4 digits, need to incorporate z information
					filenameStructure = 'Pos*_S*_Crop*_t*'+fileName[25:]
					if len(fileName) == 32 or len(fileName) == 37:
						z = fileName[27:28]
				else:
					print "error in file name", fileName
				if timeDilation == '':
					scaling = 1
				else:
					scaling = int(timeDilation[1])

				outputData.append((fullDate, Pos, fileName, orientation, filenameStructure, mry, trachea, species, temperature, scaling, z))
				outputFile.write(str(fullDate) + '\t' + 'Pos' + str(Pos) + '\t' + str(fileName) + '\t' + str(orientation) + '\t' + str(filenameStructure) + '\t' + str(mry) + '\t' + str(trachea) + '\t' + str(species) + '\t' + str(temperature) + '\t' + str(scaling)  + '\t' + str(z) + '\n')

				# outputFileAbridged.write((fullDate, Pos, fileName, orientation, filenameStructure, mry, trachea, species, temperature, scaling, z))
				# print alreadyCompletedDates[0]
				if (fullDate, Pos) not in alreadyCompletedDates:
					if species <> '':
						if mry <> '':
							if trachea <> '':
								if orientation <> '':
									outputFileNotTested.write(str(fullDate) + '\t' + 'Pos' + str(Pos) + '\t' + str(fileName) + '\t' + str(orientation) + '\t' + str(filenameStructure) + '\t' + str(mry) + '\t' + str(trachea) + '\t' + str(species) + '\t' + str(temperature) + '\t' + str(scaling)  + '\t' + str(z) + '\n')
								else:
									outputFileNotAnnotated.write(str(fullDate) + '\t' + 'Pos' + str(Pos) + '\t' + str(fileName) + '\t' + str(orientation) + '\t' + str(filenameStructure) + '\t' + str(mry) + '\t' + str(trachea) + '\t' + str(species) + '\t' + str(temperature) + '\t' + str(scaling)  + '\t' + str(z) + '\n')
							else:
								outputFileNotAnnotated.write(str(fullDate) + '\t' + 'Pos' + str(Pos) + '\t' + str(fileName) + '\t' + str(orientation) + '\t' + str(filenameStructure) + '\t' + str(mry) + '\t' + str(trachea) + '\t' + str(species) + '\t' + str(temperature) + '\t' + str(scaling)  + '\t' + str(z) + '\n')
						else:
							outputFileNotAnnotated.write(str(fullDate) + '\t' + 'Pos' + str(Pos) + '\t' + str(fileName) + '\t' + str(orientation) + '\t' + str(filenameStructure) + '\t' + str(mry) + '\t' + str(trachea) + '\t' + str(species) + '\t' + str(temperature) + '\t' + str(scaling)  + '\t' + str(z) + '\n')
				elif (fullDate, Pos) in alreadyCompletedDates:
					outputFileTested.write(str(fullDate) + '\t' + 'Pos' + str(Pos) + '\t' + str(fileName) + '\t' + str(orientation) + '\t' + str(filenameStructure) + '\t' + str(mry) + '\t' + str(trachea) + '\t' + str(species) + '\t' + str(temperature) + '\t' + str(scaling)  + '\t' + str(z) + '\n')
				else:
					print fullDate, Pos