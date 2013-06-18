# try writing a gui for image analysis (select the image that is at the right time point)
# python ImageSelectionGUI_v7.py 5 -m -e cephalic_furrow_formation
# python ImageSelectionGUI_v7.py 5 -e cephalic_furrow_formation
from Tkinter import *
from PIL import Image, ImageTk
import os, re, sys, string

# select events
eventList = []
possibleEventList = ['posterior_gap','pole_bud_appears','nuclei_at_periphery','pole_cells_form','yolk_contraction','cellularization_begins','membrane_reaches_yolk',
				'pole_cells_migrate','cephalic_furrow_formation','pole_cell_invagination','cephalic_furrow_reclines','amnioproctodeal_invagination',
				'transversal_fold_formation','anterior-midgut_primordium','stomodeal_plate_forms','stomodeum_invagination','clypeolabral_lobe_forms',
				'germband_maxima','clypeolabrum_rotates','posterior_gap_before_germband_shortens','germband_retraction_starts',
				'amnioserosa','germband_retracted','dorsal_divot','clypeolabrum_retracts','anal_plate_forms','midgut_unified','clypeolabrum_ventral_lobe_even',
				'gnathal_lobes_pinch','head_involution_done','heart-shaped_midgut','convoluted_gut','muscle_contractions','trachea_fill','hatch']
				
if '-e' in sys.argv:
	for entry in sys.argv:
		if entry in possibleEventList:
			eventList.append(entry)
			mainEvent = entry
else:
	eventList = possibleEventList
	
if mainEvent not in os.listdir('.'):
	os.mkdir(mainEvent)

# read in the TimeLapseInstructions_run_trialX.txt, this can be tailored to what T and Sp are desired
	
speciesToAnalyze = str(raw_input('species to analyze: '))
temperatureToAnalyze = str(raw_input('temperature to analyze: '))

if '-m' in sys.argv:
	timeLapseInfoFile = open(str(raw_input('Time lapse instructions file input (eg. TimeLapseInstructions_run_trialX.txt): ')), 'r')
	recordOutput = open(str(raw_input('Output file (eg. TimeLapseCorrections_trialX.txt): ')), 'w')
else:
	timeLapseInfoFile = open('TimeLapseInstructions_run_trial15.txt','r') # default, can update to other trials
	if len(temperatureToAnalyze) == 4:
		recordOutput = open(mainEvent + '/'+'TimeLapseCorrections_trial15_'+ speciesToAnalyze + '_' + temperatureToAnalyze +'.txt', 'w')
	elif len(temperatureToAnalyze) ==2:
		recordOutput = open(mainEvent + '/'+'TimeLapseCorrections_trial15_'+ speciesToAnalyze + '_' + temperatureToAnalyze +'.0.txt', 'w')

recordOutput.write('conditions' + '\t' + 'event' + '\t' + 'date/position' + '\t' + 'finalTimepoint' + '\t' + 'dilation' + '\n')
recordlist = []
numberOfImages = int(sys.argv[1])

imageDirectory = '/Volumes/Agkistrodon/Imaging/'
# imageDirectory = '/Desktop/Imaging/'
	
timelapseInfoList = []
conditionDict = {}
for i in timeLapseInfoFile:
	[date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temperature, scaling, zStack, filler] = re.split("[\s]",i)
	# line = re.split("[\s]",i)
	if speciesToAnalyze == strain:
		if temperatureToAnalyze == temperature:
			timelapseInfoList.append((date, position, firstFilename, orientation, fileNameStructure, mryTime, tfTime, strain, temperature, scaling, zStack))
			if temperature+strain not in conditionDict:
				conditionDict[temperature+strain] = []
			if (date+'_'+position, firstFilename[:25], orientation, len(firstFilename), mryTime, tfTime, strain, temperature, scaling) not in conditionDict[temperature+strain]:
				conditionDict[temperature+strain].append((date+'_'+position, firstFilename[:25], orientation, len(firstFilename), mryTime, tfTime, strain, temperature, scaling))
			conditionDict[temperature+strain] = sorted(list(set(conditionDict[temperature+strain])))

# determine the conditions present	
conditionList = conditionDict.keys()
conditionList.sort()

statusCall = 'proceed'
imageCall = 'aligned'

class App: # define the App, which should record which button is pressed
	def __init__(self, master, timepoint, animal, event, condition, date, position, nameA, nameB, timeDepth, tNameLength, dilation):
		frame = Frame(master)
		frame.grid()
		self.frame = frame
		
		self.quit_button = Button(text="No matches", fg = "red", padx = 50, pady = 15, command=self.noMatches)
		self.quit_button.grid(row=1, column=2)
		
		self.next_button = Button(text="Next Images", padx = 50, pady = 15, command=self.nextImages)
		self.next_button.grid(row=1, column = 3)
		
		self.previous_button = Button(text="Previous Images", padx = 40, pady = 15, command=self.previousImages)
		self.previous_button.grid(row=1, column = 1)

		self.jump_next = Button(text=">>", padx = 20, pady = 15, command=self.jumpNext)
		self.jump_next.grid(row=1, column = 4)

		self.jump_previous = Button(text="<<", padx = 20, pady = 15, command=self.jumpPrevious)
		self.jump_previous.grid(row=1, column = 0)
		
		if imageCall == 'aligned':
			self.switch_to_originals = Button(text="Switch to originals", padx = 40, pady = 10, command=self.switchOriginal)
			self.switch_to_originals.grid(row=0, column = 4)
		elif imageCall == 'original':
			self.switch_to_aligned = Button(text="Switch to aligned", padx = 40, pady = 10, command=self.switchAligned)
			self.switch_to_aligned.grid(row=0, column = 4)
			
		self.increaseZstack = Button(text="+ z", padx = 20, pady = 5, command=self.increaseZ)
		self.increaseZstack.grid(row=0, column = 3)
		
		self.decreaseZstack = Button(text="- z", padx = 20, pady = 5, command=self.decreaseZ)
		self.decreaseZstack.grid(row=0, column = 2)
		
		self.v = IntVar()
		self.v.set(0)
		self.timepoint = int(timepoint)
		self.animal = animal
		self.event = event
		self.condition = condition
		self.dilation = dilation

		for number in range(timeDepth):
			if self.timepoint > timeDepth:
				if imageCall == 'aligned':
					importImage = Image.open(imageDirectory+ date +'/'+ position +'/aligned/'+ nameA+'0'*(tNameLength-len(str(self.timepoint-3+number)))+str(self.timepoint-3+number)+nameB)
				elif imageCall == 'original':
					importImage = Image.open(imageDirectory+ date +'/'+ position +'/'+ nameA+'0'*(tNameLength-len(str(self.timepoint-3+number)))+str(self.timepoint-3+number)+nameB)
					importImage = importImage.resize((250,705), Image.ANTIALIAS)
			else:
				if imageCall == 'aligned':
					importImage = Image.open(imageDirectory+ date +'/'+ position +'/aligned/'+ nameA+'0'*(tNameLength-len(str(self.timepoint+number)))+str(self.timepoint+number)+nameB)
				elif imageCall == 'original':
					importImage = Image.open(imageDirectory+ date +'/'+ position +'/'+ nameA+'0'*(tNameLength-len(str(self.timepoint+number)))+str(self.timepoint+number)+nameB)
					importImage = importImage.resize((250,705), Image.ANTIALIAS)
			embryoImage = ImageTk.PhotoImage(importImage)
			label = Label(image=embryoImage)
			label.image = embryoImage
			b = Radiobutton(master, image = label.image, variable = self.v, value = number, indicatoron=0, command= self.matchFound)
			b.grid(row=2, column=number)
		
	def matchFound(self):
		self.frame.quit()
		recordlist.append((self.condition, self.event, self.animal, self.timepoint, 1+self.v.get()))
		print 'recording data:', self.condition, self.event, self.animal, self.timepoint, 1+self.v.get()
		recordOutput.write(str(self.condition) + '\t' + str(self.event) + '\t' + str(self.animal) + '\t' + str(self.timepoint + 1 + self.v.get() - 3) + '\t' + self.dilation + '\n')
		global statusCall
		statusCall = 'proceed'
	def nextImages(self):
		self.frame.quit()
		global statusCall
		statusCall = 'next'
	def previousImages(self):
		self.frame.quit()
		global statusCall
		statusCall = 'previous'
	def jumpNext(self):
		self.frame.quit()
		global statusCall
		statusCall = 'jumpnext'		
	def jumpPrevious(self):
		self.frame.quit()
		global statusCall
		statusCall = 'jumpprevious'
	def noMatches(self):
		self.frame.quit()
		recordlist.append((self.condition, self.event, self.animal, self.timepoint, 0))
		recordOutput.write(str(self.condition) + '\t' + str(self.event) + '\t' + str(self.animal) + '\t' + 'NaN' + '\t' + self.dilation + '\n')
		global statusCall
		statusCall = 'proceed'
	def switchOriginal(self):
		self.frame.quit()
		global imageCall
		imageCall = 'original'
		global statusCall
		statusCall = 'refresh'
	def switchAligned(self):
		self.frame.quit()
		global imageCall
		imageCall = 'aligned'
		global statusCall
		statusCall = 'refresh'
	def increaseZ(self):
		self.frame.quit()
		global statusCall
		statusCall = 'zUp'
	def decreaseZ(self):
		self.frame.quit()
		global statusCall
		statusCall = 'zDown'
		

root = Tk() # this generates the window

recordlist_AllConditions = []
focusFileNames = os.listdir("./FocusFiles/")
for condition in conditionList:
	recordlist_AllEvents = []
	for devEvent in eventList:
		event = possibleEventList.index(devEvent)
#		for index, x in enumerate(possibleEventList):
#			if x == devEvent:
#				event = index
		print devEvent, event
#	for event in range(35): # len(eventList)):
		recordlist_AllAnimals = []
		for (animal, firstFilename, orientation, fileStructure, mry, tf, species, temperature, scaling) in conditionDict[condition]:
			try:
				zRecord = []
				tList = []
				focusList = []
				for focusFileName in focusFileNames:
					# print focusFileName, animal[:8], animal[9:]
					if animal[:8]+'_'+animal[9:] in focusFileName:
						zRecord.append(focusFileName[-5:-4])
				if zRecord == []:
					print "no focus files", firstFilename
					break # zPos = 0
				for zPos in zRecord:
					# print 'ImageProcessing/TrialSubsets_'+animal[:8]+'_'+animal[9:]+'_'+str(zPos)+'_EventTimesTrial_avg_'+species+'_'+temperature+'_'+ scaling+'.txt'
					eventFile = open('ImageProcessing/TrialSubsets_'+animal[:8]+'_'+animal[9:]+'_'+str(zPos)+'_EventTimesTrial_avg_'+species+'_'+temperature+'_'+ scaling+'.txt', 'r')
					line = eventFile.read()
					entry = re.split("[\s]",line)
					focusFile = open('FocusFiles/FocusScoresUnaligned_'+animal[:8]+'_'+animal[9:]+'_'+zPos+'.txt','r')
					lines = focusFile.readlines()
					focusScore = re.split("[\s]",lines[int(float(entry[event]))])
					tList.append((int(float(entry[event])),int(zPos)))
					focusList.append(float(focusScore[1]))
				# find the max zScore and return the zPos and time
				[estTime,zPos] = tList[focusList.index(max(focusList))]
				if fileStructure == 28 or fileStructure == 31 or fileStructure == 33 or fileStructure == 36:
					tCallLength = 3
				elif fileStructure == 29 or fileStructure == 32 or fileStructure == 34 or fileStructure == 37:
					tCallLength = 4
				if 	fileStructure == 28 or fileStructure == 29:
					filenameB = '.jpg'
				elif fileStructure == 31 or fileStructure == 32:
					filenameB = '_z'+str(zPos)+'.jpg'
				elif fileStructure == 33 or fileStructure == 34:
					filenameB = '_ch00.jpg'
				elif fileStructure == 36 or fileStructure == 37:
					filenameB = '_z'+str(zPos)+'_ch00.jpg'
				w = Label(root, text="Pick the image that matches "+ devEvent, wraplength=275)
				w.grid(row=0, column = 0)
				app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
				print 'matlab estimate is:', estTime # 'program will run, matlab estimate is:', estTime
				recordlist = []
				root.mainloop() # start running the program
				while statusCall <> 'proceed':
					while statusCall == 'refresh':
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
					while statusCall == 'zUp':
						if len(filenameB) == 12 or len(filenameB) == 7:
							filenameC = '_z' + str(int(filenameB[2:3])+1) + filenameB[3:] # update z via filenameB, see if file exists
							zStatus = os.path.exists(imageDirectory+ animal[:8] +'/'+ animal[9:] +'/aligned/'+ firstFilename[:21]+'0'*(tCallLength)+filenameC)
						else: zStatus = False
						if zStatus == True:
							filenameB = filenameC
						else:
							print 'already at max z'
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
					while statusCall == 'zDown':
						if len(filenameB) == 12 or len(filenameB) == 7:
							filenameC = '_z' + str(int(filenameB[2:3])-1) + filenameB[3:] # update z via filenameB, see if file exists 
							zStatus = os.path.exists(imageDirectory+ animal[:8] +'/'+ animal[9:] +'/aligned/'+ firstFilename[:21]+'0'*(tCallLength)+filenameC)
						else: zStatus = False
						if zStatus == True:
							filenameB = filenameC
						else:
							print 'already at min z'
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
					while statusCall == 'next':
						estTime += numberOfImages
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
					while statusCall == 'previous':
						if estTime > numberOfImages*1.5:
							estTime -= numberOfImages
						else:
							print 'already at beginning'
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
					while statusCall == 'jumpnext':
						estTime += numberOfImages*5
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
					while statusCall == 'jumpprevious':
						if estTime > numberOfImages*6:
							estTime -= numberOfImages*5
						app = App(root, estTime, animal, devEvent, condition, animal[:8], animal[9:], firstFilename[:21], filenameB, numberOfImages, tCallLength, scaling)
						root.mainloop()
				# print 'program has run, status is: ', statusCall
				recordlist_AllAnimals.append(recordlist)
			except:
				print "no EventTimesTrial file", firstFilename
		recordlist_AllEvents.append((event, recordlist_AllAnimals)) # save times for all events	
	recordlist_AllConditions.append((condition, recordlist_AllEvents)) # save times for all conditions
# save the data or process it
