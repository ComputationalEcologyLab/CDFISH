# -------------------------------------------------------------------------------------------------
# CDFISH_Disperse.py
# Author: Erin L Landguth
# Created: December 2010
# Description: This is the function/module file for dispersal processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# CDFISH functions
try:
	from CDFISH_Modules import *
except ImportError:
	raise ImportError, "CDFISH Modules required."

# Python specific functions
import pdb, random, copy, sys, heapq

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# --------------------------------------------------------------------------
def logMsg(outf,msg):
	'''
	logMsg() --log file message handler.
	Inputs:
	outf - open file handle
	msg -- string containing formatted message
	--always outputs to log file by default.
	--using msgVerbose, can be set to "Tee" output to stdout as well
	'''
	outf.write(msg+ '\n')
	if msgVerbose:
		print("%s"%(msg))
		
	# End::logMsg()

# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(x[1] for x in lst)
	n=random.uniform(0,wtotal)
	count = 0
	for item, weight in lst:
		if n < weight:
			break
		n = n-weight
		count = count + 1
	return item,count
	
	#End::w_choice_general()
	
# ---------------------------------------------------------------------------------------------------	 
def GetOffspringProb(offspring,subpop,dispmovethresh,dispmoveno,source_subpop):
	'''
	GetOffspringProb()
	Adds subpopulation, source to subpop cost and probability to offspring list.
	'''	
	for ioff in xrange(len(offspring)):
	
		# Get subpop assignment -- append to offspring in spot 4
		subpopassignment = subpop[offspring[ioff][0]]
		offspring[ioff].append(subpopassignment)		
		
		# Get source to subpop cost (-1 for index) -- append to offspring in spot 5
		costtosource = source_subpop[int(subpopassignment)-1]
		offspring[ioff].append(costtosource)
			
		# Check if cost is greater than dispersal threshold, force 0
		if costtosource >= float(max(dispmovethresh)):
			probability = 0.0
			
		# Else assign linear probability draw if dispmoveno == 1
		elif dispmoveno == '1':
			probability = -(float(costtosource)/dispmovethresh[0]) + 1.0
		# Else assign inverse square probability draw
		elif dispmoveno == '2':
			probability = 1./(pow(float(costtosource),2)+1)
		# Else assign nearest neighbor prob draw
		elif dispmoveno == '3':
			print('Not coded yet...email Erin.')
			sys.exit(-1)
		# Else assign panmictic prob draw
		elif dispmoveno == '4':
			probability = 0.0
		# Else assign neg binomial prob draw
		elif dispmoveno == '5':
			print('Not coded yet...email Erin.')
			sys.exit(-1)
		# Else assign the min-max linear prob draw
		elif dispmoveno == '6':
			# If all source locations are the same.
			if min(source_subpop)==max(source_subpop):
				probability = 1.
			else:
				probability = (1./(min(source_subpop)-max(source_subpop)))*(costtosource-min(source_subpop))+1
		# Else assign the min-max scaled by use defined cost
		elif dispmoveno == '7':			
			# If they are the same.
			if min(dispmovethresh)==max(dispmovethresh):
				probability = 1.
			else:
				probability = (1./(float(min(dispmovethresh))-float(max(dispmovethresh))))*(costtosource-float(min(dispmovethresh)))+1
		else:
			print('Dispersal movement function is not defined.')
			sys.exit(-1)		
		
		# Then append probability to offspring to spot 6
		offspring[ioff].append(probability)
		
	return offspring
		
	#End::GetOffspringProb
	
# ---------------------------------------------------------------------------------------------------	
def Do1LocusSelection(offspring,fitvals1,gridspot):
	'''
	Do1LocusSelection()
	This function calculates offsprings differential mortality, ie,
	offspring viability selection, for the 1-locus selection model.
	'''
	
	# If L0A0|L0A0 -- loci under selection:
	if offspring[3][0] == 2:

		# The grab it's fitness values
		differentialmortality = fitvals1[gridspot][0]/100
																
	# If L0A0|L0A1 -- loci under selection:
	elif offspring[3][0] == 1:

		# The grab it's fitness values
		differentialmortality = fitvals1[gridspot][1]/100
																															
	# If L0A1|L0A1 -- loci under selection
	else:
		
		# The grab it's fitness values
		differentialmortality = fitvals1[gridspot][2]/100
	
	return differentialmortality
	
	# End::Do1LocusSelection()

# ---------------------------------------------------------------------------------------------------	
def GetStrayers(tempoffspring,straying,nosubpops,cdmatrix,dispmovethresh,\
	dispmoveno,source_subpop,OffDispMortality):
	'''
	GetStrayers()
	From offspring pool calculate the the number of strayers: First based on 
	random draw based on straying probabilities, then based on movement function
	weighted random draw. Reassign offpsring home population to new population it 
	will stray to.
	'''
	
	# Loop through offspring
	for ioff in xrange(len(tempoffspring)):
	
		# Subpopulation of offspring - careful of indexing
		homesubpop = int(tempoffspring[ioff][4])
	
		# Storage variables for straying
		probtostray = []
	
		# Flip a coin to see if this offspring strays
		randstray = random.random()
		
		# If it strays then: (BE CAREFUL of indexing here...)
		if randstray < (straying[homesubpop-1])/100.:
			
			# Get cost values from the subpopulation - subce subpops start at 1 count
			for icost in xrange(1,nosubpops+1):
				
				# Skip own subpopulation and less than dispersal thresholds
				if icost != homesubpop and \
				cdmatrix[homesubpop][icost] <= float(max(dispmovethresh)) and \
				float(dispmovethresh[0]) != 0.0:
					
					# Assign a probability: linear function
					if dispmoveno == '1':
						# Weight by the sources probability to get back (-1 into index)
						source_subpop_weight = -(float(source_subpop[icost-1])/dispmovethresh[0]) + 1.0
						# If it is less than dispmovement
						if source_subpop_weight < 0.0:
							source_subpop_weight = 0.0
						probtostray.append([icost,\
						source_subpop_weight*(-(float(cdmatrix[homesubpop][icost])/dispmovethresh[0]) + 1.0)])
					# Assign a probability: inverse square function
					elif dispmoveno == '2':
						# Weight by the sources probability to get back 
						source_subpop_weight = 1./(pow(source_subpop[icost-1],2)+1)
						# If it is less than dispmovement
						if source_subpop_weight < 0.0:
							source_subpop_weight = 0.0
						probtostray.append([icost,\
						source_subpop_weight*(1./(pow(float(cdmatrix[homesubpop][icost]),2)+1))])
					# Assign a probability: NN
					elif dispmoveno == '3':
						print('Not coded yet...email Erin.')
						sys.exit(-1)
					# Assign probability: panmixia
					elif dispmoveno == '4':
						probtostray.append([icost,0.0])
					# Assign probability: neg binomial
					elif dispmoveno == '5':
						print('Not coded yet...email Erin.')
						sys.exit(-1)
					# Else assign the min-max linear prob draw
					elif dispmoveno == '6':
						# If all source locations are the same.
						if min(source_subpop)==max(source_subpop):
							source_subpop_weight = 1.
						else:
							source_subpop_weight = (1./(min(source_subpop)-max(source_subpop)))*(source_subpop[icost-1]-min(source_subpop))+1
						cdmatarray = np.asarray(cdmatrix[homesubpop])
						minval = np.min(cdmatarray[np.nonzero(cdmatarray)[0]])
						maxval = np.min(heapq.nlargest(2,cdmatarray))
						# If all the same
						if minval == maxval:
							probtostray.append([icost,\
							source_subpop_weight*1.])									
						else:
							#If it is less than dispmovement
							if source_subpop_weight < 0.0:
								source_subpop_weight = 0.0
							probtostray.append([icost,\
							source_subpop_weight*((1./(minval-maxval))*(cdmatrix[homesubpop][icost]-minval)+1)])										
					# Else assign the min-max scaled by use defined cost
					elif dispmoveno == '7':
						# If they are the same.
						if min(dispmovethresh)==max(dispmovethresh):
							source_subpop_weight = 1.
						else:
							source_subpop_weight = (1./(float(min(dispmovethresh))-float(max(dispmovethresh))))*(source_subpop[icost-1]-float(min(dispmovethresh)))+1
						cdmatarray = np.asarray(cdmatrix[homesubpop])
						minval = np.min(cdmatarray[np.nonzero(cdmatarray)[0]])
						maxval = np.min(heapq.nlargest(2,cdmatarray))
						# If all the same
						if minval == maxval:
							probtostray.append([icost,\
							source_subpop_weight*1.])				
						else:
							#If it is less than dispmovement
							if source_subpop_weight < 0.0:
								source_subpop_weight = 0.0
							probtostray.append([icost,\
							source_subpop_weight*((1./(minval-maxval))*(cdmatrix[homesubpop][icost]-minval)+1)])
					
			# If probarray equals zero, then that means this offspring did not stray (but died)
			if len(probtostray) == 0.0 or sum(np.asarray(probtostray)[:,1]) == 0.0:
			
				# Keep track of offspring that died on dispersal route
				OffDispMortality.append(tempoffspring[ioff])
				
			# Else if values in probarray, then grab a weighted random draw
			else:
				tupitemselect = w_choice_general(probtostray)
				itemselect = tupitemselect[0]
				
				# Then reassign the subpopulation number for this offspring, appending 'from'
				tempoffspring[ioff][4] = str(itemselect)+'from'+str(homesubpop)
	
	StayerTup = tempoffspring,OffDispMortality	
	return StayerTup
	#End::GetStrayers()
	
# ---------------------------------------------------------------------------------------------------	 
def DoDisperse(offspringno,freegrid,offspring,dispmoveno,dispmovethresh,gen,\
Migrants,OpenLocations,loci,alleles,nogrids,\
xgridcopy,ygridcopy,allelefreqlst,logfHndl,\
newmortperc,subpop,source_subpop,residency,\
nosubpops,Residors,cdmatrix,Strayer,straying,DispDeaths,\
OffLeftOver,subpoptotals,cdevolveans,fitvals1,burningen,
SelectionDeaths):

	'''
	DoDisperse()
	Disperse the new offspring to empty spots on grid either through 
	residency or migration or straying
	Input: Units of dipsersal, movement function,
	offspring, freegrid, cdmatrix 
	Output: OffDisperseIN = [offspring,freegrid,name,[offspringgenes]] 	
	'''	
		
	# Recall offspring = [0-bearpairmother][1-bearpairfather]
	# [2-sex][3-genes][4-subpop][5-costtosubpop][6-probabilitytosubpop]
	# [7-freegrid dispersed to][8-newidentification]
	
	# Create variable to store offspring that return home
	# THis will be [[offspring],freegrid,uniquename]
	OffDisperseIN=[]
	
	# Create variable to store offspring that died on dispersal route
	OffDispMortality = []
	
	# Create variable to store offspring that strayed
	OffDispStray = []
	
	# Create variable to store offspring that migrate
	OffDispMigrate = []
			
	# Get unique subpopulations
	unisubpops = np.unique(subpop)
	
	# Loop through offspring to assign subpop, cost, and probability
	offspring = GetOffspringProb(offspring,subpop,dispmovethresh,dispmoveno,source_subpop)
		
	# Copy offpspring to delete from
	tempoffspring = copy.deepcopy(offspring)
	
	# Shuffle offspring for random selection
	shuffle(tempoffspring)
	shuffle(freegrid)
	
	# Copy freegrid for deleting and adding to, to avoid inf loop
	tempfreegrid = copy.deepcopy(freegrid)
	
	# Add spot to track dispersing deaths for cdevolve
	SelectionDeaths.append([])
			
	# ------------------
	# Residency
	# ------------------	
	
	# Storage variables for residency section
	tallyresidor = []	
	opensubgrids = []
	subgridswithresident = []
	
	# Create list of lists storage spots for number of subgrids
	for i in xrange(nosubpops):
		opensubgrids.append([])
	
	# Count up the unique number of subgrids appending to opensubgrids
	for i in xrange(len(freegrid)):
		# Loop through unique subpops
		for j in xrange(nosubpops):
			# Store orignal locations in variable
			if subpop[freegrid[i]] == unisubpops[j]:
				opensubgrids[int(unisubpops[j])-1].append(freegrid[i])
	
	# Loop through the number of subpops to get the number of residors
	for isubgrid in xrange(nosubpops):
	
		# Calculate the number of freegrid locations that will get a residor in each subgrid
		residornumber = round(len(opensubgrids[isubgrid])*(residency[isubgrid]/100.))
		tallyresidor.append(residornumber)
	
	# Get the random free subgrids that an offspring will reside in
	for isubgrid in xrange(nosubpops):
		subgridswithresident.append(random.sample(opensubgrids[isubgrid],int(tallyresidor[isubgrid])))
	#pdb.set_trace() # when the offspring is less than grid throughs index error....
	# Loop through each subpop
	for isubgrid in xrange(nosubpops):
	
		# Loop through each subpop residency spot
		for ires in xrange(int(tallyresidor[isubgrid])):
		
			# Grab the offspring in that subpop and loop through those....
			#which(tempoffspring==)
			
			# Index for not stepping out of bounds
			ioff = 0
			
			# Loop through offspring
			while ioff < len(tempoffspring):
				
				# If selection is on
				if cdevolveans == '1' and gen >= burningen:
			
					# If offspring is from the grid's subpopulation
					if tempoffspring[ioff][4] == subpop[subgridswithresident[isubgrid][ires]]:
						
						# Call 1-locus selection model
						differentialmortality = Do1LocusSelection(tempoffspring[ioff],fitvals1,subgridswithresident[isubgrid][ires])
												
						# Then flip the coin to see if offspring survives its location
						randcheck = rand()
						
						# If offspring did not survive: break from loop, move to next offspring
						if randcheck < differentialmortality:
							SelectionDeaths[gen].append(1)
							# Delete this offspring from pool							
							del(tempoffspring[ioff])
							ioff = ioff + 1
							continue
						
						# If offspring survived						
						else:
							
							# Assign offspring to that freegrid spot
							recd = [tempoffspring[ioff],subgridswithresident[isubgrid][ires],'Residor_'+'gen'+str(gen)+\
								'pop'+tempoffspring[ioff][4]+'F'+str(tempoffspring[ioff][0])+'M'+str(tempoffspring[ioff][1])]
							OffDisperseIN.append(recd)
							SelectionDeaths[gen].append(0)
							
							# Delete this offspring from pool
							del(tempoffspring[ioff])
							
							# Delete this freegrid spot
							delspot = np.where(freegrid==subgridswithresident[isubgrid][ires])[0]
							freegrid = np.delete(freegrid,delspot)
							tempfreegrid = np.delete(tempfreegrid,delspot)
												
							# Break from loop
							break
							
					# If offspring is not from the grid's subpopulation
					else:
							
						# Update the counter
						ioff = ioff + 1
						
				# If selection is off
				else:
			
					# If offspring is from the grid's subpopulation
					if tempoffspring[ioff][4] == subpop[subgridswithresident[isubgrid][ires]]:
						
						# Assign offspring to that freegrid spot
						recd = [tempoffspring[ioff],subgridswithresident[isubgrid][ires],'Residor_'+'gen'+str(gen)+\
							'pop'+tempoffspring[ioff][4]+'F'+str(tempoffspring[ioff][0])+'M'+str(tempoffspring[ioff][1])]
						OffDisperseIN.append(recd)
						SelectionDeaths[gen].append(0)
						
						# Delete this offspring from pool
						del(tempoffspring[ioff])
						
						# Delete this freegrid spot
						delspot = np.where(freegrid==subgridswithresident[isubgrid][ires])[0]
						freegrid = np.delete(freegrid,delspot)
						tempfreegrid = np.delete(tempfreegrid,delspot)
											
						# Break from loop
						break
							
					# If offspring is not from the grid's subpopulation
					else:
							
						# Update the counter
						ioff = ioff + 1
	
	# Sum up the Residors
	Residors.append(len(OffDisperseIN))
	
	# --------------------
	# Straying
	# --------------------
	
	StrayTup = GetStrayers(tempoffspring,straying,nosubpops,cdmatrix,dispmovethresh,\
	dispmoveno,source_subpop,OffDispMortality)
	tempoffspring = StrayTup[0]
	OffDispMortality = StrayTup[1]
					
	# --------------------
	# Migration
	# --------------------
	
	# Loop through grid spots: looping through shuffled freegrid
	while len(freegrid) != 0 and len(tempoffspring) != 0:
		
		# But grab a random freegrid 
		randgrid = random.sample(freegrid,1)[0]
				
		# Index for not stepping out of bounds for offspring loop
		offcount = 0
		
		# Loop through shuffled offspring pool
		while offcount < len(tempoffspring):
						
			# If offspring is from the grid's subpopulation: here we index into this spot for subpop
			if tempoffspring[offcount][4].split('from')[0] == subpop[randgrid]:
			
				# flip coin for probability of returning home or not
				randdisperse = random.random()
				
				# --------------------
				# RETURN
				# --------------------				
				# If coin flip <= probability: RETURN OFFSPRING HOME -- fill up the grid spot
				if randdisperse < tempoffspring[offcount][6]:
					
					# If selection is on
					if cdevolveans == '1' and gen >= burningen:
					
						# Call 1-locus selection model
						differentialmortality = Do1LocusSelection(tempoffspring[offcount],fitvals1,randgrid)
												
						# Then flip the coin to see if offspring survives its location
						randcheck = rand()
						
						# If offspring did not survive: break from loop, move to next offspring
						if randcheck < differentialmortality:
							SelectionDeaths[gen].append(1)
							# Delete this offspring from pool
							if offcount == len(tempoffspring):
								pdb.set_trace()
							del(tempoffspring[offcount])
							offcount = offcount + 1
							continue
							
						# If offspring did survive:
						else:
						
							# Assign offspring to that freegrid spot: case for Migrator
							if len(tempoffspring[offcount][4].split('from')) == 1:
								recd = [tempoffspring[offcount],randgrid,'Migrator_'+'gen'+str(gen)+\
								'pop'+tempoffspring[offcount][4]+'F'+str(tempoffspring[offcount][0])+\
								'M'+str(tempoffspring[offcount][1])]
								
								# Append to grid fills IN
								OffDisperseIN.append(recd)
								
								# Store this offspring as a Migrator
								OffDispMigrate.append(tempoffspring[offcount])
								
							# Assign offspring to that freegrid spot: case for Strayer
							elif len(tempoffspring[offcount][4].split('from')) == 2:
								recd = [tempoffspring[offcount],randgrid,'Strayer_'+'gen'\
								+str(gen)+'pop'+tempoffspring[offcount][4].partition('from')[2]+'F'\
								+str(tempoffspring[offcount][0])+'M'+str(tempoffspring[offcount][1])]
								
								# Append to grid fills IN
								OffDisperseIN.append(recd)
								
								# Store this offspring as a Strayer
								OffDispStray.append(tempoffspring[offcount])
							
							# Something wrong
							else:
								print('Something wrong with strayers and migrators...email Erin.')
								sys.exit(-1)
																						
							# Delete this freegrid location from outer loop and from tempgrid
							freegrid = np.delete(freegrid,np.where(freegrid==randgrid)[0])
							tempfreegrid = np.delete(tempfreegrid,np.where(tempfreegrid==randgrid)[0])
												
							# Delete this offspring from pool
							del(tempoffspring[offcount])
							
							# Break from loop
							break
							
					# If selection is off
					else:
					
						# Assign offspring to that freegrid spot: case for Migrator
						if len(tempoffspring[offcount][4].split('from')) == 1:
							recd = [tempoffspring[offcount],randgrid,'Migrator_'+'gen'+str(gen)+\
							'pop'+tempoffspring[offcount][4]+'F'+str(tempoffspring[offcount][0])+\
							'M'+str(tempoffspring[offcount][1])]
							
							# Append to grid fills IN
							OffDisperseIN.append(recd)
							
							# Store this offspring as a Migrator
							OffDispMigrate.append(tempoffspring[offcount])
							
						# Assign offspring to that freegrid spot: case for Strayer
						elif len(tempoffspring[offcount][4].split('from')) == 2:
							recd = [tempoffspring[offcount],randgrid,'Strayer_'+'gen'\
							+str(gen)+'pop'+tempoffspring[offcount][4].partition('from')[2]+'F'\
							+str(tempoffspring[offcount][0])+'M'+str(tempoffspring[offcount][1])]
							
							# Append to grid fills IN
							OffDisperseIN.append(recd)
							
							# Store this offspring as a Strayer
							OffDispStray.append(tempoffspring[offcount])
						
						# Something wrong
						else:
							print('Something wrong with strayers and migrators...email Erin.')
							sys.exit(-1)
																					
						# Delete this freegrid location from outer loop and from tempgrid
						freegrid = np.delete(freegrid,np.where(freegrid==randgrid)[0])
						tempfreegrid = np.delete(tempfreegrid,np.where(tempfreegrid==randgrid)[0])
											
						# Delete this offspring from pool
						del(tempoffspring[offcount])
						
						# Break from loop
						break
				
				# ----------------
				# NOT RETURN
				# ----------------
				# Else: DO NOT RETURN OFFSPRING HOME -- keep track of that offspring that did not disperse
				else:
										
					# Keep track of offspring that do not disperse home: assume mortality
					OffDispMortality.append(tempoffspring[offcount])
					
					# Delete this offspring from pool
					del(tempoffspring[offcount])
					
					# Break from offspring loop, start over with another grid
					break
															
			# If offspring is not from the grid's subpopulation
			else:
					
				# Update the counter
				offcount = offcount + 1
		
		# If there were no offspring that belonged to that subpopulation randgrid
		# 	Delete this freegrid location from outer loop
		freegrid = np.delete(freegrid,np.where(freegrid==randgrid)[0])
			
	# Grid spots that did not fill up
	freegridisolated = tempfreegrid
		
	# Storage numbers
	SelectionDeaths[gen] = sum(SelectionDeaths[gen])
	Migrants.append(len(OffDispMigrate))
	OpenLocations.append(len(tempfreegrid))
	Strayer.append(len(OffDispStray))
	DispDeaths.append(len(OffDispMortality))
	OffLeftOver.append(len(tempoffspring))
							
	# Return variables from this argument
	tupDoDisp = OffDisperseIN,freegridisolated
	return tupDoDisp
	
	# End::DoDisperse()