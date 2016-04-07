# -------------------------------------------------------------------------------------------------
# CDFISH_Modules.py
# Author: Erin L Landguth
# Created: June 2010
# Description: This is the function/module file for CDFISH v0.88
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."

# CDFISH functions
try:
	from CDFISH_PreProcess import *
except ImportError:
	raise ImportError, "CDFISH Modules required."
	
# Python specific functions
import os, random, copy, pdb, sys
from collections import Counter

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# --------------------------------------------------------------------------
def PrepTextFile(textpath):
	'''
	PrepTextFile() - Prepare the input files
	'''
	
	return textpath.strip('\n').strip('\r')
	
	# End::PrepTextFile()

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
def ReadGrid(FIDnew,idnew,agenew,xgridnew,ygridnew,genesnew,sexnew,subpopnew,\
gen,equalsexratio,nosubpops):	
	
	'''
	DoReadGrid()
	This function reads the previous generations
	grid information at the start of the time loop.  This function 
	also equalizes the sex of the new incoming individuals if chosen.
	Input: Previous individual's information.
	Output: Renamed individual's information to begin loop.
	'''
	
	# Reassign old information to new information.
	FID = FIDnew
	subpop = subpopnew
	xgrid = xgridnew
	xgridcopy = copy.deepcopy(xgridnew)
	ygrid = ygridnew
	ygridcopy = copy.deepcopy(ygridnew)
	id = idnew
	age = agenew	
	genes = genesnew
				
	# Get unique subpop names
	unique_subpops = np.unique(subpop)
	
	# Count up the unique number of subgrids appending to subgrids
	subtotal = []
	# Create list of lists storage spots for number of subgrids
	for i in xrange(nosubpops):
		subtotal.append([])
	for i in xrange(len(subpop)):
		# Loop through unique subpops
		for j in xrange(nosubpops):
			# If subpop exits append to subgrid spot
			if subpop[i] == unique_subpops[j] and genes[i][0][0] != 'NA':
				subtotal[int(unique_subpops[j])-1].append(1)
				
	# And then sum them up
	for i in xrange(nosubpops):
		subtotal[i] = sum(subtotal[i])
	
	# If sex is not equal ratio
	if equalsexratio == 'N':
		sex = sexnew
	
	# If equal sex ratio is Y, then split up sex into equal parts in each subpopulation
	if equalsexratio == 'Y':
	
		# Some temp storage variables for sex
		sexequal = []
				
		# Create list of lists storage spots for number of subpopulations
		for i in xrange(nosubpops):
			sexequal.append([])
						
		# Want to split equal males and females per subgrid totals
		for i in xrange(nosubpops):
			# Loop over subpop totals to have half 1s and half 0s
			for j in xrange(int(np.floor(subtotal[i]/2.))):
				sexequal[i].append(0)
				sexequal[i].append(1)
			# If nogrids just happens to be odd (doesn't make sense with equal sex ratio), but whatever...grab a random draw
			if np.mod(subtotal[i],2) == 1:
				sexequal[i].append(int(2*rand()))	
		
		# Then flatten sexequal to sex string values
		sex = []
		for i in xrange(nosubpops):
			for j in xrange(subtotal[i]):
				sex.append(str(sexequal[i][j]))
		
		# Now loop through xy values appending new equal sex value to spot unless it is a 'NA'
		for i in xrange(len(subpop)):
			if age[i] == 'NA':
				sex.insert(i,'NA')
									
	# Store the number of grids
	nogrids = len(FID)
	
	# Store the number of filled grids
	tempfilledgrids = []
	for i in xrange(nogrids):
		if age[i] != 'NA':
			tempfilledgrids.append(1)
	filledgrids = sum(tempfilledgrids)
	
	#Return this functions variables
	tupReadGrid = FID,sex,id,age,xgrid,xgridcopy,ygrid,ygridcopy,\
	genes,nogrids,subpop,filledgrids,subtotal
	return tupReadGrid
	
	#End::DoReadGrid()

# ---------------------------------------------------------------------------------------------------	 
def GetMetrics(newnogrids,loci,alleles,genes,gen,Ho,Alleles,He,subpop,subgridtotal,\
p1,q1,FST=None,F=None,FIS=None,FIT=None):
	'''
	GetMetrics()
	This function summarizes the genotypes and
	produces a suite of genetic metrics.
	Ho - Observed heterozygosity per generation
	He - Expected heterozygoisty per generation
	Alleles - Total number of unique alleles in genotype*individuals
	Right and Left populations should be user-specified eventually. 
	F, FST, FIT, FIS values as well for entire population per generation.	
	'''
	
	# List for total, left, and right
	unique_alleles = Alleles
	
	# Get allele location as seqence from alleles array
	allele_numbers = []
	for i in xrange(loci):
		for j in xrange(alleles[i]):
			allele_numbers.append(j)
	allele_numbers = np.asarray(allele_numbers)
	
	# Remove the 'NA' gene values
	tempgenes = []
	for i in xrange(len(genes)):
		if genes[i][0][0] != 'NA':
			tempgenes.append(genes[i])
	
	# Cast genes as an numpy array as byte type
	genes_array_woNA = np.asarray(tempgenes,dtype='float')
	genes_array_wNA = np.asarray(genes)
	
	# The total number of alleles
	total_alleles = len(allele_numbers)
			
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	
	# Remove a unique number from above if 'NA' value exists
	for i in xrange(nosubpops):
		if np.unique(subpop)[i] == 'NA':
			nosubpops = nosubpops -1
			
	# Create a list to store the subpopulation grid number information
	subgrids = []
	all_freq_sub = []
	ho_count_sub = []
	ho_sub = []
	all_freq_sq_sub = []
	homozygosity_sub = []
	he_sub = []
	sumsubpopsHo = []
	sumsubpopsHe = []
	alleles_sub = []
	
	# Create list of lists storage spots for number of subgrids
	for i in xrange(nosubpops):
		subgrids.append([])
		all_freq_sub.append([])
		ho_count_sub.append([])
		ho_sub.append([])
		all_freq_sq_sub.append([])
		homozygosity_sub.append([])
		he_sub.append([])
		alleles_sub.append([])
		
	# Count up the unique number of subgrids appending to subgrids
	for i in xrange(len(subpop)):
		# Loop through unique subpops
		for j in xrange(nosubpops):
			# If subpop exits append to subgrid spot
			if subpop[i] == unique_subpops[j] and genes[i][0][0] != 'NA':
				subgrids[int(unique_subpops[j])-1].append(i)
	
	# Get allele frequency for total
	if newnogrids != 0:
		all_freq_tot = np.asarray(np.sum(genes_array_woNA,axis=0),dtype = 'float').reshape(total_alleles)
		all_freq_tot = all_freq_tot/(2*newnogrids)
	else:
		all_freq_tot = np.zeros(total_alleles,float)
		
	# Get the sqare of the allele frequency for total
	all_freq_sq_tot = all_freq_tot**2
	
	# Get allele frequency for subpopulations
	for i in xrange(nosubpops):
		if subgridtotal[i] != 0:
			all_freq_sub[i].append(np.asarray(np.sum(np.asarray(genes_array_wNA[subgrids[i],:,:],dtype='float'),axis=0),dtype = 'float').reshape(total_alleles))
			all_freq_sub[i] = all_freq_sub[i][0]/(2*subgridtotal[i])
		else:
			all_freq_sub[i].append(np.zeros(total_alleles,float))
			all_freq_sub[i] = all_freq_sub[i][0]
		# The square of all freq.
		all_freq_sq_sub[i].append(all_freq_sub[i]**2)
	
	# Calculate the number of homogenous alleles for total
	ho_count_tot = (np.array(genes_array_woNA==2)).sum()
	
	# Calculate the number of homogenous alleles in each subpop
	for i in xrange(nosubpops):
		ho_count_sub[i].append((np.array(np.asarray(genes_array_wNA[subgrids[i],:,:],dtype='float')==2)).sum())
	
	# Calculate the observed het for total
	if newnogrids != 0:
		ho_tot = (float(newnogrids*loci - ho_count_tot)/(loci*newnogrids))
	else:
		ho_tot = 0.0
		
	# Calculate the observed het in each subpop
	for i in xrange(nosubpops):
		if subgridtotal[i] != 0:
			ho_sub[i].append((float(subgridtotal[i]*loci - ho_count_sub[i][0])/(loci*subgridtotal[i])))
		else:
			ho_sub[i].append(0)
		
	# Append Ho information (Observed Het)
	Ho.append([ho_tot])
	for i in xrange(nosubpops):
		Ho[gen].append(ho_sub[i][0])
		
	# Calculate the homozygosity for total populations
	homozygosity_tot = sum(all_freq_sq_tot)/loci
	
	# Calculate the homozygosity for subpopulations
	for i in xrange(nosubpops):
		homozygosity_sub[i].append(sum(all_freq_sq_sub[i][0])/loci)
	
	# Store He for [Total]
	if newnogrids != 0:
		he_tot = (1. - homozygosity_tot)
	else:
		he_tot = 0.0
	
	# Store He for subpopulations
	for i in xrange(nosubpops):
		if subgridtotal[i] != 0:
			he_sub[i].append(1. - homozygosity_sub[i][0])
		else:
			he_sub[i].append(0.0)
	
	# Append He information (Expected Het)
	He.append([he_tot])
	for i in xrange(nosubpops):
		He[gen].append(he_sub[i][0])
		
	# Get total number of alleles
	alleles_tot = np.array(all_freq_tot>0.).sum()
	
	# Get the total number of alleles in each subpop
	for i in xrange(nosubpops):
		alleles_sub[i].append(np.array(all_freq_sub[i]>0.).sum())
	
	# Append allele total information
	unique_alleles.append([alleles_tot])
	for i in xrange(nosubpops):
		unique_alleles[gen].append(alleles_sub[i][0])
	
	# Get allele frequency totals for selection section
	p1.append(all_freq_tot[0])
	q1.append(all_freq_tot[1])
		
	# Get hetero indices
	for i in xrange(nosubpops):
		sumsubpopsHo.append(ho_sub[i][0]*subgridtotal[i])
		sumsubpopsHe.append(he_sub[i][0]*subgridtotal[i])
	HI = sum(sumsubpopsHo)/newnogrids
	HS = sum(sumsubpopsHe)/newnogrids
	
	# Return variables from this function
	tupGetMetrics = Ho,unique_alleles,He,all_freq_tot
	return tupGetMetrics
	
	#End::GetMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def InheritGenes(gen,AllelesMutated,offspringno,offspring,genes,loci,muterate):
	'''
	InheritGenes()
	Pass along gentic information to survived offspring from parents
	following Mendal inheritance.  Future work here involves recombination. 
	The KAM mutatiaon model applied.  Future work includes more mutation
	models.
	Input: offspring, genes 
	Output: [femaleid,maleid,cdmatidofmother,cdmatidoffather,sex],[
	genetic information]		
	'''		
		
	# Storage for tracking how many alleles mutated
	noallelesmutated = []
	
	# If there are offspring
	if int(offspringno) != int(0):
	
		# Begin loop through offspring
		for i in xrange(offspringno):
			
			# Add spot in offspring array for individual i's genes
			offspring[i].append([])
			
			# Temp storage for i's mother's genes
			mothergenes=genes[offspring[i][0]]
			# Temp storage for i's father's genes
			fathergenes=genes[offspring[i][1]]
							
			# Loop through each locus
			for jspot in xrange(loci):
				
				# Temporary index storage
				tempindfather = []
				tempindmother = []
									
				# Loop through each allele-mother and father have same len of alleles
				for kspot in xrange(len(fathergenes[jspot])):
											
					# Check the homogeneous 2 case - in father genes
					if int(fathergenes[jspot][kspot])==2:
						
						# Get a random number for allele mutation
						mutationrandno = rand()
						
						# Check if random number is less than or equal to muterate
						if mutationrandno <= muterate:
						
							# Randomly choose another allele
							randallelespot = int(len(fathergenes[jspot])*rand())
							tempindfather.append(randallelespot)
							
							# Count a mutation
							noallelesmutated.append(1)
						
						# and if random number is not less than or equal to muterate
						else:
							tempindfather.append(kspot)
					
					# Check the homogeneous 2 case - in mother genes	
					if int(mothergenes[jspot][kspot])==2:
						
						# Get a random number for allele mutation
						mutationrandno = rand()
						
						# Check if random number is less than or equal to muterate
						if mutationrandno <= muterate:
						
						# Randomly choose another allele
							randallelespot = int(len(mothergenes[jspot])*rand())
							tempindmother.append(randallelespot)

							# Count a mutation
							noallelesmutated.append(1)

						# and if random number is not less than or equal to muterate
						else:
							tempindmother.append(kspot)
						
					# Check the hetero 1 case for father genes
					if int(fathergenes[jspot][kspot])==1:
						tempindfather.append(kspot)
						
					# Check the hetero 1 case for mother genes
					if int(mothergenes[jspot][kspot])==1:
						tempindmother.append(kspot)
				
				# Check if the tempindex has a length of 2 (which means it was not homo at
				#	at this locus), then randomly select one of them, and then check for mutation
				# Check from father genes
				if len(tempindfather) == 2:
				
					# Then randomly select one of the homo alleles
					temprandnofather = int(2*rand())
					
					# Delete from list
					del(tempindfather[temprandnofather])
					
					# Get a random number for allele mutation
					mutationrandno = rand()
					
					# Check if random number is less than or equal to muterate
					if mutationrandno <= muterate:
						
						# Randomly choose another allele
						randallelespot = int(len(fathergenes[jspot])*rand())
						
						# and then reassign this spot
						tempindfather[0] = randallelespot
						
						# Count a mutation
						noallelesmutated.append(1)
				
				# Check from mother genes
				if len(tempindmother) == 2:
				
					# THen randomly select on of the homo alleles
					temprandnomother = int(2*rand())
					
					# Delete from list
					del(tempindmother[temprandnomother])
					
					# Get a random number for allele mutation
					mutationrandno = rand()
					
					# Check if random number is less than or equal to muterate
					if mutationrandno <= muterate:
					
						# Randomly choose another allele
						randallelespot = int(len(mothergenes[jspot])*rand())
						
						# and then reassign this spot
						tempindmother[0] = randallelespot
						
						# Count a mutation
						noallelesmutated.append(1)
				
				# Now write to offspring genes array the selected alleles in locus j
				for kspot in xrange(len(fathergenes[jspot])):
					
					# Hetero case 1 AB
					if tempindfather[0] == kspot and tempindmother[0] != kspot:
						offspring[i][3].append(1)
					# Homo case AA or BB
					elif tempindfather[0] == kspot and tempindmother[0] == kspot:
						offspring[i][3].append(2)
					# Hetero case 2 BA
					elif tempindmother[0] == kspot and tempindfather[0] != kspot:
						offspring[i][3].append(1)
					# Or nothing there at all
					elif tempindmother[0] != kspot and tempindfather[0] != kspot:
						offspring[i][3].append(0)
												
		# Now store the total number of alleles that mutated
		AllelesMutated.append(sum(noallelesmutated))
		
		# Delete temp variables to free up space
		del(mothergenes)
		del(fathergenes)
		del(tempindmother)
		del(tempindfather)
		
	# If there are no offspring
	elif int(offspringno) == int(0):
	
		# Store the total number of alleles that mutated
		AllelesMutated.append(0.0)				
	
	# Return variables from this argument
	return offspring
	
	# End::InheritGenes()
	
# ---------------------------------------------------------------------------------------------------	 
def DoAdultMortality(filledgrids,nogrids,sex,id,\
age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,agemort):
	'''
	DoAdultMortality()
	Mortality of old generation
	Input: Adult mortality% 
	Output: Old files minus the killed off individuals:
	freegrid = [xycdmatid location of free grid spot 
	in random order]		
	'''
		
	# Get total number of deaths for each age class, be careful of NAs and strings
	uniqueages = Counter(np.asarray(np.asarray(age,dtype='|S10')[np.where(np.asarray(age,dtype='|S10') != 'NA')[0]],dtype=np.int8))
	agedeaths = []
	# Swith for if max is 0 age
	if max(uniqueages) == 0:
		agedeaths.append(round(agemort[0]*uniqueages[int(0)]))
	# Swith for if max is greater than 0
	else:
		for i in xrange(max(uniqueages)+1):
			# Switch for if age class goes beyond agemort file given, then make last mortality.
			if i <= len(agemort)-1:
				agedeaths.append(round(agemort[i]*uniqueages[int(i)]))
			else:
				agedeaths.append(round(agemort[-1]*uniqueages[int(i)]))
	
	# Store total number of Deaths information
	Deaths.append(agedeaths)
						
	# Grab locations that are open
	openindex = np.where(np.asarray(sex) == 'NA')[0]
	
	# Grab locations that are not open
	filledindex = np.where(np.asarray(sex) != 'NA')[0]
	
	# Then take a sample from the possible age class indices to delete from
	deleteoldindex = []
	for i in xrange(len(agedeaths)):
		# If agemort is just of length 1 and less than 100s, then there can be age classes appearing
		#	in which case they will all get uniform agemort mortality.
		if agedeaths[i] > 0 and len(agedeaths) != 0:
			# NA switch
			if len(openindex) == 0:
				deleteoldindex.append(random.sample(np.where(np.asarray(age) == i)[0],int(agedeaths[i])))
			else:
				deleteoldindex.append(random.sample(np.where(np.asarray(age) == str(i))[0],int(agedeaths[i])))
		elif agedeaths[i] >0 and len(agedeaths) == 0:
			deleteoldindex.append(random.sample(filledindex,int(agedeaths[i])))
		
	# Flatten and turn into array
	deleteoldindex = np.asarray([item for sublist in deleteoldindex for item in sublist])
		
	# Then add indices together
	deleteallindex = np.append(openindex,deleteoldindex)
	
	# Switch for if anyone died or not.
	if len(deleteoldindex) != 0:
		# Store freegrid locations
		freegrid = np.asarray(FID)[deleteallindex]
		# Delete all of the old generations information.
		sex = np.delete(sex,deleteallindex)
		id = np.delete(id,deleteallindex)
		age = np.delete(age,deleteallindex)
		xgrid = np.delete(xgrid,deleteallindex)
		ygrid = np.delete(ygrid,deleteallindex)
		genes = np.delete(genes,deleteallindex,axis=0)
		FID = np.delete(FID,deleteallindex)
	else:
		# Store freegrid locations
		freegrid = []
		# Delete all of the old generations information.
		sex = np.asarray(sex)
		id = np.asarray(id)
		age = np.asarray(age)
		xgrid = np.asarray(xgrid)
		ygrid = np.asarray(ygrid)
		genes = np.asarray(genes)
		FID = np.asarray(FID)
		
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)
	freegrid = list(np.asarray(freegrid,dtype="int"))
	
	# Return variables from this argument
	tupAMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID
	return tupAMort
	
	
	'''
	# Find index of random mort number from total nogrids and Delete indexes
	# Store delete xycdmatid information: these are the free grid spots
	freegrid=[]
	oldmort = []
		
	# Get total number of deaths for each generation.
	oldmort = round(oldmortperc*filledgrids)
	
	# Store total number of Deaths information
	Deaths.append(oldmort)
			
	# Create a delete index temp variables
	deleteindex = []
	
	# Add any actual free grid locations to freegrid and then delete that information
	for i in xrange(nogrids):
	
		# Just an error check, make sure it is actually a free grid spot.
		if sex[i] == 'NA':	

			# Store these indices
			deleteindex.append(i)
			
			# Add this free grid location to freegrid
			freegrid.append(int(FID[i]))
			
	# Delete all of the old generations information.
	sex = np.delete(sex,deleteindex)
	id = np.delete(id,deleteindex)
	age = np.delete(age,deleteindex)
	xgrid = np.delete(xgrid,deleteindex)
	ygrid = np.delete(ygrid,deleteindex)
	if len(deleteindex) > 0:
		for idel in xrange(len(deleteindex)):
			genes = np.delete(genes,range(deleteindex[idel]*sum(alleles),deleteindex[idel]*sum(alleles)+sum(alleles)))
	FID = np.delete(FID,deleteindex)
	
	# Loop through the total adult kill off number
	for i in xrange(int(oldmort)):
	
		# - i from nogrids so we don't step out of index list
		oldrandindtemp = (filledgrids-i)*rand()
		
		# Store deleted 'grids' locations as FID
		freegrid.append(int(FID[int(oldrandindtemp)]))
		
		# Delete all of the old generations information.
		sex = np.delete(sex,int(oldrandindtemp))
		id = np.delete(id,int(oldrandindtemp))
		age = np.delete(age,int(oldrandindtemp))
		xgrid = np.delete(xgrid,int(oldrandindtemp))
		ygrid = np.delete(ygrid,int(oldrandindtemp))
		genes = np.delete(genes,range(int(oldrandindtemp)*sum(alleles),int(oldrandindtemp)*sum(alleles)+sum(alleles)))
		FID = np.delete(FID,int(oldrandindtemp))		
	
	# Return variables from this argument
	tupAMort = freegrid,Deaths,id,sex,age,xgrid,ygrid,genes,oldmort,FID
	return tupAMort
	'''
	# End::DoAdultMortality()