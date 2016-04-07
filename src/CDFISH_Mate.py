# -------------------------------------------------------------------------------------------------
# CDFISH_Mate.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for mate processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# Python specific functions
import pdb, random, sys

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
def DoRandomMate1(nomales,females,nofemales,males,Bearpairs,subpop,count,subtally):	
	'''
	DoRandomMate1()
	This function is the random mating function for
	sexual reproduction
	females	without replacement
	males with replacement
	'''	
	
	# Check here in case there are no males in this females subpop
	subpopF = int(subpop[females[count]])
	
	# Loop through males
	for i in xrange(len(males)):
		#subpopulation of the male
		subpopM = int(subpop[males[i]])
		# Check if in same subpop
		if subpopF == subpopM:
			check = 1
			break
		else:
			check = 0
	
	# If there are no males...
	if len(males) == 0:
		check = 0
		print('There are no males.')
	
	# To ensure that mates are within the same population
	while check:
	
		# Randomly choose male
		intmale = int((nomales)*rand())
		
		# Get male's subpopulation
		subpopM = int(subpop[males[intmale]])
		
		# Make sure it is in the same subpopulations
		if subpopF == subpopM:
		
			# Store tally numbers
			subtally[int(subpop[females[count]])-1].append(1)
			
			# Append grid mates to relavant arrays
			Bearpairs.append([females[count],males[intmale]])
			break
			
		# If it is not in the same subpopulation
		else:
			continue
	
	# Return Variables from this function
	return Bearpairs,subtally
	
	#End::DoRandomMate1()
	
# ---------------------------------------------------------------------------------------------------		
def DoRandomMate2(nomales,females,nofemales,males,Bearpairs,subpop,count,\
subpoptotals,nosubpops,sex,subtally):	
	'''
	DoRandomMate2()
	This function is the random mating function for
	sexual reproduction
	females	with replacement
	males with replacement
	'''	
		
	# To ensure mates are within the same populations
	while 1:
	
		# Randomly choose male and female grid, but equal from each subpopulation
		intmale = int((nomales)*rand())
		intfemale = int((nofemales)*rand())
		
		# If mates are in the same subpopulations and if not more than subgrid total
		if subpop[females[intfemale]] == subpop[males[intmale]] and \
		sum(subtally[int(subpop[females[intfemale]])-1]) < subpoptotals[int(subpop[females[intfemale]])-1]:
			
			# Store tally numbers
			subtally[int(subpop[females[intfemale]])-1].append(1)
			
			# Append grid mates to relavant arrays
			Bearpairs.append([females[intfemale],males[intmale]])
			break

		# If not in the same subpopulation
		else:
			continue
							
	# Return Variables from this function
	return Bearpairs,subtally
	
	#End::DoRandomMate2()
	
# ---------------------------------------------------------------------------------------------------	 
def DoMate(nogrids,sex,age,freplace,mreplace,matemoveno,matemovethresh,\
xgridcopy,ygridcopy,ToTMales,ToTFemales,Population,subpop,nosubpops,\
subpoptotals,gen,reproage):
	'''
	DoMate()
	This is the mating function for choosing
	individual mate pairs.
	'''
	# ---------------------------------------------------
	# Select males and females for mating
	# ---------------------------------------------------
	
	# Storage variables males and females
	females = []		# These are the zeros
	males = []			# These are the ones
	
	# Storage variables for total numbers of males and females
	allfemales = []
	allmales = []

	# Loop through and grab each index for each catagory
	for i in xrange(nogrids):
		if sex[i] == '0':
			allfemales.append(i)
			if age[i] >= reproage:
				females.append(i)
		elif sex[i] == '1':
			allmales.append(i)
			if age[i] >= reproage:
				males.append(i)
		
	# Then get the length of each sex that are reproducing
	nomales = len(males)
	nofemales = len(females)

	# And grab this information for storing in output.csv file
	ToTMales.append([nomales])
	ToTFemales.append([nofemales])
	Population.append([nomales+nofemales])
	
	# Tally up subpopulation random grabs
	unique_subpops = np.unique(subpop)
	subtally = []
	for i in xrange(nosubpops):
		subtally.append([])
		
	# Loop through and sum up females and males within each subpop and store
	subpopfemales = []
	subpopmales = []
	# Create list of lists storage spots for number of subgrids
	for i in xrange(nosubpops):
		subpopfemales.append([])
		subpopmales.append([])
	
	# Store the number of females and males within each subpop
	for i in xrange(len(subpop)):
		# Loop through unique subpops
		for j in xrange(nosubpops):
			# If subpop exits append to subgrid spot
			if subpop[i] == unique_subpops[j] and sex[i] == '0':
				subpopfemales[int(unique_subpops[j])-1].append(1)
			if subpop[i] == unique_subpops[j] and sex[i] == '1':
				subpopmales[int(unique_subpops[j])-1].append(1)
	# Then sum them
	for i in xrange(nosubpops):
		ToTFemales[gen].append(sum(subpopfemales[i]))
		ToTMales[gen].append(sum(subpopmales[i]))
		Population[gen].append(sum(subpopfemales[i])+sum(subpopmales[i]))
	
	# Choose mate for each female or individual
	Bearpairs = []	# Empty matrix: xy indexes
	
	# Shuffle the females
	shuffle(females)
	
	# If there were no reproducing males or females
	if nomales == 0 or nofemales == 0:
		Bearpairs.append([-9999,-9999])
		
	# If there were reproducing males and females
	if nomales != 0 and nofemales != 0:
	
		# For the case of a Female without replacement and a male with replacement
		if freplace == 'N' and mreplace == 'Y':
						
			# Loop through while loop until all females paired up.		
			count = 0		# Initialize the while loop
			while count < nofemales:
							
				# Linear movement
				if (matemoveno=='1'):
					print('Linear mating not coded into CDFISH yet.')
					sys.exit(-1)
					
				# Inverse-square movement
				elif (matemoveno=='2'):
					print('Inverse-square mating not coded into CDFISH yet.')
					sys.exit(-1)
				
				# Nearest neighbor movement
				elif (matemoveno=='3'):
					print('Nearest neighbor mating not coded into CDFISH yet.')
					sys.exit(-1)
								
				# Random mixing case is outside of other cases, this case does not use the cdmatrix.
				elif (matemoveno=='4'):
					
					tupBears = DoRandomMate1(nomales,females,nofemales,males,\
					Bearpairs,subpop,count,subtally)	
					Bearpairs = tupBears[0]
					subtally = tupBears[1]
												
				# Negative exponential movement
				elif (matemoveno=='5'):
					print('Negative exponential mating not coded into CDFISH yet.')
					sys.exit(-1)	
					
				# Update count
				count = count + 1
			
		# For the case of a Female with replacement and a male with replacement
		elif freplace == 'Y' and mreplace == 'Y':

			# Loop through while loop until all females paired up, but do this nogrid times.		
			count = 0		# Initialize the while loop
			while count < nogrids:
				
				# Linear movement
				if (matemoveno=='1'):
					print('Linear mating not coded into CDFISH yet.')
					sys.exit(-1)
					
				# Inverse-square movement
				elif (matemoveno=='2'):
					print('Inverse-square mating not coded into CDFISH yet.')
					sys.exit(-1)
				
				# Nearest neighbor movement
				elif (matemoveno=='3'):
					print('Nearest neighbor mating not coded into CDFISH yet.')
					sys.exit(-1)
				
				# Random mixing case is outside of other cases, this case does not use the cdmatrix.
				elif (matemoveno=='4'):
					
					tupBears = DoRandomMate2(nomales,females,nofemales,males,\
					Bearpairs,subpop,count,subpoptotals,nosubpops,sex,subtally)	
					Bearpairs = tupBears[0]
					subtally = tupBears[1]
					
				# Negative exponential movement
				elif (matemoveno=='5'):
					print('Negative exponential mating not coded into CDFISH yet.')
					sys.exit(-1)
				
				# Update count
				count = count + 1
					
		# For the case of Female with replacement and male without replacement
		elif freplace == 'Y' and mreplace == 'N':
		
			print('Female with replacement and Male without replacement not coded yet.')
			sys.exit(-1)
			
		# For the case of Female without replacement and male without replacement
		else:
		
			print('Female without replacement and Male without replacement not coded yet.')
			sys.exit(-1)
			
	# ----------------------------------------
	# Summary Stats on Mate functions
	# ----------------------------------------
	
	# Return variables from this function
	return Bearpairs
	
	#End::DoMate()