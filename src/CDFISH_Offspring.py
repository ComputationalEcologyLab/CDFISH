# -------------------------------------------------------------------------------------------------
# CDFISH_Offspring.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for offspring processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
import pdb

# ---------------------------------------------------------------------------------------------------	
def DoOffspringSex(Bearpairs,Femalepercent,Births):
	'''
	DoOffspringSex()
	This function assigns the sex of each offspring.
	'''	
	
	# Create empty variable for storing individual offspring information
	offspring=[]
	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
		
		# And then loop through each offspring from that mate pair
		for j in xrange(Bearpairs[i][2]):
			
			# Select sex of the jth offspring - select a random number
			randsex = int(100*rand())
			
			# If that random number is less the Femalepercent, assign it to be a female
			if randsex < Femalepercent:
				offsex = 0
			
			# If the random number is greater than the Femalepercent, assign it to be a male
			else:
				offsex = 1

			# And then append all information onto a list storage variable offspring [mother,father,sex]
			offspring.append([Bearpairs[i][0],Bearpairs[i][1],offsex])
	
	# Store number of Births
	Births.append(len(offspring))
	
	# Variables returned
	tupDoOffsex = offspring,Births
	return tupDoOffsex
	
	# End::DoOffspringSex()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringRandom(Bearpairs,lmbda):
	'''
	DoOffspringRandom()
	This function chooses a random number of 
	offspring for a mated pair.
	'''	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
		
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
		
		# If females did mate up, then assign random drawn number
		else:
			
			# Randomly choose a number between 0 and 4
			randkidno = int(int((lmbda))*rand())
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(randkidno)
	
	# Variables returned
	return Bearpairs
	# End::DoOffspringRandom()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringPoisson(Bearpairs,lmbda):
	'''
	DoOffspringPoisson()
	This function chooses a number of offspring 
	from a Poisson distribution for a mated pair.
	'''		
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
		
		# If female did mate up, then assign a poisson draw with mean lmbda
		else:
			# Poisson number with mean, lambda
			poissonkidno = poisson(float(lmbda))
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(poissonkidno)
			
	# Variables returned
	return Bearpairs
	# End::DoOffspringPoisson()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringConstant(Bearpairs,lmbda):
	'''
	DoOffspringConstant()
	This function chooses a constant number of 
	offspring for each mated pair.
	'''	
	# Loop through each mate pair
	for i in xrange(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
		
		# If females did mate up, then assign lmbda constant
		else:
			# Assign a 1 [F,M,#offspring]
			Bearpairs[i].append(int(lmbda))
			
	# Variables returned
	return Bearpairs
	# End::DoOffspringConstant()

# ---------------------------------------------------------------------------------------------------	 
def DoOffspring(offno,lmbda,Bearpairs,Femalepercent,Births):
	'''
	DoOffspring()
	Choose number of Offspring for each mated pair and assign sex.
	Input: selection choice for offspring number distrubution draw
	offno, Bearpairs, lmbda.
	Output: Bear Pairs + # of offspring [Female,Male,#Offspring]
	'''
	
	# Function 1 is a uniform random draw between 0 and lmdba number	
	if (offno=='1'):
		
		Bearpairs = DoOffspringRandom(Bearpairs,lmbda)		
			
	# Function 2 is a Poisson draw
	elif (offno=='2'):
	
		Bearpairs = DoOffspringPoisson(Bearpairs,lmbda)		
				
	# Function 3 is a constant of lmbda offspring per each pairing
	elif (offno=='3'):
	
		Bearpairs = DoOffspringConstant(Bearpairs,lmbda)		
		
	# Assign sex of each offspring
	tupOffsex = DoOffspringSex(Bearpairs,Femalepercent,Births)	
	# Unpack tuple
	offspring = tupOffsex[0]
	Births = tupOffsex[1]
	
	# Return variables from this argument
	tupDoOff = offspring,len(offspring)
	return tupDoOff
		
	# End::DoOffspring()