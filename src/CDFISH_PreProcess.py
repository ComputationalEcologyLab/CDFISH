# -------------------------------------------------------------------------------------------------
# CDFISH_PreProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file pre processing.
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError, "Numpy required."
	
# Scipy function KDTree
try:
	from scipy.spatial import KDTree
	scipyAvail = True
except ImportError:
	raise ImportError, "Scipy required."
	scipyAvail = False

# CDFISH functions
try:
	from CDFISH_PostProcess import *
	from CDFISH_Mate import *
	from CDFISH_Offspring import *
	from CDFISH_GetMetrics import *
	from CDFISH_Disperse import *
except ImportError:
	raise ImportError, "CDFISH Modules required."
	
# General python imports
import os,sys

# ---------------------------------------------------------------------------------------------------	 
def GetMaxCDValue(matemovethresh,cdmatrix,minmaxdispthresh):	
	'''
	GetMaxCDValue()
	This function calculates the maximum movement thresholds.
	'''	
	# mating movement threshold if max specified
	if str(matemovethresh).endswith('max'):
		# If max
		if len(matemovethresh.strip('max')) == 0:
			matemovethresh = float(max(max(cdmatrix)))
		else:
			matemovethresh = (int(matemovethresh.strip('max'))/100.)\
			*float(max(max(cdmatrix)))
	else:
		matemovethresh = float(matemovethresh)
	
	if len(minmaxdispthresh) == 1:
		# Dispersal movement threshold if max specified
		if str(minmaxdispthresh[0]).endswith('max'):
			# If max
			if len(minmaxdispthresh[0].strip('max')) == 0:
				minmaxdispthresh = [float(max(max(cdmatrix)))]
			else:
				minmaxdispthresh = [(int(minmaxdispthresh[0].strip('max'))/100.)\
				*float(max(max(cdmatrix)))]
		else:
			minmaxdispthresh = [float(minmaxdispthresh[0])]
	else:
		minmaxdispthresh = minmaxdispthresh
	tupGetMax = matemovethresh,minmaxdispthresh
	return tupGetMax
	# End::GetMaxCDValue()
	
# ---------------------------------------------------------------------------------------------------	 
def ReadFitnessSurface(fitsurface):	
	'''
	ReadFitnessSurface()
	This function reads in the ascii fitness surface.
	'''	
	# Open file for reading
	inputfile = open(fitsurface,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	values = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		i = i.strip('\n').strip('\r').strip(' ')
		thisline = i.split(' ')
		values.append(thisline)
	
	# Grab some information from fitvalues: number of columns
	lenfirstline = len(values[0])
	ncols = int(values[0][lenfirstline-1])
		
	# Grab some information from values: number of rows
	lensecondline = len(values[1])
	nrows = int(values[1][lensecondline-1])

	# Grab some information from values: x coord value in lower left
	lenthirdline = len(values[2])
	xllcorner = float(values[2][lenthirdline-1])

	# Grab some information from values: y coord value in lower left
	lenfourthline = len(values[3])
	yllcorner = float(values[3][lenfourthline-1])

	# Grab some information from values: cell size
	lenfifthline =len(values[4])
	cellsize = float(values[4][lenfifthline-1])
	
	# Grab some information from values: Nodataval
	lensixthline =len(values[5])
	Nodataval = float(values[5][lensixthline-1])
	
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)
	
	# Error statement
	if len(values)-6 != nrows or len(values[6]) != ncols:
		print('Spatial selection surface file is not in correct format. Return error number of dimensions.')
		sys.exit(-1)
	# And turn rasterfile into a list with float values without header information
	fitnessvalues = []
	for i in xrange(len(values)-6):
		fitnessvalues.append([])
		for j in xrange(len(values[6])):
			fitnessvalues[i].append(float(values[6+i][j]))
	
	# Return tuple
	tupReadFit = fitnessvalues,ncols,nrows,xll,yll,cellsize,Nodataval
	return tupReadFit
	# End::ReadFitnessSurface()

# ---------------------------------------------------------------------------------------------------	
def CreateAlleleList(loci,alleles,xgenes):
	'''
	CreateAlleleList()
	This function creates a list for the allele storage.
	'''		
	# Store all information in a list [loci][allele#,probability]
	allelst = []
	for i in xrange(loci):
		allelst.append([])
		for k in xrange(alleles[i]):
			allelst[i].append([int(k),float(xgenes[alleles[k]*i+1+k][1])])
	
	# Return variables
	return allelst
	
	# End::CreateAlleleList()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeGenes(intgenesans,allefreqfilename,loci,alleles):	
	'''
	InitializeGenes() -- This function initializes the genetic structure 
	of the population from an allele frequency file.
	Input: genetic information
	Output: 0
	'''	
		
	# If genetic structure intialized by a file...
	if intgenesans == 'file' and allefreqfilename != 'N':
		
		# Open file for reading
		inputfile = open(allefreqfilename,'r')
		
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		xgenes = []
		
		# Split up each line in file and append to empty matrix, x
		for i in lines:
			thisline = i.split(',')
			xgenes.append(thisline)
			
		# Error check here
		if (len(xgenes)-1) != sum(alleles):
			print('Allele frequency file is not the specified number of loci and alleles as in in input file.')
			sys.exit(-1)
			
		# Delete lines from earlier
		del(lines)
		
		# Call CreateAlleleList()
		allelst = CreateAlleleList(loci,alleles,xgenes)
			
		# Delete x variable
		del(xgenes)
		
		# Return variables
		return allelst
		
		# End::InitializeGenes()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeAge(agefilename,nogrids,agemort):
	'''
	InitializeAge()
	This function initializes the age of each population
	with an age distribution list.
	'''

	# Only run this section if agefilename was specified, else skip if == 'N'
	if agefilename != 'N':
	
		# Open file for reading
		inputfile = open(agefilename,'r')
		
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		xage = []
				
		# Split up each line in file and append to empty matrix, x
		for i in lines:
			thisline = i.split(',')
			xage.append(thisline)
			
		# Delete lines from earlier
		del(lines)
		
		# Store all information in a list 
		agelst = [] # [age,probability] for age distribution
		agemort = []
		Mprob_mature = []
		Mprob_migrate = []
		Mprob_juvenile = []
		Fprob_mature = []
		Fprob_migrate = []
		Fprob_juvenile = []
		temp = 0
		for i in xrange(len(xage)-1):
			agelst.append([int(i),float(xage[i+1][1])/nogrids])
			temp = temp + float(xage[i+1][1])
			agemort.append(float(xage[i+1][2])/100.)
			Mprob_mature.append(float(xage[i+1][3]))
			Mprob_migrate.append(float(xage[i+1][4]))
			Mprob_juvenile.append(float(xage[i+1][5]))
			Fprob_mature.append(float(xage[i+1][6]))
			Fprob_migrate.append(float(xage[i+1][7]))
			Fprob_juvenile.append(float(xage[i+1][8]))
	
		# Error checks here: if distribution does not equal one
		if temp != float(nogrids):			
			print('Number in age class must total the carrying capacity.')
			sys.exit(-1)
		# Error checks here: if number of classes does not equal mortality age classes
		if len(agemort) != len(agelst):
			print('Number in age class must equal the total age structure mortality list.')
			sys.exit(-1)
		# Error probabilties in age list
		if len(Mprob_mature) != len(Mprob_migrate) != len(Mprob_juvenile)\
		!= len(Fprob_mature) != len(Fprob_migrate) != len(Fprob_juvenile):
			print('Age distribution file in the wrong format.')
			sys.exit(-1)
		
		# Deletes
		del(xage)
		del(temp)
	
	# if agedistribution == 'N'
	elif agefilename == 'N':
		
		# Create an aglst of zeros
		agelst = [[0,1.0]]
		Mprob_mature = [1.0]
		Mprob_migrate = [1.0]
		Mprob_juvenile = [0.0]
		Fprob_mature = [1.0]
		Fprob_migrate = [1.0]
		Fprob_juvenile = [0.0]
	
	# Return variables
	tupAgeFile = agelst, Mprob_mature, Mprob_migrate, Mprob_juvenile, \
	Fprob_mature, Fprob_migrate, Fprob_juvenile, agemort
	return tupAgeFile
	
	# End::InitializeAge()
	
# ---------------------------------------------------------------------------------------------------	 
def ReadCDMatrix(cdmatrixfilename):
	'''
	ReadCDMatrix()
	This function reads in the cost distance matrix.
	'''	
	
	# Check statements
	if os.path.exists(cdmatrixfilename):
		# Open file for reading
		inputfile = open(cdmatrixfilename,'rU')
	else:
		print("CDFISH ReadCDMatrix() error: open failed, could not open %s"%(cdmatrixfilename))
		sys.exit(-1)
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	# Close the file
	inputfile.close()
	
	# Create an empty matrix to append to 
	bigCD = []
	
	# Split up each line in file and append to empty matrix, x
	for spot in lines:
		thisline = spot.split(',')
		bigCD.append(thisline)
	
	# Delete lines from earlier
	del(lines)
		
	# Store number of files
	nofiles = len(bigCD)
	
	# Create a matrix of to be filled 	
	cdmatrix = []
	
	# Fill up matrix with float value of array x
	for j in xrange(nofiles):
		cdmatrix.append([])
		for k in xrange(nofiles):
			cdmatrix[j].append(float(bigCD[j][k]))
			
	# Delete x variable
	del(bigCD)
	
	# Return variables
	tupReadMat = cdmatrix,nofiles
	return tupReadMat
	
	# End::ReadMateCDMatrix

# ---------------------------------------------------------------------------------------------------	 
def ReadXY(xyfilename):
	'''
	ReadXY()
	This function reads in the xy individual point locations.
	Input: Directory location and XY filename
	Output: xy list
	'''		
	
	# Check statements
	if os.path.exists(xyfilename):
		# Open file for reading
		inputfile = open(xyfilename,'r')
	else:
		print("CDFISH ReadXY() error: open failed, could not open %s"%(xyfilename))
		sys.exit(-1)
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	#Close the file
	inputfile.close()
	
	# Create an empty matrix to append to
	xy = []
	
	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.split(',')
		xy.append(thisline)
		
	# Delete lines from earlier
	del(lines)
		
	# Return variables
	return xy
	
	# End::ReadXY()

# ---------------------------------------------------------------------------- 
def DoGridOut_cdpop0(ithmcrundir,gen,loci,alleles,nogrids,subpop,xgrid,ygrid,\
id,sex,age,agelst,genes,allelst,intgenesans):
	'''
	DoGridOut_cdpop0()
	Output grid0.csv in cdpop format.  This is the initial genotypes and 
	individual information stored to file.
	'''		
	
	# Create file to write matrix to
	outputfile = open(ithmcrundir+'grid'+str(0)+'.csv','w')
		
	# Write out the titles
	title = ['Subpopulation','XCOORD','YCOORD','ID','sex','age']
	
	# Write out the title from xy points
	for i in xrange(len(title)):
		outputfile.write(title[i]+',')
	
	# Write out the genes title informations
	# Loop through loci
	for i in xrange(loci-1):		
		
		# Loop for allele length
		for j in xrange(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	
	# To get a return character on the end of the title
	for i in xrange(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1))
	# Get return character
	outputfile.write('\n')
	
	# Write out all of the information.		
	for i in xrange(nogrids):
		outputfile.write(subpop[i]+',')
		outputfile.write(str(float(xgrid[i]))+',')
		outputfile.write(str(float(ygrid[i]))+',')
		if sex[i] == 'NA':
			outputfile.write('OPEN,')
			outputfile.write('NA,')
			outputfile.write('NA,')
			age.append('NA')
		else:
			outputfile.write(id[i]+',')
			outputfile.write(str(sex[i])+',')
			agetemp = w_choice_general(agelst)[0]
			outputfile.write(str(agetemp)+',')
			age.append(agetemp)
		
		# if known genes
		if intgenesans == 'known':
			
			# Write out gene info
			for jk in xrange(loci-1):
				for kl in xrange(alleles[jk]):
					if sex[i] == 'NA':
						outputfile.write('NA,')
					else:
						outputfile.write(str(int(genes[i][jk][kl]))+',')
			# To get return character on end
			for jk in xrange(alleles[loci-1]-1):
				if sex[i] == 'NA':
					outputfile.write('NA,')
				else:
					outputfile.write(str(int(genes[i][loci-1][jk]))+',')
			if sex[i] == 'NA':
				outputfile.write('NA\n')
			else:
				outputfile.write(str(int(genes[i][loci-1][alleles[loci-1]-1]))+'\n')
		
		# if file genes or random genes 			
		elif intgenesans == 'file' or intgenesans == 'random':
		
			#Store empty array to be appended to for gene info
			indall = []
						
			# And store genes information
			genes.append([])
			
			# For each loci:
			for j in xrange(loci):
				
				# Select 2 alleles at random with replacement if random
				if intgenesans == 'random':
					rand1 = int(alleles[j]*random.random())
					rand2 = int(alleles[j]*random.random())
					
				# Take a random draw from the w_choice function at jth locus
				elif intgenesans == 'file':
					rand1 = w_choice_general(allelst[j])[0]
					rand2 = w_choice_general(allelst[j])[0]

				# Store genes loci spot
				genes[i].append([])
				
				# Append assinment onto indall array - run through each condition for assignment of 1s or 2s or 0s
				# 	1s = heterozygous at that locus
				#	2s = homozygous at that locus
				#	0s = absence of allele
				for k in xrange(alleles[j]):
					
					# Somebody not in this spot
					if sex[i] == 'NA':
						# THen append tempindall to indall
						indall.append('NA')
						# And to genes list
						genes[i][j].append('NA')
					
					# Else if somebody is in spot, assign genes
					else:
					
						# Assignment of 2, the rest 0
						if rand1 == rand2: 
							if k < rand1 or k > rand1:
								tempindall = 0
							elif k == rand1:
								tempindall = 2
							
						# Assignment of 1s, the rest 0
						if rand1 != rand2:
							if k < min(rand1,rand2) or k > max(rand1,rand2):
								tempindall = 0
							elif k == rand1 or k == rand2:
								tempindall = 1
							else:
								tempindall = 0
								
						# THen append tempindall to indall
						indall.append(tempindall)
						
						# And to genes list
						genes[i][j].append(tempindall)
		
			# Add indall information to outputfile text
			for j in xrange(len(indall)-1):
				outputfile.write(str(indall[j])+',')
				
			# To get return character on the end
			outputfile.write(str(indall[len(indall)-1])+'\n')
			
	# Close file
	outputfile.close()
	
	# if file genes or random genes 			
	if intgenesans == 'file' or intgenesans == 'random':
		# Delete extra
		del(indall)
		
	# Return variables
	return genes
	
	# End::DoGridOut_cdpop0()	
	
# ---------------------------------------------------------------------------------------------------	 
def DoPreProcess(outdir,cdmatfilename,foldertime,ibatch,ithmcrun,\
xyfilename,loci,alleles,gen,logfHndl,matemovethresh,minmaxdispthresh,\
intgenesans,allefreqfilename,agefilename,agemort,cdevolveans,\
offspringfitsurfaceAA,offspringfitsurfaceAa,offspringfitsurfaceaa):
	'''
	DoPreProcess()
	This function does all the pre-processing work before
	CDFISH begins its time loops.
	'''
	# ----------------------------
	# Create directory
	# ----------------------------		
	ithmcrundir = outdir+'batchrun'+\
	str(ibatch)+'mcrun'+str(ithmcrun)+'/'
	os.mkdir(ithmcrundir)

	# ------------------------------------------------------------------
	# Read in cdmatrix # Eventually convert to probability here...
	# ------------------------------------------------------------------ 
	tupReadMat = ReadCDMatrix(cdmatfilename)
	
	# Unpack tuple, cdmatrix and 1 less number of subpopulations (length
	#	of cdmatrix included source as a point).
	cdmatrix = tupReadMat[0]
	nosubpops = tupReadMat[1] - 1
	
	# Extract source to subpop, this is the first row of cdmatrix
	source_subpop = cdmatrix[0][1:tupReadMat[1]]
	
	# Get maximum cdvalue to use for movethreshold if specified
	tupGetMax = GetMaxCDValue(matemovethresh,cdmatrix,minmaxdispthresh)
	
	# Unpack tuple
	matemovethresh = tupGetMax[0]
	minmaxdispthresh = tupGetMax[1]
	
	# ------------------------------------------------------------------
	# Read in xy points file and store info in list
	# ------------------------------------------------------------------ 
	xy = ReadXY(xyfilename)
	
	# Error statement for 5 column data
	if len(xy[1]) != 5 and intgenesans!='known':
		print('XY input file must be 5 columns, see example input files.')
		sys.exit(-1)
	
	# Store all information in lists by variable name
	FID = []
	subpop = []
	xgrid = []
	ygrid=[]
	id = []
	sex = []	
	age = []	
	genes = []
	for i in xrange(len(xy)-1):
		FID.append(i)
		subpop.append(xy[i+1][0])
		if xy[i+1][1] == 'NA' or xy[i+1][2] == 'NA':
			print('You need to specify the (x,y) locations for this location even if it is a NA value.')
			sys.exit(-1)
		xgrid.append(float(xy[i+1][1]))
		ygrid.append(float(xy[i+1][2]))
		id.append(xy[i+1][3])
		sex.append(xy[i+1][4].strip('\n'))

		# If sex was F or M...change to 0 and 1
		for isexspot in xrange(len(sex)):
			if sex[isexspot] == 'F':
				sex[isexspot] = 0
			elif sex[isexspot] == 'M':
				sex[isexspot] = 1
					
		# Store genetic information: genes[individual][locus][allele]
		if intgenesans=='known':
			# Error check here to make sure gene file matches specified loci and alleles
			if sum(alleles) != len(xy[i+1][6:len(xy[i+1])]):
				print('Known genes file does not match loci and alleles given.')
				sys.exit(-1)
			genes.append([])			
			for j in xrange(loci):
				genes[i].append(xy[i+1][int(6+sum(alleles[0:j])):int(6+sum(alleles[0:j+1]))])
			# Strip the last end charater
			genes[i][loci-1][alleles[loci-1]-1] = genes[i][loci-1][alleles[loci-1]-1].strip('\n')
			
	# Store the number of grids
	nogrids = len(xy)-1
	
	# Delete x variable
	del(xy)
	
	# --------------------------
	# Error Checks
	# --------------------------
	# For now, subpops need to be ordered 1 to N and not skipping, no 0s
	if len(np.where(np.unique(subpop)=='0')[0]) != 0:
		print('Subpopulation identification field can not have 0 values.')
		sys.exit(-1)
	tempcheck = []
	for i in xrange(len(np.unique(subpop))):
		tempcheck.append(int(np.unique(subpop)[i]))
	tempcheck = np.sort(tempcheck)
	if len(tempcheck) > 1:
		for i in xrange(len(tempcheck)-1):
			if tempcheck[i+1]-tempcheck[i] > 1:
				print('Subpopulation identification field must be labeled sequentially or a single value.')
				sys.exit(-1)
	
	# -------------------------------------------
	# Read in fitness surfaces
	# -------------------------------------------
	
	# 1 Locus fitness surfaces files read in
	# ---------------------------------------
	fitvals1 = []
	if cdevolveans == '1':
		
		# Read in fitness ascii grid for offspring fitness surface
		tupReadFitAA = ReadFitnessSurface(offspringfitsurfaceAA)
		# Unpack tuple
		offspringfitnessvaluesAA = tupReadFitAA[0]
		ncols = tupReadFitAA[1]
		nrows = tupReadFitAA[2]
		xll = tupReadFitAA[3]
		yll = tupReadFitAA[4]
		cellsize = tupReadFitAA[5]
		Nodataval = tupReadFitAA[6]
		
		# Read in fitness ascii grid for offspring fitness surface
		tupReadFitAa = ReadFitnessSurface(offspringfitsurfaceAa)
		# Unpack tuple
		offspringfitnessvaluesAa = tupReadFitAa[0]
		
		# Read in fitness ascii grid for offspring fitness surface
		tupReadFitaa = ReadFitnessSurface(offspringfitsurfaceaa)
		# Unpack tuple
		offspringfitnessvaluesaa = tupReadFitaa[0]
		
		# Error check on fitness surface dimensions
		if (len(offspringfitnessvaluesAA) == len(offspringfitnessvaluesAa) == len(offspringfitnessvaluesaa)) == False:
			print('Spatial gradient surfaces are not the same size. Crop to the same dimensions (ncols x nrows).')
			sys.exit(-1)
		if (len(offspringfitnessvaluesAA[0]) == len(offspringfitnessvaluesAa[0]) == len(offspringfitnessvaluesaa[0])) == False:
			print('Spatial gradient surfaces are not the same size. Crop to the same dimensions (ncols x nrows).')
			sys.exit(-1)
					
		# -------------------------------------------
		# Store fitness values for each grid spot
		# -------------------------------------------
		
		# Store x and y locations of each fitness value pixel and fitnessvalue 
		xfitval1 = []
		yfitval1 = []
		fitnessvalues1 = []
		
		# Now loop row
		for irow in xrange(nrows):
			
			# Column loop
			for icol in xrange(ncols):
			
				# Get key spot name
				yspot = yll+(cellsize*(nrows-1-irow))
				yfitval1.append(yspot)				
				
				# Get key spot name			
				xspot = xll+(cellsize*icol)
				xfitval1.append(xspot)
				
				# Append fitness values
				fitnessvalues1.append([float(offspringfitnessvaluesAA[irow][icol]),\
				float(offspringfitnessvaluesAa[irow][icol]),\
				float(offspringfitnessvaluesaa[irow][icol])])
			
		# Make numpy xgrid and ygrid
		dataxy = np.zeros((len(xgrid),2))
		dataxy[:,0] = xgrid
		dataxy[:,1] = ygrid
		
		# Make numpy xfitval and yfitval
		datafitv = np.zeros((len(xfitval1),2))
		datafitv[:,0] = xfitval1
		datafitv[:,1] = yfitval1
		
		# Call the KDTree function to get closest values
		tree = KDTree(datafitv)
		
		# Query the tree to get fixed points
		fixed_pts = tree.query(dataxy)
		
		# Then grab fitnessvalue
		for ifit in xrange(nogrids):		
			fitvals1.append([fitnessvalues1[fixed_pts[1][ifit]][0],\
			fitnessvalues1[fixed_pts[1][ifit]][1],\
			fitnessvalues1[fixed_pts[1][ifit]][2]])
		
		# THen delete finessvalues...only fitvals passes on.
		del(fitnessvalues1)
		del(xfitval1)
		del(yfitval1)
		del(dataxy)
		del(datafitv)
		
	# -----------------------------------------------------------
	# Initialize age structure and return age file information 
	# ----------------------------------------------------------- 
	tupAgeFile = InitializeAge(agefilename,nogrids,agemort)
	agelst = tupAgeFile[0]
	Mprob_mature = tupAgeFile[1] 
	Mprob_migrate = tupAgeFile[2]
	Mprob_juvenile = tupAgeFile[3]
	Fprob_mature = tupAgeFile[4] 
	Fprob_migrate = tupAgeFile[5] 
	Fprob_juvenile = tupAgeFile[6]
	agemort = tupAgeFile[7]
	
	# --------------------------------------------------------------
	# Initialize genetic structure
	# --------------------------------------------------------------
	allelst = InitializeGenes(intgenesans,allefreqfilename,loci,alleles)
	
	# --------------------------------------------------------------------
	# Create output file grid0.csv and write to it and return genes
	# -------------------------------------------------------------------- 
	genes = DoGridOut_cdpop0(ithmcrundir,0,loci,alleles,nogrids,subpop,\
	xgrid,ygrid,id,sex,age,agelst,genes,allelst,intgenesans)
	
	# --------------------------------
	# Get sub population totals
	# --------------------------------
	# Count up the unique number of subgrids appending to subgrids
	subtotal = []
	# Create list of lists storage spots for number of subgrids
	for i in xrange(nosubpops):
		subtotal.append([])
	
	# Get unique subpop names
	unique_subpops = np.unique(subpop)
	
	for i in xrange(len(id)):
		# Loop through unique subpops
		for j in xrange(nosubpops):
			# If subpop exists append to subgrid spot
			if subpop[i] == unique_subpops[j] and age[i] != 'NA':
				subtotal[int(unique_subpops[j])-1].append(1)
				
	# And then sum them up
	for i in xrange(nosubpops):
		subtotal[i] = sum(subtotal[i])
	filledgrids = sum(subtotal)
	
	# Return this functions variables
	tupPreProcess = ithmcrundir,FID,id,sex,age,xgrid,ygrid,genes,nogrids,\
	filledgrids,subpop,cdmatrix,nosubpops,source_subpop,subtotal,\
	matemovethresh,minmaxdispthresh,fitvals1,Mprob_mature,\
	Mprob_migrate,Mprob_juvenile,Fprob_mature,Fprob_migrate,\
	Fprob_juvenile,agemort
	return tupPreProcess
	
	#End::DoPreProcess()
	
# ---------------------------------------------------------------------------------------------------	 		
def DoUserInput(fileans):
	
	'''
	DoUserInput()
	This function reads in the user input and 
	stores the variables.
	'''
	
	# Open file for reading
	inputfile = open(fileans,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	inputvariables = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.split(',')
		inputvariables.append(thisline)
		
	# Delete lines
	del(lines)

	return inputvariables
	
	#End::DoUserInput()