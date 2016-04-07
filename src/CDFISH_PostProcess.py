# -------------------------------------------------------------------------------------------------
# CDFISH_PostProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file post processing.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
except ImportError:
	raise ImportError, "Numpy required."
import pdb,random
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
def DoGridOut_cdpop(ithmcrundir,gen,loci,alleles,nogrids,subpopnew,xgridnew,\
ygridnew,idnew,sexnew,agenew,genesnew,logfHndl):
	'''
	DoGridOut_cdpop()
	Output grid.csv in cdpopformat	
	'''	
	# Create file to write info to
	outputfile = open(ithmcrundir+'grid'+str(gen+1)+'.csv','w')
	
	# Write out the titles - Add Titles from xypoints
	title = ['Subpopulation','XCOORD','YCOORD','ID','sex','age']
	
	# Write out the title from xy points
	for i in xrange(len(title)):
		# Write out FID number
		outputfile.write(title[i]+',')
			
	# Write out the loci title info
	# Loop for loci length
	for i in xrange(loci-1):
		# Loop for allele length
		for j in xrange(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	# To get a return character on the end of the title
	for i in xrange(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1)+'\n')

	# Loop through each grid spot and output
	for i in xrange(len(subpopnew)):
		
		outputfile.write(subpopnew[i]+',')
		outputfile.write(str(float(xgridnew[i]))+',')
		outputfile.write(str(float(ygridnew[i]))+',')
		outputfile.write(idnew[i]+',')
		outputfile.write(sexnew[i]+',')
		outputfile.write(str(agenew[i])+',')
		
		# Write out gene info
		for jk in xrange(loci-1):
			for kl in xrange(alleles[jk]):
				outputfile.write(str(genesnew[i][jk][kl])+',')
		# TO get return character on end
		for jk in xrange(alleles[loci-1]-1):
			outputfile.write(str(genesnew[i][loci-1][jk])+',')
		outputfile.write(str(genesnew[i][loci-1][alleles[loci-1]-1])+'\n')
				
	# Logging message
	stringout = 'The file grid'+str(gen+1)+'.csv has been created'
	logMsg(logfHndl,stringout)		
	
	# Close file
	outputfile.close()
	# End::DoGridOut_cdpop()
	
# ---------------------------------------------------------------------------------------------------	 
def DoOutput(nogrids,FID,OffDisperseIN,xgridcopy,ygridcopy,gen,id,\
sex,age,xgrid,ygrid,genes,nthfile,ithmcrundir,loci,alleles,subpop,\
logfHndl,freegridisolated):
	'''
	DoOutput()
	Generate .txt file.
	Input: ithmcrundir
	Output: ithmcrundir will have .csv files of x,y coord location values of
	cost distance dispersal of the old+new generation with gene info	
	'''	
	
	# ----------------------------------------------------------------------------------------
	# Order the grids back from the random processes - 0,1,2,...nogrids
	# ----------------------------------------------------------------------------------------
	
	# Storage lists for ordering id and no
	orderofgridid = []
	orderofgridno = []
	
	# Loop through all grid points
	for i in xrange(nogrids):
	
		# Loop through the FID values if any
		for jFID in xrange(len(FID)):
			if int(FID[jFID])==i:
				orderofgridid.append('FID'+str(jFID))
				orderofgridno.append(jFID)
		
		# Loop through the dispersal values if any
		for jFG in xrange(len(OffDisperseIN)):
			if int(OffDisperseIN[jFG][1])==i:
				orderofgridid.append('FG'+str(jFG))
				orderofgridno.append(jFG)
				
		# Loop through the free grid values if any
		for jopen in xrange(len(freegridisolated)):
			if int(freegridisolated[jopen])==i:
				orderofgridid.append('OPEN'+str(jopen))
				orderofgridno.append(jopen)
	
	# ------------------------------------------------------------------------------------------------------------------
	# Store the OffDisperseIN generation information: Update grid...
	# ------------------------------------------------------------------------------------------------------------------
	
	# Store new grid values
	FIDnew = []
	subpopnew = []
	xgridnew = []
	ygridnew = []
	idnew = []
	sexnew = []	
	agenew = []	
	genesnew = []
		
	# Extract information from OffDisperseIN
	if len(OffDisperseIN) != 0:
	
		# Temp variables storage
		xtempoff = []
		ytempoff = []
		offsex=[]
		offid=[]
		
		# Extract grid x and y location, id,sex
		for dispot in xrange(len(OffDisperseIN)):
			xtempoff.append(xgridcopy[OffDisperseIN[dispot][1]])
			ytempoff.append(ygridcopy[OffDisperseIN[dispot][1]])
			offsex.append(OffDisperseIN[dispot][0][2])
			offid.append(OffDisperseIN[dispot][2])
	
	# Loop through each grid spot to find the right output index...
	for jgrid in xrange(nogrids):
				
		# Write out the old generation
		if orderofgridid[jgrid][0:3] == 'FID':
			
			# Write out the old generation
			FIDnew.append(FID[orderofgridno[jgrid]])
			idnew.append(id[orderofgridno[jgrid]])
			sexnew.append(str(sex[orderofgridno[jgrid]]))
			# This adds a year to the age class
			agenew.append(int(age[orderofgridno[jgrid]])+1)
			xgridnew.append(float(xgrid[orderofgridno[jgrid]]))
			ygridnew.append(float(ygrid[orderofgridno[jgrid]]))
			# Make sure this is the correct order writing to file!
			# Write out gene info
			genesnew.append(genes[orderofgridno[jgrid]])
			subpopnew.append(subpop[jgrid])
					
		# Get the ordered grid location for the jth FG value (the OffDisperseIN)
		elif orderofgridid[jgrid][0:2] == 'FG':
					
			# Write out the OffDisperseIN generation
			FIDnew.append(str(OffDisperseIN[orderofgridno[jgrid]][1]))
			idnew.append(offid[orderofgridno[jgrid]])
			sexnew.append(str(int(offsex[orderofgridno[jgrid]])))
			# This makes the age of the offspring 0
			agenew.append(int(0))
			xgridnew.append(xtempoff[orderofgridno[jgrid]])
			ygridnew.append(ytempoff[orderofgridno[jgrid]])
			subpopnew.append(subpop[jgrid])
			# Write out gene info
			genesnew.append([])
			for jloci in xrange(loci):
				genesnew[jgrid].append([])
				for jalleles in xrange(alleles[jloci]):
					genesnew[jgrid][jloci].append(OffDisperseIN[orderofgridno[jgrid]][0][3][jalleles+sum(alleles[0:jloci])])
			
		# Get the ordered grid location for the jth OPEN value (freegridisolated)
		elif orderofgridid[jgrid][0:4] == 'OPEN':

			# Write out the freegridisolated information
			FIDnew.append(str(freegridisolated[orderofgridno[jgrid]]))
			idnew.append('OPEN')
			sexnew.append('NA')
			agenew.append('NA')
			xgridnew.append(float(xgridcopy[freegridisolated[orderofgridno[jgrid]]]))
			ygridnew.append(float(ygridcopy[freegridisolated[orderofgridno[jgrid]]]))
			subpopnew.append(subpop[jgrid])
			# Write out gene info
			genesnew.append([])
			for jloci in xrange(loci):
				genesnew[jgrid].append([])
				for jalleles in xrange(alleles[jloci]):
					genesnew[jgrid][jloci].append('NA')
										
	# ------------------------------------------------------------
	# Write out text file for generations specified by nthfile
	# ------------------------------------------------------------
	# Check if nthfile == generation
	for inthfile in xrange(len(nthfile)):				
		if gen == nthfile[inthfile]:
			
			# Call DoGridOut_cdpop()
			DoGridOut_cdpop(ithmcrundir,gen,loci,alleles,nogrids,\
			subpopnew,xgridnew,ygridnew,idnew,sexnew,agenew,genesnew,\
			logfHndl)
			
	# Return variables from this argument
	tupDoOut = FIDnew,idnew,sexnew,agenew,xgridnew,ygridnew,\
	genesnew,subpopnew
	return tupDoOut
	
	# End::DoOutput()
	
# ---------------------------------------------------------------------------------------------------	 
def DoPostProcess(ithmcrundir,nogrids,xgridcopy,ygridcopy,\
loci,alleles,looptime,Population,ToTFemales,ToTMales,Migrants,OpenLocations,Births,\
Deaths,Alleles,He,Ho,AllelesMutated,nthfile,gen,logfHndl,Residors,\
Strayer,DispDeaths,OffLeftOver,subpop,SelectionDeaths,p1,q1):
	'''
	DoPostProcess()
	Creates the output.csv file.
	'''	
		
	# Create time array
	time = np.arange(0,gen+1,1)
	
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	
	# Create file to write info to
	outputfile = open(ithmcrundir+'output.csv','w')
	
	# Write out the titles
	outputtitle = ['Generation','Population','Females','Males',\
	'Migrants','Open','Residors','Stayers','DispersalDeaths','SelectionDeaths','OffspringDiscarded',\
	'Births','AdultDeaths','Alleles','He',\
	'Ho','Mutations','p1','q1']
	
	# Write out the title
	for i in xrange(len(outputtitle)-1):
		outputfile.write(outputtitle[i]+',')
	# To get return character on the end
	outputfile.write(str(outputtitle[len(outputtitle)-1])+'\n')		
	
	# Write to file
	for i in xrange(len(time)):
		outputfile.write(str(time[i])+',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Population[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTFemales[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(ToTMales[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(Migrants[i])+',')
		outputfile.write(str(OpenLocations[i])+',')
		outputfile.write(str(Residors[i])+',')
		outputfile.write(str(Strayer[i])+',')
		outputfile.write(str(DispDeaths[i])+',')
		outputfile.write(str(SelectionDeaths[i])+',')
		outputfile.write(str(OffLeftOver[i])+',')
		outputfile.write(str(Births[i])+',')
		for j in xrange(len(Deaths[i])):
			outputfile.write(str(Deaths[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Alleles[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(He[i][j])+'|')
		outputfile.write(',')
		for j in xrange(nosubpops+1):
			outputfile.write(str(Ho[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(AllelesMutated[i])+',')
		outputfile.write(str(p1[i])+',')
		outputfile.write(str(q1[i])+'\n')
						
	# Logging message
	stringout = 'The file outputfile.csv has been created'
	logMsg(logfHndl,stringout)	
	
	# Close file
	outputfile.close()