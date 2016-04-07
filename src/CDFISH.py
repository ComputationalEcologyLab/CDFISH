# CDFISH.py
# Author: Erin L Landguth
# Created: January 2011
# v 0.52 -- 2011 Sept 22
# ----------------------------------------------------------------------------
# General CDFISH information
appName = "CDFISH"
appVers = "0.56"
appRele = "2013.11.20"
authorNames = "Erin L. Landguth"

# ---------------
# Global symbols
#----------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False
# File absolute paths for importing functions
SRC_PATH =  "../src/"

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import datetime,time,pdb,os,sys,shutil

# Numpy functions
try:
	import numpy as np                    
except ImportError as eMsg:
	print("ImportError (%s) Numpy required."%(eMsg))
	sys.exit(-1)
	
#Import the package specific folders
CDFISH_folder = os.path.dirname(os.path.abspath(SRC_PATH+"CDFISH"))

if CDFISH_folder not in sys.path:
     sys.path.insert(0, CDFISH_folder)

# CDFISH functions
try:
	from CDFISH_Modules import * 
except ImportError:
	raise ImportError, "CDFISH_Modules required."

#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
	
	# ------------------------------------------------------	
	# Start timer, get script arguments
	# ------------------------------------------------------
	# Timing events: start
	start_time = datetime.datetime.now()
	foldertime = int(time.time())
	
	if len(sys.argv) >= 4:
		datadir = sys.argv[1]+'/'
		fileans = datadir+sys.argv[2]
		outdir = datadir+sys.argv[3]+str(foldertime)+'/'		
	
	# If user did not specify .rip file
	else:
		print "User must specify data directory, input file name, and output file directory (e.g., at command line type CDFISH.py ../CDFISH_data/ cdfish_test.csv cdfish_test_out)."
		sys.exit(-1)		
	
	# If .ip file does not exist
	if not os.path.exists(fileans):
		print("Cannot find or open runtime inputs file(%s)"%(fileans))
		sys.exit(-1)
		
	# Create output file directory - will automatically put in the data directory
	os.mkdir(outdir)
	
	# This properly names log file
	logSessionPath = outdir+"cdfish.log"
	logfHndl =open(logSessionPath,'w')
	
	msgVerbose = True
	logMsg(logfHndl,"\n%s Release %s Version %s\n"%(appName,appRele,appVers))
	logMsg(logfHndl,"Author(s): %s"%(authorNames)+'\n')
	logMsg(logfHndl,"Session runtime inputs from: %s"%(fileans)+'\n\n')    
	msgVerbose = False
				
	# ------------------------------------	
	# Call DoUserInput()
	# ------------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Call Function DoUserInput() -- Reads in file.
	inputvariables = DoUserInput(fileans)
	
	# Print to log
	stringout = 'DoUserInput(): '+str(datetime.datetime.now() -start_time1) + ''
	logMsg(logfHndl,stringout)
	print 'DoUserInput(): ',str(datetime.datetime.now() -start_time1),''

	# -------------------------------------	
	# Begin Batch Looping
	# -------------------------------------
	# This loop is defined by the number of rows in inputvariables.csv
	for ibatch in xrange(len(inputvariables)-1):
	
		# Timing events: start
		start_timeB = datetime.datetime.now()
			
		# Store all information and the type of each, also do some error checks 
		xyfilename = datadir+str(inputvariables[ibatch+1][0])
		agefilename = str(inputvariables[ibatch+1][1]).strip('\n')
		cdmatfilename = datadir+inputvariables[ibatch+1][2]
		mcruns = int(inputvariables[ibatch+1][3])
		looptime = int(inputvariables[ibatch+1][4])
		nthfile_choice = str(inputvariables[ibatch+1][5])
		nthfile_list = inputvariables[ibatch+1][6]
		nthfile_seq = inputvariables[ibatch+1][7]
		matemoveno = str(inputvariables[ibatch+1][8])
		matemovethresh = str(inputvariables[ibatch+1][9])
		freplace = str(inputvariables[ibatch+1][10])
		mreplace = str(inputvariables[ibatch+1][11])
		reproage = int(inputvariables[ibatch+1][12])
		dispmoveno = str(inputvariables[ibatch+1][13])
		dispmovethresh = inputvariables[ibatch+1][14]
		residencypercent = inputvariables[ibatch+1][15]
		straypercent = inputvariables[ibatch+1][16]
		offno = str(inputvariables[ibatch+1][17])
		lmbda = str(inputvariables[ibatch+1][18])
		Femalepercent = int(inputvariables[ibatch+1][19])
		newmortperc = float(inputvariables[ibatch+1][20])/100
		agemortperc = inputvariables[ibatch+1][21]
		equalsexratio = str(inputvariables[ibatch+1][22])
		muterate = float(inputvariables[ibatch+1][23])
		loci = int(inputvariables[ibatch+1][24])
		allelelist = inputvariables[ibatch+1][25]
		intgenesans = str(inputvariables[ibatch+1][26])
		allefreqfilename = str(inputvariables[ibatch+1][27])
		cdevolveans = str(inputvariables[ibatch+1][28])
		burningen = int(inputvariables[ibatch+1][29])
		offspringfitsurfaceAA = datadir+str(inputvariables[ibatch+1][30])
		offspringfitsurfaceAa = datadir+str(inputvariables[ibatch+1][31])
		offspringfitsurfaceaa = datadir+str(inputvariables[ibatch+1][32]).strip('\n')
		
		# Split up alleles, removing values, and appending
		alleles = []
		for iall in xrange(len(allelelist.split('|'))):
			alleles.append(int(allelelist.lstrip('"').split('|')[iall]))
		# If just one value entered, turn into an array
		if len(alleles) == 1:
			alleles = int(inputvariables[ibatch+1][25])*np.ones(loci,int)
		# If more than one value, then do an error check and turn into array
		else:
			if len(alleles) != loci:
				print('Number of loci does not match with the allele length.')
				sys.exit(-1)
			else:
				alleles = np.asarray(alleles)
				print('Code section not finished yet.')
				sys.exit(-1)
			
		# Split up residency subpopulations, removing values, and appending
		residency = []
		for isource in xrange(len(residencypercent.split('|'))):
			residency.append(int(residencypercent.lstrip('"').split('|')[isource]))
			
		# Split up straying percentages for each subpopulation, removing values, and appending
		straying = []
		for isource in xrange(len(straypercent.split('|'))):
			straying.append(int(straypercent.lstrip('"').split('|')[isource]))
				
		# Grab the nthfile list range specific to user input, list or sequence
		if nthfile_choice == 'Sequence' or nthfile_choice == 'sequence':
			# Check if mod == 0, to compute nthfile sequence
			if np.mod(looptime,nthfile_seq) == 0:		
				nthfile = range(0,looptime+int(nthfile_seq),int(nthfile_seq))
			# If mod != 0 then truncate nthfile by one
			else:
				nthfile = range(0,looptime+int(nthfile_seq),int(nthfile_seq))
				del(nthfile[-1])
		if nthfile_choice == 'List' or nthfile_choice == 'list':
			nthfile = []
			# Split up list, removing space values, and appending to nthfile
			for inum in xrange(len(nthfile_list.split('|'))):
				# Don't append 0, this gets written out automatically
				if int(nthfile_list.split('|')[inum]) != 0:
					# And then append 1 less then what user entered, due to the
					#	way I have indexed grid values
					nthfile.append(int(nthfile_list.split('|')[inum])-1)
						
		# Error check on nthfile, must be 1 less than looptime for indexing
		if max(nthfile) >= looptime:
			print 'nthfile selection maximum value must be 1 less than your looptime.'
			sys.exit(-1)
			
		# Check here on split of mortality - in file now
		'''
		agemort = []
		for isource in xrange(len(agemortperc.split('|'))):
			agemort.append(float(agemortperc.lstrip('"').split('|')[isource])/100)
		oldmortperc = agemort[-1]
		'''
		agemort = []
		# If just one returned
		if len(agemortperc.split('|')) == 1:
			oldmortperc = float(agemortperc)/100
			agemort = [oldmortperc]
		elif len(agemortperc.split('|')) != 1:
			print('For age specific mortality, use file format.')
			sys.exit(-1)
			
		# Set directories for agefilename and allelefreqency file
		if agefilename != 'N':
			agefilename = datadir+agefilename
		if intgenesans == 'file' and allefreqfilename != 'N':
			allefreqfilename = datadir+allefreqfilename
		elif intgenesans == 'file' and allefreqfilename == 'N':
			print('Allele frequency file option specified, must give name of file.')
			sys.exit(-1)
			
		# Split up dispersal movement for option 7 min and max
		minmaxdispthresh = []
		for iminmax in xrange(len(dispmovethresh.split('|'))):
			minmaxdispthresh.append(dispmovethresh.lstrip('"').split('|')[iminmax])
		# If just one returned
		if len(dispmovethresh) == 1:
			minmaxdispthresh = str(dispmovethresh)
		
		# ---------------------------------
		# Some Error checking
		# ---------------------------------
		if dispmoveno == '7' and len(minmaxdispthresh) != 2:
			print('dispmovethresh must have 2 min and max threshold bounds.')
			sys.exit(-1)
			
		# If cdevolve is turned on must have 2 alleles
		if cdevolveans != 'N' and alleles[0] != 2:
			print 'Use 2 alleles per locus when CDEVOLVE is turned on.'
			sys.exit(-1)
		# Must have more than 1 loci
		if loci <= 1:
			print('Currently, CDFISH needs more than 1 locus to run.')
			sys.exit(-1)
		
		# ---------------------------------------------	
		# Begin Monte-Carlo Looping
		# ---------------------------------------------
		
		# xrange(mcruns) is typically 10 - 50
		for ithmcrun in xrange(mcruns):	
		
			# Timing events: start
			start_timeMC = datetime.datetime.now()
		
			# -----------------------------------------
			# Create storage variables
			# ------------------------------------------	
			# These variables will be stored in output.csv at the end of the simulation
			Population = []
			Migrants = []
			Deaths = []
			Births = []
			OpenLocations = []
			ToTFemales = []
			ToTMales = []
			Alleles = []
			He = []
			Ho = []
			AllelesMutated = []
			Residors = []
			Strayer = []
			DispDeaths = []
			OffLeftOver = []
			SelectionDeaths = []
			p1 = []
			q1 = []
					
			# ------------------------------------	
			# Call DoPreProcess()
			# ------------------------------------

			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			# Prepare fitness surface grid file
			if cdevolveans == '1':
				offspringfitsurfaceAA = PrepTextFile(offspringfitsurfaceAA)
				offspringfitsurfaceAa = PrepTextFile(offspringfitsurfaceAa)
				offspringfitsurfaceaa = PrepTextFile(offspringfitsurfaceaa)
				
			# Prepare fitness surface grid file
			if cdevolveans == '2':
				print('Not coded in yet.')
				sys.exit(-1)
			
			tupPreProcess = DoPreProcess(outdir,cdmatfilename,foldertime,ibatch,ithmcrun,\
			xyfilename,loci,alleles,0,logfHndl,matemovethresh,minmaxdispthresh,\
			intgenesans,allefreqfilename,agefilename,agemort,
			cdevolveans,offspringfitsurfaceAA,\
			offspringfitsurfaceAa,offspringfitsurfaceaa)
			
			# Unpack tuple
			ithmcrundir = tupPreProcess[0]
			FID = tupPreProcess[1]
			id = tupPreProcess[2]
			sex = tupPreProcess[3]
			age = tupPreProcess[4]
			xgrid = tupPreProcess[5]
			xgridcopy = copy.deepcopy(xgrid)
			ygrid = tupPreProcess[6]
			ygridcopy = copy.deepcopy(ygrid)
			genes = tupPreProcess[7]
			nogrids = tupPreProcess[8]
			filledgrids = tupPreProcess[9]
			subpop = tupPreProcess[10]
			cdmatrix = tupPreProcess[11]
			nosubpops = tupPreProcess[12]
			source_subpop = tupPreProcess[13]
			subpoptotals = tupPreProcess[14]
			matemovethresh = tupPreProcess[15]
			minmaxdispthresh = tupPreProcess[16]
			fitvals1 = tupPreProcess[17]
			Mprob_mature = tupPreProcess[18] 
			Mprob_migrate = tupPreProcess[19]
			Mprob_juvenile = tupPreProcess[20]
			Fprob_mature = tupPreProcess[21] 
			Fprob_migrate = tupPreProcess[22] 
			Fprob_juvenile = tupPreProcess[23]
			agemort = tupPreProcess[24]
			
			# Print to log
			stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'DoPreProcess(): ',str(datetime.datetime.now() -start_time1),''
			
			# -------------------------------------------
			# Start Generation Looping 
			# -------------------------------------------
			
			# Begin generation loop
			for gen in xrange(looptime):
			
				# CDCLIMATE READ IN HERE #
				
				# Timing events: start
				start_timeGen = datetime.datetime.now()
			
				# -------------------------------	
				# Call ReadGrid0()
				# -------------------------------
				# Use information generated from PreProcess step for first 
				#	generation, else use the following updated grid information
				if gen != 0:
				
					# Timing events: start
					start_time1 = datetime.datetime.now()
					
					tupReadGrid = ReadGrid(FIDnew,idnew,agenew,xgridnew,\
					ygridnew,genesnew,sexnew,subpopnew,gen,equalsexratio,\
					nosubpops)
					FID = tupReadGrid[0]
					sex = tupReadGrid[1]
					id = tupReadGrid[2]
					age = tupReadGrid[3]
					xgrid = tupReadGrid[4]
					xgridcopy = tupReadGrid[5]
					ygrid = tupReadGrid[6]
					ygridcopy = tupReadGrid[7]
					genes = tupReadGrid[8]
					nogrids = tupReadGrid[9]
					subpop = tupReadGrid[10]
					filledgrids = tupReadGrid[11]
					subpoptotals = tupReadGrid[12]
					
					# Here we system exit if number of grids left is 0 or 1
					if filledgrids == 0 or filledgrids == 1:
						# Print to log
						stringout = 'Population went extinct, program ended.'
						logMsg(logfHndl,stringout)
						print('Population went extinct after generation '+str(gen-1)+'.\n')
						Population.append(np.zeros(len(subpoptotals)+1))
						Migrants.append(0)
						Deaths.append([0])
						Births.append(0)
						OpenLocations.append(nogrids)
						ToTFemales.append(np.zeros(len(subpoptotals)+1))
						ToTMales.append(np.zeros(len(subpoptotals)+1))
						Alleles.append(np.zeros(len(subpoptotals)+1))
						He.append(np.zeros(len(subpoptotals)+1))
						Ho.append(np.zeros(len(subpoptotals)+1))
						AllelesMutated.append(0)
						Residors.append(0)
						Strayer.append(0)
						DispDeaths.append(0)
						OffLeftOver.append(0)
						SelectionDeaths.append(0)
						p1.append(0)
						q1.append(0)
						break						
					# Here we system exit if there are only F or M left
					if len(np.where(np.unique(sex)=='0')[0]) == 0:
						# Print to log
						stringout = 'No females left in population, program ended.'
						logMsg(logfHndl,stringout)
						print('There are no more females left in population after generation '+str(gen-1)+'.\n')
						Population.append(np.zeros(len(subpoptotals)+1))
						Migrants.append(0)
						Deaths.append([0])
						Births.append(0)
						OpenLocations.append(nogrids)
						ToTFemales.append(np.zeros(len(subpoptotals)+1))
						ToTMales.append(np.zeros(len(subpoptotals)+1))
						Alleles.append(np.zeros(len(subpoptotals)+1))
						He.append(np.zeros(len(subpoptotals)+1))
						Ho.append(np.zeros(len(subpoptotals)+1))
						AllelesMutated.append(0)
						Residors.append(0)
						Strayer.append(0)
						DispDeaths.append(0)
						OffLeftOver.append(0)
						SelectionDeaths.append(0)
						p1.append(0)
						q1.append(0)
						break
					if len(np.where(np.unique(sex)=='1')[0]) == 0:
						# Print to log
						stringout = 'No males left in population, program ended.'
						logMsg(logfHndl,stringout)
						print('There are no more males left in population after generation '+str(gen-1)+'.\n')
						Population.append(np.zeros(len(subpoptotals)+1))
						Migrants.append(0)
						Deaths.append([0])
						Births.append(0)
						OpenLocations.append(nogrids)
						ToTFemales.append(np.zeros(len(subpoptotals)+1))
						ToTMales.append(np.zeros(len(subpoptotals)+1))
						Alleles.append(np.zeros(len(subpoptotals)+1))
						He.append(np.zeros(len(subpoptotals)+1))
						Ho.append(np.zeros(len(subpoptotals)+1))
						AllelesMutated.append(0)
						Residors.append(0)
						Strayer.append(0)
						DispDeaths.append(0)
						OffLeftOver.append(0)
						SelectionDeaths.append(0)
						p1.append(0)
						q1.append(0)
						break
					
					# Print to log
					stringout = 'ReadGrid(): '+str(datetime.datetime.now() -start_time1) + ''
					logMsg(logfHndl,stringout)
					print 'ReadGrid(): ',str(datetime.datetime.now() -start_time1),''
										
				# ---------------------------------
				# Call GetMetrics()
				# ---------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupGetMetrics = GetMetrics(filledgrids,loci,alleles,genes,gen,Ho,\
				Alleles,He,subpop,subpoptotals,p1,q1)
				# Unpack Tuple
				Ho = tupGetMetrics[0]
				Alleles = tupGetMetrics[1]
				He = tupGetMetrics[2]
				allelefreqlst = tupGetMetrics[3]
				
				# Print to log
				stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'GetMetrics(): ',str(datetime.datetime.now() -start_time1),''
												
				# ---------------------------------------
				# Call DoMate()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				Bearpairs = DoMate(nogrids,sex,age,freplace,mreplace,matemoveno,\
				matemovethresh,xgridcopy,ygridcopy,ToTMales,ToTFemales,\
				Population,subpop,nosubpops,subpoptotals,gen,\
				reproage)
				
				# Print to log
				stringout = 'DoMate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoMate(): ',str(datetime.datetime.now() -start_time1),''
								
				# ---------------------------------------
				# Call DoOffspring()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupDoOff = DoOffspring(offno,lmbda,Bearpairs,Femalepercent,\
				Births)				
				offspring = tupDoOff[0]	
				offspringno = tupDoOff[1]
				
				# Print to log
				stringout = 'DoOffspring(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoOffspring(): ',str(datetime.datetime.now() -start_time1),''
								
				# ---------------------------------------
				# Call InheritGenes()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				offspring = InheritGenes(gen,AllelesMutated,offspringno,\
				offspring,genes,loci,muterate)				
									
				# Print to log
				stringout = 'InheritGenes(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'InheritGenes(): ',str(datetime.datetime.now() -start_time1),''
								
				# ------------------------------------------
				# Call DoAdultMortality()
				# ------------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupAMort = DoAdultMortality(filledgrids,nogrids,sex,id,\
				age,xgrid,ygrid,gen,genes,Deaths,alleles,FID,agemort)
				
				freegrid = tupAMort[0]
				id = tupAMort[1]
				sex = tupAMort[2]
				age = tupAMort[3]
				xgrid = tupAMort[4]
				ygrid = tupAMort[5]
				genes = tupAMort[6]	
				FID = tupAMort[7]
				
				# Print to log
				stringout = 'DoAdultMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoAdultMortality(): ',str(datetime.datetime.now() -start_time1),''
								
				# ------------------------------------------
				# Call DoDisperse()
				# ------------------------------------------			
			
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupDoDisp = DoDisperse(offspringno,freegrid,offspring,dispmoveno,\
				minmaxdispthresh,gen,Migrants,OpenLocations,loci,alleles,nogrids,\
				xgridcopy,ygridcopy,allelefreqlst,logfHndl,newmortperc,subpop,\
				source_subpop,residency,nosubpops,Residors,cdmatrix,Strayer,straying,\
				DispDeaths,OffLeftOver,subpoptotals,cdevolveans,fitvals1,burningen,\
				SelectionDeaths)
				
				OffDisperseIN = tupDoDisp[0]
				freegridisolated = tupDoDisp[1]
							
				# Print to log
				stringout = 'DoDisperse(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoDisperse(): ',str(datetime.datetime.now() -start_time1),''
				
				# ------------------------------------------
				# Call DoOutput()
				# ------------------------------------------	
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupDoOut = DoOutput(nogrids,FID,OffDisperseIN,xgridcopy,ygridcopy,gen,id,sex,age,xgrid,\
				ygrid,genes,nthfile,ithmcrundir,loci,alleles,subpop,logfHndl,freegridisolated)
				# Unpack Tuple
				FIDnew = tupDoOut[0]
				idnew = tupDoOut[1]
				sexnew = tupDoOut[2]
				agenew = tupDoOut[3]
				xgridnew = tupDoOut[4]
				ygridnew = tupDoOut[5]
				genesnew = tupDoOut[6]
				subpopnew = tupDoOut[7]
				
				# Print to log
				stringout = 'DoOutput(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print 'DoOutput(): ',str(datetime.datetime.now() -start_time1),''
				
				# Print to log
				stringout = 'End Generation Loop'+str(gen)+': '+str(datetime.datetime.now() -start_timeGen) + '\n'
				logMsg(logfHndl,stringout)
				print 'End Generation Loop',str(gen),': ',str(datetime.datetime.now() -start_timeGen),'\n'
					
			# End::generation loop
						
			# ------------------------------------------
			# Call DoPostProcess()
			# ------------------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			DoPostProcess(ithmcrundir,nogrids,xgridcopy,ygridcopy,\
			loci,alleles,looptime,Population,ToTFemales,ToTMales,Migrants,OpenLocations,Births,\
			Deaths,Alleles,He,Ho,AllelesMutated,nthfile,gen,logfHndl,Residors,\
			Strayer,DispDeaths,OffLeftOver,subpop,SelectionDeaths,p1,q1)
			
			# Print to log
			stringout = 'DoPostProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print 'DoPostProcess(): ',str(datetime.datetime.now() -start_time1),''
				
			# Print to log
			stringout = 'End Monte Carlo Loop'+str(ithmcrun)+': '+str(datetime.datetime.now() -start_timeMC) + '\n'
			logMsg(logfHndl,stringout)
			print 'End Monte Carlo Loop',str(ithmcrun),': ',str(datetime.datetime.now() -start_timeMC),'\n'
						
		# End::Monte Carlo Loop
		
		# Print to log
		stringout = 'End Batch Loop'+str(ibatch)+': '+str(datetime.datetime.now() -start_timeB) + '\n'
		logMsg(logfHndl,stringout)
		print 'End Batch Loop',str(ibatch),': ',str(datetime.datetime.now() -start_timeB),'\n'
		
	#End::Batch Loop
	
# End::Main Loop
	
# Print to log
stringout = 'Total CDFISH Simulation Time: '+str(datetime.datetime.now() -start_time) + ''
logMsg(logfHndl,stringout)
logfHndl.close()
print 'Total CDFISH Simulation Time: ',str(datetime.datetime.now() -start_time),''