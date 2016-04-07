# -------------------------------------------------------------------------------------------------
# CDFISH_GetMetrics.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for get metrics processes.
# --------------------------------------------------------------------------------------------------

def calcF(F,Ho,He,gen):
	#	F is the inbreeding coefficient: 1-Hot/Het
	F.append(1. - (Ho[gen][0]/He[gen][0]))
	return F
def calcFST(FST,he_tot,HS):
	#	Fst is the effect of subpop to total: (1-((Hs)/Het)
	FST.append((he_tot - HS)/he_tot)
	return FST
def calcFIS(FIS,HS,HI):
	#	Fis is the effect of departure within subpops: (1-((1./2)*(Ho1+Ho2))/((1./2)*(He1+He2)))
	FIS.append((HS-HI)/HS)
	return FIS
def calcFIT(FIT,he_tot,HI):
	#	Fit is the effect of departure amoung subpops: (1-((1./2)*(Ho1+Ho2))/(Het))
	FIT.append((he_tot-HI)/he_tot)
	return FIT