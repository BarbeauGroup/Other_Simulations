#                                                                                            #
# ========================================================================================== #
# ========================================================================================== #
# =================== Python Based Inverse Beta Decay Event Generator ====================== #
# ========================================================================================== #
# ============================= Written by C. Awe - May 2018 =============================== #
# ========================================================================================== #
# ========================================================================================== #
#                                                                                            #
#                      Code to generate Monte Carlo IBD events.                              #
#                                                                                            #
# Based on "Precise quasielastic neutrino/nucleon cross-section" - A. Strumia and F. Vissani #
#                                                                                            #
#   Theta is defined as the angle between the positron and neutrino in the lab frame.        #
#                                Energy/Mass in MeV                                          #
#                                      c = 1                                                 #
#                                    hbar = 1                                                #
#                              cross sections in MeV^-2                                      #

import ROOT
import random
from ROOT import TGraph, TH1D, TH1F, TH2F, TCanvas, TFile, TTree
import numpy as np
import scipy as sp
import matplotlib as mpl
import pandas as pds
import sklearn as skl
from array import array
from tqdm import tqdm

cos_theta_c = 0.9742915 #Cosine of the Cabibo angle.
g_f = 1.16637e-11 #Fermi coupling constant in MeV^-2.
m_p = 938.27203 #Proton mass in MeV/c^2
m_n = 939.56536 #Neutron mass in MeV/c^2
m_e = 0.5109989 #Electron mass in MeV/c^2
m_A = 1e3 #Constant term in form factor expressions.
m_V = 842.615 #Constant term in form factor expressions.
m = 938.9 #Average nucleon mass in MeV/c^2.
m_pi = 134.976 #Charged pion mass.
#Constant term defined in Strumia.
delta = float( ( np.square( m_n ) - np.square( m_p ) - np.square( m_e ) ) / ( 2 * m_p ) ) 
D = 1.293 #Mass difference between proton and neutron.
#Threshold set 1 eV higher to allow leptons a full angular distribution.
E_thr = 1.807 #IBD threshold in MeV.
ksi = 3.706 * 1.0 / (2.0 * m_p ) #Difference between proton and neutron anomalous MM"s (nuclear magnetons).
numSamples = 1000 #Number of values to take for numerical integration.

################################################
############# Basic Math Functions #############
################################################

#Computes epsilon (we go to NLO in epsilon).
def epsilon( E_v ):
	return float( E_v / m_p )
	
#Computes kappa, just a number that shows up in the equations.	
def kappa( e, cos_theta_e ):
	return float( np.square( 1 + e ) - np.square( e * cos_theta_e ) )

#Computes positron energy given E_nu and theta_e	
def posEnergy( E_v, cos_theta_e ):
	e = epsilon( E_v )
	k = kappa( e, cos_theta_e )
	return float( ( ( E_v - delta ) * ( 1 + e ) + e * cos_theta_e * np.sqrt(
	np.square( E_v - delta ) - np.square( m_e ) * k ) ) / k )

#Computes positron momentum given E_e	
def posMomentum( E_e ):
	return float( np.sqrt( np.square( E_e ) - np.square( m_e ) ) )

#Computes neutron energy given E_nu and E_e	
def ntronEnergy( E_v, E_e ):
	global m_p
	return float( E_v + m_p - E_e )

#Computes neutron momentum given neutron energy.
def ntronMomentum( E_n ):
	#print E_n - m_n
	return float( np.sqrt( np.square( E_n ) - np.square( m_n ) ) )
	
#Computes neutron angle given neutrino energy and positron angle.
def ntronAngle( E_v, cos_theta_e, P_e, P_n ):
	#Calculate neutron scattering angle.
	ntronAngl = float( ( E_v - P_e * cos_theta_e ) / P_n )
	#Check that cos_theta_n <= 1. This can fail in extreme cases where rounding becomes important.
	if ( ntronAngl > 1.0 ):
		#print "Warning: Numerical Rounding Error, cos_theta_n > 1"
		#print "This may indicate unphysical neutrino energies. Setting to 1.0.""
		ntronAngl = 1.0
	return ntronAngl

#Computes the opening angle between the daughter particles.	
def openingAngle( cos_theta_e, cos_theta_n ):
	return np.arccos( cos_theta_e ) * 180.0 / np.pi +  np.arccos( cos_theta_n ) * 180.0 / np.pi
	
#############################################
############## Build 4-Vectors ##############
#############################################

#All neutrinos assumed to be moving along the x-axis.
#The positron and neutron define the x-y plane.
#z-axis is neglected (azimuthal symmetry).
#Neutrino set to be massless.
#This section not strictly necessary, but may be useful for debugging.

def get_v_vec( E_v ):
	v_vec = np.array[E_v, E_v, 0, 0]
	return v_vec

#Proton assumed at rest in the lab frame.
def get_P_vec():
	P_vec = np.array[m_p, 0, 0, 0]
	return P_vec

def get_pos_vec( E_e, cos_theta_e, P_e ):
	pos_vec = np.array[E_e, P_e * cos_theta_e, P_e * np.sqrt( 1 - cos_theta_e
	* cos_theta_e ), 0]
	return pos_vec

def get_ntron_vec( E_n, cos_theta_n, P_n ):
	ntron_vec = np.array[E_n, P_n * cos_theta_n, P_n * np.sqrt( 1 - cos_theta_n
	* cos_theta_n ), 0]
	return ntron_vec

##############################################
########### Build Mandelstrom Vars ###########
##############################################

def getS( E_v ):
	return float( m_p * m_p + 2 * E_v * m_p )

def getT( E_v, E_e, cos_theta_e, P_e ):
	sin_theta_e = float( np.sqrt( 1 - np.square( cos_theta_e ) ) )
	return float( np.square( E_v - E_e ) - np.square( E_v - P_e * cos_theta_e ) 
	- np.square( P_e * sin_theta_e ) )

def getU( E_v, E_n, cos_theta_n, P_n ):
	sin_theta_n = np.sqrt( 1 - np.square( cos_theta_n ) )
	return float( np.square( E_v - E_n ) - np.square( E_v - P_n * cos_theta_n ) 
	- np.square( P_n * sin_theta_n ) )

#############################################
############# Build Form Factors ############
#############################################

def getf1( t ):
	return float( ( 1 - ( 1 + ksi ) * t / ( 4 * np.square( m ) ) ) / ( ( 1 - t / ( 4 * np.square( m ) ) )
	* np.square( 1 - t / np.square( m_V ) ) ) )

def getf2( t ):
	return float( (ksi ) / ( ( 1 - t / ( 4 * np.square( m ) ) )
	* np.square( 1 - t / np.square( m_V ) ) ) )

def getg1( t ):
	g0 = -1.270
	return float( g0 / np.square( 1 - t / np.square( m_A ) ) )
	
#Not needed for NLO approximation, but appears in the full expression.	
def getg2( t, g1 ):
	return float( 2 * np.square( m ) * g1 / ( np.square( m_pi ) - t ) )
	
def getA( t, f1, g1, f2 ):
	return float( np.square( m ) * ( np.square( f1 ) - np.square( g1 ) )
	* ( t - np.square( m_e ) ) - np.square( m ) * np.square( D ) * 
	( np.square( f1 ) - np.square( g1 ) ) - 2 * np.square( m_e ) *
	m * D * g1 * ( f1 + f2 ) )

def getB( t, f1, g1, f2 ):
	return float( t * g1 * ( f1 + f2 ) )

def getC( f1, g1 ):
	return float( ( np.square( f1 ) + np.square( g1 ) ) / 4 )

#############################################
########### Build Matrix Elements ###########
#############################################

def getMSquared( s, t, u, A, B, C ):
	return float( A - ( s - u ) * B + np.square( s - u ) * C )

#############################################
########## Compute Differential CC ##########
#############################################

def getdCC( cos_theta_e, E_e, P_e, e, s, M ):
	return float( 2 * m_p * P_e * e / ( 1 + e * ( 1 - E_e / P_e * cos_theta_e ) )
	* np.square( g_f ) * np.square( cos_theta_c ) / ( 2 * np.pi * 
	np.square( s - np.square( m_p ) ) ) * M )

#############################################	
###### Sample the reactor nu spectrum #######
#############################################

def getEnergy( filename ):
	f = ROOT.TFile( filename )
	hist = f.Get("specHist")
	return hist.GetRandom()

##############################################	
#### Sample Positron Angular Distribution ####
##############################################

#Create a random positron angle (cos_theta_e between -1 and 1).
#May need to add a line seeding this to be fully rigorous.
def getAngle():
	return random.uniform(-1,1)
	
#############################################
######## Compute Total Cross Section ########
#############################################

def getCC( E_v ):
	totCC = 0
	for i in range( -numSamples / 2, numSamples / 2 ):
		cos_theta_e = float( 2 * i / numSamples )
		e = epsilon( E_v )
		k = kappa( e, cos_theta_e )
		E_e = posEnergy( E_v, cos_theta_e )
		P_e = posMomentum( E_e )
		E_n = ntronEnergy( E_v, E_e )
		P_n = ntronMomentum( E_n )
		cos_theta_n = ntronAngle( E_v, cos_theta_e, P_e, P_n )
		angl = openingAngle( cos_theta_e, cos_theta_n )
		s = getS( E_v )
		t = getT( E_v, E_e, cos_theta_e, P_e )
		u = getU( E_v, E_n, cos_theta_n, P_n )
		f1 = getf1( t )
		f2 = getf2( t )
		g1 = getg1( t )
		A = getA( t, f1, g1, f2 )
		B = getB( t, f1, g1, f2 )
		C = getC( f1, g1 )
		M = getMSquared( s, t, u, A, B, C )
		dCC = getdCC( cos_theta_e, E_e, P_e, e, s, M )
		totCC += float( dCC * 3.89105e-22 * 2 / numSamples ) #differential cc * step size
	return totCC

#############################################	
############### Main Function ###############
#############################################

#Supply this with the name of the root file containing the spectrum and 
#the number of IBD events to generate.
def main():
	#Initialize everything.
	E_v = array( 'd', [0] )
	E_e = array( 'd', [0] )
	E_n = array( 'd', [0] )
	cos_theta_e = array( 'd', [0] )
	cos_theta_n = array( 'd', [0] )
	dCC = array( 'd', [0] )
	angl = array( 'd', [0] )

	#Get run parameters from the user.
	filename = input( "Enter the path to the neutrino spectrum .root file: " )
	numEvents = int( input( "Enter the number of IBD events you wish to generate: " ) )
	output = input( "Enter the path and name of the output file you wish to produce: " )
	
	#Set up root stuff.
	eventFile = ROOT.TFile( output, "recreate" )
	eventTree = ROOT.TTree( "eventTree", "IBD Event Tree" )
	eventTree.Branch( "E_v", E_v, "E_v/D" )
	eventTree.Branch( "E_e", E_e, "E_e/D" )
	eventTree.Branch( "E_n", E_n, "E_n/D" )
	eventTree.Branch( "cos_theta_e", cos_theta_e, "cos_theta_e/D" )
	eventTree.Branch( "cos_theta_n", cos_theta_n, "cos_theta_n/D" )
	eventTree.Branch( "dCC", dCC, "dCC/D" )
	eventTree.Branch( "openingAngle", angl, "angl/D" )
	posSpecHist = ROOT.TH1F("posSpecHist","Positron Spectrum",1000,0,100)
	ntronSpecHist = ROOT.TH1F("ntronSpecHist","Neutron Spectrum",1000,0,5)
	posAnglHist = ROOT.TH1F("posAnglHist","Positron Angular Distribution",100,-1,1)
	ntronAnglHist = ROOT.TH1F("ntronAnglHist","Neutron Angular Distribution",50,0,1)
	anglHist = ROOT.TH1F("anglHist","Opening Angle Distribution",1800,0,180)
	ntronVposHist = ROOT.TH2F("ntronVposHist","Positron Energy vs. Neutron Energy",1000,0,2,1000,0,90)
	
	#Loop through generating events. tqdm gives a progress bar.
	ibdCount = 0
	for i in tqdm( range( 0, numEvents ) ):
		#E_v[0] = getEnergy( filename )
		E_v[0] = 5.0;
		#Check that we"re above IBD threshold.
		if E_v[0] > E_thr:
			cos_theta_e[0] = getAngle()
			e = epsilon( E_v[0] )
			k = kappa( e, cos_theta_e[0] )
			E_e[0] = posEnergy( E_v[0], cos_theta_e[0] )
			P_e = posMomentum( E_e[0] )
			E_n[0] = ntronEnergy( E_v[0], E_e[0] )
			P_n = ntronMomentum( E_n[0] )
			cos_theta_n[0] = ntronAngle( E_v[0], cos_theta_e[0], P_e, P_n )
			angl[0] = openingAngle( cos_theta_e[0], cos_theta_n[0] )
			s = getS( E_v[0] )
			t = getT( E_v[0], E_e[0], cos_theta_e[0], P_e )
			u = getU( E_v[0], E_n[0], cos_theta_n[0], P_n )
			f1 = getf1( t )
			f2 = getf2( t )
			g1 = getg1( t )
			A = getA( t, f1, g1, f2 )
			B = getB( t, f1, g1, f2 )
			C = getC( f1, g1 )
			M = getMSquared( s, t, u, A, B, C )
			dCC[0] = getdCC( cos_theta_e[0], E_e[0], P_e, e, s, M )
			eventTree.Fill()
			posSpecHist.Fill( E_e[0] - m_e, dCC[0] )
			ntronSpecHist.Fill( E_n[0] - m_n, dCC[0] )
			posAnglHist.Fill( cos_theta_e[0], dCC[0] )
			ntronAnglHist.Fill( cos_theta_n[0], dCC[0] )
			anglHist.Fill( angl[0], dCC[0] )
			ntronVposHist.Fill( E_n[0] - m_n, E_e[0] - m_e, dCC[0] )
			ibdCount += 1
			
	#Report the number of neutrinos above threshold and write our file.
	print str( ibdCount ) + " neutrinos out of " + str( numEvents ) + " above threshold."
	eventFile.Write()
	eventFile.Close()

#Execute main function 	
if __name__== "__main__":
  main()		
		
		
		
		
		
		
		
		
		
