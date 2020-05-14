# Process PyBD output for panel detector concept.
# C. Awe  - 9/18/2018

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

m_e = 0.5109989 # Electron mass in MeV/c^2
m_n = 939.56536 # Neutron mass in MeV/c^2
m_n_kg = 1.6749286E-27 # Neutron mass in kg
MeV_to_J = 1.6021773e-13 # Converts MeV to Joules
mps_to_cmpns = 1e-7 # Converts m/s to cm/ns

# Compute neutron velocity in cm/ns
def getNtronVelocity( E_n ):
	# Get neutron kinetic energy.
	k_n = E_n - m_n
	# Convert to joules.
	k_n *= MeV_to_J
	# Get velocity (assumes non-relativistic energies, good for most IBD scenarios)
	return float( mps_to_cmpns * np.sqrt( 2 * k_n / m_n_kg ) )

# Full 3D case
def getPhi():
	return random.uniform( 0, 2.0 * np.pi )

# Randomly generate IBD position (horizontal)	
def getX( h ):
 return random.uniform( 0, h )

# Randomly generate IBD position (vertical)	
def getY( h ):
 return random.uniform( 0, h )

# Compute neutron TOF
# d = panel seperation in cm
# h = panel height in cm 
# theta = neutron angle in radians
# y = IBD vertical position in cm
# x = IBD horizontal position in cm
# Neglecting panel thickness for now
def getNtronTOF( d, h, theta, x, y, E_n):
	# Check if the neutron is going to miss the next panel.
	phi = getPhi()
	v_n = getNtronVelocity( E_n )
	y_n = float( d * np.tan( theta ) )
	# Dummy variable to record neutron displacement off neutrino axis.
	h_n = y_n
	# Compute horizontal displacement.
	# Adjust neutron position by a phi rotation.
	x_n = float( y_n * np.sin( phi ) ) 
	y_n *= np.cos( phi )
	y_n += y # Neutron vertical position.
	x_n += x # Neutron horizontal position.
	if ( y_n > h ):
		# If so, assign an unphysical value.
		return -1.0
	elif ( y_n < 0 ):
		# If so, assign an unphysical value.
		return -1.0
	elif ( x_n > h ):
		# If so, assign an unphysical value.
		return -1.0
	elif ( x_n < 0 ):
		# If so, assign an unphysical value.
		return -1.0
	else: 
		# Good event, compute TOF.
		return float( np.sqrt( np.square( d ) + np.square( h_n ) ) / v_n )
	
def main():
	# Initialize dummy arrays
	TOF = array( 'd', [0] )
	E_v = array( 'd', [0] )
	E_e = array( 'd', [0] )
	E_n = array( 'd', [0] )
	cos_theta_n = array( 'd', [0] )
	cos_theta_e = array( 'd', [0] )
	dCC = array( 'd', [0] )
	angl = array( 'd', [0] )
	theta_n = array( 'd', [0] )
	
	# Prompt user for input file and run parameters
	#filename = input( "Enter the path to the PyBD output file: " )
	#d = input( "Enter panel seperation in cm: " )
	#h = input( "Enter panel height in cm: " )
	#output = input( "Enter the path and name of the output file you wish to produce: " )
	filename = "~/Desktop/PyBD/RootSpectra/reactor_1e6.root"
	d = 30.0
	h = 30.0
	output = "~/Desktop/PyBD/PostProcessed/reactorPanels.root"

	# Load in root file.
	f = ROOT.TFile( filename )
	ibdTree = f.Get( "eventTree" )
	ibdTree.SetBranchAddress( "E_v", E_v )
	ibdTree.SetBranchAddress( "E_e", E_e )
	ibdTree.SetBranchAddress( "E_n", E_n )
	ibdTree.SetBranchAddress( "cos_theta_n", cos_theta_n )
	ibdTree.SetBranchAddress( "cos_theta_e", cos_theta_e )
	angl = array( 'd', [0] )
	ibdTree.SetBranchAddress( "dCC", dCC )
		
	# Declare a new TTree with branches to hold our processed values
	processedFile = ROOT.TFile( output, "recreate" )
	processedTree = ibdTree.CloneTree(0)
	processedTree.Branch( "TOF", TOF, "TOF/D" )
	
	# Declare histograms to hold CC weighted events
	En_v_TOF = ROOT.TH2F("En_v_TOF","Time of Flight vs. Neutron Energy",2000,0,200,1000,0,1000)
	engyHist_1 = ROOT.TH2F("engyHist_1","Neutron Energy vs. Positron Energy for TOF = 200 +/- 2 ns",50,0,5,400,0,40)
	engyHist_2 = ROOT.TH2F("engyHist_2","Neutron Energy vs. Positron Energy for TOF = 150 +/- 2 ns",50,0,5,400,0,40)
	engyHist_3 = ROOT.TH2F("engyHist_3","Neutron Energy vs. Positron Energy for TOF = 300 +/- 2 ns",50,0,5,400,0,40)
	engyHist_4 = ROOT.TH2F("engyHist_3","Neutron Energy vs. Positron Energy for TOF = 250 +/- 2 ns",50,0,5,400,0,40)
			
	# Step through each event and compute TOF
	numEvents = ibdTree.GetEntries()
	for i in tqdm( range( 0, 500000 ) ):
		ibdTree.GetEntry( i )
		phi = getPhi()
		x = getX( h )
		y = getY( h )
		theta_n = np.arccos( cos_theta_n[0] )
		TOF[0] = getNtronTOF( d, h, theta_n, x, y, E_n[0] )
		processedTree.Fill()
		En_v_TOF.Fill( 1000 * ( E_n[0] - m_n ), TOF[0], dCC[0] )
		if ( TOF[0] > 198 and TOF[0] < 202 ):
			engyHist_1.Fill( E_e[0] - m_e, 1000 * ( E_n[0] - m_n ), dCC[0] )
			En_v_TOF.Fill( 1000 * ( E_n[0] - m_n ), TOF[0], dCC[0] )
		if ( TOF[0] > 148 and TOF[0] < 152 ):
			engyHist_2.Fill( E_e[0] - m_e, 1000 * ( E_n[0] - m_n ), dCC[0] )
			En_v_TOF.Fill( 1000 * ( E_n[0] - m_n ), TOF[0], dCC[0] )
		if ( TOF[0] > 298 and TOF[0] < 302 ):
			engyHist_3.Fill( E_e[0] - m_e, 1000 * ( E_n[0] - m_n ), dCC[0] )
		if ( TOF[0] > 248 and TOF[0] < 252 ):
			engyHist_4.Fill( E_e[0] - m_e, 1000 * ( E_n[0] - m_n ), dCC[0] )
		#print "TOF = " + str( TOF[0] ) + " ns."
		
	# Write to an output file
	processedFile.Write()
	processedFile.Close()
	
#Execute main function 	
if __name__== "__main__":
  main()
	
	
		
