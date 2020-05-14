import ROOT
import random
from ROOT import TGraph, TH1D, TH1F, TH2F, TCanvas, TFile, TTree
import numpy as np
import scipy as sp
import matplotlib as mpl
import pandas as pds
import sklearn as skl
from array import array

#Write a crude reactor antineutrino spectrum to a root file. Based on Klapdor 1982.
def writeReactorSpectrum():
	#Initialize output file and histogram.
	specFile = ROOT.TFile("reactorNuSpec.root","recreate")
	specHist = ROOT.TH1D("specHist","Reactor Neutrino Spectrum",100,1,9)
	
	#Populate two arrays with values sampled from the spectrum.
	numPoints = 17
	x, y = array( "d" ), array( "d" )
	x.append( 1.0 );
	y.append( 2.36 );
	x.append( 1.5 );
	y.append( 1.71 );
	x.append( 2.0 );
	y.append( 1.31 );
	x.append( 2.5 );
	y.append( 0.888 );
	x.append( 3.0 );
	y.append( 0.613 );
	x.append( 3.5 );
	y.append( 0.412 );
	x.append( 4.0 );
	y.append( 0.268 );
	x.append( 4.5 );
	y.append( 0.16 );
	x.append( 5.0 );
	y.append( 0.097 );
	x.append( 5.5 );
	y.append( 0.0596 );
	x.append( 6.0 );
	y.append( 0.0346 );
	x.append( 6.5 );
	y.append( 0.0189 );
	x.append( 7.0 );
	y.append( 0.01 );
	x.append( 7.5 );
	y.append( 0.00399 );
	x.append( 8.0 );
	y.append( 0.00131 );
	x.append( 8.5 );
	y.append( 0.00031 );
	x.append( 9.0 );
	y.append( 0.000141 );
	#Make a TGraph to interpolate these values.
	specGraph = ROOT.TGraph(numPoints,x,y)
	
	#Generate 1000 interpolated points based on our raw spectrum.
	for i in range(0,10000):
		energy = float(8.0) * i / 10000 + 1
		#print energy
		weight = float(specGraph.Eval(energy)) * 100 / 10000
		#print weight
		specHist.Fill(energy,weight)
	
	#Write the output file.
	specFile.Write()
	specFile.Close()
	
#Parse a datatheif csv file and generate a SN spectrum in root.
def writeSNSpectrum( csvFile, outputFilename ):	
	#Open our csv file and read the values into memory.
	csv = open( csvFile, "r" )
	values = csv.read()
	csv.close()
	
	#Split our csv file into separate lines.
	lines = values.split( '\n' )
	
	#Buil arrays to hold a TGraph.
	x, y = array( "d" ), array( "d" )
	
	#Initialize output file and histogram.
	numPoints = len(lines)
	specFile = ROOT.TFile( outputFilename, "recreate" )
	specHist = ROOT.TH1D( "specHist", "Neutrino Spectrum", 100, 0.3, 54.3 )
	
	for index, line in enumerate( lines[ 0: numPoints - 1 ] ):
		args = line.split( ',' )
		x.append( float( args[0] ) )
		y.append( float( args[1] ) )
		
	#Fill TGraph
	specGraph = ROOT.TGraph(numPoints,x,y)
	
	#Generate 1000 interpolated points based on our raw spectrum.
	for i in range(0,10000):
		energy = float( 54.0 ) * i / 10000 + 1
		#print energy
		weight = float(specGraph.Eval(energy)) * 100 / 10000
		#print weight
		specHist.Fill(energy,weight)
	
	#Write the output file.
	specFile.Write()
	specFile.Close()
