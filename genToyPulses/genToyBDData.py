# Generates toy liquid scintillator pulses. Uses the functional form
# of a liquid scintillator approach from here: https://arxiv.org/abs/1912.07682
#
# I'm using parameters I fit to get, I've left in Forrests as well if you want to try others.
#
# Functions:
#   plotPulse(pulseList) - Plots a pulse on a TH1 based on a list of samples
#   generateNoiseHistogram(filename, channel) - Looks at baseline samples from
#     a specific channel to generate a noise histogram. Returns mean and stdev.
#     This is used to identify a threshold for determining whether a baseline
#     region is empty.
#   generateNoiseSamples(filename,channel,mean,sigma) - Steps through a TTree
#     and harvests as many empty baseline samples as possible. Returns a list
#     of empty baselines.
#   combineNoiseToMakeWaveforms(noiseList) - Takes a list of empty waveform
#     regions, and sticks them together to make continuous noise waveforms
#   generatePulses(nPulses) - Generates toy BD pulses, add them to a noise
#     waveform. Returns a list of fake pulses.
#   generateToyDataTree(pulses,Rs,onsets,noiseTraces,channel) - Makes a fake
#     data tree based on a list of fake pulses
#
# Usage:
#   python genToyBDData.py <root file with sis3316tree>
#
#Notes:
#   - Expects two root files of neutron energies (neutronHist.root) and
#     gammas (gammaHist.root) with histograms named "htemp" to draw amplitude
#     distributions from.
#   - Set up for a sis3316tree, you'll have to change a lot of things if
#     trying to read/generate data for a different digitizer
#   - Generates an equal number of neutron and gamma pulses, this can be changed
#     in the code

# Used to generate toy
import ROOT
import sys
import random
import numpy
import math
import array

#
fitNeutrons=1

#Name of file to store fake pulses in
outputName="fakeBDPulses.root"

#Samples to use for harvesting noise. Note that we actually need baselineSamples+1
#element available because we'll compare the baselineSamples+1 element of one
#trace to the 0th sample of the next to stitch them together continuously
baselineSamples=50

#Combine 6 baseline regions to form a full trace of noise. Should be an integer
#multiple of baselineSamples
fullTraceSamples=6*baselineSamples

#Integration lengths, there must be this many samples after the onset in trace
integrationLength=100 #Samples, starts integrating at onset
tailIntegralDelay=7 #Start tail integral this many samples after integral starts

#Which channel to use
channel=24

#Pulse onset location. Not really gaussian, but approximating it as such
#Draw onset from this distribution
onsetMean=106 #samples
onsetSigma=0.5#samples

#How many fake pulses to generate
nPulses=5030



#Plots a pulse based on a list of samples passed in.
def plotPulse(pulseList):
  c1=ROOT.TCanvas("c1","c1")
  hist=ROOT.TH1D("hist","hist",len(pulseList),0,len(pulseList))
  for i in range(0,len(pulseList)):
    hist.SetBinContent(i+1,pulseList[i])
  hist.Draw()
  c1.Modified()
  c1.Update()
  try:
    input("Press enter to continue")
  except SyntaxError:
    pass
  hist.Delete()

#Read in waveforms, analyze baseline regions to get noise TH1D
def generateNoiseHistogram(filename,channel):
  
  #Read in sis3316tree, use first x samples to generate a histogram of noise.
  inpFile=ROOT.TFile(filename,"READ")
  tree=inpFile.Get("sis3316tree")
  nEntries=tree.GetEntries()
  
  #Optimizations
  channelIDBranch=tree.GetBranch("channelID")
  waveformBranch=tree.GetBranch("waveform")
  cacheSize=500000000
  tree.SetCacheSize(cacheSize)
  tree.AddBranchToCache(channelIDBranch,1)
  tree.AddBranchToCache(waveformBranch,1)
  tree.StopCacheLearningPhase()
  
  #Make Histogram
  noiseHist = ROOT.TH1D("noiseHist","noiseHist;ADC;Counts",16384,0,16384)
  
  #Loop through TTree
  for entry in range(0,nEntries):
    tree.LoadTree(entry);
    channelIDBranch.GetEntry(entry)
    if tree.channelID==channel:
      waveformBranch.GetEntry(entry)
      
      for sample in range(0,baselineSamples):
        noiseHist.Fill(tree.waveform[sample])

    if entry%100000==0:
      print("On entry "+str(entry)+" of "+str(nEntries))
  
  fit=ROOT.TF1("fit","gaus",500,2500)
  noiseHist.Fit(fit,"QR0","",500,2500)
  mean=fit.GetParameter(1)
  sigma=fit.GetParameter(2)
  return mean,sigma

def generateNoiseSamples(filename,channel,mean,sigma):
  #Samples to use for harvesting noise
  threshold=mean+2*sigma

  #Read in sis3316tree, use first x samples to generate a histogram of noise.
  inpFile=ROOT.TFile(filename,"READ")
  tree=inpFile.Get("sis3316tree")
  nEntries=tree.GetEntries()

  #Optimizations
  channelIDBranch=tree.GetBranch("channelID")
  waveformBranch=tree.GetBranch("waveform")
  cacheSize=500000000
  tree.SetCacheSize(cacheSize)
  tree.AddBranchToCache(channelIDBranch,1)
  tree.AddBranchToCache(waveformBranch,1)
  tree.StopCacheLearningPhase()

  noiseList=[]
  
  #Loop through TTree
  for entry in range(0,nEntries):
    tree.LoadTree(entry);
    channelIDBranch.GetEntry(entry)
    if tree.channelID==channel:
      waveformBranch.GetEntry(entry)
      
      waveformList=[]
      for sample in range(0,baselineSamples+1):
        waveformList.append(tree.waveform[sample])

      #print(sum(waveformList)/baselineSamples)
      if sum(waveformList)*1./(baselineSamples+1) < threshold:
        noiseList.append(waveformList)

    if entry%100000==0:
      print("On entry "+str(entry)+" of "+str(nEntries))
      
  return noiseList

def combineNoiseToMakeWaveforms(noiseList):
  
  noiseTraces=[]
  currentTrace=[]
  
  while len(noiseList)>0:
  
    if (len(noiseList)%1000==0):
      print(str(len(noiseList))+" partial noise waveform samples remaining")
  
    #Grab the first entry if there is not a current entry
    if len(currentTrace)==0:
      currentTrace=noiseList.pop(0)
      lastADC=currentTrace.pop(len(currentTrace)-1)
    else:
      #Get list of possible waveforms to add to it
      possibleWaveforms = [i for i in noiseList if i[0]==lastADC]
      
      #If there's at least one possible waveform to add, choose a
      #random valid waveform, remove from list, and add to our current trace
      if len(possibleWaveforms)!=0:
        nextElement=possibleWaveforms.pop(random.randrange(len(possibleWaveforms)))
        noiseList.remove(nextElement)
        currentTrace.extend(nextElement)
        lastADC=currentTrace.pop(len(currentTrace)-1)
      else:
        #No valid succesors found, delete the trace we've built so far
        currentTrace=[]
        
    if len(currentTrace)==fullTraceSamples:
      noiseTraces.append(currentTrace)
      currentTrace=[]
      
  return noiseTraces
  
#Generates pulses
def generatePulses(nPulses):

  nPulsesToGenerate=nPulses
  
  pulseLength=fullTraceSamples #samples

  '''
  #Forrest settings
  riseTime=2.1 #ns
  fastDecayComponent=5.45 #ns
  slowDecayComponent=52.0#ns
  
  t_r = riseTime/nsPerSample #rise time, samples
  t_f = t_f_n = t_f_g = fastDecayComponent/nsPerSample #fast decay comp., samples
  t_s = t_f_n = t_f_h = slowDecayComponent/nsPerSample #slow decay comp., samples
  
  #Use this for pulse shape https://arxiv.org/pdf/1912.07682.pdf
  R_neutrons = 0.949
  R_neutrons_sigma = 0.008
  R_gammas = 0.999
  R_gammas_sigma = 0.003
  '''
  
  #Use this for pulse shape https://arxiv.org/pdf/1912.07682.pdf
  #Data-based settings settings
  t_r = 0.289 #+-0.044
  t_f = 1.343 #+-0.378
  t_s = 10.831 #+-4.036
  
  R_neutrons = 0.920
  R_neutrons_sigma = 0.022
  R_gammas = 0.984
  R_gammas_sigma = 0.007
  
  #Use data-pulled distributions of energy of pulses
  neutronHistFile=ROOT.TFile("neutronHist.root","READ")
  neutronHist=neutronHistFile.Get("htemp")
  gammaHistFile=ROOT.TFile("gammaHist.root","READ")
  gammaHist=gammaHistFile.Get("htemp")
  
  pulses=[]
  onsets=[]
  Rs=[]
  
  for ipulse in range(0,nPulsesToGenerate):

    pulse=[]
    
    #Get onset, assume normally distributed independent of shape, amplitude
    t0 = numpy.random.normal(onsetMean,onsetSigma)
    
    #Determine whether this is a neutron or gamma with equal probability
    isNeutron = numpy.random.uniform(0,2)
    
    #Sample from appropriate energy distribution
    if isNeutron>=1:
      A = neutronHist.GetRandom()
      R = numpy.random.normal(R_neutrons,R_neutrons_sigma)
    else:
      A = gammaHist.GetRandom()
      R = numpy.random.normal(R_gammas,R_gammas_sigma)
      
    #Sample from timing distributions

    #Generate shape
    for sample in range(0,pulseLength):
      f = 1./(math.exp(-(sample-t0)/t_r)+1)
      g = 1./(math.exp((sample-t0)/t_f)+1)
      h = 1./(math.exp((sample-t0)/t_s)+1)
      
      pulse.append(int(A * f * (R*g+(1-R)*h)))
      onsets.append(int(t0))

    Rs.append(R)
    pulses.append(pulse)
    
  return pulses,Rs,onsets
      
def generateToyDataTree(pulses,Rs,onsets,noiseTraces,channel):
  
  numSamples=fullTraceSamples
  
  #Output file
  outFile=ROOT.TFile(outputName,"RECREATE")
  sis3316tree=ROOT.TTree("sis3316tree","Unsorted events")
  channelID=array.array('H',[0])
  timestamp=array.array('L',[0])
  peakHighIndex=array.array('H',[0])
  peakHighValue=array.array('H',[0])
  pileupFlag=array.array('H',[0])
  nSamples=array.array('i',[0])
  accumulatorSum=array.array('d',8*[0])
  waveform=array.array('H',numSamples*[0])
  trueIntegral=array.array('d',[0])
  R=array.array('d',[0])
  psd=array.array('d',[0])
  
  sis3316tree.Branch('channelID',channelID,'channelID/s')
  sis3316tree.Branch('timestamp',timestamp,'timestamp/l')
  sis3316tree.Branch('peakHighIndex',peakHighIndex,'peakHighIndex/s')
  sis3316tree.Branch('peakHighValue',peakHighValue,'peakHighValue/s')
  sis3316tree.Branch('pileupFlag',pileupFlag,'pileupFlag/O')
  sis3316tree.Branch('nSamples',nSamples,'nSamples/i')
  sis3316tree.Branch('waveform',waveform,'waveform[300]/s')
  sis3316tree.Branch('accumulatorSum',accumulatorSum,'accumulatorSum[8]/i')
  
  sis3316tree.Branch('trueIntegral',trueIntegral,'trueIntegral/D')
  sis3316tree.Branch('R',R,'R/D')
  sis3316tree.Branch('psd',psd,'psd/D')
  
  channelID[0]=channel
  timestamp[0]=0
  peakHighIndex[0]=0
  peakHighValue[0]=0
  pileupFlag[0]=0
  nSamples[0]=fullTraceSamples
  
  for pulseNum in range(0,len(pulses)):
    pulse=pulses[pulseNum]
    R[0]=Rs[pulseNum]
    noiseTrace=noiseTraces[random.randrange(len(noiseTraces))]
    realPulse=[(a + b) for a, b in zip(pulse, noiseTrace)]
    #plotPulse(noiseTrace)
    #plotPulse(pulse)
    #plotPulse(realPulse)
    saturatedEvents=[i for i in realPulse if i >= 16384]
    if len(saturatedEvents)==0:
      for i in range(0,len(pulse)):
        waveform[i]=realPulse[i]
      
      if onsets[i]+integrationLength < fullTraceSamples:
        integratedSection=[pulse[i] for i in range(onsets[i],onsets[i]+integrationLength)]
        tailIntegralSection=[pulse[i] for i in range(onsets[i]+tailIntegralDelay,onsets[i]+integrationLength)]
        trueIntegral[0]=sum(integratedSection)
        tailIntegral=sum(tailIntegralSection)
        if tailIntegral>0:
          psd[0]=tailIntegral/trueIntegral[0]
          sis3316tree.Fill()
      
  sis3316tree.Write()
  outFile.Close()



#Make histogram of baselines, generate mean and sigma
print("Calculating Baseline...")
mean,sigma = generateNoiseHistogram(sys.argv[1],channel)
print("Baseline is "+str(mean)+" +- "+str(sigma)+"\n")

#Make list of noise trace samples
print("Generating noise samples...")
noiseList=generateNoiseSamples(sys.argv[1],channel,mean,sigma)
print("Found "+str(len(noiseList))+" valid baseline windows\n")

#make noise traces
print("Generating noise traces...")
noiseTraces=combineNoiseToMakeWaveforms(noiseList)
print("Generated "+str(len(noiseTraces))+" noise traces\n")

#Make raw pulses
print("Generating fake pulse shapes...")
pulses,Rs,onsets=generatePulses(nPulses)
print("Generated "+str(nPulses)+" fake pulse shapes\n")

#Make fake pulses
print("Making fake pulses...")
generateToyDataTree(pulses,Rs,onsets,noiseTraces,channel)
print("Done!")
