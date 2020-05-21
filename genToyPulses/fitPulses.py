# Procedure to fit LS pulses with a double exponential corresponding to a slow and fast decay time.
# Based on the paper: https://arxiv.org/abs/1912.07682
#
# The procedure there consisted of multiple parts:
#   1. Fit neutron population to get t_f, t_s components
#   2. Fit combined gamma/neutron population to get distributions of t_r, R, etc.
#
# Usage:
#       python fitPulsesNonRooFit.py
#
# Results will be stored in "fitResults.root"
#
# Notes:
#   - Both this and the secondary code use samples rather than ns
#   - The neutron cut, and all event cut (probably just an integral cut since low-E events
#     are hard to fit) need to be set in the main loop. That's also where you'll set the
#     t_f and t_s guesses and ranges, and if you're fixing them their fixed values.
#
import ROOT
import sys
import array

def plotWaveform(wf):
  try:
    c1
  except NameError:
    c1=ROOT.TCanvas("c1","c1")
    
  wflen=len(wf)
  ahist=ROOT.TH1D("ahist","ahist",wflen,0,wflen)
  for i in range(0,len(wf)):
    ahist.SetBinContent(i+1,wf[i])
  ahist.Draw()
  c1.Modified()
  c1.Update()
  try:
    input("Press enter to continue.")
  except SyntaxError:
    pass
  ahist.Delete()
  
fitNeutrons=1 #Set to 1 first, then look at fit results to get mean of t_f and t_s values
              #Then set to 0 for full population
baselineBins=40 #How many sample for baseline
waveformLength=300 #Full waveform length in samples
channel=24 #Channel in tree to fit

#Limit for integration, assumed data was collected with internal triggers so pulses
# will always occur in the same location
integral_startSample=100
integral_endSample=300

#Same but for the tail integral for PSD
tailIntegral_startSample=113
tailIntegral_endSample=300

#Parameter limits for fitting, all in samples
onsetTime_guess=105
onsetTime_min=95
onsetTime_max=115

riseTime_guess=1.0
riseTime_min=0.0
riseTime_max=1.5

plotRange_min=60
plotRange_max=180

inpFile=ROOT.TFile(sys.argv[1],"READ")
tree=inpFile.Get("sis3316tree")
nEntries=tree.GetEntries()
print("Found "+str(nEntries)+" entries")
 
A=array.array('d',[0])
onset=array.array('d',[0])
riseTime=array.array('d',[0])
R=array.array('d',[0])
S=array.array('d',[0])
fastTime=array.array('d',[0])
slowTime=array.array('d',[0])
slowestTime=array.array('d',[0])
baseline=array.array('d',[0])
psd=array.array('d',[0])

outFile=ROOT.TFile("fitResults.root","RECREATE")
fitTree=ROOT.TTree("fitTree","")
fitTree.Branch('A',A,'A/D')
fitTree.Branch('onset',onset,'onset/D')
fitTree.Branch('R',R,'R/D')
fitTree.Branch('S',S,'S/D')
fitTree.Branch('riseTime',riseTime,'riseTime/D')
fitTree.Branch('fastTime',fastTime,'fastTime/D')
fitTree.Branch('slowTime',slowTime,'slowTime/D')
fitTree.Branch('slowestTime',slowestTime,'slowestTime/D')
fitTree.Branch('baseline',baseline,'baseline/D')
fitTree.Branch('psd',psd,'psd/D')
 
for entry in range(0,int(nEntries)):

  if (entry%1000==0):
    print("On entry "+str(entry)+" of "+str(nEntries))

  #Get entry
  tree.GetEntry(entry)
  
  #If it's ch. 0, we'll fit
  if tree.channelID==channel:
  
    ########################
    ##Make data set to fit##
    ########################
    #Make a histogram to fit it
    hist=ROOT.TH1D("hist","hist",waveformLength,0,waveformLength)
    #Fill histogram
    for i in range(0,waveformLength):
      hist.SetBinContent(i+1,tree.waveform[i])
      
    ####################################
    ##Get integral to guess amplitudes##
    ####################################
    base=0
    for i in range(0,baselineBins):
      base+=tree.waveform[i]
    base=base*1./baselineBins
    
    integral=0
    for i in range(integral_startSample,integral_endSample):
      integral+=(tree.waveform[i]-base)
    tailIntegral=0
    for i in range(tailIntegral_startSample,tailIntegral_endSample):
      tailIntegral+=(tree.waveform[i]-base)
    
    if integral>0:
      psd[0]=tailIntegral*1./integral
    else:
      psd[0]=0
      
    ############
    ##Make TF1##
    ############
    neutronCut = integral>=1000 and integral < 4500 and psd[0]>=0.3 and psd[0]<0.6
    allEventCut = integral>=1000
    if fitNeutrons==1:
      cut=neutronCut
    else:
      cut=allEventCut
      
    if allEventCut==1:
      fit=ROOT.TF1("fit","[0]/(exp(([1]-x)/[2])+1) * ([3]/(exp((x-[1])/[4])+1) + (1-[3])/(exp((x-[1])/[5])+1)) + [6]",0,waveformLength)
      
      fit.SetParameter(0,integral)
      fit.SetParLimits(0,0,integral*2.0)
      fit.SetParName(0,"A")
      
      fit.SetParameter(1,onsetTime_guess)
      fit.SetParLimits(1,onsetTime_min,onsetTime_max)
      fit.SetParName(1,"t0")
      
      fit.SetParameter(2,riseTime_guess)
      fit.SetParLimits(2,riseTime_min,riseTime_max)
      fit.SetParName(2,"t_r")
      
      #Shouldn't need to change
      fit.SetParameter(3,0.90)
      fit.SetParLimits(3,0,1)
      fit.SetParName(3,"R")
      
      if fitNeutrons==1:
        fit.SetParName(4,"t_f")
        fit.SetParameter(4,1.2)
        fit.SetParLimits(4,0,8.0)
        
        fit.SetParName(5,"t_s")
        fit.SetParameter(5,20.)
        fit.SetParLimits(5,1.5,60)
      else:
        fit.SetParName(4,"t_f")
        fit.FixParameter(4,1.343)
      
        fit.SetParName(5,"t_s")
        fit.FixParameter(5,10.831)
      
      fit.SetParameter(6,base)
      fit.SetParLimits(6,base*0.8,base*1.2)
      fit.SetParName(6,"baseline")
      
      hist.Fit(fit,"QM0","",0,waveformLength)
      hist.GetXaxis().SetRangeUser(plotRange_min,plotRange_max)
      
      A[0]=fit.GetParameter(0)
      onset[0]=fit.GetParameter(1)
      riseTime[0]=fit.GetParameter(2)
      R[0]=fit.GetParameter(3)
      fastTime[0]=fit.GetParameter(4)
      slowTime[0]=fit.GetParameter(5)
      baseline[0]=fit.GetParameter(6)
      fitTree.Fill()
      
    ##Memory Management##
    hist.Delete()

outFile.cd()
fitTree.Write()
outFile.Close()
inpFile.Close()
