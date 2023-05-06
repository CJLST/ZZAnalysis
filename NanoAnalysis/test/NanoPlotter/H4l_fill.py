### Example macro for filling standard histograms from H4l nanoAODs.
### Histograms are stored on a file and can then be plotted with

from __future__ import print_function
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


pathMC = "/eos/user/n/namapane/H4lnano/220420/" # FIXME: Use 2018 MC for the time being
pathDATA = "/eos/user/n/namapane/H4lnano/230506/Data2022/"

ZmassValue = 91.1876

maxEntriesPerSample = 1e12 # Use only up to this number of events in each MC sample, for quick tests.



ROOT.TH1.SetDefaultSumw2()

####################################
def fillHistos(samplename, filename) :

    h_ZZMass2 = ROOT.TH1F("ZZMass_2GeV_"+samplename,"ZZMass_2GeV_"+samplename,65,70.,200.)
    h_ZZMass2.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass2.GetYaxis().SetTitle("Events / 2 GeV")
    
    h_ZZMass4 = ROOT.TH1F("ZZMass_4GeV_"+samplename,"ZZMass_4GeV_"+samplename,233,70.,1002.)
    h_ZZMass4.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass4.GetYaxis().SetTitle("Events / 4 GeV")

    h_ZZMass10 = ROOT.TH1F("ZZMass_10GeV_"+samplename,"ZZMass_10GeV_"+samplename,93,70.,1000.)
    h_ZZMass10.GetXaxis().SetTitle("m_{#it{4l}} (GeV)")
    h_ZZMass10.GetYaxis().SetTitle("Events / 10 GeV")


    f = ROOT.TFile.Open(filename)

    events = f.Events
    events.SetBranchStatus("*", 0)
    events.SetBranchStatus("run", 1)
    events.SetBranchStatus("luminosityBlock", 1)
    events.SetBranchStatus("event", 1)
    events.SetBranchStatus("ZZCand_*", 1)
    events.SetBranchStatus("bestCandIdx", 1)
    events.SetBranchStatus("HLT_passZZ4l", 1)
    nEntries = events.GetEntries() 

    isMC = False
    if(samplename == "Data"):
        print("Data: sel=", nEntries)
    else:
        isMC = True
        events.SetBranchStatus("overallEventWeight",1)

        # Get sum of weights
        runs = f.Runs
        nRuns = runs.GetEntries()
        iRun = 0
        genEventCount = 0
        genEventSumw = 0.
        while iRun < nRuns and runs.GetEntry(iRun) :
            genEventCount += runs.genEventCount
            genEventSumw += runs.genEventSumw
            iRun +=1
        print (samplename, ": gen=", genEventCount, "sel=",nEntries, "sumw=", genEventSumw)
        if nEntries>maxEntriesPerSample :
            genEventSumw = genEventSumw*maxEntriesPerSample/nEntries
            nEntries=maxEntriesPerSample
            print("   scaling to:", nEntries, "sumw=", genEventSumw )


    iEntry=0
    printEntries=max(5000,nEntries/10)
    while iEntry<nEntries and events.GetEntry(iEntry):
        iEntry+=1
        if iEntry%printEntries == 0 : print("Processing", iEntry)

        bestCandIdx = events.bestCandIdx
        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        if(bestCandIdx != -1 and events.HLT_passZZ4l): 
            weight = 1.
            if isMC : weight = (events.overallEventWeight*events.ZZCand_dataMCWeight[bestCandIdx])/genEventSumw
            m4l=events.ZZCand_mass[bestCandIdx]
            h_ZZMass2.Fill(m4l,weight)
            h_ZZMass4.Fill(m4l,weight)
            h_ZZMass10.Fill(m4l,weight)
        
    f.Close()
    
    return h_ZZMass2,h_ZZMass4,h_ZZMass10

def runMC():
    outFile = "H4l_MC.root" 

    samples = [
        dict(name = "WWZ",filename = pathMC+"WWZ/ZZ4lAnalysis.root"),
        dict(name = "WZZ",filename = pathMC+"WZZ/ZZ4lAnalysis.root"),
        dict(name = "ZZZ",filename = pathMC+"ZZZ/ZZ4lAnalysis.root"),
        
        dict(name = "VBFToZZTo4l",filename = pathMC + "VBFToContinToZZTo4l/ZZ4lAnalysis.root"),
        dict(name = "TTZToLLNuNu",filename = pathMC + "TTZToLLNuNu_M10ext1/ZZ4lAnalysis.root"),
        dict(name = "TTZJets",filename = pathMC + "TTZJets_M10_MLMext1/ZZ4lAnalysis.root"),

        dict(name = "ggTo4mu",filename = pathMC+"ggTo4mu_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo4e",filename = pathMC+"ggTo4e_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo4tau",filename = pathMC+"ggTo4tau_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo2e2mu",filename = pathMC+"ggTo2e2mu_Contin_MCFM701/ZZ4lAnalysis.root"),       
        dict(name = "ggTo2e2tau",filename = pathMC+"ggTo2e2tau_Contin_MCFM701/ZZ4lAnalysis.root"),
        dict(name = "ggTo2mu2tau",filename = pathMC+"ggTo2mu2tau_Contin_MCFM701/ZZ4lAnalysis.root"),

        dict(name = "ZZTo4l",filename = pathMC+"ZZTo4lext1/ZZ4lAnalysis.root"),
        
        dict(name = "VBF125",filename = pathMC+"VBFH125/ZZ4lAnalysis.root"),
        dict(name = "ggH",filename = pathMC+"ggH125/ZZ4lAnalysis.root"),
        dict(name = "WplusH125",filename = pathMC+"WplusH125/ZZ4lAnalysis.root"),
        dict(name = "WminusH125",filename = pathMC+"WminusH125/ZZ4lAnalysis.root"),
        dict(name = "ZH125",filename = pathMC+"ZH125/ZZ4lAnalysis.root"),
        dict(name = "ttH125",filename = pathMC+"ttH125/ZZ4lAnalysis.root"),
    ]


    of = ROOT.TFile.Open(outFile,"recreate") 
    
    for s in samples:
         histos = fillHistos(s["name"], s["filename"])
         for h in histos:
             of.WriteObject(h,h.GetName())

    of.Close()

def runData():
    outFile = "H4l_Data.root" 

    of = ROOT.TFile.Open(outFile,"recreate") 
                
    hs_data = fillHistos("Data", pathDATA+ "/ZZ4lAnalysis.root")
    for h in hs_data:
        h.SetBinErrorOption(ROOT.TH1.kPoisson)
        of.WriteObject(h,h.GetName())

    of.Close()

if __name__ == "__main__" :
#    runMC()
    runData()
