### Draw and decorate plots produced with H4l_fill.py.
from __future__ import print_function
import math
import ROOT
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True

inFilenameMC   = "H4l_MC.root"
inFilenameData = "H4l_Data.root"
outFilename = "Plots.root"

### 2018 plots
#Lum = 59.74 # 1/fb
#pathMC = "/eos/user/n/namapane/H4lnano/220420/"
#pathDATA = "/eos/user/n/namapane/H4lnano/220420/Data2018/"

blindPlots=True
blindHLow = 105.
blindHHi  = 140.
blindHM   = 500.


### 2022 plots
#Lum = 7.55 # 1/fb #2022 C-D
Lum = 14.88 # 1/fb #2022 C-F

# Set style matching the one used for HZZ plots
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPadTopMargin(0.05)  
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadLeftMargin(0.16) 
ROOT.gStyle.SetPadRightMargin(0.03)
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetAxisColor(1, "XYZ")
ROOT.gStyle.SetStripDecimals(True)
ROOT.gStyle.SetTickLength(0.03, "XYZ")
ROOT.gStyle.SetNdivisions(510, "XYZ")
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetTitleSize(0.05, "XYZ")
ROOT.gStyle.SetTitleOffset(1.00, "X")
ROOT.gStyle.SetTitleOffset(1.25, "Y")
ROOT.gStyle.SetLabelOffset(0.008, "XYZ")
ROOT.gStyle.SetLabelSize(0.04, "XYZ")

canvasSizeX=910
canvasSizeY=700

#####################
def printCanvases(type="png", path=".") :
    canvases = ROOT.gROOT.GetListOfCanvases()
    for c in canvases :
        c.Print(path+"/"+c.GetTitle()+"."+type)

def printCanvas(c, type="png", name=None, path="." ) :
    if name == None : name = c.GetTitle()
    name=name.replace(">","")
    name=name.replace("<","")
    name=name.replace(" ","_")
    c.Print(path+"/"+name+"."+type)


######################
def Stack (f, version = "_4GeV_"):
    name = "ZZMass" + version

    #------------EW------------------#
    WWZ = f.Get(name+"WWZ")
    WZZ = f.Get(name+"WZZ")
    ZZZ = f.Get(name+"ZZZ")
    VBFToZZTo4l = f.Get(name+"VBFToZZTo4l")
    TTZToLLNuNu = f.Get(name+"TTZToLLNuNu")
    TTZJets = f.Get(name+"TTZJets")

   
    EWSamples = [WZZ, ZZZ, VBFToZZTo4l, TTZToLLNuNu, TTZJets]
    EW = WWZ.Clone("h_EW")
    for i in EWSamples:
        EW.Add(i,1.)
    EW.Scale(Lum*1000.)        
    EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))
    

    #------------ggTo-----------------#
    ggTo4mu = f.Get(name+"ggTo4mu") 
    ggTo4e = f.Get(name+"ggTo4e")
    ggTo4tau = f.Get(name+"ggTo4tau")
    ggTo2e2mu = f.Get(name+"ggTo2e2mu")
    ggTo2e2tau = f.Get(name+"ggTo2e2tau")
    ggTo2mu2tau = f.Get(name+"ggTo2mu2tau")

    ggZZSamples = [ ggTo4e, ggTo4tau, ggTo2e2mu, ggTo2e2tau, ggTo2mu2tau]
    ggToZZ = ggTo4mu.Clone("h_ggTo")
    for i in ggZZSamples:
        ggToZZ.Add(i,1.)
    ggToZZ.Scale(Lum*1000.)    
    ggToZZ.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ggToZZ.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))  

    #-----------qqZZ---------------#
    ZZTo4l = f.Get(name+"ZZTo4l")
    ZZTo4l.Scale(Lum*1000.)    
    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    
    #-----------signal------------#
    VBF125 = f.Get(name+"VBF125")
    ggH = f.Get(name+"ggH")
    WplusH125 = f.Get(name+"WplusH125")
    WminusH125 = f.Get(name+"WminusH125")
    ZH125 = f.Get(name+"ZH125")
    ttH125 = f.Get(name+"ttH125")
    
    signalSamples = [ggH, WplusH125, WminusH125, ZH125, ttH125]
    signal = VBF125.Clone("h_signal")
    
    for i in signalSamples:
        signal.Add(i, 1.)
    signal.Scale(Lum*1000.)    
    signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))

    
    #------------------Stack----------#
    if version=="_4GeV_" :
        hs = ROOT.THStack("Stack_4GeV", "; m_{#it{4l}} (GeV) ; Events / 4 GeV" )
    elif version=="_10GeV_" :
        hs = ROOT.THStack("Stack_10GeV", "; m_{#it{4l}} (GeV) ; Events / 10 GeV" )
    else:
        hs = ROOT.THStack("Stack_2GeV", "; m_{#it{4l}} (GeV) ; Events / 2 GeV" )

    
    hs.Add(EW,"HISTO")
    hs.Add(ggToZZ,"HISTO")
    hs.Add(ZZTo4l,"HISTO")
    hs.Add(signal,"HISTO")
    
    return hs

### Get a TGraph for data, blinded if required
def dataGraph (f, version = "_4GeV_", blind = True):
    name = "ZZMass"+ version
    hd = f.Get(name+"Data")

    nbinsIn = hd.GetNbinsX()
    nbins = 0
    
    x = np.array([0.]*nbinsIn, dtype='double')
    y = np.array([0.]*nbinsIn, dtype='double')
    errX = np.array([0.]*nbinsIn, dtype='double')  
    UpErr = np.array([0.]*nbinsIn, dtype='double')
    LowErr = np.array([0.]*nbinsIn, dtype='double')

    for i in range (1, nbinsIn):
        if blind and ((hd.GetBinCenter(i)>=blindHLow and hd.GetBinCenter(i)<=blindHHi) or hd.GetBinCenter(i)>=blindHM) : continue
        x[nbins] = hd.GetBinCenter(i)
        y[nbins] = hd.GetBinContent(i)
        UpErr[nbins] = hd.GetBinErrorUp(i)
        LowErr[nbins] = hd.GetBinErrorLow(i)
        nbins += 1
    
    Data = ROOT.TGraphAsymmErrors(nbins,x,y,errX,errX,LowErr,UpErr)
    Data.SetMarkerStyle(20)
    Data.SetLineColor(ROOT.kBlack)
    Data.SetMarkerSize(0.9)
    return Data


fMC = ROOT.TFile.Open(inFilenameMC,"READ")
fData = ROOT.TFile.Open(inFilenameData,"READ")
of = ROOT.TFile.Open(outFilename,"recreate")




# Labels for log plots
xlabelsv = [80, 100, 200, 300, 400, 500]
label_margin = -0.4
xlabels=[None]*len(xlabelsv)
for i, label in enumerate(xlabelsv): 
    xlabels[i] = ROOT.TLatex(label, label_margin , str(label));
    xlabels[i].SetTextAlign(23)
    xlabels[i].SetTextFont(42)
    xlabels[i].SetTextSize(0.04)


HStack = Stack(fMC)
HData = dataGraph(fData)
HStack_hm = HStack.Clone()
HData_hm = HData.Clone()

Canvas = ROOT.TCanvas("M4l","M4l",canvasSizeX,canvasSizeY)
Canvas.SetTicks()
Canvas.SetLogx()
#ymaxd=HData.GetMaximum()
xmin=ROOT.Double(0.)
ymin=ROOT.Double(0.)
xmax=ROOT.Double(0.)
ymax=ROOT.Double(0.)
HData.ComputeRange(xmin,ymin,xmax,ymax)
HStack.SetMaximum(int(max(HStack.GetMaximum(), ymax+1)))
HStack.Draw("histo")
HStack.GetXaxis().SetRangeUser(70., 500.)
HData.Draw("samePE1")
# Hide labels and rewrite them
HStack.GetXaxis().SetLabelSize(0)
for label in xlabels :
    label.Draw()
Canvas.Write()
print(HStack.GetMaximum(), HData.GetMaximum()) 
    
### Zoomed m4l
HStack_z = Stack(fMC, "_2GeV_")
HData_z = dataGraph(fData, "_2GeV_")
Canvas_z = ROOT.TCanvas("M4l_z","M4l_z",canvasSizeX,canvasSizeY)
Canvas_z.SetTicks()
HData_z.ComputeRange(xmin,ymin,xmax,ymax)
HStack_z.SetMaximum(int(max(HStack_z.GetMaximum(), ymax+1)))
HStack_z.Draw("histo")
HStack_z.GetXaxis().SetRangeUser(70., 170.)
HData_z.Draw("samePE1")

### Zoomed high mass
Canvas_hm = ROOT.TCanvas("M4l_hm","M4l_hm",canvasSizeX,canvasSizeY)
Canvas_hm.SetTicks()
Canvas_hm.SetLogy()
HStack_hm.Draw("histo")
HStack_hm.GetXaxis().SetRangeUser(170.,1000.)
HData_hm.Draw("samePE1")

### High mass, 10 Ge
HStack10 = Stack(fMC, "_10GeV_")
HData10 = dataGraph(fData, "_10GeV_")
Canvas10 = ROOT.TCanvas("M4l_hm10","M4l_hm10",canvasSizeX,canvasSizeY)
Canvas10.SetTicks()
Canvas10.SetLogy()
HStack10.Draw("histo")
HStack10.GetXaxis().SetRangeUser(170.,1000.)
HData10.Draw("samePE0E1")

printCanvases()
