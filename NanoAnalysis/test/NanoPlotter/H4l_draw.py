### Draw and decorate plots produced with H4l_fill.py.
from __future__ import print_function
import math
import ROOT
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True

inFilename = "H4l_220420.root"
outFilename = "Plots_220420.root"

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
    
    
    EW = [WZZ, ZZZ, VBFToZZTo4l, TTZToLLNuNu, TTZJets]
    h_EW = WWZ.Clone("h_EW")
    for k in EW:
        h_EW.Add(k,1.)
        
    h_EW.SetLineColor(ROOT.TColor.GetColor("#000099"))
    h_EW.SetFillColor(ROOT.TColor.GetColor("#0331B9"))
    

    #------------ggTo-----------------#
    ggTo4mu = f.Get(name+"ggTo4mu") 
    ggTo4e = f.Get(name+"ggTo4e")
    ggTo4tau = f.Get(name+"ggTo4tau")
    ggTo2e2mu = f.Get(name+"ggTo2e2mu")
    ggTo2e2tau = f.Get(name+"ggTo2e2tau")
    ggTo2mu2tau = f.Get(name+"ggTo2mu2tau")

    ggTo = [ ggTo4e, ggTo4tau, ggTo2e2mu, ggTo2e2tau, ggTo2mu2tau]
    h_ggTo = ggTo4mu.Clone("h_ggTo")

    for i in ggTo:
        h_ggTo.Add(i,1.)
    
    h_ggTo.SetLineColor(ROOT.TColor.GetColor("#000099"))
    h_ggTo.SetFillColor(ROOT.TColor.GetColor("#4b78ff"))  

    #-----------qqZZ---------------#
    ZZTo4l = f.Get(name+"ZZTo4l")
    
    ZZTo4l.SetLineColor(ROOT.TColor.GetColor("#000099"))
    ZZTo4l.SetFillColor(ROOT.TColor.GetColor("#99ccff"))
    
    #-----------signal------------#
    VBF125 = f.Get(name+"VBF125")
    ggH = f.Get(name+"ggH")
    WplusH125 = f.Get(name+"WplusH125")
    WminusH125 = f.Get(name+"WminusH125")
    ZH125 = f.Get(name+"ZH125")
    ttH125 = f.Get(name+"ttH125")
    
    signal = [ggH, WplusH125, WminusH125, ZH125, ttH125]
    h_signal = VBF125.Clone("h_signal")
    
    h_signal.SetLineColor(ROOT.TColor.GetColor("#cc0000"))
    h_signal.SetFillColor(ROOT.TColor.GetColor("#ff9b9b"))
    
    for j in signal:
        h_signal.Add(j, 1.)
    
    #------------------Stack----------#
    if version=="_4GeV_" :
        hs = ROOT.THStack("Stack_4GeV", "; m_{#it{4l}} (GeV) ; Events / 4 GeV" )
    elif version=="_10GeV_" :
        hs = ROOT.THStack("Stack_10GeV", "; m_{#it{4l}} (GeV) ; Events / 10 GeV" )
    else:
        hs = ROOT.THStack("Stack_2GeV", "; m_{#it{4l}} (GeV) ; Events / 2 GeV" )
        
    hs.Add(h_EW,"HISTO")
    hs.Add(h_ggTo,"HISTO")
    hs.Add(ZZTo4l,"HISTO")
    hs.Add(h_signal,"HISTO")
    #hs.SetMaximum(160.)
    
    return hs

##################
def StackData (f, version = "_4GeV_"):
    name = "ZZMass"+ version
    Data = f.Get(name+"Data")    
    Data.SetMarkerStyle(20)
    Data.SetLineColor(ROOT.kBlack)
    Data.SetMarkerSize(0.9)
    return Data


f = ROOT.TFile.Open(inFilename,"READ")
of = ROOT.TFile.Open(outFilename,"recreate")




# Labels for log plots
xlabelsv = [80, 100, 200, 300, 400, 500]
label_margin = -1.5
xlabels=[None]*len(xlabelsv)
for i, label in enumerate(xlabelsv): 
    xlabels[i] = ROOT.TLatex(label, label_margin , str(label));
    xlabels[i].SetTextAlign(23)
    xlabels[i].SetTextFont(42)
    xlabels[i].SetTextSize(0.04)


HStack = Stack(f)
HData = StackData(f)
HStack_hm = HStack.Clone()
HData_hm = HData.Clone()

Canvas = ROOT.TCanvas("M4l","M4l",canvasSizeX,canvasSizeY)
Canvas.SetTicks()
Canvas.SetLogx()
#HStacks.SetMaximum(140./1.05)
HStack.Draw("histo")
HData.Draw("samePE1")
# Hide labels and rewrite them
HStack.GetXaxis().SetLabelSize(0)
for label in xlabels :
    label.Draw()

Canvas.Write()
#ROOT.gPad.Modified()
#ROOT.gPad.Update()
    
### Zoomed m4l
HStack_z = Stack(f, "_2GeV_")
HData_z = StackData(f, "_2GeV_")
Canvas_z = ROOT.TCanvas("M4l_z","M4l_z",canvasSizeX,canvasSizeY)
Canvas_z.SetTicks()
#HStack_z.SetMaximum(100./1.05)
HStack_z.Draw("histo")
#HStack_z.GetYaxis().SetRangeUser(0., 100.)
HStack_z.GetXaxis().SetRangeUser(70., 170.)
HData_z.Draw("samePE1")

### Zoomed high mass
Canvas_z = ROOT.TCanvas("M4l_hm","M4l_hm",canvasSizeX,canvasSizeY)
Canvas_z.SetTicks()
Canvas_z.SetLogy()
HStack_hm.Draw("histo")
HStack_hm.GetXaxis().SetRangeUser(170.,1000.)
HData_hm.Draw("samePE1")

### High mass, 10 Ge
HStack10 = Stack(f, "_10GeV_")
HData10 = StackData(f, "_10GeV_")
Canvas10 = ROOT.TCanvas("M4l_hm10","M4l_hm10",canvasSizeX,canvasSizeY)
Canvas10.SetTicks()
Canvas10.SetLogy()
HStack10.Draw("histo")
HStack10.GetXaxis().SetRangeUser(170.,1000.)
HData10.Draw("samePE0E1")

printCanvases()
