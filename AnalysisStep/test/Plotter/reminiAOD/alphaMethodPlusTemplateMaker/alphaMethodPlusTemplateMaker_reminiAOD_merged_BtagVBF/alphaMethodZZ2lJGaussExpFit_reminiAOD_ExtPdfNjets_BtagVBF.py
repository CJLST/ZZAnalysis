
import ROOT, os, string
import math
from math import *

from ROOT import *
from array import array

import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re

from tdrStyle import *
setTDRStyle()

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)

def parseOptions():
    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    # fit opt
    parser.add_option("--fitOpt", dest="fitOpt",default='MaxLikelihood',help="ML or Chi2 fit")
    #parser.add_option("--doVariableBinning", action="store_true", dest="doVariableBinning", default=False, help="whether use variable binning")
    # sample opt
    parser.add_option("--dyMC", dest="dyMC",default='njetLoDY',help="Type of DY sample")
    parser.add_option("--isRefitSR",action="store_true", dest="isRefitSR", default=True,help="whether use refitted mZZ for SR data/mc")

    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)

    #(option, args) = parser.parse_args()
    return parser 

def ZJetsIntegral(minMass,maxMass):
    return ZJetsNorm(minMass,maxMass)

def ZJetsGaussExpIntegral(minMass,maxMass,mean,sigma,l0):

    # integral of Zjets
    prescale = str(pow(math.atan(1)*2,0.5)) #TMath::Sqrt(PI/2)=str(pow(math.atan(1)*2,0.5))
    prescale = prescale+"*sigma*TMath::Exp(p0-0.5*l0*l0)"

    if(minMass<mean and maxMass<mean)  :

      left = str(pow(0.5,0.5))+"*(mean-"+str(minMass)+")/sigma"
      right = str(pow(0.5,0.5))+"*(mean-"+str(maxMass)+")/sigma"
      gauss_left=prescale+"*TMath::Erf("+left+")"
      gauss_right=prescale+"*TMath::Erf("+right+")"
      return gauss_left+'-'+gauss_right

    elif(minMass<mean and maxMass>=mean and maxMass<(mean+sigma*l0)) :

      left = str(pow(0.5,0.5))+"*(mean-"+str(minMass)+")/sigma"
      right = str(pow(0.5,0.5))+"*("+str(maxMass)+"-mean)/sigma"
      gauss_left=prescale+"*TMath::Erf("+left+")"
      gauss_right=prescale+"*TMath::Erf("+right+")"
      return gauss_left+'+'+gauss_right

    elif(minMass<mean and maxMass>(mean+sigma*l0)) :

      left = str(pow(0.5,0.5))+"*(mean-"+str(minMass)+")/sigma"
      right = str(pow(0.5,0.5))+"*"+str(l0)
      gauss_left=prescale+"*TMath::Erf("+left+")"
      gauss_right=prescale+"*TMath::Erf("+right+")"
      exp = "(sigma/l0)*(TMath::Exp(p0-l0*l0)-TMath::Exp(p0-l0*("+str(maxMass)+"-mean)/sigma))"
      return gauss_left+'+'+gauss_right+'+'+exp

    elif(minMass>=mean and minMass<(mean+sigma*l0) and maxMass>mean and maxMass<(mean+sigma*l0)) :

      left = str(pow(0.5,0.5))+"*("+str(minMass)+"-mean)/sigma"
      right = str(pow(0.5,0.5))+"*("+str(maxMass)+"-mean)/sigma"
      gauss_left=prescale+"*TMath::Erf("+left+")"
      gauss_right=prescale+"*TMath::Erf("+right+")"

      return gauss_right+'-'+gauss_left

    elif(minMass>=mean and minMass<(mean+sigma*l0) and maxMass>=(mean+sigma*l0)) :

      left = str(pow(0.5,0.5))+"*("+str(minMass)+"-mean)/sigma"
      right = str(pow(0.5,0.5))+"*"+str(l0)
      gauss_left=prescale+"*TMath::Erf("+left+")"
      gauss_right=prescale+"*TMath::Erf("+right+")"
      exp = "(sigma/l0)*(TMath::Exp(p0-l0*l0)-TMath::Exp(p0-l0*("+str(maxMass)+"-mean)/sigma))"
      return gauss_right+'-'+gauss_left+'+'+exp

    elif(minMass>=(mean+sigma*l0) and maxMass>=(mean+sigma*l0)) :

      return "(sigma/l0)*(TMath::Exp(p0-l0*("+str(minMass)+"-mean)/sigma)-TMath::Exp(p0-l0*("+str(maxMass)+"-mean)/sigma))"

#alphaMethod(fs,jet,cat, parser)
def alphaMethod(fs,jet,cat, minMass, parser):
    
    (options,args) = parser.parse_args()

    print "Option setting"
    fitOpt = options.fitOpt
    
    dyMC = 'njetLoDY'
    print 'Use ',dyMC
    isRefitSR = False
    doVariableBinning = False # ML needs equal binning

    if(fitOpt=="chi2" or fitOpt=="Chi2") :
      doVariableBinning = True # chi2 needs varaible binning
    if(doVariableBinning) :
      fitOpt+="VariableBinning"
    else :
      fitOpt+="EqualBinning"

    print 'Fit options is ',fitOpt

    inputfile = TFile('TMVAAndRoofitInputs_all_'+dyMC+'.root',"READ")

    bins=[]
    bins_alpha=[]
    bins_alpha_ALT=[]

    # default binning 5GeV
    nRebin=1
    if(cat!='btag' and cat!='vbf') :
      nRebin=10
      for i in range(0,8):
        bins_alpha.append(500+i*50)
        
      bins_alpha.append(900)
      bins_alpha.append(1000)
      bins_alpha.append(1100)
      bins_alpha.append(1300)
      bins_alpha.append(1500)
      bins_alpha.append(2000)
      bins_alpha.append(3500)

      bins = bins_alpha
      bins_alpha_ALT = bins_alpha

    if(cat=='btag') :
      nRebin=10
      bins_alpha=[400,500,600,700,800,1000,1200,3500]

      bins = bins_alpha
      bins_alpha_ALT = bins_alpha

    if(cat=='vbf') :
      nRebin=10
      bins_alpha=[400,500,600,700,800,1000,1200,3500]

      bins = bins_alpha
      bins_alpha_ALT = bins_alpha

    nbins = len(bins)-1

    print "mVV spectrum"
    Shape_SB = {}
    Shape_SR = {}

    jet = "merged"

    processes = ['Data','DY','TT','VV']
    for process in processes :

       Shape_SB[process]=(inputfile.Get("hmass_"+jet+"SB"+cat+"_"+process)).Clone()

       if(isRefitSR):
         Shape_SR[process]=(inputfile.Get("hmassRefit_"+jet+"SR"+cat+"_"+process)).Clone()
       else:
         Shape_SR[process]=(inputfile.Get("hmass_"+jet+"SR"+cat+"_"+process)).Clone()

       Shape_SB[process].Rebin(nRebin)
       Shape_SR[process].Rebin(nRebin)

    vzSB_norm = Shape_SB['VV'].Integral()
    vzSR_norm = Shape_SR['VV'].Integral()
    ttbarSB_norm = Shape_SB['TT'].Integral()
    ttbarSR_norm = Shape_SR['TT'].Integral()

    print "alpha factor"
    alpha=alphaFactor(Shape_SR['DY'],Shape_SB['DY'],bins_alpha)
    alpha.SetName("alpha")
    alpha_ALT=alphaFactor(Shape_SR['DY'],Shape_SB['DY'],bins_alpha_ALT)
    alpha_ALT.SetName("alpha_ALT")

    print "alpha uncertainty"
    alphaUp = alpha.Clone()
    alphaUp.SetName('alphaUp'+cat)
    alphaDn = alpha.Clone()
    alphaDn.SetName('alphaDn'+cat)

    for i in range(1,alpha.GetXaxis().GetNbins()+1):
        alphaUp.SetBinContent(i,alpha.GetBinContent(i)+alpha.GetBinError(i))
        alphaDn.SetBinContent(i,max(0.0001,alpha.GetBinContent(i)-alpha.GetBinError(i)) )
        '''
        errorJEC = 0.5*abs(alpha_JECup.GetBinContent(i)-alpha_JECdn.GetBinContent(i))
        errorBinning = abs(alpha_ALT.GetBinContent(i)-alpha.GetBinContent(i))
        binerror = pow(alpha.GetBinError(i)*alpha.GetBinError(i)+errorJEC*errorJEC+errorBinning*errorBinning,0.5)
        alpha.SetBinError(i,binerror)
        '''

    c = TCanvas("c","c",10,10,1000,800)
    c.cd()
    alpha.SetFillColor(kGreen+3)
    alpha.SetMaximum(1.5)
    alpha.SetMinimum(0.0)
    alpha.Draw("ep2hist")
    alpha_ALT.SetLineColor(kGreen+2)
    alphaDn.SetLineColor(kRed)
    alphaUp.SetLineColor(kRed)
    alphaDn.SetFillColor(0)
    alphaUp.SetFillColor(0)
    alphaDn.Draw("histsame")
    alphaUp.Draw("histsame")
    alpha_ALT.Draw("histsame")
    if(isRefitSR):
      c.SaveAs("alpha"+jet+cat+"_"+dyMC+"_RefitSR.png")
    else :
      c.SaveAs("alpha"+jet+cat+"_"+dyMC+"_RecoSR.png")

    # SR data
    dataSR = Shape_SR['Data'].Clone()
    dataSR.SetName("dataSR")
    # CR data
    dataCR = Shape_SB['Data'].Clone()
    dataCR.Multiply(alpha)
    dataCR.SetName("alphaXdataCR")
    #############
    dataCRup = Shape_SB['Data'].Clone()
    dataCRup.Multiply(alphaUp)
    dataCRup.SetName("alphaUpXdataCR")
    ########
    dataCRdn = Shape_SB['Data'].Clone()
    dataCRdn.Multiply(alphaDn)
    dataCRdn.SetName("alphaDnXdataCR")
    alphaUnc = 1-dataCRdn.Integral()/dataCR.Integral()
    print 'SB data X alpha up/dn is ',alphaUnc

    # SR ttbar/vz
    ttbarSR = Shape_SR['TT'].Clone()
    vzSR = Shape_SR['VV'].Clone()
    vzSR.Scale(1.1)

    # dataCR-TTCR as dySR estimation from raw hist
    ttbarCR = Shape_SB['TT'].Clone()
    ttbarCR.Multiply(alpha)
    ttbarCR.SetName("alphaXttbarCR")
    vzCR = Shape_SB['VV'].Clone()
    vzCR.Multiply(alpha)
    vzCR.SetName("alphaXvzCR")
    vzCR.Scale(1.1)

    dataSR_variableBinning = rebinMerge(dataSR,bins)
    dataCR_variableBinning = rebinMerge(dataCR,bins)
    ttbarCR_variableBinning = rebinMerge(ttbarCR,bins)
    ttbarSR_variableBinning = rebinMerge(ttbarSR,bins)
    vzCR_variableBinning = rebinMerge(vzCR,bins)
    vzSR_variableBinning = rebinMerge(vzSR,bins)

    print 'Debug dataCR ',dataCR.Integral(),' dataCR_variableBinning ',dataCR_variableBinning.Integral()

    ## Z+jets as alpha X (dataSB-ttbarVZSB)
    dataCR_variableBinning.Add(ttbarCR_variableBinning,-1)
    dataCR_variableBinning.Add(vzCR_variableBinning,-1)
    dataCR.Add(ttbarCR,-1)
    dataCR.Add(vzCR,-1)
    if(isRefitSR) :
      plotDataVSbkg(dataSR_variableBinning,dataCR_variableBinning,ttbarSR_variableBinning,vzSR_variableBinning,"SRNoFitting_"+jet+cat+fs+"_"+dyMC+"_Refit",500)
    else :
      plotDataVSbkg(dataSR_variableBinning,dataCR_variableBinning,ttbarSR_variableBinning,vzSR_variableBinning,"SRNoFitting_"+jet+cat+fs+"_"+dyMC+"_Reco",500)

    ## Add ttbarVZSB back to dataCR for fitting dataset
    dataCR_variableBinning.Add(ttbarCR_variableBinning)
    dataCR_variableBinning.Add(vzCR_variableBinning)
    dataCR.Add(ttbarCR)
    dataCR.Add(vzCR)

    #obserble
    mass=RooRealVar("mass","mass", minMass,3500)
    #parameters
    mean=RooRealVar("mean", "mean", 700, 600, 800)
    sigma=RooRealVar("sigma", "sigma", 50, 30, 250)
    l0=RooRealVar("l0", "l0", 0.1, 0.05, 1.0)

    arg="((mass-mean)/sigma)";
    mathExpression1 = "TMath::Exp(-0.5*l0*l0-0.5*"+arg+"*"+arg+")*(0.5+TMath::Sign(0.5,l0-"+arg+"))";
    mathExpression2 = "TMath::Exp(-l0*"+arg+")*(0.5-TMath::Sign(0.5,l0-"+arg+"))";
    mathExpression = mathExpression1+"+"+mathExpression2;

    zjets = RooGenericPdf("zjets",mathExpression,RooArgList(mass,l0,mean,sigma))
    Nexp = dataCR.Integral(dataCR.FindBin(minMass),dataCR.GetXaxis().GetNbins())
    Nexp = Nexp - ttbarCR.Integral(dataCR.FindBin(minMass),dataCR.GetXaxis().GetNbins())
    Nexp = Nexp - vzCR.Integral(dataCR.FindBin(minMass),dataCR.GetXaxis().GetNbins())
    print "Nexp, ",Nexp
    Nzjets = RooRealVar("Nzjets","Nzjets",Nexp,0.5*Nexp, 5*Nexp)

    #p0=RooRealVar("p0", "p0", log(Nexp/20),-50,50)
    #rfvSigRate_zjets = ROOT.RooFormulaVar("zjets_norm",ZJetsNorm(500,3500),ROOT.RooArgList(mean,sigma,l0,p0) )
    zjets_ext = RooExtendPdf("zjets_ext","zjets_ext",zjets,Nzjets)

    print "TT+VV Pdf"
    ttbarVZCR=ttbarCR.Clone()
    ttbarVZCR.SetName("ttbarVZCR")
    ttbarVZCR.Add(vzCR)
    print "TT+VV Pdf dataHist"
    ttbarVZCR_templates = ROOT.RooDataHist("ttbarVZCRdataHist","ttbarVZCRdataHist",RooArgList(mass),RooFit.Import(Shape_SB['TT'],kFALSE))
    print "TT+VV Pdf dataHistPdf"
    ttbarVZCRPdf = ROOT.RooHistPdf("ttbarVZCRPdf","ttbarVZCRPdf",RooArgSet(mass),ttbarVZCR_templates)
    print "TT+VV Pdf NttbarVZCR"
    NttbarVZCR = RooRealVar("NttbarVZCR","NttbarVZCR",ttbarVZCR.Integral())
    print "TT+VV Pdf EXT dataHistPdf"
    ttbarVZCR_ext = RooExtendPdf("ttbarVZCR_ext","ttbarVZCR_ext",ttbarVZCRPdf,NttbarVZCR)

    print "sum of DY+(TT+VV) extended PDF"
    sum_ext = RooAddPdf("sum","DY+ttbarVZ",RooArgList(zjets,ttbarVZCRPdf),RooArgList(Nzjets,NttbarVZCR))

    print 'load the dataset'
    dataCRHist = RooDataHist("dataCRdataHist","alpha-times-sbData",RooArgList(mass),RooFit.Import(dataCR,kFALSE))
    dataCRHist_variableBinning = RooDataHist("dataCRdataHist_variableBinning","alpha-times-sbData",RooArgList(mass),RooFit.Import(dataCR_variableBinning,kFALSE))

    print "set fit options in Roofit"
    currentlist = RooLinkedList()
    cmd = ROOT.RooFit.Save()
    currentlist.Add(cmd)
    cmd1 = ROOT.RooFit.Extended()
    currentlist.Add(cmd1)
    cmd2 = ROOT.RooFit.PrintLevel(-1)
    currentlist.Add(cmd2)

    print "perform fit with ",fitOpt
    r = RooFitResult()
    #print 'fitOpt.find("MaxLikelihood") ',fitOpt.find("MaxLikelihood")
    if(fitOpt.find("MaxLikelihood")==0):
       print 'fit with ML equalbinning'
       cmd3 = ROOT.RooFit.SumW2Error(kTRUE)
       currentlist.Add(cmd3)
       r = sum_ext.fitTo(dataCRHist,currentlist)
       r = sum_ext.fitTo(dataCRHist,currentlist)
       r = sum_ext.fitTo(dataCRHist,currentlist)
    elif(fitOpt.find('Chi2')==0 or fitOpt.find('chi2')==0) :
       print 'fit with Chi2 variablebinning'
       cmd3 = ROOT.RooFit.DataError(ROOT.RooDataHist.Poisson)
       currentlist.Add(cmd3)
       r = sum_ext.chi2FitTo(dataCRHist_variableBinning,currentlist)
       r = sum_ext.chi2FitTo(dataCRHist_variableBinning,currentlist)
       r = sum_ext.chi2FitTo(dataCRHist_variableBinning,currentlist)

    #cmd3 = ROOT.RooFit.Strategy(2)
    #currentlist.Add(cmd3)
    #cmd4 = ROOT.RooFit.Minimizer("Minuit2","migrad")
    #currentlist.Add(cmd4)

    print 'Nzjets ',Nzjets.getVal(), ' vs ',Nexp

    czz = TCanvas("c", "c", 800, 800 );
    czz.cd()
   
    czz.SetLogy()

    frame = mass.frame(RooFit.Title("alpha method for Zjets"))#,RooFit.Bins(30))
    dataCRHist.plotOn(frame)
    sum_ext.paramOn(frame)
    sum_ext.plotOn(frame,RooFit.VisualizeError(r,1,kTRUE),RooFit.FillColor(kYellow+2))
    sum_ext.plotOn(frame,RooFit.LineColor(kGreen))
    dataCRHist.plotOn(frame)
    #chi2 = frame.chiSquare(3) # 3 is number of DOF
    #print 'chi2 ',chi2
    frame.SetMinimum(0.01)
    frame.Draw()
    czz.SaveAs("RooFit_"+jet+cat+"_"+fitOpt+"_"+dyMC+".png")


    k = Nexp/Nzjets.getVal()
    print 'k, correction to normalization in the fit is ',k
    k = 1.0

    zjetsCR_variableBinning = dataCR_variableBinning.Clone()
    argset = RooArgSet(mass)
    for i in range(1,zjetsCR_variableBinning.GetXaxis().GetNbins()+1):

        low = zjetsCR_variableBinning.GetXaxis().GetBinLowEdge(i)
        high = zjetsCR_variableBinning.GetXaxis().GetBinUpEdge(i)
        mass.setRange(str(low)+'_'+str(high),low,high)
        DY_integral_i = zjets.createIntegral(argset,RooFit.NormSet(argset),RooFit.Range(str(low)+'_'+str(high)))
        rfvSigRate_zjets_i = ROOT.RooFormulaVar("zjets_norm_"+str(i),"@0*@1",ROOT.RooArgList(Nzjets,DY_integral_i) )

        zjetsCR_variableBinning.SetBinContent(i,k*rfvSigRate_zjets_i.getVal())
        zjetsCR_variableBinning.SetBinError(i,k*rfvSigRate_zjets_i.getPropagatedError(r))

    print 'plot fitting results'
    plotDataVSbkg(dataCR_variableBinning,zjetsCR_variableBinning,ttbarCR_variableBinning, vzCR_variableBinning,"CRFitting"+jet+cat+fs+"_"+fitOpt+"_"+dyMC,minMass)
    plotDataVSbkg(dataSR_variableBinning,zjetsCR_variableBinning,ttbarSR_variableBinning,vzSR_variableBinning,"SRFitting"+jet+cat+fs+"_"+fitOpt+"_"+dyMC,minMass)

    print "Nzjets ",zjetsCR_variableBinning.Integral()," compare with raw ",dataSR.Integral()-Shape_SR['TT'].Integral()-Shape_SR['VV'].Integral()

    ##############################################################

    print 'Save mZZ templates and ALT templates'

    l0_nominal_str = str(l0.getVal())
    m_nominal_str = str(mean.getVal())
    s_nominal_str = str(sigma.getVal())

    l0l0=r.correlation("l0","l0")*l0.getError()*l0.getError()
    l0m=r.correlation("l0","mean")*l0.getError()*mean.getError()
    l0s=r.correlation("l0","sigma")*l0.getError()*sigma.getError()

    ml0=r.correlation("mean","l0")*mean.getError()*l0.getError()
    mm=r.correlation("mean","mean")*mean.getError()*mean.getError()
    ms=r.correlation("mean","sigma")*mean.getError()*sigma.getError()

    sl0=r.correlation("sigma","l0")*sigma.getError()*l0.getError()
    sm=r.correlation("sigma","mean")*sigma.getError()*mean.getError()
    ss=r.correlation("sigma","sigma")*sigma.getError()*sigma.getError()

    cov_elements = [l0l0,l0m,l0s,
                    ml0,mm,ms,
                    sl0,sm,ss]

    n = 3
    cov = TMatrixDSym(n,array('d',cov_elements))
    eigen = TMatrixDSymEigen(cov)
    vecs = eigen.GetEigenVectors()
    vals  = eigen.GetEigenValues()

    eigCoeff0 = []
    eigCoeff1 = []
    eigCoeff2 = []

    eig0Name="eig0"
    eig1Name="eig1"
    eig2Name="eig2"

    print 'eigen vals ',vals(0),' ',vals(1),' ',vals(2)

    eigCoeff0 = [vecs(0,0)*pow(vals(0),0.5),vecs(0,1)*pow(vals(1),0.5),vecs(0,2)*pow(vals(2),0.5)]
    eigCoeff1 = [vecs(1,0)*pow(vals(0),0.5),vecs(1,1)*pow(vals(1),0.5),vecs(1,2)*pow(vals(2),0.5)]
    eigCoeff2 = [vecs(2,0)*pow(vals(0),0.5),vecs(2,1)*pow(vals(1),0.5),vecs(2,2)*pow(vals(2),0.5)]

    l0_unc_str = "( ("+str(eigCoeff0[0])+")*"+eig0Name+"+("+str(eigCoeff0[1])+")*"+eig1Name+"+("+str(eigCoeff0[2])+")*"+eig2Name+")"
    ######
    m_unc_str = "( ("+str(eigCoeff1[0])+")*"+eig0Name+"+("+str(eigCoeff1[1])+")*"+eig1Name+"+("+str(eigCoeff1[2])+")*"+eig2Name+")"
    ######
    s_unc_str = "( ("+str(eigCoeff2[0])+")*"+eig0Name+"+("+str(eigCoeff2[1])+")*"+eig1Name+"+("+str(eigCoeff2[2])+")*"+eig2Name+")"
    ######

    l0_str = "("+l0_nominal_str+"+"+l0_unc_str+")"
    ######
    m_str = "("+m_nominal_str+"+"+m_unc_str+")"
    ######
    s_str = "("+s_nominal_str+"+"+s_unc_str+")"

    #### original template to start with
    print l0_str
    print m_str
    print s_str


    minMass = 500
    zjetsAlphaSR = TH1D("zjetsAlphaSR","zjetsAlphaSR",300,500,3500)
    zjetsTemplates = {}
    zjetsTemplates['nominal'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,0,0)

    zjetsTemplates['eig0up'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,0,1)
    zjetsTemplates['eig0dn'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,0,-1)

    zjetsTemplates['eig1up'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,1,1)
    zjetsTemplates['eig1dn'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,1,-1)

    zjetsTemplates['eig2up'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,2,1)
    zjetsTemplates['eig2dn'] = templatesMaker_forTail(minMass,zjetsAlphaSR,mass,Nzjets,l0_str,m_str,s_str,2,-1)

    print 'Calculate branching ratios'
    Shape_SB_ee = {}
    Shape_SB_mumu = {}

    inputfile_fs = {} 
    inputfile_fs['all'] = TFile('TMVAAndRoofitInputs_all_'+dyMC+'.root',"READ")
    inputfile_fs['ee'] = TFile('TMVAAndRoofitInputs_ee_'+dyMC+'.root',"READ")
    inputfile_fs['mumu'] = TFile('TMVAAndRoofitInputs_mumu_'+dyMC+'.root',"READ")
    processes = ['Data','DY','TT','VV']
    for process in processes :
       Shape_SB_ee[process]=(inputfile_fs['ee'].Get("hmass_"+jet+"SB"+cat+"_"+process)).Clone()
       Shape_SB_ee[process].SetName("hmass_"+jet+"ee"+"SB"+cat+"_"+process)
       Shape_SB_mumu[process]=(inputfile_fs['mumu'].Get("hmass_"+jet+"SB"+cat+"_"+process)).Clone()    
       Shape_SB_mumu[process].SetName("hmass_"+jet+"mumu"+"SB"+cat+"_"+process)
       Shape_SB_ee[process].Rebin(10)
       Shape_SB_mumu[process].Rebin(10)

    dataCR_ee = Shape_SB_ee['Data'].Clone()
    dataCR_ee.SetName("alphaXdataCR_ee")
    dataCR_ee.Multiply(alpha)
    dataCR_mumu = Shape_SB_mumu['Data'].Clone()
    dataCR_mumu.SetName("alphaXdataCR_mumu")
    dataCR_mumu.Multiply(alpha)

    ttbarVZCR_ee = Shape_SB_ee['TT'].Clone()
    ttbarVZCR_ee.Add(Shape_SB_ee['VV'])
    ttbarVZCR_ee.SetName("alphaXdataCR_ee")
    ttbarVZCR_ee.Multiply(alpha)
    ttbarVZCR_mumu = Shape_SB_mumu['TT'].Clone()
    ttbarVZCR_mumu.Add(Shape_SB_mumu['VV'])
    ttbarVZCR_mumu.SetName("alphaXdataCR_mumu")
    ttbarVZCR_mumu.Multiply(alpha)

    Total = (dataCR_mumu.Integral()-ttbarVZCR_mumu.Integral()+dataCR_ee.Integral()-ttbarVZCR_ee.Integral())
    BR_mumu = (dataCR_mumu.Integral()-ttbarVZCR_mumu.Integral())/Total
    BR_ee = (dataCR_ee.Integral()-ttbarVZCR_ee.Integral())/Total

    BR = 1.0
    if(fs=="ee") :
      BR = BR_ee
    if(fs=="mumu") :
      BR = BR_mumu

    print 'BR ',BR,' ',cat,' ',fs
    zjetsTemplates['nominal'].Scale(BR)
    zjetsTemplates['eig0up'].Scale(BR)
    zjetsTemplates['eig0dn'].Scale(BR)
    zjetsTemplates['eig1up'].Scale(BR)
    zjetsTemplates['eig1dn'].Scale(BR)
    zjetsTemplates['eig2up'].Scale(BR)
    zjetsTemplates['eig2dn'].Scale(BR)

    path="../alphaMethodPlusTemplateMaker_reminiAOD_merged/"
    inputfile_fs['all'] = TFile(path+'TMVAAndRoofitInputs_all_'+dyMC+'.root',"READ")
    inputfile_fs['ee'] = TFile(path+'TMVAAndRoofitInputs_ee_'+dyMC+'.root',"READ")
    inputfile_fs['mumu'] = TFile(path+'TMVAAndRoofitInputs_mumu_'+dyMC+'.root',"READ")

    ttbarSR_ws=(inputfile_fs['all'].Get("hmass_"+jet+"SR_TT")).Clone()
    ttbarSR_ws.SetName("ttbar_50GeV")
    ttbarSR_norm_fs=inputfile_fs[fs].Get("hmass_"+jet+"SR"+cat+"_TT").Integral()
    ttbarSR_ws.Rebin(10)
    ttbarSR_ws.Smooth(25)
    ttbarSR_ws.Scale(ttbarSR_norm_fs/ttbarSR_ws.Integral())
    vzSR_ws=(inputfile_fs['all'].Get("hmass_"+jet+"SR_VV")).Clone()
    vzSR_ws.SetName("vz_50GeV")
    vzSR_norm_fs=inputfile_fs[fs].Get("hmass_"+jet+"SR"+cat+"_VV").Integral()
    vzSR_ws.Rebin(10)
    vzSR_ws.Smooth(10)
    vzSR_ws.Scale(ttbarSR_norm_fs/vzSR_ws.Integral())
    vzSR_ws.Scale(1.1)

    for i in range(1,ttbarSR_ws.GetXaxis().GetNbins()+1):
        bincontent = ttbarSR_ws.GetBinContent(i)
        if(bincontent<=0) :
          ttbarSR_ws.SetBinContent(i,1e-10)
    for i in range(1,vzSR_ws.GetXaxis().GetNbins()+1):
        bincontent = vzSR_ws.GetBinContent(i)
        if(bincontent<=0) :
          vzSR_ws.SetBinContent(i,1e-10)

    vznew=TH1F("vz","vz",300,500,3500)
    ttbarnew=TH1F("ttbar","ttbar",300,500,3500)

    for i in range(0,vznew.GetNbinsX()) :
        for j in range(0,5) :
            vznew.SetBinContent(i*5+1+j, vzSR_ws.GetBinContent(i+1)/5.)
            ttbarnew.SetBinContent(i*5+1+j, ttbarSR_ws.GetBinContent(i+1)/5.)

    print 'ttbarSR_ws ',ttbarSR_ws.Integral(),' ',fs,' ',cat
    print 'vzSR_ws ',vzSR_ws.Integral(),' ',fs,' ',cat

    dataSR_obs=(inputfile_fs[fs].Get("hmass_"+jet+"SR"+cat+"_Data")).Clone()
    dataSR_obs.SetName("data_obs")
    dataSR_obs.Rebin(2)

    cat_prime=cat
    if(cat==''):
      cat_prime='untagged'
    if(cat=='btag'):
      cat_prime='btagged'
    if(cat=='vbf'):
      cat_prime='vbftagged'
    print 'Write templates to root file'
    f=TFile('mzzTemplate_'+dyMC+'_'+fs+'qq_'+jet+'_'+cat_prime+'.root','RECREATE')
    f.cd()

    zjetsTemplates['nominal'].Write()

    zjetsTemplates['eig0up'].Write()
    zjetsTemplates['eig0dn'].Write()
    zjetsTemplates['eig1up'].Write()
    zjetsTemplates['eig1dn'].Write()
    zjetsTemplates['eig2up'].Write()
    zjetsTemplates['eig2dn'].Write()

    print 'zjets ',zjetsTemplates['nominal'].Integral()
    zjetsTemplates['eig3up']=zjetsTemplates['nominal'].Clone()
    zjetsTemplates['eig3up'].SetName("zjets_eig3Up")
    zjetsTemplates['eig3up'].Scale(1+pow(alphaUnc*alphaUnc+pow(Nzjets.getError()/Nzjets.getVal(),2),0.5))
    zjetsTemplates['eig3up'].Write()
    print 'zjetseig3up ',zjetsTemplates['eig3up'].Integral()

    zjetsTemplates['eig3dn']=zjetsTemplates['nominal'].Clone()
    zjetsTemplates['eig3dn'].SetName("zjets_eig3Dn")
    zjetsTemplates['eig3dn'].Scale(1-pow(alphaUnc*alphaUnc+pow(Nzjets.getError()/Nzjets.getVal(),2),0.5))
    zjetsTemplates['eig3dn'].Write()
    print 'zjetseig3dn ',zjetsTemplates['eig3dn'].Integral()

    print 'alphaUnc ',alphaUnc,' Nzjets.getError()/Nzjets.getVal() ',Nzjets.getError()/Nzjets.getVal()
    print '1+pow(alphaUnc*alphaUnc+pow(Nzjets.getError()/Nzjets.getVal(),2),0.5) ',1+pow(alphaUnc*alphaUnc+pow(Nzjets.getError()/Nzjets.getVal(),2),0.5)
    print '1-pow(alphaUnc*alphaUnc+pow(Nzjets.getError()/Nzjets.getVal(),2),0.5) ',1-pow(alphaUnc*alphaUnc+pow(Nzjets.getError()/Nzjets.getVal(),2),0.5)

    ttbarnew.Write()
    vznew.Write()

    dataSR_obs.Write()

    f.Close()

    print 'making final plots'
    dataSR_plot=dataSR_obs.Clone()
    dataSR_plot.SetName("dataSR_plot")
    dataSR_plot.Rebin(5)

    zjetsTemplates_nominal=zjetsTemplates['nominal'].Clone()
    zjetsTemplates_nominal.SetName("zjetsTemplates_nominal")
    zjetsTemplates_nominal.Rebin(5)

    if(isRefitSR) :
      plotDataVSbkg(dataSR_plot,zjetsTemplates_nominal,ttbarSR_ws,vzSR_ws,"SRTemplate_"+jet+cat+fs+"_"+dyMC+"_Refit",500)
    else :
      plotDataVSbkg(dataSR_plot,zjetsTemplates_nominal,ttbarSR_ws,vzSR_ws,"SRTemplate_"+jet+cat+fs+"_"+dyMC+"_Reco",500)


def templatesMaker_forTail(minMass,zjetsCR,mass,Nzjets,l0_str,m_str,s_str,index,nSigma) :

    zjetsCRtmp = zjetsCR.Clone()
    postfix = ""
    if(nSigma>0) :
      postfix = "Up"
    if(nSigma<0) :
      postfix = "Dn"
    zjetsCRtmp.SetName("zjets_eig"+str(index)+postfix)
    if(nSigma==0) :
      zjetsCRtmp.SetName("zjets")

    eig0Name="eig0"
    eig1Name="eig1"
    eig2Name="eig2"

    eig0 = RooRealVar(eig0Name,eig0Name,0,-3,3)
    eig1 = RooRealVar(eig1Name,eig1Name,0,-3,3)
    eig2 = RooRealVar(eig2Name,eig2Name,0,-3,3)

    if(index==0) :
      eig0.setVal(nSigma)
      eig1.setVal(0)
      eig2.setVal(0)
    if(index==1) :
      eig0.setVal(0)
      eig1.setVal(nSigma)
      eig2.setVal(0)
    if(index==2) :
      eig0.setVal(0)
      eig1.setVal(0)
      eig2.setVal(nSigma)

    ####
    arg="((mass-"+m_str+")/"+s_str+")";
    mathExpression1 = "TMath::Exp(-0.5*"+l0_str+"*"+l0_str+"-0.5*"+arg+"*"+arg+")*(0.5+TMath::Sign(0.5,"+l0_str+"-"+arg+"))";
    mathExpression2 = "TMath::Exp(-"+l0_str+"*"+arg+")*(0.5-TMath::Sign(0.5,"+l0_str+"-"+arg+"))";
    mathExpressionDiagonal = mathExpression1+"+"+mathExpression2;

    zjetsDiagnol = RooGenericPdf("zjetsDiagonal",mathExpressionDiagonal,RooArgList(eig0,eig1,eig2,mass))

    argset = RooArgSet(mass)
    for i in range(1,zjetsCRtmp.GetXaxis().GetNbins()+1):

        low = zjetsCR.GetXaxis().GetBinLowEdge(i)
        high = zjetsCR.GetXaxis().GetBinUpEdge(i)

        mass.setRange(str(low)+'_'+str(high),low,high)
        DY_integral_i = zjetsDiagnol.createIntegral(argset,RooFit.NormSet(argset),RooFit.Range(str(low)+'_'+str(high)))
        rfvSigRateZjetsDiagnol_i = ROOT.RooFormulaVar("zjets_norm_"+str(i),"@0*@1",ROOT.RooArgList(Nzjets,DY_integral_i) )

        zjetsCRtmp.SetBinContent(i,rfvSigRateZjetsDiagnol_i.getVal())
        zjetsCRtmp.SetBinError(i,0)

    return zjetsCRtmp

def rebinMerge(hist,bins) :

    Nbins = hist.GetXaxis().GetNbins()
    TH1_tmp = ROOT.TH1D(hist.GetName()+"_rebinMerge", "tmp", len(bins)-1, array('d',bins))

    print 'mergebin for hist ',hist.GetName()
    for i in range(1,len(bins)):
      lower_edge = TH1_tmp.GetXaxis().GetBinLowEdge(i)
      higher_edege = TH1_tmp.GetXaxis().GetBinUpEdge(i)
      bincontent = 0.0
      binerror = 0.0

      for j in range(1,Nbins+1):
        bincenter = hist.GetXaxis().GetBinCenter(j)

        if(bincenter>=lower_edge and bincenter<higher_edege) :
          #print 'low ',lower_edge,' center ',bincenter,' high ',higher_edege,' content ',hist.GetBinContent(j)
          bincontent = bincontent + hist.GetBinContent(j)
          binerror = binerror + hist.GetBinError(j)*hist.GetBinError(j)
        else :
          continue

      binerror = pow(binerror,0.5)
      TH1_tmp.SetBinContent(i, bincontent)
      TH1_tmp.SetBinError(i, binerror)

    '''
    c = TCanvas("c","c",10,10,1000,800)
    c.cd()
    c.SetLogy()
    TH1_tmp.SetMarkerStyle(20)
    TH1_tmp.Draw("ep3")
    c.SaveAs(TH1_tmp.GetName()+".png")
    '''

    return TH1_tmp

def splitAlpha(refHist,alpha) :

    TH1_tmp = refHist.Clone()
    TH1_tmp.SetName(alpha.GetName()+"_splitAlpha")
    TH1_tmp.SetTitle(alpha.GetTitle())
    
    for i in range(1,refHist.GetXaxis().GetNbins()+1):
      bincenter = refHist.GetXaxis().GetBinCenter(i)

      for j in range(1,alpha.GetXaxis().GetNbins()+1):
        lower_edge = alpha.GetXaxis().GetBinLowEdge(j)
        higher_edege = alpha.GetXaxis().GetBinUpEdge(j)
   
        if(bincenter>=lower_edge and bincenter<higher_edege) :
          bincontent = alpha.GetBinContent(j)
          binerror = alpha.GetBinError(j)

          TH1_tmp.SetBinContent(i,bincontent)
          TH1_tmp.SetBinError(i,binerror)

          break

        else :
          continue

    return TH1_tmp   

def alphaFactor(Shape_SR,Shape_SB,bins) :
    print bins

    alpha_num=rebinMerge(Shape_SR,bins)
    alpha_den=rebinMerge(Shape_SB,bins)
    
    alpha_num.Divide(alpha_den)
    alpha_num.SetName("alpha_"+alpha_num.GetName())

    '''
    c = TCanvas("c","c",10,10,1000,800)
    c.cd()
    alpha_num.Draw("ep3")
    c.SaveAs(alpha_num.GetName()+".png")
    '''

    alpha=splitAlpha(Shape_SR,alpha_num)

    #alpha=Shape_SR.Clone()
    #alpha.Divide(Shape_SB)

    return alpha

def density50GeV(hist) :

    hist_density=hist.Clone()
    hist_density.SetName(hist.GetName()+"_density50GeV")
    for i in range(1,hist.GetXaxis().GetNbins()+1):
        low = hist.GetXaxis().GetBinLowEdge(i)
        high = hist.GetXaxis().GetBinUpEdge(i)

        bincontent = hist.GetBinContent(i)
        binerror = hist.GetBinError(i)

        hist_density.SetBinContent(i,bincontent*50.0/(high-low))
        hist_density.SetBinError(i,binerror*50.0/(high-low))

    return hist_density

def plotDataVSbkg(dataSR,zjetsSR,ttbarSR,VVSR,title,minMass):

    dataSR_density=density50GeV(dataSR)
    zjetsSR_density=density50GeV(zjetsSR)

    VVSR_density=density50GeV(VVSR)
    ttbarVVSR_density=density50GeV(ttbarSR)
    ttbarVVSR_density.Add(VVSR_density)

    bkgSR_density=zjetsSR_density.Clone()
    bkgSR_density.Add(ttbarVVSR_density)
    bkgSR_density.SetName("bkgSR_density")

    bkgSR_noUnc_density=bkgSR_density.Clone()
    bkgSR_noUnc_density.SetName("bkgSR_noUnc_density")
    for i in range(1,bkgSR_noUnc_density.GetXaxis().GetNbins()+1):
        bkgSR_noUnc_density.SetBinError(i,0)

    bkgSR_density.SetLineColor(kGreen+2)
    bkgSR_density.SetFillColor(kGreen+2)

    bkg_unc_density = bkgSR_density.Clone()
    bkg_unc_density.SetName("bkgSR_unc_density")
    bkg_unc_density.Divide(bkgSR_noUnc_density)

    dataSR_density.SetMarkerStyle(20)
    dataSR_density.SetMarkerColor(kBlack)

    dataToMC_density = dataSR_density.Clone()
    dataToMC_density.SetName("dataToMC")
    dataToMC_density.Divide(bkgSR_noUnc_density)

    ttbarVVSR_density.SetLineColor(kYellow+2)
    ttbarVVSR_density.SetFillColor(kYellow+2)

    VVSR_density.SetLineColor(kMagenta-3)
    VVSR_density.SetFillColor(kMagenta-3)

    #plotting data/MC
    c1 = TCanvas("c1","c1", 800, 800)
    c1.SetLogy()
    c1.SetBottomMargin(0.3)
    c1.SetRightMargin(0.03)

    dummy = ROOT.TH1D("dummy","dummy", 1, minMass,3500)
    dummy.SetLabelSize(0)
    dummy.GetYaxis().SetTitle('Nevents/50GeV')
    dummy.SetMaximum(dataSR_density.GetBinContent(1)*2)
    dummy.SetMinimum(0.01)
    dummy.Draw()
    bkgSR_density.Draw("histsame")
    bkgSR_unc_density=bkgSR_density.Clone()
    bkgSR_unc_density.SetName("bkgSR_unc_density")
    bkgSR_unc_density.SetFillColor(kBlue+2)
    bkgSR_unc_density.SetFillStyle(3114)
    bkgSR_unc_density.SetMarkerSize(0.001)
    bkgSR_unc_density.SetMarkerColor(kBlue+2)
    bkgSR_unc_density.SetMarkerStyle(0)
    #bkgSR_unc_density.Draw("e2same")

    ttbarVVSR_density.SetLabelSize(0)
    ttbarVVSR_density.Draw("histsame")
    VVSR_density.SetLabelSize(0)
    VVSR_density.Draw("histsame")
    dataSR_density.SetLabelSize(0)
    dataSR_density.Draw("ep3same")

    legend = TLegend(.7,.7,.90,.90)
    if("CR" in title):
      print title, 'find CR'
      legend.AddEntry(dataSR_density, 'alpha#times sideband data', "p")
    else:
      print title, 'not find CR'
      legend.AddEntry(dataSR_density, 'SR data', "p")
    legend.AddEntry(bkgSR_density, 'z+jets', "f")
    legend.AddEntry(ttbarVVSR_density, 'WW/TT', "f")
    legend.AddEntry(VVSR_density, 'VV', "f")
    legend.AddEntry(bkgSR_unc_density, 'bkg uncertainty', "f")
    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.6*c1.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.92, 0.94," 35.6 fb^{-1} (13 TeV)")
    latex2.SetTextSize(0.4*c1.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.2, 0.94, "CMS")
    latex2.SetTextSize(0.4*c1.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    latex2.DrawLatex(0.3, 0.94, "Work in progress")

    gPad.RedrawAxis()

    pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
    pad.SetTopMargin(0.70);
    pad.SetRightMargin(0.03);
    pad.SetFillColor(0);
    #pad.SetGridy(1);
    pad.SetFillStyle(0);
    pad.Draw();
    pad.cd(0);

    dummyR = ROOT.TH1D("dummyR","dummyR", 1, minMass,3500)
    dummyR.SetBinContent(1,0.0)
    dummyR.GetXaxis().SetTitle('m_{ZZ} [GeV]')
    dummyR.SetMinimum(0.5)
    dummyR.SetMaximum(2.0)
    dummyR.GetXaxis().SetMoreLogLabels(kTRUE)
    dummyR.GetXaxis().SetNoExponent(kTRUE)
    dummyR.GetYaxis().SetTitleSize(0.02);
    dummyR.GetYaxis().SetTitleOffset(1.8);
    dummyR.GetYaxis().SetTitle("Data/MC")
    dummyR.GetYaxis().CenterTitle();
    dummyR.GetYaxis().SetLabelSize(0.02);
    dummyR.GetXaxis().SetMoreLogLabels(kTRUE)

    dummyR.GetXaxis().SetNoExponent(kTRUE)
    dummyR.GetXaxis().SetTitle("m_{ZZ} [GeV]")
    dummyR.Draw()

    line2 = TLine(dummyR.GetXaxis().GetBinLowEdge(1),1.,
                  dummyR.GetXaxis().GetBinUpEdge(dummyR.GetNbinsX()),1.)
    line2.SetLineColor(kRed)
    line2.SetLineStyle(kDashed)
    line2.Draw("same")

    dataToMC_density.Draw("ep3same")
    bkg_unc_density.SetFillColor(kBlue+2)
    bkg_unc_density.SetFillStyle(3114)
    bkg_unc_density.SetMarkerSize(0.001)
    bkg_unc_density.SetMarkerColor(kBlue+2)
    bkg_unc_density.SetMarkerStyle(0)
    bkg_unc_density.Draw("e2same")

    #c1.SaveAs('DataVSbkg_'+title+'.pdf')
    c1.SaveAs('DataVSbkg_'+title+'.png')

### Define function for processing of os command
def processCmd(cmd, quiet = 0):
    #print cmd
    #status, output = commands.getstatusoutput(cmd)
    #output = subprocess.check_output(cmd, shell=True)
    output = '\n'
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=STDOUT,bufsize=-1)
    for line in iter(p.stdout.readline, ''):
        output=output+str(line)
        print line
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError("%r failed, exit status: %d" % (cmd, p.returncode))
    #if (status !=0 and not quiet):
    #    print 'Error in processing command:\n   ['+cmd+']'
    #    print 'Output:\n   ['+output+'] \n'
    if (not quiet):
        print 'Output:\n   ['+output+'] \n'

    return output


def Run():
    parser=parseOptions()
    (options,args) = parser.parse_args()
    #alphaMethod(fs,jet,cat, parser)

    #alphaMethod('all','merged','btag', 400,parser)
    #alphaMethod('all','merged','vbf', 400,parser)

    alphaMethod('ee','merged','btag', 400,parser)
    alphaMethod('mumu','merged','btag', 400,parser)

    alphaMethod('ee','merged','vbf', 400,parser)
    alphaMethod('mumu','merged','vbf', 400,parser)

if __name__ == "__main__":
    Run()