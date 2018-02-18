# Due to PU problem in 2017 MC per-sample reweighting is needed. This script reweights CJLST Ntuple to data profile changing PUweights, overallEventWeight and sumof weights in hCounters
# Run with python reweight_sample -i /SampleName/ZZ4lAnalysis.root
# Script produces /SampleName/ZZ4lAnalysis_corrPU.root file with correct PU weights

import ROOT,json,array,math,sys
from ROOT import *
from optparse import OptionParser


def main(options):
	f = ROOT.TFile.Open(options.input, 'READ')
	
	MCPU = f.Get("PlotsZZ/hNTrueIntNoWeight")
	MCPU.Scale(1./MCPU.Integral())
	MCPU_rw = TH1F("MCPU_rw","MCPU_rw",100,0,100)
	
	data    = ROOT.TFile.Open("DataPUProfile/DataPileupHistogram2017_69200_100bins.root",'READ')
	data_up = ROOT.TFile.Open("DataPUProfile/DataPileupHistogram2017_72383_100bins.root",'READ')
	data_dn = ROOT.TFile.Open("DataPUProfile/DataPileupHistogram2017_66017_100bins.root",'READ')

	dataPU    = data.Get("pileup")
	dataPU.Scale(1./dataPU.Integral())
	dataPU_up = data_up.Get("pileup")
	dataPU_up.Scale(1./dataPU_up.Integral())
	dataPU_dn = data_dn.Get("pileup")
	dataPU_dn.Scale(1./dataPU_dn.Integral())
	
	PUweights    = []
	PUweights_up = []
	PUweights_dn = []
	
	if (dataPU.GetSize() != MCPU.GetSize()):
		print "ERROR!! Data and MC do not have same number of bins!!! Can not calculate PU weights."
		sys.exit()
	else:
		for bin in range(1,dataPU.GetSize()-1):
			if( dataPU.GetBinContent(bin) == 0 or MCPU.GetBinContent(bin) == 0):
				PUweights.append(0.)
			else:
				PUweights.append(dataPU.GetBinContent(bin)/MCPU.GetBinContent(bin))
			if( dataPU_up.GetBinContent(bin) == 0 or MCPU.GetBinContent(bin) == 0):
				PUweights_up.append(0.)
			else:
				PUweights_up.append(dataPU_up.GetBinContent(bin)/MCPU.GetBinContent(bin))
			if( dataPU_dn.GetBinContent(bin) == 0 or MCPU.GetBinContent(bin) == 0):
				PUweights_dn.append(0.)
			else:
				PUweights_dn.append(dataPU_dn.GetBinContent(bin)/MCPU.GetBinContent(bin))

	for bin in range(1,dataPU.GetSize()-1):
		MCPU_rw.SetBinContent(bin, MCPU.GetBinContent(bin)*PUweights[bin - 1])

#	c = TCanvas("c","c",900,900)
#	c.cd()
#	MCPU_rw.SetLineColor(ROOT.kBlack)
#	MCPU_rw.Draw("HIST")
#	MCPU.SetLineColor(ROOT.kRed)
#	MCPU.Draw("HIST SAME")
#	dataPU.SetLineColor(ROOT.kGreen)
#	dataPU.Draw("SAME")
#	c.SaveAs("controlPlot.pdf")


	outFile = options.input.split(".root")[0]+"_corrPU.root"
	newFile = ROOT.TFile.Open(outFile, "RECREATE")

   #Copy control plots first
	directory_new = newFile.mkdir("PlotsZZ")
	directory_new.cd()
	hNvtxNoWeight = f.Get("PlotsZZ/hNvtxNoWeight").Clone()
	hNvtxWeight = f.Get("PlotsZZ/hNvtxWeight").Clone()
	hNTrueIntNoWeight = f.Get("PlotsZZ/hNTrueIntNoWeight").Clone()
	hNTrueIntWeight = f.Get("PlotsZZ/hNTrueIntWeight").Clone()
	hRhoNoWeight = f.Get("PlotsZZ/hRhoNoWeight").Clone()
	hRhoWeight = f.Get("PlotsZZ/hRhoWeight").Clone()
	hGenZZMass = f.Get("PlotsZZ/hGenZZMass").Clone()
	hGenZZMass_4e = f.Get("PlotsZZ/hGenZZMass_4e").Clone()
	hGenZZMass_4mu = f.Get("PlotsZZ/hGenZZMass_4mu").Clone()
	hGenZZMass_2e2mu = f.Get("PlotsZZ/hGenZZMass_2e2mu").Clone()
	hGenZZMass_2l2tau = f.Get("PlotsZZ/hGenZZMass_2l2tau").Clone()
	MCPU_rw.Write()

	directories = {"ZZTree", "ZTree", "CRZLTree", "CRZLLTree"}

	for directory in directories:
		if (not f.GetListOfKeys().Contains(directory)):
			continue
		print directory
		gen_sumWeights = 0.
		gen_sumGenMCWeight = 0.
		gen_sumPUWeight = 0.
		tree = f.Get(directory+"/candTree")
		tree.SetBranchStatus("*", 1);
		tree.SetBranchStatus("PUWeight", 0);
		tree.SetBranchStatus("overallEventWeight",0);
		PUWeight = array.array('f',[0])
		overallEventWeight = array.array('f',[0])

		entries = tree.GetEntries()
		hCounter = f.Get(directory+"/Counters")
		#--- Write to new file
		directory_new = newFile.mkdir(directory)
		directory_new.cd()
		tree_new = tree.CloneTree(0)
		tree_new.Branch("PUWeight", PUWeight, "PUWeight/F")
		tree_new.Branch("overallEventWeight", overallEventWeight, "overallEventWeight/F")
		hCounter_new = hCounter.Clone()



		for z in range(entries):
			tree.GetEntry(z)
#			if (z%100000 == 0):
#				percent = z/float(entries)*100
#				print "%.2f %% processed..." %percent
			if (tree.NTrueInt > 99. or tree.NTrueInt < 0.):#protect agains crazy values of NTrueInt
				PUWeight[0] = 0.
				overallEventWeight[0] = 0.
			else:
				PUWeight[0] = PUweights[int(tree.NTrueInt)]
				overallEventWeight[0] = PUWeight[0]*tree.genHEPMCweight*tree.dataMCWeight
			gen_sumWeights += overallEventWeight[0]
			gen_sumGenMCWeight += tree.genHEPMCweight
			gen_sumPUWeight += PUWeight[0]
			tree_new.Fill()

		if(directory == "ZTree"):
			hCounter_new.SetBinContent(0,gen_sumWeights);
			hCounter_new.SetBinContent(1,gen_sumWeights);
			hCounter_new.SetBinContent(2,gen_sumGenMCWeight);
			hCounter_new.SetBinContent(3,gen_sumPUWeight);
		else:
			hCounter_new.SetBinContent(0,gen_sumWeights);
			hCounter_new.SetBinContent(40,gen_sumWeights);
			hCounter_new.SetBinContent(41,gen_sumGenMCWeight);
			hCounter_new.SetBinContent(42,gen_sumPUWeight);

		tree_new.GetCurrentFile().Write("",ROOT.TObject.kOverwrite)
		hCounter_new.Write("",ROOT.TObject.kOverwrite)

		if (directory != "ZZTree" or (not options.failedTree)):
			continue
		print directory + "/candTree_failed"
		ftree = f.Get(directory+"/candTree_failed")
		ftree.SetBranchStatus("*", 1);
		#--- Write to new file
		directory_new.cd()
		ftree_new = ftree.CloneTree()
		ftree_new.GetCurrentFile().Write("",ROOT.TObject.kOverwrite)




if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-i", "--input",     default="~/Moriond2016/CMSSW_7_6_3/src/PhysicsTools/TagAndProbe/test/TnPTree_data.root",           help="Input filename")
    parser.add_option("--failedTree", action="store_true", default=False, dest="failedTree")

    (options, arg) = parser.parse_args()
     
    main(options)
