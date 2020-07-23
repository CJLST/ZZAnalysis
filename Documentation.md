ZZAnalysis
==========

See README for instructions on how to install the framework.

## Samples corresponding to this tag
### Data
Data samples can be found at

	/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200430_LegacyRun2

These samples are in full sync with BBF group. 

### Simulation
MC samples can be found at

	/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased

Recall that for 2016 MC, you should use the samples under `MC_2016_CorrectBTag`.

## Correction of xsec and BR on the fly
During the sync exercise between HIG-19-001 and HIG-19-009 we found different
values of xsec and BR to be off w.r.t. the reference ones (i.e. Yellow Report 4).
The `xsec` branch containted in the `TTrees` represents the `xsec times BR`. This 
information is used to weight properly the samples and compute the yields.
Recall that the total weight to be applied to the samples is:

	weight= L1prefiringWeight * xsec* overallEventWeight/gen_sum_weights*_lumi*1000

To update `xsec` with the correct values, you can use [this patch](https://mbonanom.web.cern.ch/mbonanom/ZZ4l/CutBased_130220/updated_xsec.cc).

Note that the csv files present in this tag have been updated with the correct values.
Hence this patch won't be needed in case of a new production.

## VBS samples
The `xsec` value of VBS samples is off by a factor 1E3. Hence the weight to be applied
to these samples (`VBFToContinToZZ4l`) is

	weight= L1prefiringWeight * xsec* overallEventWeight/gen_sum_weights*_lumi*1000

## Muon and electron SF application
The samples were not generated with the latest muon and electron scale factors.
The updated scale factors have been included in the [LeptonSFHelper](https://github.com/CJLST/ZZAnalysis/blob/200717_LegacyRun2_PAPER/AnalysisStep/src/LeptonSFHelper.cc) and can be applied *on the fly* using this tool.
To do so:

	LeptonSFHelper *lepSFHelper = new LeptonSFHelper();
	_event_weight= L1prefiringWeight * xsec* overallEventWeight/gen_sum_weights*_lumi*1000
	float _updatedSF = 1.0;
	for (int lepl = 0; lepl < 4; ++lepl)
	{
	      float lid = LepLepId->at(lepl);
	      float lpt = LepPt->at(lepl);
	      float leta = LepEta->at(lepl);
	      if(abs(lid) == 11) leta = LepSCEta->at(lepl);
	      bool isCrack = LepisCrack->at(lepl);
	      _updatedSF *= lepSFHelper->getSF(year, lid, lpt, leta, leta, isCrack);
	}
	// check very rare cases in which _updatedSF = 0 and use the old one
	if(_updatedSF == 0) _updatedSF = dataMCWeight;

	_event_weight *= _updatedSF/dataMCWeight; 

If you are doing this on a piece of code external to the CJLST framework, 
you have to set `cmsenv` in a `CMSSW` working area containing the CJLST framework 
and import the correct module:

	#include "ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h"
	#include "ZZAnalysis/AnalysisStep/src/LeptonSFHelper.cc"

# Fake rate files and ZX yields
The fake rate files corresponding to this production can be found in the [`data/FakeRates` folder](https://github.com/CJLST/ZZAnalysis/tree/200717_LegacyRun2_PAPER/AnalysisStep/data/FakeRates), under the names

	newData_FakeRates_*S_*.root

The latest ZX yields, and the corresponding rations used for the combined estimation, can be found:

 * [Here](https://elfontan.web.cern.ch/elfontan/CutBased_analysis/ZX_ESTIMATE/2016_yields.txt) for **2016**, 
    both in the `m4l>70 GeV` and `118<m4l<140 GeV` mass ranges;
 * [Here](https://elfontan.web.cern.ch/elfontan/CutBased_analysis/ZX_ESTIMATE/2017_yields.txt) for **2017**, 
    both in the `m4l>70 GeV` and `118<m4l<140 GeV` mass ranges;
 * [Here](https://elfontan.web.cern.ch/elfontan/CutBased_analysis/ZX_ESTIMATE/2016_yields.txt) for **2018**, 
    both in the `m4l>70 GeV` and `118<m4l<140 GeV` mass ranges;
 * [Here for all three years in `105<m4l<140 GeV` mass range](https://elfontan.web.cern.ch/elfontan/CutBased_analysis/ZX_ESTIMATE/105_140_ZXyields.txt)

The systematic uncertainties associated to the z+jets yields are given by three contributions:

 * Statistical composition;
 * FR uncertainty;
 * Background composition;

The most up to date values can be found [here](https://elfontan.web.cern.ch/elfontan/CutBased_analysis/ZX_ESTIMATE/ZX_SystematicUncertainties.txt).




