#include "root_lib/TFileServiceWrapper.h"
#include "../../interface/Histograms.h"
#include "ZZTreeReader.h"
#include <algorithm>




void ZZTreeReader::Loop(bool useWeight) {
  if (fChain == 0) return;

  HCand* hCandZZCandM70     = new HCand("hCandZZCandM70");
  HCand* hCandZZCandM100    = new HCand("hCandZZCandM100");
  HCand* hCandZZCandM100_zz = new HCand("hCandZZCandM100_zz");
  HCand* hCandZZCandM100_Mela05 = new HCand("hCandZZCandM100_Mela05");
  HCand* hCandZZCandM121_131    = new HCand("hCandZZCandM121_131");


  //--- Loop over events
  const Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = 0;
  for (Long64_t jentry=0;jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);

    // Find passing candidates
    if (iBC>=0 ) { //iBC includes: goodLeptons, best Z1, Z1mass, Z2 choice

      float weight = 1.;
      if (useWeight) weight = PUWeight12;


      float ZZMassBC = ZZMass->at(iBC);
      float Z1MassBC = Z1Mass->at(iBC);
      float Z2MassBC = Z2Mass->at(iBC);
      float ZZPtBC   = ZZPt->at(iBC);
	  float LDBC     = p0plus_VAJHU->at(iBC)/(p0plus_VAJHU->at(iBC) + bkg_VAMCFM->at(iBC));
	  float pseudoLDBC     = p0plus_VAJHU->at(iBC)/(p0plus_VAJHU->at(iBC) + p0minus_VAJHU->at(iBC));
	  float cosThetaStarBC = costhetastar->at(iBC);
      float cosTheta1BC    = helcosthetaZ1->at(iBC);
      float cosTheta2BC    = helcosthetaZ2->at(iBC);
      float phiBC          = helphi->at(iBC);
      float phistar1BC     = phistarZ1->at(iBC);
      

      vector<float> pt(4);
      vector<float> isoPF(4);
      vector<float> SIP(4);
      pt[0]    = Lep1Pt->at(iBC);
      pt[1]    = Lep2Pt->at(iBC);
      pt[2]    = Lep3Pt->at(iBC);
      pt[3]    = Lep4Pt->at(iBC);
      isoPF[0] = Lep1combRelIsoPF->at(iBC);
      isoPF[1] = Lep2combRelIsoPF->at(iBC);
      isoPF[2] = Lep3combRelIsoPF->at(iBC);
      isoPF[3] = Lep4combRelIsoPF->at(iBC);

      vector< float> ptS(pt);
      vector<float> isoPFS(isoPF);
      vector<float> SIPS(SIP);
      sort(ptS.begin(),ptS.end());
      sort(isoPFS.begin(),isoPFS.end());
      sort(SIPS.begin(),SIPS.end());
      float pt1 =ptS[3]; // leading-pT
      float pt2 =ptS[2]; // sub-leading pT
      
	  bool FullSel70  = ZZsel->at(iBC)>=70;
      bool FullSel100 = ZZsel->at(iBC)>=100;
      bool passMz_zz = (Z1MassBC>60. && Z1MassBC<120. && Z2MassBC>60. && Z2MassBC<120.);

      if (FullSel70) {	
	hCandZZCandM70->Fill(ZZMassBC, Z1MassBC, Z2MassBC, 1., 1., ptS, isoPFS, SIPS, weight,
						 cosThetaStarBC, cosTheta1BC, cosTheta2BC, phiBC, phistar1BC, LDBC, pseudoLDBC);
      }

      if (FullSel100) {
	hCandZZCandM100->Fill(ZZMassBC, Z1MassBC, Z2MassBC, 1., 1., ptS, isoPFS, SIPS, weight, 
			      cosThetaStarBC, cosTheta1BC, cosTheta2BC, phiBC, phistar1BC, LDBC, pseudoLDBC);

	if (FullSel70 && passMz_zz) { // for ZZ x-section measurement
	  hCandZZCandM100_zz->Fill(ZZMassBC, Z1MassBC, Z2MassBC, 1., 1., ptS, isoPFS, SIPS, weight,
							cosThetaStarBC, cosTheta1BC, cosTheta2BC, phiBC, phistar1BC, LDBC, pseudoLDBC);
	}

	if (FullSel70 && LDBC > 0.5) { // for high-mela events
	  hCandZZCandM100_Mela05->Fill(ZZMassBC, Z1MassBC, Z2MassBC, 1., 1., ptS, isoPFS, SIPS, weight,
				       cosThetaStarBC, cosTheta1BC, cosTheta2BC, phiBC, phistar1BC, LDBC, pseudoLDBC);
	}
	
	if (FullSel70 && ZZMassBC>121.5 && ZZMassBC<130.5) { //three central signal bins
	  hCandZZCandM121_131->Fill(ZZMassBC, Z1MassBC, Z2MassBC, 1., 1., ptS, isoPFS, SIPS, weight,
				    cosThetaStarBC, cosTheta1BC, cosTheta2BC, phiBC, phistar1BC, LDBC, pseudoLDBC);
	}
	}
    }
  }
  delete hCandZZCandM70;
  delete hCandZZCandM100;
  delete hCandZZCandM100_zz;
  delete hCandZZCandM100_Mela05;
  delete hCandZZCandM121_131;

}

