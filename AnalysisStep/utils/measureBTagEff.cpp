// A utility script to measure the b tag efficiency / mistag rate for MC to be used with a dedicated csv file to applu b-tagg scale factor
// To run just do root measureBTagEff.cpp++()

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cctype>
#include <algorithm>
#include "TString.h"
#include "TChain.h"
#include "TCut.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGraph2DErrors.h"
#include "TFile.h"
#include "TEfficiency.h"

using namespace std;

void measureBTagEff()
{
    // read MC sample
    vector<TString> samples;
    samples.push_back("");
	

    // -------------------------------------------------------------------------------------------
    // histos
    float PtBins[]  =  {30, 40, 50, 70, 100, 150, 200, 300, 600} ;
    float EtaBins[] =  {0, 0.6, 1.2, 2.1, 2.5} ;
    int nPtBins  = sizeof(PtBins)/sizeof(float) - 1;
    int nEtaBins = sizeof(EtaBins)/sizeof(float) - 1;

    //float WPtag[3] = {0.1522, 0.4941, 0.8001}; // L, M, T -- 90X DeepCSV for Moriond 2018, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    
    float WPtag[3] = {0.1241, 0.4184, 0.7527};// L, M, T -- 90X DeepCSV for Moriond 2019, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    string WPname[3] = {"L", "M", "T"};


    // FIXME: do for the 3 WP tag x the number of selections
    vector<TH2F*> h2_BTaggingEff_Denom_b(3); // .at(WP).at(selection)
    vector<TH2F*> h2_BTaggingEff_Denom_c(3);
    vector<TH2F*> h2_BTaggingEff_Denom_udsg(3);
    vector<TH2F*> h2_BTaggingEff_Num_b(3);
    vector<TH2F*> h2_BTaggingEff_Num_c(3);
    vector<TH2F*> h2_BTaggingEff_Num_udsg(3);

    string outputFile = "bTagEfficiencies_102X_Moriond19.root";
    cout << "Saving efficiency weights into: " << outputFile << endl;
    TFile* fOut = new TFile (outputFile.c_str(), "recreate");
    fOut->cd();

    for (int iWP = 0; iWP < 3; iWP++)
    {
            h2_BTaggingEff_Denom_b    .at(iWP) = new TH2F (Form("h2_BTaggingEff_Denom_b_%s_ALL"   , WPname[iWP].c_str()), ";p_{T} [GeV];#eta", nPtBins, PtBins, nEtaBins, EtaBins);
            h2_BTaggingEff_Denom_c    .at(iWP) = new TH2F (Form("h2_BTaggingEff_Denom_c_%s_ALL"   , WPname[iWP].c_str()), ";p_{T} [GeV];#eta", nPtBins, PtBins, nEtaBins, EtaBins);
            h2_BTaggingEff_Denom_udsg .at(iWP) = new TH2F (Form("h2_BTaggingEff_Denom_udsg_%s_ALL", WPname[iWP].c_str()), ";p_{T} [GeV];#eta", nPtBins, PtBins, nEtaBins, EtaBins);
            h2_BTaggingEff_Num_b      .at(iWP) = new TH2F (Form("h2_BTaggingEff_Num_b_%s_ALL"     , WPname[iWP].c_str()), ";p_{T} [GeV];#eta", nPtBins, PtBins, nEtaBins, EtaBins);
            h2_BTaggingEff_Num_c      .at(iWP) = new TH2F (Form("h2_BTaggingEff_Num_c_%s_ALL"     , WPname[iWP].c_str()), ";p_{T} [GeV];#eta", nPtBins, PtBins, nEtaBins, EtaBins);
            h2_BTaggingEff_Num_udsg   .at(iWP) = new TH2F (Form("h2_BTaggingEff_Num_udsg_%s_ALL"  , WPname[iWP].c_str()), ";p_{T} [GeV];#eta", nPtBins, PtBins, nEtaBins, EtaBins);
    }

    // control histos
    TH1F* h_pT_b = new TH1F ("h_pT_b", "b jets pT distribution; pT; events", 200, 0, 2000);
    TH1F* h_pT_c = new TH1F ("h_pT_c", "c jets pT distribution; pT; events", 200, 0, 2000);
    TH1F* h_pT_udsg = new TH1F ("h_pT_udsg", "udsg jets pT distribution; pT; events", 200, 0, 2000);

    TH1F* h_eta_b = new TH1F ("h_eta_b", "b jets #eta distribution; #eta; events", 50, -2.5, 2.5);
    TH1F* h_eta_c = new TH1F ("h_eta_c", "c jets #eta distribution; #eta; events", 50, -2.5, 2.5);
    TH1F* h_eta_udsg = new TH1F ("h_eta_udsg", "udsg jets #eta distribution; #eta; events", 50, -2.5, 2.5);

    TH1F* h_CSV_b = new TH1F ("h_CSV_b", "b jets CSV distribution; CSV; events", 120, 0, 1.2);
    TH1F* h_CSV_c = new TH1F ("h_CSV_c", "c jets CSV distribution; CSV; events", 120, 0, 1.2);
    TH1F* h_CSV_udsg = new TH1F ("h_CSV_udsg", "udsg jets CSV distribution; CSV; events", 120, 0, 1.2);

    // ------------------------------------------------------------------------------------------

    // additional jets
    std::vector<float>* JetPt   = 0;
    std::vector<float>* JetEta  = 0;
    std::vector<float>* JetBTagger = 0;
    std::vector<int>*   JetHadronFlavour = 0;
    float genHEPMCweight;

    // loop on samples and fill num / denom histos
    for (unsigned int iSample = 0 ; iSample < samples.size () ; ++iSample)
    {
		  TFile *f = TFile::Open(samples.at(iSample));
        TTree *tree = (TTree*)f->Get("ZTree/candTree") ;

		
        tree->SetBranchStatus ("JetPt"   , 1);
        tree->SetBranchStatus ("JetEta"  , 1);
        tree->SetBranchStatus ("JetBTagger" , 1);
        tree->SetBranchStatus ("JetHadronFlavour" , 1);
        tree->SetBranchStatus ("genHEPMCweight" , 1);
		 
        // set addresses
        tree->SetBranchAddress ("JetPt", &JetPt);
        tree->SetBranchAddress ("JetEta", &JetEta);
        tree->SetBranchAddress ("JetBTagger", &JetBTagger);
        tree->SetBranchAddress ("JetHadronFlavour", &JetHadronFlavour);
        tree->SetBranchAddress ("genHEPMCweight", &genHEPMCweight);
		 
        // event loop
        int nEvts = tree->GetEntries();
        //int nEvts = 1000;

        cout << "Will process " << nEvts << " events" << endl;
        for (int iEvent = 0 ; iEvent < nEvts ; ++iEvent)
        {
            if (iEvent % 100000 == 0) cout << iEvent << " / " << nEvts << endl;
            //cout << iEvent << " / " << nEvts << endl;
            tree->GetEntry (iEvent) ;
				bool bTag [3] ;

				for (unsigned int ijet = 0; ijet < JetPt->size(); ijet++)
				{
					bTag[0] = (JetBTagger->at(ijet) > WPtag[0]) ;
					bTag[1] = (JetBTagger->at(ijet) > WPtag[1]) ;
					bTag[2] = (JetBTagger->at(ijet) > WPtag[2]) ;
					int flav =  JetHadronFlavour->at(ijet);
					float this_pt =  JetPt->at(ijet);
					float this_eta = TMath::Abs(JetEta->at(ijet));
					float this_eta_noAbs = JetEta->at(ijet);

					if (abs(flav) == 5) // b jets
					{
						h2_BTaggingEff_Denom_b.at(0) -> Fill (this_pt, this_eta);
						h2_BTaggingEff_Denom_b.at(1) -> Fill (this_pt, this_eta);
						h2_BTaggingEff_Denom_b.at(2) -> Fill (this_pt, this_eta);
						if (bTag[0]) h2_BTaggingEff_Num_b.at(0) -> Fill (this_pt, this_eta);
						if (bTag[1]) h2_BTaggingEff_Num_b.at(1) -> Fill (this_pt, this_eta);
						if (bTag[2]) h2_BTaggingEff_Num_b.at(2) -> Fill (this_pt, this_eta);

						h_pT_b->Fill (this_pt);
						h_eta_b->Fill (this_eta_noAbs);
						h_CSV_b->Fill (JetBTagger->at(ijet));
					}
					else if (abs(flav) == 4) // c jets
					{
						h2_BTaggingEff_Denom_c.at(0) -> Fill (this_pt, this_eta);
						h2_BTaggingEff_Denom_c.at(1) -> Fill (this_pt, this_eta);
						h2_BTaggingEff_Denom_c.at(2) -> Fill (this_pt, this_eta);
						if (bTag[0]) h2_BTaggingEff_Num_c.at(0) -> Fill (this_pt, this_eta);
						if (bTag[1]) h2_BTaggingEff_Num_c.at(1) -> Fill (this_pt, this_eta);
						if (bTag[2]) h2_BTaggingEff_Num_c.at(2) -> Fill (this_pt, this_eta);

						h_pT_c->Fill (this_pt);
						h_eta_c->Fill (this_eta_noAbs);
						h_CSV_c->Fill (JetBTagger->at(ijet));

					}
					else // udsg jets
					{
						h2_BTaggingEff_Denom_udsg.at(0) -> Fill (this_pt, this_eta);
						h2_BTaggingEff_Denom_udsg.at(1) -> Fill (this_pt, this_eta);
						h2_BTaggingEff_Denom_udsg.at(2) -> Fill (this_pt, this_eta);
						if (bTag[0]) h2_BTaggingEff_Num_udsg.at(0) -> Fill (this_pt, this_eta);
						if (bTag[1]) h2_BTaggingEff_Num_udsg.at(1) -> Fill (this_pt, this_eta);
						if (bTag[2]) h2_BTaggingEff_Num_udsg.at(2) -> Fill (this_pt, this_eta);
     
						h_pT_udsg->Fill (this_pt);
						h_eta_udsg->Fill (this_eta_noAbs);
						h_CSV_udsg->Fill (JetBTagger->at(ijet));
					}
			}
	}

    }

    cout << "...Finished histo filling, going to compute efficiencies" << endl;
    ///////////////////////////////////////////////////////////
    // num / denom are ready, compute efficiencies and store them

    fOut->cd();
    vector<TH2F*> eff_b    (3) ; // WP indexed
    vector<TH2F*> eff_c    (3) ; // WP indexed
    vector<TH2F*> eff_udsg (3) ; // WP indexed

    for (int iWP = 0; iWP < 3; iWP++)
    {
		TEfficiency* pEff_b = 0;
		TEfficiency* pEff_c = 0;
		TEfficiency* pEff_udsg = 0;

		eff_b.at(iWP) = new TH2F    (Form("eff_b_%s_ALL"   , WPname[iWP].c_str() ) , "eff_b; pT; |#eta|", nPtBins, PtBins, nEtaBins, EtaBins);
		eff_c.at(iWP) = new TH2F    (Form("eff_c_%s_ALL"   , WPname[iWP].c_str()) , "eff_c; pT; |#eta|", nPtBins, PtBins, nEtaBins, EtaBins);
		eff_udsg.at(iWP) = new TH2F (Form("eff_udsg_%s_ALL", WPname[iWP].c_str()) , "eff_udsg; pT; |#eta|", nPtBins, PtBins, nEtaBins, EtaBins);

		// b jets
		if(TEfficiency::CheckConsistency( *(h2_BTaggingEff_Num_b.at(iWP)), *(h2_BTaggingEff_Denom_b.at(iWP)) ) )
		{
			pEff_b = new TEfficiency(*(h2_BTaggingEff_Num_b.at(iWP)), *(h2_BTaggingEff_Denom_b.at(iWP)));
			for (int gBin = 0; gBin < (h2_BTaggingEff_Num_b.at(iWP))->GetSize(); gBin++)
			{
				(eff_b.at(iWP))->SetBinContent(gBin, pEff_b->GetEfficiency(gBin));
			}

		}
		// c jets
		if(TEfficiency::CheckConsistency( *(h2_BTaggingEff_Num_c.at(iWP)), *(h2_BTaggingEff_Denom_c.at(iWP)) ) )
		{
			pEff_c = new TEfficiency(*(h2_BTaggingEff_Num_c.at(iWP)), *(h2_BTaggingEff_Denom_c.at(iWP)));
			for (int gBin = 0; gBin < (h2_BTaggingEff_Num_c.at(iWP))->GetSize(); gBin++)
			{
				(eff_c.at(iWP))->SetBinContent(gBin, pEff_c->GetEfficiency(gBin));
			}
		
		}
		// udsg jets
		if(TEfficiency::CheckConsistency( *(h2_BTaggingEff_Num_udsg.at(iWP)), *(h2_BTaggingEff_Denom_udsg.at(iWP)) ) )
		{
			pEff_udsg = new TEfficiency(*(h2_BTaggingEff_Num_udsg.at(iWP)), *(h2_BTaggingEff_Denom_udsg.at(iWP)));
			for (int gBin = 0; gBin < (h2_BTaggingEff_Num_udsg.at(iWP))->GetSize(); gBin++)
			{
				(eff_udsg.at(iWP))->SetBinContent(gBin, pEff_udsg->GetEfficiency(gBin));
			}
			
		}
		 
    }
    // save
    fOut->Write();
}
