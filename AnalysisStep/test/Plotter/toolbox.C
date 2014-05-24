
// Compute ratio of a step plot bin w.r.t. Nevt_Gen and compare to standard selection

void ratios(TString finalState="2e2mu") {
  TString base = "root://cmsphys05//data1/b/botta/HZZ_root/";
  TString dir[5] ={"PRODSTD/", // The reference selection
		   "PRODNo12/", 
		   "PROD4Over4/",
		   "PROD4Over4No12/",
		   "PRODLOWm/"
  };

  vector<float> ratio120(5);
  vector<float> ratio125(5);

    for (int i=0; i<5; ++i){
      TFile* fh120 = TFile::Open(base+dir[i]+"ZZ4lAnalysis_H120.root");
      TString rootdir = TString("Plots") + finalState + "/";    
      TH1F*  f120 = (TH1F*) fh120->Get(rootdir + "nEventComplete");     
      TFile* fh125 = TFile::Open(base+dir[i]+"ZZ4lAnalysis_H125.root");
      TString rootdir = TString("Plots") + finalState + "/";    
      TH1F*  f125 = (TH1F*) fh125->Get(rootdir + "nEventComplete");     

      ratio120[i]=f120->GetBinContent(11)/f120->GetBinContent(0);
      ratio125[i]=f125->GetBinContent(11)/f125->GetBinContent(0);
      
      
      cout <<  base+dir[i] << endl
	   << " h120: " << ratio120[i] << " R/Rstd: " << ratio120[i]/ratio120[0] << endl
	   << " h125: " << ratio125[i] << " R/Rstd: " << ratio125[i]/ratio125[0] << endl; 
      
    }
}
