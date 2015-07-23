void rootFileIntegrity(TString filename) {
  bool isZombie=false;
  bool isRecovered=false;
  {
    TFile f(filename, "read");
    isZombie=f.IsZombie();
    isRecovered=f.TestBit(TFile::kRecovered);    

    // Test for empty file that still pass the above checks
    if (!(isZombie || isRecovered)) {
      TList* l = f.GetListOfKeys();
      if(l->IsEmpty()){
	isRecovered=true;
      }
    }
  }

  if (isZombie ){
    cout << "File " << filename << " corrupted; renaming" << endl;
    TString command = "mv " + filename + " " + filename + ".corrupted";
    gSystem->Exec(command);
    exit(EXIT_FAILURE);
  } else if (isRecovered ){
    cout << "File " << filename << " was not closed correctly ; renaming" << endl;
    TString command = "mv " + filename + " " + filename + ".recovered";
    gSystem->Exec(command);
    exit(EXIT_FAILURE);
  } else {
    cout << "Integrity check succeeded on " << filename << endl;
  }
}
