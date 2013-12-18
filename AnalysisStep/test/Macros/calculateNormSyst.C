{
  gROOT->SetStyle("Plain");

  ifstream inputFile1;
  inputFile1.open("Nevt_4mu.txt");

  if(!inputFile1){
    cout << "Error: an input file is missing!" << endl;
  }


  //const int Npoints = 45;
  //double mHvec[Npoints] = {115.,116.,117.,118.,119.,120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,130.,135.,140.,145.,150.,160.,170.,180.,190.,200.,220.,250.,275.,300.,325.,350.,375.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.,900.,950.,1000.};
  //double err_mHvec[Npoints] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  const int Npoints = 34;
  double mHvec[Npoints] = {120.,124.,125.,126.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,250.,275.,300.,325.,350.,400.,425.,450.,475.,525.,550.,575.,600.,650.,700.,750.,800.,900.,950.,1000.};
  double err_mHvec[Npoints] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  double parMatrix1[Npoints];

  for(int i1=0; i1<Npoints; i1++){ 

    double number1 = 0.;
    inputFile1 >> number1;

    parMatrix1[i1] = number1;

    //cout << number1 << endl;
  }  

  TGraph *graph_N_diff = new TGraph(Npoints, mHvec, parMatrix1);
  graph_N_diff->SetTitle("");
  graph_N_diff->GetXaxis()->SetTitle("m_{H}");
  graph_N_diff->GetYaxis()->SetTitle("Normalization systematic error for 4#mu");
  //graph_N_diff->GetYaxis()->SetRangeUser(0.,2.);

  TCanvas c1;
  c1.cd();
  graph_N_diff->Draw("A*");
  c1.SaveAs("NormSyst.png");





}
