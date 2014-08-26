/** 
 * Print x-section files for existing signal samples.
 * Must be run from HZZ4L_Combination/CombinationPy/CreateDatacards
 * run as:
 * root -b -q computeCSW.cc
 *
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <iomanip>

#include <TSystem.h>
#include <TROOT.h>

using namespace std;

// // Higgs branching ratios
/***********************IDs************************/
/*                       Total = 0                */
/*                       H->bb = 1                */
/*                   H->tautau = 2                */
/*                     H->mumu = 3                */
/*                       H->ss = 4                */
/*                       H->cc = 5                */
/*                       H->tt = 6                */
/*                       H->gg = 7                */
/*                   H->gamgam = 8                */
/*                     H->gamZ = 9                */
/*                       H->WW = 10               */
/*                       H->ZZ = 11               */
/*                       H->4e = 12               */
/*                    H->2e2mu = 13               */
/*              H->4lep (e,mu) = 14               */
/*          H->4lep (e,mu,tau) = 15               */
/*                H->e+nu e-nu = 16               */
/*               H->e+nu mu-nu = 17               */
/*    H->2l2nu(l=e,mu)(nu=any) = 18               */
/* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
/*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
/* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
/* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
/*  H->2nu2q (nu=any)(q=udcsb) = 23               */
/*            H->4q (q=udcsb) = 24                */
/*      H->4f (f=any fermion) = 25                */
/**************************************************/

// // Higgs inclusive cross sections 
/**********IDs*************/
/*     ggToH = 1          */
/*       VBF = 2          */
/*        WH = 3          */
/*        ZH = 4          */
/*       ttH = 5          */
/**************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include <TROOT.h>
#include "../../../../Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h"
#endif

int computeCSW()
{
  gSystem->Load("libHiggsHiggs_CS_and_Width.so");
  gROOT->LoadMacro("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/include/HiggsCSandWidth.h+");

  double sqrts;
  cout.precision(10);

 
  HiggsCSandWidth *myCSW = new HiggsCSandWidth("YR3",gSystem->ExpandPathName("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/txtFiles/"));
  //  HiggsCSandWidth *myCSW = new HiggsCSandWidth(gSystem->ExpandPathName("$CMSSW_BASE/src/Higgs/Higgs_CS_and_Width/txtFiles/"));

  /*//////////
  8 TeV
  *///////////
  sqrts = 8;
  double total8TeV[] = {110, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 135, 140, 145, 150, 160, 170, 180, 190, 200, 220, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
  int l_total8TeV = sizeof(total8TeV)/8;

  double gg8TeV[] = {110, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 135, 140, 145, 150, 160, 170, 175, 180, 185, 190, 200, 220, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
  int l_gg8TeV = sizeof(gg8TeV)/8;

  double gg8TeV_4l[] = {91.2, 110, 120, 125, 125.6, 126, 130, 140};
  int l_gg8TeV_4l = sizeof(gg8TeV_4l)/8;


  double VBF8TeV[] = {115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 125.6, 126, 127, 128, 129, 130, 135, 140, 145, 150, 160, 170, 180, 190, 200, 220, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
  int l_VBF8TeV = sizeof(VBF8TeV)/8;

  // WH_ZH_TTH
  double assoc8TeV_split[] = {110, 115, 120, 125, 125.6, 126, 130, 140, 150, 160, 180, 200};
  int l_assoc8TeV_split = sizeof(assoc8TeV_split)/8;

  // minlo samples 
  double minlo8TeV[] = {90, 95, 100, 105, 110, 115, 120, 124, 125, 126, 130, 135, 140, 145, 150, 155, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
  int l_minlo8TeV = sizeof(minlo8TeV)/8;


  //Create and open the output file
  ofstream fileOut8TeV;
  fileOut8TeV.open("Xsection8TeV_tmp.txt");

  //Total 8 TeV
  for( int i = 0; i < l_total8TeV; i++)
    {
      double mH = total8TeV[i];      
      double BR = myCSW->HiggsBR(15,total8TeV[i]);
      double ggToH      = myCSW->HiggsCS(1,mH,sqrts);
      double VBF        = myCSW->HiggsCS(2,mH,sqrts);
      double WH         = myCSW->HiggsCS(3,mH,sqrts);
      double ZH         = myCSW->HiggsCS(4,mH,sqrts);
      double ttH        = myCSW->HiggsCS(5,mH,sqrts);
      double total      = myCSW->HiggsCS(0,mH,sqrts);
      double totalBySum = ggToH + VBF + WH + ZH + ttH;
      
      //Write the output file
      fileOut8TeV << "all         totH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << totalBySum*BR  <<  "  1" << 

	" " << total << " " << totalBySum << " " << BR << endl; 
  
    }

  fileOut8TeV << endl << endl << endl;


  //ggH 8 TeV
  for( int i = 0; i < l_gg8TeV; i++)
    {
      double mH = gg8TeV[i];      
      double BR = myCSW->HiggsBR(15,gg8TeV[i]);
      double ggToH      = myCSW->HiggsCS(1,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut8TeV << "Only ggH" << endl;
      fileOut8TeV << "all          ggH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ggToH*BR  <<  "  1" << endl; 
  
    }

  fileOut8TeV << endl << endl << endl;

  //ggH 8 TeV (only l=e,mu)
  for( int i = 0; i < l_gg8TeV_4l; i++)
    {
      double mH = gg8TeV_4l[i];      
      double BR = myCSW->HiggsBR(14,gg8TeV_4l[i]);
      double ggToH      = myCSW->HiggsCS(1,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut8TeV << "Only ggH, only l=e,mu" << endl;
      fileOut8TeV << "all            H"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ggToH*BR  <<  "  1" << endl; 
  
    }

  fileOut8TeV << endl << endl << endl;
    

  //VBF 8 TeV
  for( int i = 0; i < l_VBF8TeV; i++)
    {
      double mH = VBF8TeV[i];      
      double BR = myCSW->HiggsBR(15,VBF8TeV[i]);
      double VBF        = myCSW->HiggsCS(2,mH,sqrts);

      //Write the output file
      fileOut8TeV << "all            VBFH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << VBF*BR  <<  "  1" << endl; 
      
    }
  

  fileOut8TeV << endl << endl << endl;


  //ggH 8 TeV (only l=e,mu)
  for( int i = 0; i < l_gg8TeV_4l; i++)
    {
      double mH = gg8TeV_4l[i];      
      double BR = myCSW->HiggsBR(14,gg8TeV_4l[i]);
      double VBF      = myCSW->HiggsCS(2,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut8TeV << "Only VBF, only l=e,mu" << endl;
      fileOut8TeV << "all           VBFH_emu"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << VBF*BR  <<  "  1" << endl; 
  
    }

  fileOut8TeV << endl << endl << endl;


  //These samples have W/Z inclusive decays and ZZ from Higgs inclusive as well, but require 4 leptons in the final state coming from W or Z
  //WH 8TeV
  for( int i = 0; i < l_assoc8TeV_split; i++)
    {

      double mH = assoc8TeV_split[i];      
      double BRHZZ = myCSW->HiggsBR(11,mH);
      double BRH4l_emu = myCSW->HiggsBR(14,mH);
      double WH    = myCSW->HiggsCS(3,mH,sqrts);

      //Write the output file
      if (i==0) fileOut8TeV << "WH, H->ZZ |  H->4l, l=e,mu" << endl;
      fileOut8TeV << "all            WH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << WH*BRHZZ  <<  "  1 " << WH*BRH4l_emu << endl; 

    }

  fileOut8TeV << endl << endl << endl;

  //ZH 8TeV
  for( int i = 0; i < l_assoc8TeV_split; i++)
    {
      
      double mH = assoc8TeV_split[i];      
      double BRHZZ = myCSW->HiggsBR(11,mH);
      double BRH4l_emu = myCSW->HiggsBR(14,mH);
      double ZH    = myCSW->HiggsCS(4,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut8TeV << "ZH, H->ZZ |  H->4l, l=e,mu" << endl;
      fileOut8TeV << "all            ZH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ZH*BRHZZ  <<  "  1 " << ZH*BRH4l_emu << endl; 
      
    }

  fileOut8TeV << endl << endl << endl;

  //ttH 8TeV
  for( int i = 0; i < l_assoc8TeV_split; i++)
    {
      
      double mH = assoc8TeV_split[i];      
      double BRHZZ = myCSW->HiggsBR(11,mH);
      double BRH4l_emu = myCSW->HiggsBR(14,mH);
      double ttH    = myCSW->HiggsCS(5,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut8TeV << "ttH, H->ZZ |  H->4l, l=e,mu" << endl;
      fileOut8TeV << "all            ttH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ttH*BRHZZ  <<  "  1 " << ttH*BRH4l_emu << endl; 
      
    }

  fileOut8TeV << endl << endl << endl;


  //minlo samples, x-sections come from elsewhere, we just need to add BRs
  for( int i = 0; i < l_minlo8TeV; i++) {

    double mH = minlo8TeV[i];      
    double BR = myCSW->HiggsBR(15,mH);
    double ggToH      = myCSW->HiggsCS(1,mH,sqrts);

    fileOut8TeV << "all            minloH"  << (int) mH <<"                      1   " << ggToH << "  "  << std::fixed << setprecision(10) << BR << endl; 

  }
  

  fileOut8TeV.close();
  
  /*//////////
  7 TeV
  *///////////
  sqrts = 7;
  double total7TeV[] = {110, 115, 120, 124, 125, 126, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
  int l_total7TeV = sizeof(total7TeV)/8;


  double gg7TeV[] = {110, 115, 120, 122, 124, 125, 125.6, 126, 128, 130, 135, 140, 145, 150, 160, 170, 175, 180, 185, 190, 200, 210, 220, 225, 230, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
  int l_gg7TeV = sizeof(gg7TeV)/8;

  double gg7TeV_4l[] = {91.2, 110, 120, 125, 125.6, 126, 130, 140};
  int l_gg7TeV_4l = sizeof(gg7TeV_4l)/8;


  double VBF7TeV[] = {115, 120, 125, 125.6, 126, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 225, 230, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 650, 700, 800, 900, 950, 1000};
  int l_VBF7TeV = sizeof(VBF7TeV)/8;

  // WH_ZH_TTH
  double assoc7TeV_split[] = {110, 115, 120, 125, 125.6, 126, 130, 140, 150, 160, 180, 200};
  int l_assoc7TeV_split = sizeof(assoc7TeV_split)/8;


  //Create and open the output file
  ofstream fileOut7TeV;
  fileOut7TeV.open("Xsection7TeV_tmp.txt");

  //Total 7 TeV
  for( int i = 0; i < l_total7TeV; i++)
    {

      double mH = total7TeV[i];      
      double BR = myCSW->HiggsBR(15,mH);
      double ggToH      = myCSW->HiggsCS(1,mH,sqrts);
      double VBF        = myCSW->HiggsCS(2,mH,sqrts);
      double WH         = myCSW->HiggsCS(3,mH,sqrts);
      double ZH         = myCSW->HiggsCS(4,mH,sqrts);
      double ttH        = myCSW->HiggsCS(5,mH,sqrts);
      double total      = myCSW->HiggsCS(0,mH,sqrts);
      double totalBySum = ggToH + VBF + WH + ZH + ttH;

      //Write the output file
      fileOut7TeV << "all            H"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << totalBySum*BR  <<  "  1" << endl; 

    }

  fileOut7TeV << endl << endl << endl;

  //ggH 7 TeV
  for( int i = 0; i < l_gg7TeV; i++)
    {

      double mH = gg7TeV[i];      
      double BR = myCSW->HiggsBR(15,mH);
      double ggToH      = myCSW->HiggsCS(1,mH,sqrts);

      //Write the output file
      if (i==0) fileOut7TeV << "Only ggH" << endl;
      fileOut7TeV << "all            H"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ggToH*BR  <<  "  1" << endl; 

    }

  fileOut7TeV << endl << endl << endl;

  //ggH 7 TeV (only l=e,mu)
  for( int i = 0; i < l_gg7TeV_4l; i++)
    {
      double mH = gg7TeV_4l[i];      
      double BR = myCSW->HiggsBR(14,gg7TeV_4l[i]);
      double ggToH      = myCSW->HiggsCS(1,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut7TeV << "Only ggH, only l=e,mu" << endl;
      fileOut7TeV << "all            H"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ggToH*BR  <<  "  1" << endl; 
  
    }

  fileOut7TeV << endl << endl << endl;


  //VBF 7 TeV
  for( int i = 0; i < l_VBF7TeV; i++)
    {
      double mH = VBF7TeV[i];      
      double BR = myCSW->HiggsBR(15,VBF7TeV[i]);
      double VBF        = myCSW->HiggsCS(2,mH,sqrts);

      //Write the output file
      fileOut7TeV << "all            VBFH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << VBF*BR  <<  "  1" << endl; 
      
    }

  fileOut7TeV << endl << endl << endl;



  //VBF 7 TeV (only l=e,mu)
  for( int i = 0; i < l_gg7TeV_4l; i++)
    {
      double mH = gg7TeV_4l[i];      
      double BR = myCSW->HiggsBR(14,gg7TeV_4l[i]);
      double VBF        = myCSW->HiggsCS(2,mH,sqrts);
      
      //Write the output file
      if (i==0) fileOut7TeV << "Only VBF, only l=e,mu" << endl;
      fileOut7TeV << "all           VBFH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << VBF*BR  <<  "  1" << endl; 
  
    }

  fileOut7TeV << endl << endl << endl;



  //These samples have W/Z inclusive decays and ZZ from Higgs inclusive as well, but require 4 leptons in the final state coming from W or Z
  //WH 7TeV
  for( int i = 0; i < l_assoc7TeV_split; i++)
    {

      double mH = assoc7TeV_split[i];      
      double BRHZZ = myCSW->HiggsBR(11,mH);
      double BRH4l_emu = myCSW->HiggsBR(14,mH);
      double WH         = myCSW->HiggsCS(3,mH,sqrts);

      //Write the output file
      if (i==0) fileOut8TeV << "WH, H->ZZ |  H->4l, l=e,mu" << endl;
      fileOut7TeV << "all            WH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << WH*BRHZZ  <<  "  1 " << WH*BRH4l_emu << endl; 

    }

  fileOut7TeV << endl << endl << endl;

  //ZH 7TeV
  for( int i = 0; i < l_assoc7TeV_split; i++)
    {
      
      double mH = assoc7TeV_split[i];      
      double BRHZZ = myCSW->HiggsBR(11,mH);
      double BRH4l_emu = myCSW->HiggsBR(14,mH);
      double ZH    = myCSW->HiggsCS(4,mH,sqrts);
      
      //Write the output file
      fileOut7TeV << "all            ZH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ZH*BRHZZ  <<  "  1 " << ZH*BRH4l_emu << endl; 
      
    }

  fileOut7TeV << endl << endl << endl;

  //ttH 7TeV
  for( int i = 0; i < l_assoc7TeV_split; i++)
    {
      
      double mH = assoc7TeV_split[i];      
      double BRHZZ = myCSW->HiggsBR(11,mH);
      double ttH    = myCSW->HiggsCS(5,mH,sqrts);
      
      //Write the output file
      fileOut7TeV << "all            ttH"  << (int) mH <<"                      1   "  << std::fixed << setprecision(10) << ttH*BRHZZ <<  "  1" << endl; 
      
    }





  fileOut7TeV.close();

  return 0;

}
