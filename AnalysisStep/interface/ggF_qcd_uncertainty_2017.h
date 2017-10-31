#ifndef ggF_qcd_uncertainty_2017_h
#define ggF_qcd_uncertainty_2017_h

#include <vector>

typedef std::vector<double> NumV;

//
// Fractional uncertainty amplitudes of the "WG1 scheme", the "STXS scheme" and the merged "2017 scheme"
// The six first numbers are the same from each method below, namely the uncertainty amplitude of the jet bins:
// mu, res, mig01, mig12, vbf2j, vbf3j
// The last numbers are pT dependent uncertainies
NumV qcd_ggF_uncert_wg1(int Njets30, double pTH, int STXS);  // 7 nuisances, 5 x jetbin, pTH, qm_t
NumV qcd_ggF_uncert_stxs(int Njets30, double pTH, int STXS); // 8 nuisances, 5 x jetbin, D60, D120, D200
NumV qcd_ggF_uncert_2017(int Njets30, double pTH, int STXS); // 8 nuisances, 5 x jetbin, pT60, pT120, qm_t
NumV qcd_ggF_uncert_jve(int Njets30, double pTH, int STXS);  // 7 nuisances, 4 x jetbin, pT60, pT120, qm_t

// Scale factors defined as "1+uncert", where uncert is the fractional uncertainty amplitude
// This can be treated as an event weight to propagate the uncertainty to any observable/distribution.
NumV qcd_ggF_uncertSF_wg1(int Njets30, double pTH, int STXS_Stage1, double Nsigma=1.0);
NumV qcd_ggF_uncertSF_stxs(int Njets30, double pTH, int STXS_Stage1, double Nsigma=1.0);
NumV qcd_ggF_uncertSF_2017(int Njets30, double pTH, int STXS_Stage1, double Nsigma=1.0);
NumV qcd_ggF_uncertSF_jve(int Njets30, double pTH, int STXS_Stage1, double Nsigma=1.0);

NumV blptw(int Njets30);
double vbf_2j(int STXS);
double vbf_3j(int STXS);
double interpol(double x, double x1, double y1, double x2, double y2);
double qm_t(double pT);
double pT120(double pT, int Njets30);
double pT60(double pT, int Njets30);
NumV jetBinUnc(int Njets30, int STXS);

#endif
