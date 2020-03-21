#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Histograms.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/Variables.h>
#include <ZZAnalysis/AnalysisStep/test/Plotter_v2/include/CMS_lumi.h>
#include <TROOT.h>
#include <sstream>
#include <stdlib.h>
#include <iomanip>

using namespace std;

// Constructor
//====================================================
Histograms::Histograms( string blinding )
{
   _blinding = blinding;
   _merge_2e2mu = true;

   _s_process.push_back("Data");
   _s_process.push_back("H125");
   _s_process.push_back("H125ggH");
   _s_process.push_back("H125VBF");
   _s_process.push_back("H125VH");
   _s_process.push_back("H125ttH");
   _s_process.push_back("H125bbH");
   _s_process.push_back("H125tqH");
   _s_process.push_back("qqZZ");
   _s_process.push_back("ggZZ");
   _s_process.push_back("VVV");
   _s_process.push_back("Zjets");
   _s_process.push_back("other");

   _s_final_state.push_back("4e");
   _s_final_state.push_back("4mu");
   _s_final_state.push_back("2e2mu");
   _s_final_state.push_back("2mu2e");
   _s_final_state.push_back("4l");

   _s_category.push_back("UnTagged");
   _s_category.push_back("VBF1jTagged");
   _s_category.push_back("VBF2jTagged");
   _s_category.push_back("VHLeptTagged");
   _s_category.push_back("VHHadrTagged");
   _s_category.push_back("ttHLeptTagged");
   _s_category.push_back("ttHHadrTagged");
   _s_category.push_back("VHMETTagged");
   _s_category.push_back("Inclusive");

   _s_STXS_category.push_back("ggH-0j/pT[0,10]");
   _s_STXS_category.push_back("ggH-0j/pT[10-200]");
   _s_STXS_category.push_back("ggH-1j/pT[0-60]");
   _s_STXS_category.push_back("ggH-1j/pT[60-120]");
   _s_STXS_category.push_back("ggH-1j/pT[120-200]");
   _s_STXS_category.push_back("ggH-2j/pT[0-60]");
   _s_STXS_category.push_back("ggH-2j/pT[60-120]");
   _s_STXS_category.push_back("ggH-2j/pT[120-200]");
   _s_STXS_category.push_back("ggH/pT>200");
   _s_STXS_category.push_back("ggH-2j/mJJ>350");
   _s_STXS_category.push_back("VBF-1j");
   _s_STXS_category.push_back("VBF-rest");
   _s_STXS_category.push_back("VBF-2j/mJJ[350,700]");
   _s_STXS_category.push_back("VBF-2j/mJJ>700 ");
   _s_STXS_category.push_back("VBF-3j/mJJ>350 ");
   _s_STXS_category.push_back("VBF-2j/pT>200");
   _s_STXS_category.push_back("VH-had/mJJ[60-120]");
   _s_STXS_category.push_back("VH-rest");
   _s_STXS_category.push_back("VH-lep/pTV[0-150]");
   _s_STXS_category.push_back("VH-lep/pTV>150");
   _s_STXS_category.push_back("ttH-lep");
   _s_STXS_category.push_back("ttH-had");

   _s_STXS_bins.push_back("ggH-0j/pT[0,10]");
   _s_STXS_bins.push_back("ggH-0j/pT[10-200]");
   _s_STXS_bins.push_back("ggH-1j/pT[0-60]");
   _s_STXS_bins.push_back("ggH-1j/pT[60-120]");
   _s_STXS_bins.push_back("ggH-1j/pT[120-200]");
   _s_STXS_bins.push_back("ggH-2j/pT[0-60]");
   _s_STXS_bins.push_back("ggH-2j/pT[60-120]");
   _s_STXS_bins.push_back("ggH-2j/pT[120-200]");
   _s_STXS_bins.push_back("ggH/pT>200");
   _s_STXS_bins.push_back("ggH-2j/mJJ>350");
   _s_STXS_bins.push_back("qqH-2j/mJJ[60-120]");
   _s_STXS_bins.push_back("qqH-2j/pT>200 ");
   _s_STXS_bins.push_back("qqH-2j/mJJ[350,700]");
   _s_STXS_bins.push_back("qqH-2j/mJJ>700");
   _s_STXS_bins.push_back("qqH-3j/mJJ>350");
   _s_STXS_bins.push_back("qqH-rest");
   _s_STXS_bins.push_back("VH/pTV[0-150]");
   _s_STXS_bins.push_back("VH/pTV>150");
   _s_STXS_bins.push_back("bbH");
   _s_STXS_bins.push_back("ttH");
   _s_STXS_bins.push_back("tH");

    _s_STXS_bins_histoName.push_back("ggH_0j_pT_0_10");
    _s_STXS_bins_histoName.push_back("ggH_0j_pT_10_200");
    _s_STXS_bins_histoName.push_back("ggH_1j_pT_0_60");
    _s_STXS_bins_histoName.push_back("ggH_1j_pT_60_120");
    _s_STXS_bins_histoName.push_back("ggH_1j_pT_120_200");
    _s_STXS_bins_histoName.push_back("ggH_2j_pT_0_60");
    _s_STXS_bins_histoName.push_back("ggH_2j_pT_60_120");
    _s_STXS_bins_histoName.push_back("ggH_2j_pT_120_200");
    _s_STXS_bins_histoName.push_back("ggH_pT_200");
    _s_STXS_bins_histoName.push_back("ggH_2j_mJJ_350");
    _s_STXS_bins_histoName.push_back("qqH_2j_mJJ_60_120");
    _s_STXS_bins_histoName.push_back("qqH_2j_pT_200 ");
    _s_STXS_bins_histoName.push_back("qqH_2j_mJJ_350_700");
    _s_STXS_bins_histoName.push_back("qqH_2j_mJJ_700");
    _s_STXS_bins_histoName.push_back("qqH_3j_mJJ_350");
    _s_STXS_bins_histoName.push_back("qqH_rest");
    _s_STXS_bins_histoName.push_back("VH_pTV_0_150");
    _s_STXS_bins_histoName.push_back("VH_pTV_150");
    _s_STXS_bins_histoName.push_back("bbH");
    _s_STXS_bins_histoName.push_back("ttH");
    _s_STXS_bins_histoName.push_back("tH");

   _s_category_label.push_back("Untagged category");
   _s_category_label.push_back("VBF-1jet tagged category");
   _s_category_label.push_back("VBF-2jet tagged category");
   _s_category_label.push_back("VH-leptonic tagged category");
   _s_category_label.push_back("VH-hadronic tagged category");
   _s_category_label.push_back("t#bar{t}H-leptonic tagged category");
   _s_category_label.push_back("t#bar{t}H-hadronic tagged category");
   _s_category_label.push_back("VH-MET tagged category");
   _s_category_label.push_back("Inclusive category");


   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
         {
            //=====
            // M4l
            //=====
            _histo_name = "M4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::M4lMain().var_X_label + ";" + Variables::M4lMain().var_Y_label;
            histos_1D[Settings::M4lMain][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::M4lMain().var_N_bin,
                                                                               Variables::M4lMain().var_min, Variables::M4lMain().var_max);

            _histo_name = "M4l_zoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::M4lMainZoomed().var_X_label + ";" + Variables::M4lMainZoomed().var_Y_label;
            histos_1D[Settings::M4lMainZoomed][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::M4lMainZoomed().var_N_bin,
                                                                                     Variables::M4lMainZoomed().var_min, Variables::M4lMainZoomed().var_max);

            _histo_name = "M4l_HighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::M4lMainHighMass().var_X_label + ";" + Variables::M4lMainHighMass().var_Y_label;
            histos_1D[Settings::M4lMainHighMass][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::M4lMainHighMass().var_N_bin,
                                                                                       Variables::M4lMainHighMass().var_min, Variables::M4lMainHighMass().var_max);

            //=====
            // MZ1
            //=====
            _histo_name = "MZ1" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::MZ1().var_X_label + ";" + Variables::MZ1().var_Y_label;
            histos_1D[Settings::MZ1][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::MZ1().var_N_bin,
                                                                           Variables::MZ1().var_min, Variables::MZ1().var_max);

            _histo_name = "MZ1_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::MZ1_M4L118130().var_X_label + ";" + Variables::MZ1_M4L118130().var_Y_label;
            histos_1D[Settings::MZ1_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::MZ1_M4L118130().var_N_bin,
                                                                                        Variables::MZ1_M4L118130().var_min, Variables::MZ1_M4L118130().var_max);

            //=====
            // MZ2
            //=====
            _histo_name = "MZ2" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::MZ2().var_X_label + ";" + Variables::MZ2().var_Y_label;
            histos_1D[Settings::MZ2][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::MZ2().var_N_bin,
                                                                              Variables::MZ2().var_min, Variables::MZ2().var_max);

            _histo_name = "MZ2_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::MZ2_M4L118130().var_X_label + ";" + Variables::MZ2_M4L118130().var_Y_label;
            histos_1D[Settings::MZ2_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::MZ2_M4L118130().var_N_bin,
                                                                                        Variables::MZ2_M4L118130().var_min, Variables::MZ2_M4L118130().var_max);

            //====
            // KD
            //====
            _histo_name = "KD" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::KD().var_X_label + ";" + Variables::KD().var_Y_label;
            histos_1D[Settings::KD][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::KD().var_N_bin, Variables::KD().var_min,
                                                                             Variables::KD().var_max);

            _histo_name = "KD_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::KD_M4L118130().var_X_label + ";" + Variables::KD_M4L118130().var_Y_label;
            histos_1D[Settings::KD_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::KD_M4L118130().var_N_bin,
                                                                                       Variables::KD_M4L118130().var_min, Variables::KD_M4L118130().var_max);

            _histo_name = "DVBFDEC" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVBFDEC().var_X_label + ";" + Variables::DVBFDEC().var_Y_label;
            histos_1D[Settings::DVBFDEC][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVBFDEC().var_N_bin, Variables::DVBFDEC().var_min,
                                                                             Variables::DVBFDEC().var_max);

            _histo_name = "DVBFDEC_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVBFDEC_M4L118130().var_X_label + ";" + Variables::DVBFDEC_M4L118130().var_Y_label;
            histos_1D[Settings::DVBFDEC_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVBFDEC_M4L118130().var_N_bin,
                                                                                       Variables::DVBFDEC_M4L118130().var_min, Variables::DVBFDEC_M4L118130().var_max);

            _histo_name = "DVHDEC" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVHDEC().var_X_label + ";" + Variables::DVHDEC().var_Y_label;
            histos_1D[Settings::DVHDEC][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVHDEC().var_N_bin, Variables::DVHDEC().var_min,
                                                                             Variables::DVHDEC().var_max);

            _histo_name = "DVHDEC_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVHDEC_M4L118130().var_X_label + ";" + Variables::DVHDEC_M4L118130().var_Y_label;
            histos_1D[Settings::DVHDEC_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVHDEC_M4L118130().var_N_bin,
                                                                                       Variables::DVHDEC_M4L118130().var_min, Variables::DVHDEC_M4L118130().var_max);

            _histo_name = "D1jet" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::D1jet().var_X_label + ";" + Variables::D1jet().var_Y_label;
            histos_1D[Settings::D1jet][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::D1jet().var_N_bin,
                                                                                Variables::D1jet().var_min, Variables::D1jet().var_max);

            _histo_name = "D1jet_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::D1jet_M4L118130().var_X_label + ";" + Variables::D1jet_M4L118130().var_Y_label;
            histos_1D[Settings::D1jet_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::D1jet_M4L118130().var_N_bin,
                                                                                          Variables::D1jet_M4L118130().var_min, Variables::D1jet_M4L118130().var_max);

            _histo_name = "D2jet" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::D2jet().var_X_label + ";" + Variables::D2jet().var_Y_label;
            histos_1D[Settings::D2jet][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::D2jet().var_N_bin,
                                                                                Variables::D2jet().var_min, Variables::D2jet().var_max);

            _histo_name = "D2jet_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::D2jet_M4L118130().var_X_label + ";" + Variables::D2jet_M4L118130().var_Y_label;
            histos_1D[Settings::D2jet_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::D2jet_M4L118130().var_N_bin,
                                                                                          Variables::D2jet_M4L118130().var_min, Variables::D2jet_M4L118130().var_max);

            _histo_name = "DWH" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DWH().var_X_label + ";" + Variables::DWH().var_Y_label;
            histos_1D[Settings::DWH][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DWH().var_N_bin,
                                                                              Variables::DWH().var_min, Variables::DWH().var_max);

            _histo_name = "DWH_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DWH_M4L118130().var_X_label + ";" + Variables::DWH_M4L118130().var_Y_label;
            histos_1D[Settings::DWH_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DWH_M4L118130().var_N_bin,
                                                                                     Variables::DWH_M4L118130().var_min, Variables::DWH_M4L118130().var_max);

            _histo_name = "DZH" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DZH().var_X_label + ";" + Variables::DZH().var_Y_label;
            histos_1D[Settings::DZH][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DZH().var_N_bin,
                                                                     Variables::DZH().var_min, Variables::DZH().var_max);

            _histo_name = "DZH_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DZH_M4L118130().var_X_label + ";" + Variables::DZH_M4L118130().var_Y_label;
            histos_1D[Settings::DZH_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DZH_M4L118130().var_N_bin,
                                                                               Variables::DZH_M4L118130().var_min, Variables::DZH_M4L118130().var_max);

            _histo_name = "DVH" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVH().var_X_label + ";" + Variables::DVH().var_Y_label;
            histos_1D[Settings::DVH][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVH().var_N_bin,
                                                                     Variables::DVH().var_min, Variables::DVH().var_max);

            _histo_name = "DVH_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVH_M4L118130().var_X_label + ";" + Variables::DVH_M4L118130().var_Y_label;
            histos_1D[Settings::DVH_M4L118130][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVH_M4L118130().var_N_bin,
                                                                               Variables::DVH_M4L118130().var_min, Variables::DVH_M4L118130().var_max);



            //===========
            // MZ1 vs MZ2
            //===========
            _histo_name = "MZ1vsMZ2" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::MZ1vsMZ2().var_X_label + ";" + Variables::MZ1vsMZ2().var_Y_label;
            histos_2D[Settings::MZ1vsMZ2][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::MZ1vsMZ2().var_X_N_bin,
                                                                                Variables::MZ1vsMZ2().var_X_min, Variables::MZ1vsMZ2().var_X_max, Variables::MZ1vsMZ2().var_Y_N_bin,
                                                                                Variables::MZ1vsMZ2().var_Y_min, Variables::MZ1vsMZ2().var_Y_max);

            _histo_name = "MZ1vsMZ2_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::MZ1vsMZ2_M4L118130().var_X_label + ";" + Variables::MZ1vsMZ2_M4L118130().var_Y_label;
            histos_2D[Settings::MZ1vsMZ2_M4L118130][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::MZ1vsMZ2_M4L118130().var_X_N_bin,
                                                                                          Variables::MZ1vsMZ2_M4L118130().var_X_min, Variables::MZ1vsMZ2_M4L118130().var_X_max,
                                                                                          Variables::MZ1vsMZ2_M4L118130().var_Y_N_bin, Variables::MZ1vsMZ2_M4L118130().var_Y_min,
                                                                                          Variables::MZ1vsMZ2_M4L118130().var_Y_max);

            //===========
            // KD vs M4l
            //===========
            _histo_name = "KDvsM4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::KDvsM4l().var_X_label + ";" + Variables::KDvsM4l().var_Y_label;
            histos_2DError[Settings::KDvsM4l][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::KDvsM4l().var_X_N_bin,
                                                                              Variables::KDvsM4l().var_X_min, Variables::KDvsM4l().var_X_max,
                                                                              Variables::KDvsM4l().var_Y_N_bin, Variables::KDvsM4l().var_Y_min,
                                                                              Variables::KDvsM4l().var_Y_max);

            _histo_name = "KDvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::KDvsM4lZoomed().var_X_label + ";" + Variables::KDvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::KDvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::KDvsM4lZoomed().var_X_N_bin,
                                                                                    Variables::KDvsM4lZoomed().var_X_min, Variables::KDvsM4lZoomed().var_X_max,
                                                                                    Variables::KDvsM4lZoomed().var_Y_N_bin, Variables::KDvsM4lZoomed().var_Y_min,
                                                                                    Variables::KDvsM4lZoomed().var_Y_max);

            _histo_name = "KDvsM4lHighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::KDvsM4lHighMass().var_X_label + ";" + Variables::KDvsM4lHighMass().var_Y_label;
            histos_2DError[Settings::KDvsM4lHighMass][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::KDvsM4lHighMass().var_X_N_bin,
                                                                                      Variables::KDvsM4lHighMass().var_X_min, Variables::KDvsM4lHighMass().var_X_max,
                                                                                      Variables::KDvsM4lHighMass().var_Y_N_bin, Variables::KDvsM4lHighMass().var_Y_min,
                                                                                      Variables::KDvsM4l().var_Y_max);

            _histo_name = "DVBFDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVBFDECvsM4lZoomed().var_X_label + ";" + Variables::DVBFDECvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVBFDECvsM4lZoomed().var_X_N_bin,
                                                                                    Variables::DVBFDECvsM4lZoomed().var_X_min, Variables::DVBFDECvsM4lZoomed().var_X_max,
                                                                                    Variables::DVBFDECvsM4lZoomed().var_Y_N_bin, Variables::DVBFDECvsM4lZoomed().var_Y_min,
                                                                                    Variables::DVBFDECvsM4lZoomed().var_Y_max);

            _histo_name = "DVHDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVHDECvsM4lZoomed().var_X_label + ";" + Variables::DVHDECvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::DVHDECvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVHDECvsM4lZoomed().var_X_N_bin,
                                                                                    Variables::DVHDECvsM4lZoomed().var_X_min, Variables::DVHDECvsM4lZoomed().var_X_max,
                                                                                    Variables::DVHDECvsM4lZoomed().var_Y_N_bin, Variables::DVHDECvsM4lZoomed().var_Y_min,
                                                                                    Variables::DVHDECvsM4lZoomed().var_Y_max);

            _histo_name = "D1jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::D1jetvsM4lZoomed().var_X_label + ";" + Variables::D1jetvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::D1jetvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::D1jetvsM4lZoomed().var_X_N_bin,
                                                                                       Variables::D1jetvsM4lZoomed().var_X_min, Variables::D1jetvsM4lZoomed().var_X_max,
                                                                                       Variables::D1jetvsM4lZoomed().var_Y_N_bin, Variables::D1jetvsM4lZoomed().var_Y_min,
                                                                                       Variables::D1jetvsM4lZoomed().var_Y_max);

            _histo_name = "D2jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::D2jetvsM4lZoomed().var_X_label + ";" + Variables::D2jetvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::D2jetvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::D2jetvsM4lZoomed().var_X_N_bin,
                                                                                       Variables::D2jetvsM4lZoomed().var_X_min, Variables::D2jetvsM4lZoomed().var_X_max,
                                                                                       Variables::D2jetvsM4lZoomed().var_Y_N_bin, Variables::D2jetvsM4lZoomed().var_Y_min,
                                                                                       Variables::D2jetvsM4lZoomed().var_Y_max);

            _histo_name = "DWHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DWHvsM4lZoomed().var_X_label + ";" + Variables::DWHvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::DWHvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DWHvsM4lZoomed().var_X_N_bin,
                                                                                     Variables::DWHvsM4lZoomed().var_X_min, Variables::DWHvsM4lZoomed().var_X_max,
                                                                                     Variables::DWHvsM4lZoomed().var_Y_N_bin, Variables::DWHvsM4lZoomed().var_Y_min,
                                                                                     Variables::DWHvsM4lZoomed().var_Y_max);

            _histo_name = "DZHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DZHvsM4lZoomed().var_X_label + ";" + Variables::DZHvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::DZHvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DZHvsM4lZoomed().var_X_N_bin,
                                                                                     Variables::DZHvsM4lZoomed().var_X_min, Variables::DZHvsM4lZoomed().var_X_max,
                                                                                     Variables::DZHvsM4lZoomed().var_Y_N_bin, Variables::DZHvsM4lZoomed().var_Y_min,
                                                                                     Variables::DZHvsM4lZoomed().var_Y_max);

            _histo_name = "DVHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::DVHvsM4lZoomed().var_X_label + ";" + Variables::DVHvsM4lZoomed().var_Y_label;
            histos_2DError[Settings::DVHvsM4lZoomed][i_fs][i_cat][i_proc] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVHvsM4lZoomed().var_X_N_bin,
                                                                                     Variables::DVHvsM4lZoomed().var_X_min, Variables::DVHvsM4lZoomed().var_X_max,
                                                                                     Variables::DVHvsM4lZoomed().var_Y_N_bin, Variables::DVHvsM4lZoomed().var_Y_min,
                                                                                     Variables::DVHvsM4lZoomed().var_Y_max);

            //===========
            // Others
            //===========
            _histo_name = "PFMET" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::PFMET().var_X_label + ";" + Variables::PFMET().var_Y_label;
            histos_1D[Settings::PFMET][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::PFMET().var_N_bin,
                                                                                        Variables::PFMET().var_min, Variables::PFMET().var_max);

            _histo_name = "Pt_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::Pt_leading().var_X_label + ";" + Variables::Pt_leading().var_Y_label;
            histos_1D[Settings::Pt_leading][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::Pt_leading().var_N_bin,
                                                                                        Variables::Pt_leading().var_min, Variables::Pt_leading().var_max);

            _histo_name = "Pt_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::Pt_trailing().var_X_label + ";" + Variables::Pt_trailing().var_Y_label;
            histos_1D[Settings::Pt_trailing][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::Pt_trailing().var_N_bin,
                                                                                        Variables::Pt_trailing().var_min, Variables::Pt_trailing().var_max);

            _histo_name = "Eta_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::Eta_leading().var_X_label + ";" + Variables::Eta_leading().var_Y_label;
            histos_1D[Settings::Eta_leading][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::Eta_leading().var_N_bin,
                                                                                        Variables::Eta_leading().var_min, Variables::Eta_leading().var_max);

            _histo_name = "Eta_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::Eta_trailing().var_X_label + ";" + Variables::Eta_trailing().var_Y_label;
            histos_1D[Settings::Eta_trailing][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::Eta_trailing().var_N_bin,
                                                                                        Variables::Eta_trailing().var_min, Variables::Eta_trailing().var_max);

            _histo_name = "SIP_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::SIP_leading().var_X_label + ";" + Variables::SIP_leading().var_Y_label;
            histos_1D[Settings::SIP_leading][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::SIP_leading().var_N_bin,
                                                                                        Variables::SIP_leading().var_min, Variables::SIP_leading().var_max);

            _histo_name = "SIP_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::SIP_trailing().var_X_label + ";" + Variables::SIP_trailing().var_Y_label;
            histos_1D[Settings::SIP_trailing][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::SIP_trailing().var_N_bin,
                                                                                        Variables::SIP_trailing().var_min, Variables::SIP_trailing().var_max);

            _histo_name = "ISO_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::ISO_leading().var_X_label + ";" + Variables::ISO_leading().var_Y_label;
            histos_1D[Settings::ISO_leading][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::ISO_leading().var_N_bin,
                                                                                        Variables::ISO_leading().var_min, Variables::ISO_leading().var_max);

            _histo_name = "ISO_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::ISO_trailing().var_X_label + ";" + Variables::ISO_trailing().var_Y_label;
            histos_1D[Settings::ISO_trailing][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::ISO_trailing().var_N_bin,
                                                                                        Variables::ISO_trailing().var_min, Variables::ISO_trailing().var_max);

            _histo_name = "Pt4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::Pt4l().var_X_label + ";" + Variables::Pt4l().var_Y_label;
            histos_1D[Settings::Pt4l][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::Pt4l().var_N_bin,
                                                                                        Variables::Pt4l().var_min, Variables::Pt4l().var_max);

            _histo_name = "NJetsBTagged" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::NJetsBTagged().var_X_label + ";" + Variables::NJetsBTagged().var_Y_label;
            histos_1D[Settings::NJetsBTagged][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::NJetsBTagged().var_N_bin,
                                                                                        Variables::NJetsBTagged().var_min, Variables::NJetsBTagged().var_max);

            _histo_name = "Eta4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::Eta4l().var_X_label + ";" + Variables::Eta4l().var_Y_label;
            histos_1D[Settings::Eta4l][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::Eta4l().var_N_bin,
                                                                                        Variables::Eta4l().var_min, Variables::Eta4l().var_max);

            _histo_name = "NExtraLep" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::NExtraLep().var_X_label + ";" + Variables::NExtraLep().var_Y_label;
            histos_1D[Settings::NExtraLep][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::NExtraLep().var_N_bin,
                                                                                        Variables::NExtraLep().var_min, Variables::NExtraLep().var_max);

            _histo_name = "NJets" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::NJets().var_X_label + ";" + Variables::NJets().var_Y_label;
            histos_1D[Settings::NJets][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::NJets().var_N_bin,
                                                                                        Variables::NJets().var_min, Variables::NJets().var_max);

            _histo_name = "M4l_110150_HighKD" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            _histo_labels = ";" + Variables::M4l_110150_HighKD().var_X_label + ";" + Variables::M4l_110150_HighKD().var_Y_label;
            histos_1D[Settings::M4l_110150_HighKD][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::M4l_110150_HighKD().var_N_bin,
                                                                                        Variables::M4l_110150_HighKD().var_min, Variables::M4l_110150_HighKD().var_max);
         }
      }
   }


   // Z+X
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         //=====
         // M4l
         //=====
         _histo_name = "M4l_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::M4lMain][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4lMain().var_N_bin, Variables::M4lMain().var_min, Variables::M4lMain().var_max);

         _histo_name = "M4l_zoomed_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4lMainZoomed().var_N_bin,
                                                                       Variables::M4lMainZoomed().var_min, Variables::M4lMainZoomed().var_max);

         _histo_name = "M4l_HighMass_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4lMainHighMass().var_N_bin,
                                                                         Variables::M4lMainHighMass().var_min, Variables::M4lMainHighMass().var_max);

         _histo_name = "M4l_ZX_shape_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX_shape[Settings::M4lMain][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4lMain().var_N_bin,
                                                                       Variables::M4lMain().var_min, Variables::M4lMain().var_max);

         _histo_name = "M4l_zoomed_ZX_shape_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX_shape[Settings::M4lMainZoomed][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4lMainZoomed().var_N_bin,
                                                                             Variables::M4lMainZoomed().var_min, Variables::M4lMainZoomed().var_max);

         _histo_name = "M4l_HighMass_ZX_shape_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX_shape[Settings::M4lMainHighMass][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4lMainHighMass().var_N_bin,
                                                                               Variables::M4lMainHighMass().var_min, Variables::M4lMainHighMass().var_max);


         //=====
         // MZ1
         //=====
         _histo_name = "MZ1_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ1][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::MZ1().var_N_bin, Variables::MZ1().var_min, Variables::MZ1().var_max);

         _histo_name = "MZ1_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::MZ1_M4L118130().var_N_bin,
                                                                       Variables::MZ1_M4L118130().var_min, Variables::MZ1_M4L118130().var_max);

         //=====
         // MZ2
         //=====
         _histo_name = "MZ2_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ2][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::MZ2().var_N_bin, Variables::MZ2().var_min, Variables::MZ2().var_max);

         _histo_name = "MZ2_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::MZ2_M4L118130().var_N_bin,
                                                                       Variables::MZ2_M4L118130().var_min, Variables::MZ2_M4L118130().var_max);

         //====
         // KD
         //====
         _histo_name = "KD_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::KD][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::KD().var_N_bin, Variables::KD().var_min, Variables::KD().var_max);

         _histo_name = "KD_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::KD_M4L118130().var_N_bin,
                                                                      Variables::KD_M4L118130().var_min, Variables::KD_M4L118130().var_max);

         _histo_name = "DVBFDEC_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVBFDEC][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DVBFDEC().var_N_bin, Variables::DVBFDEC().var_min,
                                                                             Variables::DVBFDEC().var_max);

      _histo_name = "DVBFDEC_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVBFDEC_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DVBFDEC_M4L118130().var_N_bin,
                                                                                       Variables::DVBFDEC_M4L118130().var_min, Variables::DVBFDEC_M4L118130().var_max);

      _histo_name = "DVHDEC_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVHDEC][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DVHDEC().var_N_bin, Variables::DVHDEC().var_min,
                                                                             Variables::DVHDEC().var_max);

      _histo_name = "DVHDEC_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVHDEC_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DVHDEC_M4L118130().var_N_bin,
                                                                                       Variables::DVHDEC_M4L118130().var_min, Variables::DVHDEC_M4L118130().var_max);

         _histo_name = "D1jet_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D1jet][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::D1jet().var_N_bin,
                                                               Variables::D1jet().var_min, Variables::D1jet().var_max);

         _histo_name = "D1jet_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::D1jet_M4L118130().var_N_bin,
                                                                         Variables::D1jet_M4L118130().var_min, Variables::D1jet_M4L118130().var_max);

         _histo_name = "D2jet_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D2jet][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::D2jet().var_N_bin,
                                                               Variables::D2jet().var_min, Variables::D2jet().var_max);

         _histo_name = "D2jet_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::D2jet_M4L118130().var_N_bin,
                                                                         Variables::D2jet_M4L118130().var_min, Variables::D2jet_M4L118130().var_max);

         _histo_name = "DWH_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DWH][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DWH().var_N_bin,
                                                             Variables::DWH().var_min, Variables::DWH().var_max);

         _histo_name = "DWH_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DWH_M4L118130().var_N_bin,
                                                                       Variables::DWH_M4L118130().var_min, Variables::DWH_M4L118130().var_max);

         _histo_name = "DZH_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DZH][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DZH().var_N_bin,
                                                             Variables::DZH().var_min, Variables::DZH().var_max);

         _histo_name = "DZH_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DZH_M4L118130().var_N_bin,
                                                                       Variables::DZH_M4L118130().var_min, Variables::DZH_M4L118130().var_max);

         _histo_name = "DVH_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DVH][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DVH().var_N_bin,
                                                             Variables::DVH().var_min, Variables::DVH().var_max);

         _histo_name = "DVH_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DVH_M4L118130][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::DVH_M4L118130().var_N_bin,
                                                                       Variables::DVH_M4L118130().var_min, Variables::DVH_M4L118130().var_max);

      //===========
      // KD vs M4l
      //===========
      _histo_name = "KDvsM4lZoomed_ZX" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _blinding;
      _histo_labels = ";" + Variables::KDvsM4lZoomed().var_X_label + ";" + Variables::KDvsM4lZoomed().var_Y_label;
      histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][i_cat] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::KDvsM4lZoomed().var_X_N_bin,
                                                      Variables::KDvsM4lZoomed().var_X_min, Variables::KDvsM4lZoomed().var_X_max,
                                                      Variables::KDvsM4lZoomed().var_Y_N_bin, Variables::KDvsM4lZoomed().var_Y_min,
                                                      Variables::KDvsM4lZoomed().var_Y_max);

      _histo_name = "DVBFDECvsM4lZoomed_ZX" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _blinding;
      _histo_labels = ";" + Variables::DVBFDECvsM4lZoomed().var_X_label + ";" + Variables::DVBFDECvsM4lZoomed().var_Y_label;
      histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVBFDECvsM4lZoomed().var_X_N_bin,
                                                      Variables::DVBFDECvsM4lZoomed().var_X_min, Variables::DVBFDECvsM4lZoomed().var_X_max,
                                                      Variables::DVBFDECvsM4lZoomed().var_Y_N_bin, Variables::DVBFDECvsM4lZoomed().var_Y_min,
                                                      Variables::DVBFDECvsM4lZoomed().var_Y_max);

      _histo_name = "DVHDECvsM4lZoomed_ZX" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _blinding;
      _histo_labels = ";" + Variables::DVHDECvsM4lZoomed().var_X_label + ";" + Variables::DVHDECvsM4lZoomed().var_Y_label;
      histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat] = new TH2F(_histo_name.c_str(), _histo_labels.c_str(), Variables::DVHDECvsM4lZoomed().var_X_N_bin,
                                                      Variables::DVHDECvsM4lZoomed().var_X_min, Variables::DVHDECvsM4lZoomed().var_X_max,
                                                      Variables::DVHDECvsM4lZoomed().var_Y_N_bin, Variables::DVHDECvsM4lZoomed().var_Y_min,
                                                      Variables::DVHDECvsM4lZoomed().var_Y_max);

      //===========
      // Others
      //===========
      _histo_name = "PFMET_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::PFMET][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::PFMET().var_N_bin,
                                                Variables::PFMET().var_min, Variables::PFMET().var_max);

      _histo_name = "Pt_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Pt_leading][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::Pt_leading().var_N_bin,
                                                Variables::Pt_leading().var_min, Variables::Pt_leading().var_max);

      _histo_name = "Pt_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Pt_trailing][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::Pt_trailing().var_N_bin,
                                                Variables::Pt_trailing().var_min, Variables::Pt_trailing().var_max);

      _histo_name = "Eta_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Eta_leading][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::Eta_leading().var_N_bin,
                                                Variables::Eta_leading().var_min, Variables::Eta_leading().var_max);

      _histo_name = "Eta_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Eta_trailing][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::Eta_trailing().var_N_bin,
                                                Variables::Eta_trailing().var_min, Variables::Eta_trailing().var_max);

      _histo_name = "SIP_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::SIP_leading][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::SIP_leading().var_N_bin,
                                                Variables::SIP_leading().var_min, Variables::SIP_leading().var_max);

      _histo_name = "SIP_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::SIP_trailing][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::SIP_trailing().var_N_bin,
                                                Variables::SIP_trailing().var_min, Variables::SIP_trailing().var_max);

      _histo_name = "ISO_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::ISO_leading][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::ISO_leading().var_N_bin,
                                                Variables::ISO_leading().var_min, Variables::ISO_leading().var_max);

      _histo_name = "ISO_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::ISO_trailing][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::ISO_trailing().var_N_bin,
                                                Variables::ISO_trailing().var_min, Variables::ISO_trailing().var_max);

      _histo_name = "Pt4l_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Pt4l][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::Pt4l().var_N_bin,
                                                Variables::Pt4l().var_min, Variables::Pt4l().var_max);

      _histo_name = "NJetsBTagged_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::NJetsBTagged][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::NJetsBTagged().var_N_bin,
                                                Variables::NJetsBTagged().var_min, Variables::NJetsBTagged().var_max);

      _histo_name = "Eta4l_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Eta4l][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::Eta4l().var_N_bin,
                                                Variables::Eta4l().var_min, Variables::Eta4l().var_max);

      _histo_name = "NExtraLep_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::NExtraLep][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::NExtraLep().var_N_bin,
                                                Variables::NExtraLep().var_min, Variables::NExtraLep().var_max);

      _histo_name = "NJets_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::NJets][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::NJets().var_N_bin,
                                                Variables::NJets().var_min, Variables::NJets().var_max);

      _histo_name = "M4l_110150_HighKD_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::M4l_110150_HighKD][i_fs][i_cat] = new TH1F(_histo_name.c_str(), "Z+X", Variables::M4l_110150_HighKD().var_N_bin,
                                                Variables::M4l_110150_HighKD().var_min, Variables::M4l_110150_HighKD().var_max);
      }
   }

    // STXS category data vs MC
    for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
    {
        _histo_name = "STXS_Categories_" + _s_process.at(i_proc) + "_" + _blinding;
        _histo_labels = ";" + Variables::STXS_Categories().var_X_label + ";" + Variables::STXS_Categories().var_Y_label;
        STXS_Categories[i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::STXS_Categories().var_N_bin,Variables::STXS_Categories().var_min, Variables::STXS_Categories().var_max);
        for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ ) STXS_Categories[i_proc]->GetXaxis()->SetBinLabel(i_cat + 1,_s_STXS_category.at(i_cat));

    }

    for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ )
    {
        for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
        {
            _histo_name = "Yields_" + _s_process.at(i_proc) + "_" + _s_STXS_category.at(i_cat) + "_" + _blinding;
            STXS_Yields[i_cat][i_proc] = new TH1F(_histo_name.c_str(),_histo_name.c_str(), Variables::STXS_Categories().M4lBins,Variables::STXS_Categories().M4lMin, Variables::STXS_Categories().M4lMax);
        }

        _histo_name = "YieldsZX_" + _s_STXS_category.at(i_cat) + "_" + _blinding;
        STXS_Yields_ZX[i_cat] = new TH1F(_histo_name.c_str(),_histo_name.c_str(), Variables::STXS_Categories().M4lBins,Variables::STXS_Categories().M4lMin, Variables::STXS_Categories().M4lMax);
    }

    // STXS purity
    for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
    {
        _histo_name = "STXS_Purity_" + _s_STXS_bins_histoName.at(i_proc) + "_" + _blinding;
        _histo_labels = ";" + Variables::STXS_Categories().var_X_label + ";" + Variables::STXS_Categories().var_Y_label;
        Purity_Categories[i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::STXS_Categories().var_N_bin,Variables::STXS_Categories().var_min, Variables::STXS_Categories().var_max);

    }

    for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ )
    {
        for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
        {
            _histo_name = "Purity_Yields_" + _s_STXS_bins.at(i_proc) + "_" + _s_STXS_category.at(i_cat) + "_" + _blinding;
            Purity_Yields[i_cat][i_proc] = new TH1F(_histo_name.c_str(),_histo_name.c_str(), Variables::STXS_Categories().M4lBins,Variables::STXS_Categories().M4lMin, Variables::STXS_Categories().M4lMax);
        }

        Purity_Normalization[i_cat] = 0.;
    }


}
//==================================


//==================================
//Constructor
//==================================
Histograms::Histograms( double lumi)
{
   _lumi = lumi;
   _blinding = "Unblinded";
   _merge_2e2mu = true;

   _s_process.push_back("Data");
   _s_process.push_back("H120");
   _s_process.push_back("H124");
   _s_process.push_back("H125");
   _s_process.push_back("H126");
   _s_process.push_back("H130");

   _s_process.push_back("H120ggH");
   _s_process.push_back("H124ggH");
   _s_process.push_back("H125ggH");
   _s_process.push_back("H126ggH");
   _s_process.push_back("H130ggH");

   _s_process.push_back("H120VBF");
   _s_process.push_back("H124VBF");
   _s_process.push_back("H125VBF");
   _s_process.push_back("H126VBF");
   _s_process.push_back("H130VBF");

   _s_process.push_back("H120WHlep");
   _s_process.push_back("H124WHlep");
   _s_process.push_back("H125WHlep");
   _s_process.push_back("H126WHlep");
   _s_process.push_back("H130WHlep");

  _s_process.push_back("H120WHhad");
  _s_process.push_back("H124WHhad");
  _s_process.push_back("H125WHhad");
  _s_process.push_back("H126WHhad");
  _s_process.push_back("H130WHhad");

   _s_process.push_back("H120ZHlep");
   _s_process.push_back("H124ZHlep");
   _s_process.push_back("H125ZHlep");
   _s_process.push_back("H126ZHlep");
   _s_process.push_back("H130ZHlep");

  _s_process.push_back("H120ZHhad");
  _s_process.push_back("H124ZHhad");
  _s_process.push_back("H125ZHhad");
  _s_process.push_back("H126ZHhad");
  _s_process.push_back("H130ZHhad");

   _s_process.push_back("H120ttHlep");
   _s_process.push_back("H124ttHlep");
   _s_process.push_back("H125ttHlep");
   _s_process.push_back("H126ttHlep");
   _s_process.push_back("H130ttHlep");

  _s_process.push_back("H120ttHhad");
   _s_process.push_back("H124ttHhad");
   _s_process.push_back("H125ttHhad");
   _s_process.push_back("H126ttHhad");
   _s_process.push_back("H130ttHhad");

  _s_process.push_back("H120bbH");
   _s_process.push_back("H124bbH");
   _s_process.push_back("H125bbH");
   _s_process.push_back("H126bbH");
   _s_process.push_back("H130bbH");

  _s_process.push_back("H120tqH");
   _s_process.push_back("H124tqH");
   _s_process.push_back("H125tqH");
   _s_process.push_back("H126tqH");
   _s_process.push_back("H130tqH");

   _s_process.push_back("qqZZ");
   _s_process.push_back("ggZZ");
   _s_process.push_back("VVV");
   _s_process.push_back("DY");
   _s_process.push_back("ttbar");

   _s_final_state.push_back("4e");
   _s_final_state.push_back("4mu");
   _s_final_state.push_back("2e2mu");
   _s_final_state.push_back("2mu2e");
   _s_final_state.push_back("4l");

   _s_category.push_back("UnTagged");
   _s_category.push_back("VBF1jTagged");
   _s_category.push_back("VBF2jTagged");
   _s_category.push_back("VHLeptTagged");
   _s_category.push_back("VHHadrTagged");
   _s_category.push_back("ttHLeptTagged");
   _s_category.push_back("ttHHadrTagged");
   _s_category.push_back("VHMETTagged");
   _s_category.push_back("Inclusive");

   _s_production_mode.push_back("ggH");
   _s_production_mode.push_back("qqH");
   _s_production_mode.push_back("WH_lep");
   _s_production_mode.push_back("WH_had");
   _s_production_mode.push_back("ZH_lep");
   _s_production_mode.push_back("ZH_had");
   _s_production_mode.push_back("ttH_lep");
   _s_production_mode.push_back("ttH_had");
   _s_production_mode.push_back("bbH");
   _s_production_mode.push_back("tqH");
   _s_production_mode.push_back("qqZZ");
   _s_production_mode.push_back("ggZZ");
   _s_production_mode.push_back("ttZZ");
   _s_production_mode.push_back("ttWW");
   _s_production_mode.push_back("WWZ");
   _s_production_mode.push_back("WZZ");
   _s_production_mode.push_back("ZZZ");
   _s_production_mode.push_back("TTZJets");   


   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
         {
            _histo_name = "M4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + "_" + _blinding;
            _histo_labels = ";" + Variables::M4lYields().var_X_label + ";" + Variables::M4lYields().var_Y_label;
            histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc] = new TH1F(_histo_name.c_str(), _histo_labels.c_str(), Variables::M4lYields().var_N_bin,
                                                                                    Variables::M4lYields().var_min, Variables::M4lYields().var_max);
         }
      }
   }
}
//==================================



//=======================
Histograms::~Histograms()
{
}
//=======================



//====
// M4l
//====================================================================================
void Histograms::FillM4l( float M4l, float weight, int fs, int cat, int proc )
{
  histos_1D[Settings::M4lMain][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::M4lMainZoomed][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::M4lMainHighMass][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
}
//====================================================================================


//====================================================================
void Histograms::FillM4lZX( float M4l, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::M4lMain][fs][cat]->Fill(M4l, weight);
   histos_1D_ZX[Settings::M4lMainZoomed][fs][cat]->Fill(M4l, weight);
   histos_1D_ZX[Settings::M4lMainHighMass][fs][cat]->Fill(M4l, weight);
}
//====================================================================

//====
// STXS
//====================================================================================
void Histograms::FillSTXS( float M4l, float weight, int cat_stxs, int proc )
{
    STXS_Yields[cat_stxs][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
}
//====================================================================================


//====================================================================================
void Histograms::FillSTXSZX( float M4l, float weight, int cat_stxs )
{
    STXS_Yields_ZX[cat_stxs]->Fill(M4l, weight);
}
//====================================================================================

//====================================================================================
void Histograms::FillSTXSPurity( float M4l, float weight, int cat_stxs, int bin_stxs )
{
    Purity_Yields[cat_stxs][bin_stxs]->Fill(M4l, weight);
}
//====================================================================================





//====
// MZ1
//====================================================================================
void Histograms::FillMZ1( float M4l, float MZ1, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::MZ1][fs][cat][proc]->Fill(MZ1, (proc == Settings::Data) ? 1. : weight);

   if ( M4l >= Variables::MZ1_M4L118130().cut_d && M4l <= Variables::MZ1_M4L118130().cut_u )
   {
      histos_1D[Settings::MZ1_M4L118130][fs][cat][proc]->Fill(MZ1, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillMZ1ZX( float M4l, float MZ1, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::MZ1][fs][cat]->Fill(MZ1, weight);

   if( M4l >= Variables::MZ1_M4L118130().cut_d && M4l <= Variables::MZ1_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::MZ1_M4L118130][fs][cat]->Fill(MZ1, weight);
   }
}
//====================================================================



//====
// MZ2
//===============================================================================================
void Histograms::FillMZ2( float M4l, float MZ2, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::MZ2][fs][cat][proc]->Fill(MZ2, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::MZ2_M4L118130().cut_d && M4l <= Variables::MZ2_M4L118130().cut_u)
   {
      histos_1D[Settings::MZ2_M4L118130][fs][cat][proc]->Fill(MZ2, (proc == Settings::Data) ? 1. : weight);
   }
}
//===============================================================================================

//===============================================================================
void Histograms::FillMZ2ZX( float M4l, float MZ2, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::MZ2][fs][cat]->Fill(MZ2, weight);

   if( M4l >= Variables::MZ2_M4L118130().cut_d && M4l <= Variables::MZ2_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::MZ2_M4L118130][fs][cat]->Fill(MZ2, weight);
   }
}
//===============================================================================



//===
// KD
//===
//====================================================================================
void Histograms::FillKD( float M4l, float KD, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::KD][fs][cat][proc]->Fill(KD, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::KD_M4L118130().cut_d && M4l <= Variables::KD_M4L118130().cut_u)
   {
      histos_1D[Settings::KD_M4L118130][fs][cat][proc]->Fill(KD, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillKDZX( float M4l, float KD, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::KD][fs][cat]->Fill(KD, weight);

   if( M4l >= Variables::KD_M4L118130().cut_d && M4l <= Variables::KD_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::KD_M4L118130][fs][cat]->Fill(KD, weight);
   }
}
//====================================================================

//====================================================================================
void Histograms::FillDVBFDEC( float M4l, float DVBFDEC, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::DVBFDEC][fs][cat][proc]->Fill(DVBFDEC, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::DVBFDEC_M4L118130().cut_d && M4l <= Variables::DVBFDEC_M4L118130().cut_u)
   {
      histos_1D[Settings::DVBFDEC_M4L118130][fs][cat][proc]->Fill(DVBFDEC, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillDVBFDECZX( float M4l, float DVBFDEC, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::DVBFDEC][fs][cat]->Fill(DVBFDEC, weight);

   if( M4l >= Variables::DVBFDEC_M4L118130().cut_d && M4l <= Variables::DVBFDEC_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::DVBFDEC_M4L118130][fs][cat]->Fill(DVBFDEC, weight);
   }
}
//====================================================================

//====================================================================================
void Histograms::FillDVHDEC( float M4l, float DVHDEC, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::DVHDEC][fs][cat][proc]->Fill(DVHDEC, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::DVHDEC_M4L118130().cut_d && M4l <= Variables::DVHDEC_M4L118130().cut_u)
   {
      histos_1D[Settings::DVHDEC_M4L118130][fs][cat][proc]->Fill(DVHDEC, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillDVHDECZX( float M4l, float DVHDEC, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::DVHDEC][fs][cat]->Fill(DVHDEC, weight);

   if( M4l >= Variables::DVHDEC_M4L118130().cut_d && M4l <= Variables::DVHDEC_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::DVHDEC_M4L118130][fs][cat]->Fill(DVHDEC, weight);
   }
}
//====================================================================



//======
// D1jet
//====================================================================================
void Histograms::FillD1jet( float M4l, float D1jet, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::D1jet][fs][cat][proc]->Fill(D1jet, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::D1jet_M4L118130().cut_d && M4l <= Variables::D1jet_M4L118130().cut_u)
   {
      histos_1D[Settings::D1jet_M4L118130][fs][cat][proc]->Fill(D1jet, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillD1jetZX( float M4l, float D1jet, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::D1jet][fs][cat]->Fill(D1jet, weight);

   if( M4l >= Variables::D1jet_M4L118130().cut_d && M4l <= Variables::D1jet_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::D1jet_M4L118130][fs][cat]->Fill(D1jet, weight);
   }
}
//====================================================================



//======
// D2jet
//====================================================================================
void Histograms::FillD2jet( float M4l, float D2jet, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::D2jet][fs][cat][proc]->Fill(D2jet, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::D2jet_M4L118130().cut_d && M4l <= Variables::D2jet_M4L118130().cut_u)
   {
      histos_1D[Settings::D2jet_M4L118130][fs][cat][proc]->Fill(D2jet, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillD2jetZX( float M4l, float D2jet, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::D2jet][fs][cat]->Fill(D2jet, weight);

   if( M4l >= Variables::D2jet_M4L118130().cut_d && M4l <= Variables::D2jet_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::D2jet_M4L118130][fs][cat]->Fill(D2jet, weight);
   }
}
//====================================================================



//=====
// DWH
//====================================================================================
void Histograms::FillDWH( float M4l, float DWH, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::DWH][fs][cat][proc]->Fill(DWH, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::DWH_M4L118130().cut_d && M4l <= Variables::DWH_M4L118130().cut_u)
   {
      histos_1D[Settings::DWH_M4L118130][fs][cat][proc]->Fill(DWH, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillDWHZX( float M4l, float DWH, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::DWH][fs][cat]->Fill(DWH, weight);

   if( M4l >= Variables::DWH_M4L118130().cut_d && M4l <= Variables::DWH_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::DWH_M4L118130][fs][cat]->Fill(DWH, weight);
   }
}
//====================================================================



//====
// DZH
//====================================================================================
void Histograms::FillDZH( float M4l, float DZH, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::DZH][fs][cat][proc]->Fill(DZH, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::DZH_M4L118130().cut_d && M4l <= Variables::DZH_M4L118130().cut_u)
   {
      histos_1D[Settings::DZH_M4L118130][fs][cat][proc]->Fill(DZH, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillDZHZX( float M4l, float DZH, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::DZH][fs][cat]->Fill(DZH, weight);

   if( M4l >= Variables::DZH_M4L118130().cut_d && M4l <= Variables::DZH_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::DZH_M4L118130][fs][cat]->Fill(DZH, weight);
   }
}
//====================================================================



//====
// DVH
//====================================================================================
void Histograms::FillDVH( float M4l, float DVH, float weight, int fs, int cat, int proc )
{

   histos_1D[Settings::DVH][fs][cat][proc]->Fill(DVH, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::DVH_M4L118130().cut_d && M4l <= Variables::DVH_M4L118130().cut_u)
   {
      histos_1D[Settings::DVH_M4L118130][fs][cat][proc]->Fill(DVH, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================

//====================================================================
void Histograms::FillDVHZX( float M4l, float DVH, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::DVH][fs][cat]->Fill(DVH, weight);

   if( M4l >= Variables::DVH_M4L118130().cut_d && M4l <= Variables::DVH_M4L118130().cut_u)
   {
      histos_1D_ZX[Settings::DVH_M4L118130][fs][cat]->Fill(DVH, weight);
   }
}
//====================================================================



//===========
// MZ1 vs MZ2
//====================================================================================
void Histograms::FillMZ1vsMZ2( float M4l, float MZ1, float MZ2, float weight, int fs, int cat, int proc )
{
   histos_2D[Settings::MZ1vsMZ2][fs][cat][proc]->Fill(MZ1, MZ2, (proc == Settings::Data) ? 1. : weight);

   if( M4l >= Variables::MZ1vsMZ2_M4L118130().cut_d && M4l <= Variables::MZ1vsMZ2_M4L118130().cut_u)
   {
      histos_2D[Settings::MZ1vsMZ2_M4L118130][fs][cat][proc]->Fill(MZ1, MZ2, (proc == Settings::Data) ? 1. : weight);
   }
}
//====================================================================================



//====================================================================================
void Histograms::FillVectors( float M4l, float ZZMassErrCorr, float KD, float DVBFDEC, float DVHDEC,int nCleanedJetsPt30, float D1jet, float D2jet, float DWH, float DZH, float DVH, int fs, int cat)
{
  if ( cat == Settings::VBF_2j_tagged )
  {
    vector_X[Settings::DVBFDECvsM4lZoomed][fs][cat].push_back(M4l);
    vector_Y[Settings::DVBFDECvsM4lZoomed][fs][cat].push_back(DVBFDEC);
    vector_EX[Settings::DVBFDECvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
    vector_EY[Settings::DVBFDECvsM4lZoomed][fs][cat].push_back(0.);
  }
  else if ( cat == Settings::VH_hadron_tagged )
  {
    vector_X[Settings::DVHDECvsM4lZoomed][fs][cat].push_back(M4l);
    vector_Y[Settings::DVHDECvsM4lZoomed][fs][cat].push_back(DVHDEC);
    vector_EX[Settings::DVHDECvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
    vector_EY[Settings::DVHDECvsM4lZoomed][fs][cat].push_back(0.);
  }
  else
  {
    vector_X[Settings::KDvsM4l][fs][cat].push_back(M4l);
    vector_Y[Settings::KDvsM4l][fs][cat].push_back(KD);
    vector_EX[Settings::KDvsM4l][fs][cat].push_back(ZZMassErrCorr);
    vector_EY[Settings::KDvsM4l][fs][cat].push_back(0.);

    vector_X[Settings::KDvsM4lZoomed][fs][cat].push_back(M4l);
    vector_Y[Settings::KDvsM4lZoomed][fs][cat].push_back(KD);
    vector_EX[Settings::KDvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
    vector_EY[Settings::KDvsM4lZoomed][fs][cat].push_back(0.);

    vector_X[Settings::KDvsM4lHighMass][fs][cat].push_back(M4l);
    vector_Y[Settings::KDvsM4lHighMass][fs][cat].push_back(KD);
    vector_EX[Settings::KDvsM4lHighMass][fs][cat].push_back(ZZMassErrCorr);
    vector_EY[Settings::KDvsM4lHighMass][fs][cat].push_back(0.);
  }



   if (nCleanedJetsPt30 == 1)
   {
      vector_X[Settings::D1jetvsM4lZoomed][fs][cat].push_back(M4l);
      vector_Y[Settings::D1jetvsM4lZoomed][fs][cat].push_back(D1jet);
      vector_EX[Settings::D1jetvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
      vector_EY[Settings::D1jetvsM4lZoomed][fs][cat].push_back(0.);
   }

   if (nCleanedJetsPt30 >= 2)
   {
      vector_X[Settings::D2jetvsM4lZoomed][fs][cat].push_back(M4l);
      vector_Y[Settings::D2jetvsM4lZoomed][fs][cat].push_back(D2jet);
      vector_EX[Settings::D2jetvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
      vector_EY[Settings::D2jetvsM4lZoomed][fs][cat].push_back(0.);

      vector_X[Settings::DWHvsM4lZoomed][fs][cat].push_back(M4l);
      vector_Y[Settings::DWHvsM4lZoomed][fs][cat].push_back(DWH);
      vector_EX[Settings::DWHvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
      vector_EY[Settings::DWHvsM4lZoomed][fs][cat].push_back(0.);

      vector_X[Settings::DZHvsM4lZoomed][fs][cat].push_back(M4l);
      vector_Y[Settings::DZHvsM4lZoomed][fs][cat].push_back(DZH);
      vector_EX[Settings::DZHvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
      vector_EY[Settings::DZHvsM4lZoomed][fs][cat].push_back(0.);

      vector_X[Settings::DVHvsM4lZoomed][fs][cat].push_back(M4l);
      vector_Y[Settings::DVHvsM4lZoomed][fs][cat].push_back(DVH);
      vector_EX[Settings::DVHvsM4lZoomed][fs][cat].push_back(ZZMassErrCorr);
      vector_EY[Settings::DVHvsM4lZoomed][fs][cat].push_back(0.);
   }
}
//====================================================================================



//====================================================================================
void Histograms::FillDvsM4l( float M4l, float KD, float DVBFDEC, float DVHDEC, int nCleanedJetsPt30, float D1jet, float D2jet, float DWH, float DZH, float DVH, float weight, int fs, int cat, int proc )
{
  if ( cat == Settings::VBF_2j_tagged)
  {
    histos_2DError[Settings::DVBFDECvsM4lZoomed][fs][cat][proc]  ->Fill(M4l, DVBFDEC, (proc == Settings::Data) ? 1. : weight);
  }
  else if ( cat  == Settings::VH_hadron_tagged )
  {
    histos_2DError[Settings::DVHDECvsM4lZoomed][fs][cat][proc]  ->Fill(M4l, DVHDEC, (proc == Settings::Data) ? 1. : weight);
  }
  else
  {
    histos_2DError[Settings::KDvsM4l][fs][cat][proc]        ->Fill(M4l, KD, (proc == Settings::Data) ? 1. : weight);
    histos_2DError[Settings::KDvsM4lZoomed][fs][cat][proc]  ->Fill(M4l, KD, (proc == Settings::Data) ? 1. : weight);
    histos_2DError[Settings::KDvsM4lHighMass][fs][cat][proc]->Fill(M4l, KD, (proc == Settings::Data) ? 1. : weight);
  }
   if (nCleanedJetsPt30 == 1) histos_2DError[Settings::D1jetvsM4lZoomed][fs][cat][proc]->Fill(M4l, D1jet, (proc == Settings::Data) ? 1. : weight);
   if (nCleanedJetsPt30 >= 2) histos_2DError[Settings::D2jetvsM4lZoomed][fs][cat][proc]->Fill(M4l, D2jet, (proc == Settings::Data) ? 1. : weight);
   if (nCleanedJetsPt30 >= 2) histos_2DError[Settings::DWHvsM4lZoomed][fs][cat][proc]  ->Fill(M4l, DWH,   (proc == Settings::Data) ? 1. : weight);
   if (nCleanedJetsPt30 >= 2) histos_2DError[Settings::DZHvsM4lZoomed][fs][cat][proc]  ->Fill(M4l, DZH,   (proc == Settings::Data) ? 1. : weight);
   if (nCleanedJetsPt30 >= 2) histos_2DError[Settings::DVHvsM4lZoomed][fs][cat][proc]  ->Fill(M4l, DVH,   (proc == Settings::Data) ? 1. : weight);
}
//====================================================================================


//====================================================================================
void Histograms::FillDvsM4l_ZX( float M4l, float KD, float DVBFDEC, float DVHDEC, int nCleanedJetsPt30, float weight, int fs, int cat )
{
  if ( cat == Settings::VBF_2j_tagged)
  {
    histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][fs][cat]->Fill(M4l, DVBFDEC, weight);
  }
  else if ( cat  == Settings::VH_hadron_tagged )
  {
    histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][fs][cat]->Fill(M4l, DVHDEC, weight);
  }
  else
  {
    //histos_2DError_ZX[Settings::KDvsM4l][fs][cat]      ->Fill(M4l, KD, weight);
    histos_2DError_ZX[Settings::KDvsM4lZoomed][fs][cat]  ->Fill(M4l, KD, weight);
    //histos_2DError_ZX[Settings::KDvsM4lHighMass][fs][cat]->Fill(M4l, KD, weight);
  }
}
//====================================================================================

//====================================================================================
void Histograms::FillYields( float M4l, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::M4lYields][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
}
//====================================================================================


//========
// Others
//====================================================================================
void Histograms::FillOthers( float M4l, float ZZPt, float ZZEta, float PFMET, float Pt_leading, float Pt_trailing, float Eta_leading, float Eta_trailing, float SIP_leading, float SIP_trailing, float ISO_leading, float ISO_trailing, int NExtraLep, int NJets, int NJetsBTagged, float KD, float DVBFDEC, float DVHDEC, float weight, int fs, int cat, int proc )
{
   histos_1D[Settings::PFMET][fs][cat][proc]       ->Fill(PFMET, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::Pt4l][fs][cat][proc]        ->Fill(ZZPt, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::Eta4l][fs][cat][proc]       ->Fill(ZZEta, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::Pt_leading][fs][cat][proc]  ->Fill(Pt_leading, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::Pt_trailing][fs][cat][proc] ->Fill(Pt_trailing, (proc == Settings::Data) ? 1. : weight);
  histos_1D[Settings::Eta_leading][fs][cat][proc]  ->Fill(Eta_leading, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::Eta_trailing][fs][cat][proc] ->Fill(Eta_trailing, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::SIP_leading][fs][cat][proc] ->Fill(SIP_leading, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::SIP_trailing][fs][cat][proc]->Fill(SIP_trailing, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::ISO_leading][fs][cat][proc] ->Fill(ISO_leading, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::ISO_trailing][fs][cat][proc]->Fill(ISO_trailing, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::NExtraLep][fs][cat][proc]   ->Fill(NExtraLep, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::NJets][fs][cat][proc]       ->Fill(NJets, (proc == Settings::Data) ? 1. : weight);
   histos_1D[Settings::NJetsBTagged][fs][cat][proc]->Fill(NJetsBTagged, (proc == Settings::Data) ? 1. : weight);
   if (DVBFDEC > 0.5 && cat == Settings::VBF_2j_tagged) histos_1D[Settings::M4l_110150_HighKD][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
   else if (DVHDEC > 0.5 && cat == Settings::VH_hadron_tagged) histos_1D[Settings::M4l_110150_HighKD][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
   else if (KD > 0.5 && (cat != Settings::VBF_2j_tagged && cat != Settings::VH_hadron_tagged)) histos_1D[Settings::M4l_110150_HighKD][fs][cat][proc]->Fill(M4l, (proc == Settings::Data) ? 1. : weight);
}
//====================================================================================

//====================================================================
void Histograms::FillOthersZX( float M4l, float ZZPt, float ZZEta, float PFMET, float Pt_leading, float Pt_trailing, float Eta_leading, float Eta_trailing, float SIP_leading, float SIP_trailing, float ISO_leading, float ISO_trailing, int NExtraLep, int NJets, int NJetsBTagged, float KD, float DVBFDEC, float DVHDEC, float weight, int fs, int cat )
{
   histos_1D_ZX[Settings::PFMET][fs][cat]       ->Fill(PFMET, weight);
  histos_1D_ZX[Settings::Pt4l][fs][cat]        ->Fill(ZZPt, weight);
   histos_1D_ZX[Settings::Eta4l][fs][cat]       ->Fill(ZZEta, weight);
   histos_1D_ZX[Settings::Pt_leading][fs][cat]  ->Fill(Pt_leading, weight);
   histos_1D_ZX[Settings::Pt_trailing][fs][cat] ->Fill(Pt_trailing, weight);
   histos_1D_ZX[Settings::Eta_leading][fs][cat]  ->Fill(Eta_leading, weight);
   histos_1D_ZX[Settings::Eta_trailing][fs][cat] ->Fill(Eta_trailing, weight);
   histos_1D_ZX[Settings::SIP_leading][fs][cat] ->Fill(SIP_leading, weight);
   histos_1D_ZX[Settings::SIP_trailing][fs][cat]->Fill(SIP_trailing, weight);
   histos_1D_ZX[Settings::ISO_leading][fs][cat] ->Fill(ISO_leading, weight);
   histos_1D_ZX[Settings::ISO_trailing][fs][cat]->Fill(ISO_trailing, weight);
   histos_1D_ZX[Settings::NExtraLep][fs][cat]   ->Fill(NExtraLep, weight);
   histos_1D_ZX[Settings::NJets][fs][cat]       ->Fill(NJets, weight);
   histos_1D_ZX[Settings::NJetsBTagged][fs][cat]->Fill(NJetsBTagged, weight);
   if (DVBFDEC > 0.5 && cat == Settings::VBF_2j_tagged) histos_1D_ZX[Settings::M4l_110150_HighKD][fs][cat]->Fill(M4l, weight);
   else if (DVHDEC > 0.5 && cat == Settings::VH_hadron_tagged) histos_1D_ZX[Settings::M4l_110150_HighKD][fs][cat]->Fill(M4l, weight);
   else if (KD > 0.5 && (cat != Settings::VBF_2j_tagged && cat != Settings::VH_hadron_tagged)) histos_1D_ZX[Settings::M4l_110150_HighKD][fs][cat]->Fill(M4l, weight);

}
//====================================================================


//=======================================================================================
void Histograms::MakeZXShape( int current_category , vector< vector <float> > _expected_yield_SR)
{

   M4lZX *ZXShape = new M4lZX();

   ZXShape->GetM4lZX(Variables::M4lMain().var_N_bin, Variables::M4lMain().var_min, Variables::M4lMain().var_max, current_category, _expected_yield_SR,
                     histos_1D_ZX_shape[Settings::M4lMain][Settings::fs4e][current_category],
                     histos_1D_ZX_shape[Settings::M4lMain][Settings::fs4mu][current_category],
                     histos_1D_ZX_shape[Settings::M4lMain][Settings::fs2e2mu][current_category],
                     histos_1D_ZX_shape[Settings::M4lMain][Settings::fs4l][current_category]);


   M4lZX *ZXShape_zoomed = new M4lZX();

   ZXShape_zoomed->GetM4lZX(Variables::M4lMainZoomed().var_N_bin, Variables::M4lMainZoomed().var_min, Variables::M4lMainZoomed().var_max, current_category, _expected_yield_SR,
                            histos_1D_ZX_shape[Settings::M4lMainZoomed][Settings::fs4e][current_category],
                            histos_1D_ZX_shape[Settings::M4lMainZoomed][Settings::fs4mu][current_category],
                            histos_1D_ZX_shape[Settings::M4lMainZoomed][Settings::fs2e2mu][current_category],
                            histos_1D_ZX_shape[Settings::M4lMainZoomed][Settings::fs4l][current_category]);


   M4lZX *ZXShape_HighMass = new M4lZX();

   ZXShape_HighMass->GetM4lZX(Variables::M4lMainHighMass().var_N_bin, Variables::M4lMainHighMass().var_min, Variables::M4lMainHighMass().var_max, current_category, _expected_yield_SR,
                              histos_1D_ZX_shape[Settings::M4lMainHighMass][Settings::fs4e][current_category],
                              histos_1D_ZX_shape[Settings::M4lMainHighMass][Settings::fs4mu][current_category],
                              histos_1D_ZX_shape[Settings::M4lMainHighMass][Settings::fs2e2mu][current_category],
                              histos_1D_ZX_shape[Settings::M4lMainHighMass][Settings::fs4l][current_category]);

}
//=======================================================================================

//==============================
void Histograms::FillInclusive()
{
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ ) // Change to num_of_final_states - 1
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ ) // Change to num_of_categories - 1
      {

        for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125ggH]);
          histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125VBF]);
          histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125VH]);
          histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125ttH]);
          histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125bbH]);
          histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][Settings::H125tqH]);

      }

      for ( int i_2Dhistos = 0; i_2Dhistos < num_of_2D_plot_names; i_2Dhistos++ )

        {
          histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125ggH]);
          histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125VBF]);
          histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125VH]);
          histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125ttH]);
          histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125bbH]);
          histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][Settings::H125tqH]);

      }

      for ( int i_2D_errorhistos = 0; i_2D_errorhistos < num_of_2D_error_plot_names; i_2D_errorhistos++ )

        {
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125ggH]);
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125VBF]);
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125VH]);
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125ttH]);
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125bbH]);
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][Settings::H125tqH]);

      }

      }
   }

   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
         {

        for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D[i_1Dhistos][Settings::fs4l][i_cat][i_proc]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][i_proc]);
        }

        for ( int i_2Dhistos = 0; i_2Dhistos < num_of_2D_plot_names; i_2Dhistos++ )

        {
          histos_2D[i_2Dhistos][Settings::fs4l][i_cat][i_proc]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][i_proc]);
        }

        for ( int i_2D_errorhistos = 0; i_2D_errorhistos < num_of_2D_error_plot_names; i_2D_errorhistos++ )

        {
          histos_2DError[i_2D_errorhistos][Settings::fs4l][i_cat][i_proc]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][i_proc]);
        }

         }
      }
   }

  for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {

      histos_2DError_ZX[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat]->Add(histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][i_cat]);
      histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][Settings::fs4l][i_cat]->Add(histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat]);
      histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][Settings::fs4l][i_cat]->Add(histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat]);

      }
   }


   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
         {
        for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D[i_1Dhistos][i_fs][Settings::inclusive][i_proc]->Add(histos_1D[i_1Dhistos][i_fs][i_cat][i_proc]);
        }

        for ( int i_2Dhistos = 0; i_2Dhistos < num_of_2D_plot_names; i_2Dhistos++ )

        {
          histos_2D[i_2Dhistos][i_fs][Settings::inclusive][i_proc]->Add(histos_2D[i_2Dhistos][i_fs][i_cat][i_proc]);
        }

        for ( int i_2D_errorhistos = 0; i_2D_errorhistos < num_of_2D_error_plot_names; i_2D_errorhistos++ )

        {
          histos_2DError[i_2D_errorhistos][i_fs][Settings::inclusive][i_proc]->Add(histos_2DError[i_2D_errorhistos][i_fs][i_cat][i_proc]);
        }

         }
      }
   }

  for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {

      histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive]->Add(histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][i_cat]);
      histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][Settings::inclusive]->Add(histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat]);
      histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][Settings::inclusive]->Add(histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat]);

      }
   }


   //====================
   // Sum Z+X histograms
   //====================
   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {

      for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D_ZX[i_1Dhistos][Settings::fs4l][i_cat]->Add(histos_1D_ZX[i_1Dhistos][i_fs][i_cat]);
        }

         //==========================
         // 2D plots with mass error
         //==========================
         vector_X[Settings::KDvsM4l][Settings::fs4l][i_cat].insert(vector_X[Settings::KDvsM4l][Settings::fs4l][i_cat].end(),vector_X[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_X[Settings::KDvsM4l][i_fs][i_cat].end());
         vector_Y[Settings::KDvsM4l][Settings::fs4l][i_cat].insert(vector_Y[Settings::KDvsM4l][Settings::fs4l][i_cat].end(),vector_Y[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_Y[Settings::KDvsM4l][i_fs][i_cat].end());
         vector_EX[Settings::KDvsM4l][Settings::fs4l][i_cat].insert(vector_EX[Settings::KDvsM4l][Settings::fs4l][i_cat].end(),vector_EX[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_EX[Settings::KDvsM4l][i_fs][i_cat].end());
         vector_EY[Settings::KDvsM4l][Settings::fs4l][i_cat].insert(vector_EY[Settings::KDvsM4l][Settings::fs4l][i_cat].end(),vector_EY[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_EY[Settings::KDvsM4l][i_fs][i_cat].end());

         vector_X[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_X[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_X[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::KDvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_Y[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_Y[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::KDvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EX[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EX[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::KDvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EY[Settings::KDvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EY[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::KDvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].insert(vector_X[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].end(),vector_X[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_X[Settings::KDvsM4lHighMass][i_fs][i_cat].end());
         vector_Y[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].insert(vector_Y[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].end(),vector_Y[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_Y[Settings::KDvsM4lHighMass][i_fs][i_cat].end());
         vector_EX[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].insert(vector_EX[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].end(),vector_EX[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_EX[Settings::KDvsM4lHighMass][i_fs][i_cat].end());
         vector_EY[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].insert(vector_EY[Settings::KDvsM4lHighMass][Settings::fs4l][i_cat].end(),vector_EY[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_EY[Settings::KDvsM4lHighMass][i_fs][i_cat].end());

         vector_X[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_X[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_X[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_Y[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_Y[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EX[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EX[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EY[Settings::D1jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EY[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_X[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_X[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_Y[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_Y[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EX[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EX[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EY[Settings::D2jetvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EY[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_X[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_X[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_Y[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_Y[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EX[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EX[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EY[Settings::DWHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EY[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_X[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_Y[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_Y[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EX[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EX[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EY[Settings::DZHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EY[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_X[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_X[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_Y[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_Y[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EX[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EX[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].insert(vector_EY[Settings::DVHvsM4lZoomed][Settings::fs4l][i_cat].end(),vector_EY[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());
      }
   }

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
      for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D_ZX[i_1Dhistos][i_fs][Settings::inclusive]->Add(histos_1D_ZX[i_1Dhistos][i_fs][i_cat]);
        }

         //==========================
         // 2D plots with mass error
         //==========================

         vector_X[Settings::KDvsM4l][i_fs][Settings::inclusive].insert(vector_X[Settings::KDvsM4l][i_fs][Settings::inclusive].end(),vector_X[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_X[Settings::KDvsM4l][i_fs][i_cat].end());
         vector_Y[Settings::KDvsM4l][i_fs][Settings::inclusive].insert(vector_Y[Settings::KDvsM4l][i_fs][Settings::inclusive].end(),vector_Y[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_Y[Settings::KDvsM4l][i_fs][i_cat].end());
         vector_EX[Settings::KDvsM4l][i_fs][Settings::inclusive].insert(vector_EX[Settings::KDvsM4l][i_fs][Settings::inclusive].end(),vector_EX[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_EX[Settings::KDvsM4l][i_fs][i_cat].end());
         vector_EY[Settings::KDvsM4l][i_fs][Settings::inclusive].insert(vector_EY[Settings::KDvsM4l][i_fs][Settings::inclusive].end(),vector_EY[Settings::KDvsM4l][i_fs][i_cat].begin(),vector_EY[Settings::KDvsM4l][i_fs][i_cat].end());

         vector_X[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_X[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_X[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::KDvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_Y[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_Y[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::KDvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EX[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_EX[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::KDvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EY[Settings::KDvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_EY[Settings::KDvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::KDvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].insert(vector_X[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].end(),vector_X[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_X[Settings::KDvsM4lHighMass][i_fs][i_cat].end());
         vector_Y[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].insert(vector_Y[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].end(),vector_Y[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_Y[Settings::KDvsM4lHighMass][i_fs][i_cat].end());
         vector_EX[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].insert(vector_EX[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].end(),vector_EX[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_EX[Settings::KDvsM4lHighMass][i_fs][i_cat].end());
         vector_EY[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].insert(vector_EY[Settings::KDvsM4lHighMass][i_fs][Settings::inclusive].end(),vector_EY[Settings::KDvsM4lHighMass][i_fs][i_cat].begin(),vector_EY[Settings::KDvsM4lHighMass][i_fs][i_cat].end());

         vector_X[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_X[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_X[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_Y[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_Y[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EX[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_EX[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EY[Settings::D1jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_EY[Settings::D1jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::D1jetvsM4lZoomed][i_fs][i_cat].end());

         vector_X[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_X[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_X[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_X[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());
         vector_Y[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_Y[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_Y[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_Y[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EX[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EX[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_EX[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EX[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());
         vector_EY[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EY[Settings::D2jetvsM4lZoomed][i_fs][Settings::inclusive].end(),vector_EY[Settings::D2jetvsM4lZoomed][i_fs][i_cat].begin(),vector_EY[Settings::D2jetvsM4lZoomed][i_fs][i_cat].end());

         // DWHvsM4lZoomed
         vector_X[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_X[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                              vector_X[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                              vector_X[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());

         vector_Y[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_Y[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                              vector_Y[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                              vector_Y[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());

         vector_EX[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EX[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                               vector_EX[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                               vector_EX[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());

         vector_EY[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EY[Settings::DWHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                               vector_EY[Settings::DWHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                               vector_EY[Settings::DWHvsM4lZoomed][i_fs][i_cat].end());

         // DZHvsM4lZoomed
         vector_X[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_X[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                              vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                              vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());

         vector_Y[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_Y[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                              vector_Y[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                              vector_Y[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());

         vector_EX[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EX[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                               vector_EX[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                               vector_EX[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());

         vector_EY[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EY[Settings::DZHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                               vector_EY[Settings::DZHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                               vector_EY[Settings::DZHvsM4lZoomed][i_fs][i_cat].end());

         // DVHvsM4lZoomed
         vector_X[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_X[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                              vector_X[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                              vector_X[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());

         vector_Y[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_Y[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                              vector_Y[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                              vector_Y[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());

         vector_EX[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EX[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                               vector_EX[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                               vector_EX[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());

         vector_EY[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].insert(vector_EY[Settings::DVHvsM4lZoomed][i_fs][Settings::inclusive].end(),
                                                                               vector_EY[Settings::DVHvsM4lZoomed][i_fs][i_cat].begin(),
                                                                               vector_EY[Settings::DVHvsM4lZoomed][i_fs][i_cat].end());
      }
   }

    for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ )
    {

        STXS_Categories[Settings::Zjets]->SetBinContent(i_cat + 1,STXS_Yields_ZX[i_cat]->Integral());

        for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
        {
            if(i_proc == Settings::Zjets) continue;
            STXS_Categories[i_proc]->SetBinContent(i_cat + 1,STXS_Yields[i_cat][i_proc]->Integral());
        }
    }

    for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ )// first calculate total yield in each bin so you can normalize each bin to 1
    {
        for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
        {
            Purity_Normalization[i_cat] += Purity_Yields[i_cat][i_proc]->Integral();
        }
    }

    for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ )
    {
        for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
        {
            Purity_Categories[i_proc]->SetBinContent(num_of_STXS_categories - i_cat,Purity_Yields[i_cat][i_proc]->Integral()/Purity_Normalization[i_cat]);
        }
    }
}
//==============================


//==============================
void Histograms::FillInclusiveYields()
{
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120ggH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120VBF]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120WHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120WHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120ZHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120ZHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120ttHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120ttHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120bbH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH120tqH]);

         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124ggH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124VBF]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124WHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124WHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124ZHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124ZHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124ttHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124ttHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124bbH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH124tqH]);

         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125ggH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125VBF]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125WHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125WHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125ZHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125ZHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125ttHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125ttHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125bbH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125tqH]);

         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126ggH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126VBF]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126WHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126WHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126ZHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126ZHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126ttHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126ttHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126bbH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH126tqH]);

         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130ggH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130VBF]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130WHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130WHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130ZHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130ZHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130ttHhad]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130ttHlep]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130bbH]);
         histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH130tqH]);
      }
   }

   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
         {
            histos_1D[Settings::M4lYields][Settings::fs4l][i_cat][i_proc]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]);
         }
      }
   }

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
         {
            histos_1D[Settings::M4lYields][i_fs][Settings::inclusive][i_proc]->Add(histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]);
         }
      }
   }
}
//==============================



//=================================
void Histograms::SmoothHistograms()
{
   float integral = 0;
   TH1F* CheckSmoothing;

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         //=============
         // M4l
         //=============
         integral = histos_1D_ZX[Settings::M4lMain][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::M4lMain][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::M4lMain][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::M4lMain][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::M4lMain][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat]->Integral() );
         }


         //=============
         // MZ1
         //=============
         integral = histos_1D_ZX[Settings::MZ1][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::MZ1][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::MZ1][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::MZ1][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::MZ1][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat]->Integral() );
         }

         //=============
         // MZ2
         //=============
         integral = histos_1D_ZX[Settings::MZ2][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::MZ2][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {

            histos_1D_ZX[Settings::MZ2][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::MZ2][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::MZ2][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat]->Integral() );
         }

         //=============
         // KD
         //=============
         integral = histos_1D_ZX[Settings::KD][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::KD][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {

            histos_1D_ZX[Settings::KD][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::KD][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::KD][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat]->Integral() );
         }

         //=============
         // D1jet
         //=============
         integral = histos_1D_ZX[Settings::D1jet][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::D1jet][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {

            histos_1D_ZX[Settings::D1jet][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::D1jet][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::D1jet][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat]->Integral() );
         }

         //=============
         // D2jet
         //=============
         integral = histos_1D_ZX[Settings::D2jet][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::D2jet][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {

            histos_1D_ZX[Settings::D2jet][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::D2jet][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::D2jet][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat]->Integral() );
         }

         //=============
         // DWH
         //=============
         integral = histos_1D_ZX[Settings::DWH][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::DWH][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {

            histos_1D_ZX[Settings::DWH][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::DWH][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::DWH][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat]->Integral() );
         }

         //=============
         // DZH
         //=============
         integral = histos_1D_ZX[Settings::DZH][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::DZH][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {

            histos_1D_ZX[Settings::DZH][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::DZH][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::DZH][i_fs][i_cat]->Integral() );
         }

         integral = histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat]->Integral();
         CheckSmoothing = histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat];
         CheckSmoothing->Smooth(1);
         if(integral > 0. && CheckSmoothing->Integral() > 0.)
         {
            histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat]->Smooth(1);
            histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat]->Scale( integral / histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat]->Integral() );
         }
      }
   }
}
//=================================



//==============================
void Histograms::RenormalizeZX( int year )
{
   M4lZX *ZX = new M4lZX();

   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

      {
         if ( i_1Dhistos == Settings::M4lYields ) continue;

         ZX->RenormalizeZX(i_cat, year,
                           histos_1D_ZX[i_1Dhistos][Settings::fs4e][i_cat],
                           histos_1D_ZX[i_1Dhistos][Settings::fs4mu][i_cat],
                           histos_1D_ZX[i_1Dhistos][Settings::fs2e2mu][i_cat]);
      }

   }
   ZX->~M4lZX();
}
//==============================



//=============================================
void Histograms::SaveHistos( string file_name )
{
   TFile *fOutHistos = new TFile(file_name.c_str(), "recreate");
   fOutHistos->cd();

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {

      for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D_ZX[i_1Dhistos][i_fs][i_cat]->Write();
        }

      histos_1D_ZX_shape[Settings::M4lMain][i_fs][i_cat]->Write();
      histos_1D_ZX_shape[Settings::M4lMainZoomed][i_fs][i_cat]->Write();
      histos_1D_ZX_shape[Settings::M4lMainHighMass][i_fs][i_cat]->Write();

         //==========================
         // 2D plots with mass error
         //==========================
         histos_2DError_data[Settings::KDvsM4l][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::KDvsM4l][i_fs][i_cat].size(),
                                                                                 &(vector_X[Settings::KDvsM4l][i_fs][i_cat][0]),
                                                                                 &(vector_Y[Settings::KDvsM4l][i_fs][i_cat][0]),
                                                                                 &(vector_EX[Settings::KDvsM4l][i_fs][i_cat][0]),
                                                                                 &(vector_EY[Settings::KDvsM4l][i_fs][i_cat][0]) );
         _histo_name = "KDvsM4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::KDvsM4l][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::KDvsM4l][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::KDvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::KDvsM4lZoomed][i_fs][i_cat].size(),
                                                                                       &(vector_X[Settings::KDvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_Y[Settings::KDvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_EX[Settings::KDvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_EY[Settings::KDvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "KDvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::KDvsM4lZoomed][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::KDvsM4lZoomed][i_fs][i_cat]->Write();


      histos_2DError_data[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat].size(),
                                                                                       &(vector_X[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_Y[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_EX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_EY[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "DVBFDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat]->Write();


          histos_2DError_data[Settings::DVHDECvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::DVHDECvsM4lZoomed][i_fs][i_cat].size(),
                                                                                       &(vector_X[Settings::DVHDECvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_Y[Settings::DVHDECvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_EX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat][0]),
                                                                                       &(vector_EY[Settings::DVHDECvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "DVHDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DVHDECvsM4lZoomed][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::DVHDECvsM4lZoomed][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::KDvsM4lHighMass][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::KDvsM4lHighMass][i_fs][i_cat].size(),
                                                                                         &(vector_X[Settings::KDvsM4lHighMass][i_fs][i_cat][0]),
                                                                                         &(vector_Y[Settings::KDvsM4lHighMass][i_fs][i_cat][0]),
                                                                                         &(vector_EX[Settings::KDvsM4lHighMass][i_fs][i_cat][0]),
                                                                                         &(vector_EY[Settings::KDvsM4lHighMass][i_fs][i_cat][0]) );
         _histo_name = "KDvsM4lHighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::KDvsM4lHighMass][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::KDvsM4lHighMass][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::D1jetvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::D1jetvsM4lZoomed][i_fs][i_cat].size(),
                                                                                          &(vector_X[Settings::D1jetvsM4lZoomed][i_fs][i_cat][0]),
                                                                                          &(vector_Y[Settings::D1jetvsM4lZoomed][i_fs][i_cat][0]),
                                                                                          &(vector_EX[Settings::D1jetvsM4lZoomed][i_fs][i_cat][0]),
                                                                                          &(vector_EY[Settings::D1jetvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "D1jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::D1jetvsM4lZoomed][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::D1jetvsM4lZoomed][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::D2jetvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::D2jetvsM4lZoomed][i_fs][i_cat].size(),
                                                                                          &(vector_X[Settings::D2jetvsM4lZoomed][i_fs][i_cat][0]),
                                                                                          &(vector_Y[Settings::D2jetvsM4lZoomed][i_fs][i_cat][0]),
                                                                                          &(vector_EX[Settings::D2jetvsM4lZoomed][i_fs][i_cat][0]),
                                                                                          &(vector_EY[Settings::D2jetvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "D2jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::D2jetvsM4lZoomed][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::D2jetvsM4lZoomed][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::DWHvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::DWHvsM4lZoomed][i_fs][i_cat].size(),
                                                                                        &(vector_X[Settings::DWHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_Y[Settings::DWHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_EX[Settings::DWHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_EY[Settings::DWHvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "DWHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DWHvsM4lZoomed][i_fs][i_cat]->SetName(_histo_name.c_str());
         histos_2DError_data[Settings::DWHvsM4lZoomed][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::DZHvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat].size(),
                                                                                        &(vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_Y[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_EX[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_EY[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "DZHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DZHvsM4lZoomed][i_fs][i_cat]->SetName( _histo_name.c_str() );
         histos_2DError_data[Settings::DZHvsM4lZoomed][i_fs][i_cat]->Write();


         histos_2DError_data[Settings::DVHvsM4lZoomed][i_fs][i_cat] = new TGraphErrors( vector_X[Settings::DVHvsM4lZoomed][i_fs][i_cat].size(),
                                                                                        &(vector_X[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_Y[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_EX[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]),
                                                                                        &(vector_EY[Settings::DZHvsM4lZoomed][i_fs][i_cat][0]) );
         _histo_name = "DVHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DVHvsM4lZoomed][i_fs][i_cat]->SetName( _histo_name.c_str() );
         histos_2DError_data[Settings::DVHvsM4lZoomed][i_fs][i_cat]->Write();


         histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][i_cat]->Write();
         histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat]->Write();
         histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat]->Write();

         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
         {
            for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          histos_1D[i_1Dhistos][i_fs][i_cat][i_proc]->Write();
        }

        for ( int i_2Dhistos = 0; i_2Dhistos < num_of_2D_plot_names; i_2Dhistos++ )

        {
          histos_2D[i_2Dhistos][i_fs][i_cat][i_proc]->Write();
        }

        for ( int i_2D_errorhistos = 0; i_2D_errorhistos < num_of_2D_error_plot_names; i_2D_errorhistos++ )

        {
          histos_2DError[i_2D_errorhistos][i_fs][i_cat][i_proc]->Write();
        }
         }
      }
   }

    for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
    {
        STXS_Categories[i_proc]->Write();
    }

    for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
    {
        Purity_Categories[i_proc]->Write();
    }

   fOutHistos->Close();
   delete fOutHistos;
}
//=============================================



//=============================================
void Histograms::SaveYieldHistos( string file_name )
{

   cout << "[INFO] Saving yield histograms to ROOT file." << endl;

   TFile* fOutHistos = new TFile(file_name.c_str(), "recreate");
   fOutHistos->cd();

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
         {
            histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]->Write();
         }
      }
   }

   cout << "[INFO] Closing ROOT file." << endl;

   fOutHistos->Close();
   delete fOutHistos;
}
//=============================================



//=============================
void Histograms::DeleteHistos()
{
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
       for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

        {
          if ( i_1Dhistos == Settings::M4lYields ) continue;
          delete histos_1D_ZX[i_1Dhistos][i_fs][i_cat];
        }


         //==========================
         // 2D plots with mass error
         //==========================
         delete histos_2DError_data[Settings::KDvsM4l][i_fs][i_cat];
         delete histos_2DError_data[Settings::KDvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::DVHDECvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::KDvsM4lHighMass][i_fs][i_cat];
         delete histos_2DError_data[Settings::D1jetvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::D2jetvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::DWHvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::DZHvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_data[Settings::DVHvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat];
         delete histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat];

         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
         {

        for ( int i_1Dhistos = 0; i_1Dhistos < num_of_1D_plot_names; i_1Dhistos++ )

          {
            if ( i_1Dhistos == Settings::M4lYields ) continue;
            delete histos_1D[i_1Dhistos][i_fs][i_cat][i_proc];
          }

          for ( int i_2Dhistos = 0; i_2Dhistos < num_of_2D_plot_names; i_2Dhistos++ )

          {
            delete histos_2D[i_2Dhistos][i_fs][i_cat][i_proc];
          }

          for ( int i_2D_errorhistos = 0; i_2D_errorhistos < num_of_2D_error_plot_names; i_2D_errorhistos++ )

          {
            delete histos_2DError[i_2D_errorhistos][i_fs][i_cat][i_proc];
          }

         }
      }
   }

}
//=============================



//===================================
void Histograms::DeleteYieldsHistos()
{
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
         {
            delete histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc];
         }
      }
   }
}
//===================================



//=============================================
void Histograms::GetHistos( TString file_name )
{
   TFile* histo_file = TFile::Open(file_name);

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
         {
            //=============
            // M4l
            //=============
            _histo_name = "M4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::M4lMain][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "M4l_zoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::M4lMainZoomed][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "M4l_HighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::M4lMainHighMass][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            //=============
            // MZ1
            //=============
            _histo_name = "MZ1" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::MZ1][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "MZ1_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::MZ1_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            //=============
            // MZ2
            //=============
            _histo_name = "MZ2" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::MZ2][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "MZ2_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::MZ2_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            //=============
            // KD
            //=============
            _histo_name = "KD" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::KD][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "KD_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::KD_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVBFDEC" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DVBFDEC][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVBFDEC_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DVBFDEC_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVHDEC" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DVHDEC][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVHDEC_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DVHDEC_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "D1jet" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::D1jet][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "D1jet_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::D1jet_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "D2jet" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::D2jet][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "D2jet_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::D2jet_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DWH" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DWH][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DWH_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DWH_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DZH" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DZH][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DZH_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DZH_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVH" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DVH][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVH_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::DVH_M4L118130][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            //=============
            // MZ1vsMZ2
            //=============
            _histo_name = "MZ1vsMZ2" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2D[Settings::MZ1vsMZ2][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "MZ1vsMZ2_M4L118130" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2D[Settings::MZ1vsMZ2_M4L118130][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            //=============
            // KDvsM4l
            //=============
            _histo_name = "KDvsM4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::KDvsM4l][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "KDvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::KDvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "KDvsM4lHighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::KDvsM4lHighMass][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVBFDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVHDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::DVHDECvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "D1jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::D1jetvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "D2jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::D2jetvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DWHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::DWHvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DZHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::DZHvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "DVHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_2DError[Settings::DVHvsM4lZoomed][i_fs][i_cat][i_proc] = (TH2F*)histo_file->Get(_histo_name.c_str());

            //===========
            // Others
            //===========
            _histo_name = "PFMET" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::PFMET][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "Pt_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::Pt_leading][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "Pt_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::Pt_trailing][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "Eta_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::Eta_leading][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "Eta_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::Eta_trailing][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "SIP_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::SIP_leading][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "SIP_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::SIP_trailing][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "ISO_leading" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::ISO_leading][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "ISO_trailing" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::ISO_trailing][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "Pt4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::Pt4l][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "NJetsBTagged" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::NJetsBTagged][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "Eta4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::Eta4l][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "NExtraLep" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::NExtraLep][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "NJets" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::NJets][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());

            _histo_name = "M4l_110150_HighKD" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + _blinding;
            histos_1D[Settings::M4l_110150_HighKD][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
         }
      }
   }

   // Z+X
   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         //=============
         // M4l
         //=============
         _histo_name = "M4l_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::M4lMain][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "M4l_ZX_shape_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX_shape[Settings::M4lMain][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "M4l_zoomed_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::M4lMainZoomed][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "M4l_zoomed_ZX_shape_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX_shape[Settings::M4lMainZoomed][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "M4l_HighMass_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::M4lMainHighMass][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "M4l_HighMass_ZX_shape_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX_shape[Settings::M4lMainHighMass][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());


         //=============
         // MZ1
         //=============
         _histo_name = "MZ1_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ1][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "MZ1_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ1_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         //=============
         // MZ2
         //=============
         _histo_name = "MZ2_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ2][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "MZ2_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::MZ2_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         //=============
         // KD
         //=============
         _histo_name = "KD_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::KD][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "KD_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::KD_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVBFDEC_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVBFDEC][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "DVBFDEC_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVBFDEC_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "DVHDEC_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVHDEC][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "DVHDEC_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::DVHDEC_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "D1jet_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D1jet][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "D1jet_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D1jet_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "D2jet_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D2jet][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "D2jet_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::D2jet_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DWH_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DWH][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DWH_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DWH_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DZH_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DZH][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DZH_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DZH_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVH_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DVH][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVH_M4L118130_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_1D_ZX[Settings::DVH_M4L118130][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

         //==========================
         // 2D plots with mass error
         //==========================
         _histo_name = "KDvsM4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::KDvsM4l][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "KDvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::KDvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "KDvsM4lHighMass" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::KDvsM4lHighMass][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVBFDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVHDECvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DVHDECvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "D1jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::D1jetvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "D2jetvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::D2jetvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DWHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DWHvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DZHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DZHvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVHvsM4lZoomed" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
         histos_2DError_data[Settings::DVHvsM4lZoomed][i_fs][i_cat] = (TGraphErrors*)histo_file->Get(_histo_name.c_str());

         _histo_name = "KDvsM4lZoomed_ZX" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _blinding;
         histos_2DError_ZX[Settings::KDvsM4lZoomed][i_fs][i_cat] = (TH2F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVBFDECvsM4lZoomed_ZX" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _blinding;
         histos_2DError_ZX[Settings::DVBFDECvsM4lZoomed][i_fs][i_cat] = (TH2F*)histo_file->Get(_histo_name.c_str());

         _histo_name = "DVHDECvsM4lZoomed_ZX" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _blinding;
         histos_2DError_ZX[Settings::DVHDECvsM4lZoomed][i_fs][i_cat] = (TH2F*)histo_file->Get(_histo_name.c_str());

      //===========
      // Others
      //===========
      _histo_name = "PFMET_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::PFMET][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "Pt_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Pt_leading][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "Pt_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Pt_trailing][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "Eta_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Eta_leading][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "Eta_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Eta_trailing][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "SIP_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::SIP_leading][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "SIP_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::SIP_trailing][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "ISO_leading_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::ISO_leading][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "ISO_trailing_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::ISO_trailing][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "Pt4l_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Pt4l][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "NJetsBTagged_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::NJetsBTagged][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "Eta4l_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::Eta4l][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "NExtraLep_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::NExtraLep][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "NJets_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::NJets][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());

      _histo_name = "M4l_110150_HighKD_ZX_SS_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + _blinding;
      histos_1D_ZX[Settings::M4l_110150_HighKD][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name.c_str());
      }
   }

    for ( int i_proc = 0; i_proc < num_of_processes; i_proc++ )
    {
            _histo_name = "STXS_Categories_" + _s_process.at(i_proc) + "_" + _blinding;
            STXS_Categories[i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
    }

    for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
    {
        _histo_name = "STXS_Purity_" + _s_STXS_bins_histoName.at(i_proc) + "_" + _blinding;
        Purity_Categories[i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
    }
}
//=============================================

//===================================================
void Histograms::GetYieldsHistos( TString file_name )
{
   TFile* histo_file = TFile::Open(file_name);

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
         {
            _histo_name = "M4l" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat) + "_" + _s_process.at(i_proc) + "_" + _blinding;
            histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc] = (TH1F*)histo_file->Get(_histo_name.c_str());
         }
      }
   }

}
//===================================================



//=========================================================================================================
void Histograms::plot_1D_single( TString filename, TString variable_name, TString folder, int fs, int cat )
{
   int plot_index = SetPlotName( variable_name);

   TCanvas *c;
   if(variable_name == "M4lMain") c = new TCanvas(variable_name, variable_name, 650, 500);
   else c = new TCanvas(variable_name, variable_name, 600, 600);

   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();


   is_Djet_ = plot_index == Settings::D1jet_M4L118130 || plot_index == Settings::D2jet_M4L118130 || plot_index == Settings::D1jet || plot_index == Settings::D2jet || plot_index == Settings::DVBFDEC || plot_index == Settings::DVBFDEC_M4L118130;
   is_DVH_  = plot_index == Settings::DWH_M4L118130 || plot_index == Settings::DWH || plot_index == Settings::DZH_M4L118130 || plot_index == Settings::DZH || plot_index == Settings::DVH_M4L118130 || plot_index == Settings::DVH || plot_index == Settings::DVHDEC || plot_index == Settings::DVHDEC_M4L118130;

   if ( is_Djet_ )
   {
      histos_1D[plot_index][fs][cat][Settings::H125VBF]->SetFillColor(Cosmetics::VBF().fill_color);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      histos_1D[plot_index][fs][cat][Settings::H125VH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      histos_1D[plot_index][fs][cat][Settings::H125ttH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
   }
   else if ( is_DVH_ )
   {
      histos_1D[plot_index][fs][cat][Settings::H125VH]->SetFillColor(Cosmetics::VH().fill_color);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      histos_1D[plot_index][fs][cat][Settings::H125VBF]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      histos_1D[plot_index][fs][cat][Settings::H125ttH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
   }
   else
   {
      histos_1D[plot_index][fs][cat][Settings::H125]->SetFillColor(Cosmetics::Higgs_all().fill_color);
   }

   histos_1D[plot_index][fs][cat][Settings::qqZZ]->SetFillColor(Cosmetics::qqZZ().fill_color);
   histos_1D[plot_index][fs][cat][Settings::ggZZ]->SetFillColor(Cosmetics::ggZZ().fill_color);
   histos_1D[plot_index][fs][cat][Settings::VVV]->SetFillColor(Cosmetics::VVV().fill_color);
   histos_1D[plot_index][fs][cat][Settings::qqZZ]->SetLineColor(Cosmetics::qqZZ().line_color);
   histos_1D[plot_index][fs][cat][Settings::ggZZ]->SetLineColor(Cosmetics::ggZZ().line_color);
   histos_1D[plot_index][fs][cat][Settings::VVV]->SetLineColor(Cosmetics::VVV().line_color);

   if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
   {
      histos_1D_ZX_shape[plot_index][fs][cat]->SetFillColor(Cosmetics::ZX().fill_color);
      histos_1D_ZX_shape[plot_index][fs][cat]->SetLineColor(Cosmetics::ZX().line_color);
   }
   else
   {
      histos_1D_ZX[plot_index][fs][cat]->SetFillColor(Cosmetics::ZX().fill_color);
      histos_1D_ZX[plot_index][fs][cat]->SetLineColor(Cosmetics::ZX().line_color);
   }

   histos_1D[plot_index][fs][cat][Settings::H125]->SetLineColor(Cosmetics::Higgs_all().line_color);
   histos_1D[plot_index][fs][cat][Settings::H125ggH]->SetLineColor(Cosmetics::Higgs_all().line_color);
   histos_1D[plot_index][fs][cat][Settings::H125VBF]->SetLineColor(Cosmetics::Higgs_all().line_color);
   histos_1D[plot_index][fs][cat][Settings::H125VH]->SetLineColor(Cosmetics::Higgs_all().line_color);
   histos_1D[plot_index][fs][cat][Settings::H125ttH]->SetLineColor(Cosmetics::Higgs_all().line_color);

   histos_1D[plot_index][fs][cat][Settings::Data]->SetBinErrorOption(TH1::kPoisson);
   histos_1D[plot_index][fs][cat][Settings::Data]->SetLineColor(kBlack);


   // THStack
   THStack *stack = new THStack( "stack", "stack" );

   if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
   {
      stack->Add(histos_1D_ZX_shape[plot_index][fs][cat]);
   }
   else
   {
      stack->Add(histos_1D_ZX[plot_index][fs][cat]);
   }

   stack->Add(histos_1D[plot_index][fs][cat][Settings::ggZZ]);
   stack->Add(histos_1D[plot_index][fs][cat][Settings::qqZZ]);
   stack->Add(histos_1D[plot_index][fs][cat][Settings::VVV]);

   if ( is_Djet_ )
   {
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125VH]);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125ttH]);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125bbH]);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125tqH]);
      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125ggH]);
      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125VBF]);
   }
   else if ( is_DVH_ )
   {
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125VBF]);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125ttH]);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125bbH]);
      histos_1D[plot_index][fs][cat][Settings::H125ggH]->Add(histos_1D[plot_index][fs][cat][Settings::H125tqH]);
      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125ggH]);
      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125VH]);
   }
   else
   {
      stack->Add(histos_1D[plot_index][fs][cat][Settings::H125]);
   }

   stack->Draw("HIST");

   if ( plot_index == Settings::M4lMain ) stack->GetXaxis()->SetRangeUser(70.,500.);

   float data_max = histos_1D[plot_index][fs][cat][Settings::Data]->GetBinContent(histos_1D[plot_index][fs][cat][Settings::Data]->GetMaximumBin());
   float data_max_error = histos_1D[plot_index][fs][cat][Settings::Data]->GetBinErrorUp(histos_1D[plot_index][fs][cat][Settings::Data]->GetMaximumBin());

   if ( GetVarLogY(variable_name) )
   {
      stack->SetMinimum(0.2);
      stack->SetMaximum((data_max + data_max_error)*100);
   }
   else if ( plot_index == Settings::D1jet_M4L118130 || plot_index == Settings::D2jet_M4L118130  || plot_index == Settings::DVH_M4L118130 || plot_index == Settings::DVBFDEC_M4L118130 || plot_index == Settings::DVBFDEC || plot_index == Settings::DVHDEC_M4L118130 || plot_index == Settings::DVHDEC)
   {
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.5);
   }
   else
  {
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.1);
   }

   if ( plot_index == Settings::M4lMain || plot_index == Settings::M4lMainZoomed )
   {
      if      ( fs == Settings::fs4e )    stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_4e);
      else if ( fs == Settings::fs4mu )   stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_4mu);
      else if ( fs == Settings::fs2e2mu ) stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_2e2mu);
      else if ( fs == Settings::fs2mu2e ) stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_2mu2e);
      else                                stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label);
   }
   else
   {
      stack->GetXaxis()->SetTitle(histos_1D[plot_index][fs][cat][Settings::Data]->GetXaxis()->GetTitle());
   }

   stack->GetYaxis()->SetTitle(histos_1D[plot_index][fs][cat][Settings::Data]->GetYaxis()->GetTitle());

   if ( (plot_index == Settings::M4lMainZoomed) || (plot_index == Settings::M4lMainHighMass) ) stack->GetXaxis()->SetNdivisions(1005);

   histos_1D[plot_index][fs][cat][Settings::Data]->SetMarkerSize(0.9);
   histos_1D[plot_index][fs][cat][Settings::Data]->Draw("SAME p E1 X0");


//=============
// L E G E N D
//=============

   TLegend *legend;

   if(variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
   {
      legend  = CreateLegend("right", histos_1D[plot_index][fs][cat][Settings::Data],
                                      histos_1D[plot_index][fs][cat][Settings::H125],
                                      histos_1D[plot_index][fs][cat][Settings::qqZZ],
                                      histos_1D[plot_index][fs][cat][Settings::ggZZ],
                                      histos_1D_ZX_shape[plot_index][fs][cat],
                                      histos_1D[plot_index][fs][cat][Settings::VVV]);
   }
   else if ( plot_index == Settings::D1jet_M4L118130 || plot_index == Settings::D1jet )
   {
      legend  = CreateLegendVBF("left", histos_1D[plot_index][fs][cat][Settings::Data],
                                        histos_1D[plot_index][fs][cat][Settings::H125VBF],
                                        histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VH + ttH
                                        histos_1D[plot_index][fs][cat][Settings::qqZZ],
                                        histos_1D[plot_index][fs][cat][Settings::ggZZ],
                                        histos_1D_ZX[plot_index][fs][cat], 
                                        histos_1D[plot_index][fs][cat][Settings::VVV], false);
   }
   else if ( plot_index == Settings::D2jet_M4L118130 || plot_index == Settings::D2jet  || plot_index == Settings::DVBFDEC )
   {
      legend  = CreateLegendVBF("right", histos_1D[plot_index][fs][cat][Settings::Data],
                                         histos_1D[plot_index][fs][cat][Settings::H125VBF],
                                         histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VH + ttH
                                         histos_1D[plot_index][fs][cat][Settings::qqZZ],
                                         histos_1D[plot_index][fs][cat][Settings::ggZZ],
                                         histos_1D_ZX[plot_index][fs][cat], 
                                         histos_1D[plot_index][fs][cat][Settings::VVV], false);
   }
  else if ( plot_index == Settings::DVBFDEC_M4L118130)
   {
      legend  = CreateLegendVBF("right", histos_1D[plot_index][fs][cat][Settings::Data],
                                         histos_1D[plot_index][fs][cat][Settings::H125VBF],
                                         histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VH + ttH
                                         histos_1D[plot_index][fs][cat][Settings::qqZZ],
                                         histos_1D[plot_index][fs][cat][Settings::ggZZ],
                                         histos_1D_ZX[plot_index][fs][cat], 
                                         histos_1D[plot_index][fs][cat][Settings::VVV], true);
   }
   else if ( plot_index == Settings::DWH_M4L118130 || plot_index == Settings::DWH || plot_index == Settings::DZH_M4L118130 || plot_index == Settings::DZH || plot_index == Settings::DVH_M4L118130 || plot_index == Settings::DVH || plot_index == Settings::DVHDEC )
   {
      legend  = CreateLegendVH("right", histos_1D[plot_index][fs][cat][Settings::Data],
                                        histos_1D[plot_index][fs][cat][Settings::H125VH],
                                        histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VBF + ttH
                                        histos_1D[plot_index][fs][cat][Settings::qqZZ],
                                        histos_1D[plot_index][fs][cat][Settings::ggZZ],
                                        histos_1D_ZX[plot_index][fs][cat], 
                                        histos_1D[plot_index][fs][cat][Settings::VVV], false);
   }
   else if ( plot_index == Settings::DVHDEC_M4L118130)
   {
      legend  = CreateLegendVH("right", histos_1D[plot_index][fs][cat][Settings::Data],
                                        histos_1D[plot_index][fs][cat][Settings::H125VH],
                                        histos_1D[plot_index][fs][cat][Settings::H125ggH], // ggH = ggH + VBF + ttH
                                        histos_1D[plot_index][fs][cat][Settings::qqZZ],
                                        histos_1D[plot_index][fs][cat][Settings::ggZZ],
                                        histos_1D_ZX[plot_index][fs][cat], 
                                        histos_1D[plot_index][fs][cat][Settings::VVV], true);
   }
   else
   {
      legend = CreateLegend((plot_index == Settings::MZ1 || plot_index == Settings::MZ1_M4L118130 || plot_index == Settings::MZ2 || plot_index == Settings::KD_M4L118130) ? "left" : "right",
                             histos_1D[plot_index][fs][cat][Settings::Data],
                             histos_1D[plot_index][fs][cat][Settings::H125],
                             histos_1D[plot_index][fs][cat][Settings::qqZZ],
                             histos_1D[plot_index][fs][cat][Settings::ggZZ],
                             histos_1D_ZX[plot_index][fs][cat],
                             histos_1D[plot_index][fs][cat][Settings::VVV]);
   }

   legend->Draw();

//===========
// PLOT TEXT
//===========

   TPaveText *text;

   if (  plot_index == Settings::KD_M4L118130 || plot_index == Settings::MZ1_M4L118130 )
   {
      text = CreateCutText("right top", "118 < m_{4#font[12]{l}} < 130 GeV");
      text->Draw();
   }
  else if ( plot_index == Settings::D1jet_M4L118130 )
   {
      text = CreateCutText("right top", "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) = 1}");
      text->Draw();
   }
   else if ( plot_index == Settings::MZ2_M4L118130 || plot_index == Settings::DVBFDEC_M4L118130 || plot_index == Settings::DVHDEC_M4L118130)
   {
      text = CreateCutText("left top", "118 < m_{4#font[12]{l}} < 130 GeV");
      text->Draw();
   }
    else if ( plot_index == Settings::D2jet_M4L118130 ||  plot_index == Settings::DWH_M4L118130 || plot_index == Settings::DZH_M4L118130 ||
             plot_index == Settings::DVH_M4L118130 )
   {
      text = CreateCutText("left top", "#splitline{118 < m_{4#font[12]{l}} < 130 GeV}{N(jets) #geq 2}");
      text->Draw();
   }
  else if ( plot_index == Settings::M4l_110150_HighKD)
   {
      if (cat == Settings::VBF_2j_tagged ) text = CreateCutText("right top", "D_{bkg}^{VBF+dec} > 0.5");
      else if (cat == Settings::VH_hadron_tagged ) text = CreateCutText("right top", "D_{bkg}^{VH+dec} > 0.5");
      else text = CreateCutText("left top", "#splitline{Post fit}{D_{bkg}^{kin} > 0.5}");
      text->Draw();
   }

  if ( plot_index == Settings::DVBFDEC_M4L118130 || plot_index == Settings::DVBFDEC || plot_index == Settings::DVHDEC_M4L118130 || plot_index == Settings::DVHDEC)
  {
    text = CreateCatText("left under cut text", _s_category_label.at(cat));
    text->Draw();
  }


   if ( plot_index == Settings::D2jet_M4L118130  || plot_index == Settings::D2jet || plot_index == Settings::DWH_M4L118130 || plot_index == Settings::DWH || plot_index == Settings::DZH_M4L118130 || plot_index == Settings::DZH || plot_index == Settings::DVH_M4L118130 || plot_index == Settings::DVH)
   {
      TLine *wp_line;
      wp_line = CreateDashedLine(0.5, 0.0, 0.5, (data_max + data_max_error)*1.3);
      wp_line->Draw();
   }

   if ( plot_index == Settings::D1jet_M4L118130  || plot_index == Settings::D1jet )
    {
        TLine *wp_line;
        wp_line = CreateDashedLine(0.7, 0.0, 0.7, (data_max + data_max_error)*1.3);
        wp_line->Draw();
    }

//=================
// CMS TEXT & LUMI
//=================

   CMS_lumi *lumi = new CMS_lumi;
   lumi->set_lumi_combination(c);

   // Draw X-axis log scale
   if ( plot_index == Settings::M4lMain )
   {
      stack->GetXaxis()->SetNdivisions(10);
      stack->GetXaxis()->SetLabelSize(0);
      DrawLogX(c, cat, fs);
   }

   _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + _s_final_state.at(fs) + "_" + _s_category.at(cat);
   SavePlots(c, _out_file_name, folder);
}
//=========================================================================================================





//=========================================================================================================
void Histograms::plot_STXS( TString folder )
{
    gStyle->SetPadBottomMargin(0.35);

    TCanvas *c;
    c = new TCanvas("c", "c", 650, 500);

    c->SetLogy();

    STXS_Categories[Settings::H125ggH]->SetFillColor(kRed-7);
    STXS_Categories[Settings::H125bbH]->SetFillColor(kRed-10);
    STXS_Categories[Settings::H125VBF]->SetFillColor(kMagenta+1);
    STXS_Categories[Settings::H125VH]->SetFillColor(kYellow-7);
    STXS_Categories[Settings::H125ttH]->SetFillColor(kOrange+7);
    STXS_Categories[Settings::H125tqH]->SetFillColor(kOrange-3);
    STXS_Categories[Settings::qqZZ]->SetFillColor(Cosmetics::qqZZ().fill_color);
    STXS_Categories[Settings::ggZZ]->SetFillColor(Cosmetics::ggZZ().fill_color);
    STXS_Categories[Settings::Zjets]->SetFillColor(Cosmetics::ZX().fill_color);
    STXS_Categories[Settings::VVV]->SetFillColor(Cosmetics::VVV().fill_color);

    STXS_Categories[Settings::Data]->SetBinErrorOption(TH1::kPoisson);
    STXS_Categories[Settings::Data]->SetLineColor(kBlack);

    // THStack
    THStack *stack = new THStack( "stack", "stack" );

    stack->Add(STXS_Categories[Settings::VVV]);
    stack->Add(STXS_Categories[Settings::Zjets]);
    stack->Add(STXS_Categories[Settings::ggZZ]);
    stack->Add(STXS_Categories[Settings::qqZZ]);
    stack->Add(STXS_Categories[Settings::H125tqH]);
    stack->Add(STXS_Categories[Settings::H125ttH]);
    stack->Add(STXS_Categories[Settings::H125VH]);
    stack->Add(STXS_Categories[Settings::H125VBF]);
    stack->Add(STXS_Categories[Settings::H125bbH]);
    stack->Add(STXS_Categories[Settings::H125ggH]);

    stack->Draw("HIST");

    stack->GetYaxis()->SetTitle(STXS_Categories[Settings::Data]->GetYaxis()->GetTitle());
    stack->GetXaxis()->SetLabelSize(.05);
    stack->GetXaxis()->LabelsOption("v");

    stack->SetMaximum(1e+4);

    STXS_Categories[Settings::Data]->SetMarkerSize(0.9);
    STXS_Categories[Settings::Data]->Draw("SAME p E0 X0");


    //=============
    // L E G E N D
    //=============

    TLegend *legend  = CreateLegendSTXS("right", STXS_Categories[Settings::Data],
                                                 STXS_Categories[Settings::H125ggH],
                                                 STXS_Categories[Settings::H125bbH],
                                                 STXS_Categories[Settings::H125VBF],
                                                 STXS_Categories[Settings::H125VH],
                                                 STXS_Categories[Settings::H125ttH],
                                                 STXS_Categories[Settings::H125tqH],
                                                 STXS_Categories[Settings::qqZZ],
                                                 STXS_Categories[Settings::ggZZ],
                                                 STXS_Categories[Settings::Zjets],
                                                 STXS_Categories[Settings::VVV],
                                                 false);

    legend->Draw();


    //=================
    // CMS TEXT & LUMI
    //=================

    CMS_lumi *lumi = new CMS_lumi;
    lumi->set_lumi_combination(c);

    _out_file_name = folder + "/STXS_Categorization";
    SavePlots(c, _out_file_name, folder);
}
//=========================================================================================================


//=========================================================================================================
void Histograms::plot_Purity( TString folder )
{
    gStyle->SetPadLeftMargin(0.20);
    gStyle->SetHistTopMargin(0.);

    TCanvas *c;
    c = new TCanvas("c", "c", 650, 500);

    TPad* plot_pad = new TPad("plot_pad", "plot_pad", 0.0, 0.05, 0.88, 0.85);
    TPad* legend_pad = new TPad("legend_pad", "legend_pad", 0.1762, 0.81, 0.854, 0.95);
    TPad* text_pad = new TPad("text_pad", "text_pad", 0.86, 0.0, 0.96, 0.95);

    plot_pad->Draw();
    legend_pad->Draw();
    text_pad->Draw();

    plot_pad->cd();

    Purity_Categories[Settings::bin_ggH_0J_PTH_0_10]->SetFillColor(kRed-7);
    Purity_Categories[Settings::bin_ggH_0J_PTH_10_200]->SetFillColor(kRed);
    Purity_Categories[Settings::bin_ggH_1J_PTH_0_60]->SetFillColor(kRed-9);
    Purity_Categories[Settings::bin_ggH_1J_PTH_60_120]->SetFillColor(kRed-10);
    Purity_Categories[Settings::bin_ggH_1J_PTH_120_200]->SetFillColor(kRed-8);
    Purity_Categories[Settings::bin_ggH_2J_PTH_0_60]->SetFillColor(kRed-5);
    Purity_Categories[Settings::bin_ggH_2J_PTH_60_120]->SetFillColor(kRed-6);
    Purity_Categories[Settings::bin_ggH_2J_PTH_120_200]->SetFillColor(kRed-2);
    Purity_Categories[Settings::bin_ggH_PTH_200]->SetFillColor(kRed-3);
    Purity_Categories[Settings::bin_ggH_VBF]->SetFillColor(kRed+2);

    Purity_Categories[Settings::bin_VH_had]->SetFillColor(kMagenta-4);
    Purity_Categories[Settings::bin_VBF_GT200]->SetFillColor(kMagenta-3);
    Purity_Categories[Settings::bin_VBF_2j_mjj_350_700_2j]->SetFillColor(kMagenta-2);
    Purity_Categories[Settings::bin_VBF_2j_mjj_GT700_2j]->SetFillColor(kMagenta);
    Purity_Categories[Settings::bin_VBF_2j_mjj_GT350_3j]->SetFillColor(kMagenta+1);
    Purity_Categories[Settings::bin_VBF_Rest]->SetFillColor(kMagenta-1);

    Purity_Categories[Settings::bin_VH_Lep_0_150]->SetFillColor(kYellow-7);
    Purity_Categories[Settings::bin_VH_Lep_GT150]->SetFillColor(kYellow-6);

    Purity_Categories[Settings::bin_bbH]->SetFillColor(kRed+1);
    Purity_Categories[Settings::bin_ttH]->SetFillColor(kOrange+7);
    Purity_Categories[Settings::bin_tH]->SetFillColor(kOrange-3);


    // THStack
    THStack *stack = new THStack( "stack", "stack" );

    for ( int i_proc = 0; i_proc < num_of_STXS_bins; i_proc++ )
    {
        Purity_Categories[i_proc]->SetBarWidth(0.8);
        Purity_Categories[i_proc]->SetBarOffset(0.1);
        for ( int i_cat = 0; i_cat < num_of_STXS_categories; i_cat++ ) Purity_Categories[i_proc]->GetXaxis()->SetBinLabel(i_cat + 1,_s_STXS_category.at(num_of_STXS_categories - (i_cat+1)));
        stack->Add(Purity_Categories[i_proc]);
    }

    stack->Draw("hbar");

    stack->GetYaxis()->SetTickLength(0.01);
    stack->SetMaximum(1.0);

    stack->GetXaxis()->SetTitle("Reconstructed category");
    stack->GetXaxis()->SetTitleOffset(2.3);
    stack->GetYaxis()->SetTitle("Signal fraction");

    gPad->RedrawAxis();

    text_pad->cd();

    TText *text = new TText(0.,0.,"");
    TString text_label;
    // 'ggH-0j/pT[0,10]', 'ggH-0j/pT[10-200]', 'ggH-1j/pT[0-60]', 'ggH-1j/pT[60-120]', 'ggH-1j/pT[120-200]', 'ggH-2j/pT[0-60]', 
    // 'ggH-2j/pT[60-120]', 'ggH-2j/pT[120-200]', 'ggH/pT$>$200', 'ggH-2j/mJJ$>$350', 'VBF-1j', 'VBF-rest', 'VBF-2j/mJJ[350,700]', 
    // 'VBF-2j/mJJ$>$700', 'VBF-3j/mJJ$>$350', 'VBF-2j/pT$>$200', 'VH-had/mJJ[60-120]', 'VH-rest', 'VH-lep/pTV[0-150]', 'VH-lep/pTV$>$150', 'ttH-lep', 'ttH-had'
    float exp_events[] = {0.88, 0.60, 1.38, 0.80, 6.80, 0.40, 2.23, 1.21, 1.89, 3.44, 4.74, 17.81, 1.07, 3.81, 3.69, 5.96, 3.84, 4.35, 14.78, 27.33, 91.97, 25.79};//[FIXME] Don't hard code these numbers

    //2016
    // float exp_events[] = {0.02, 0.09, 0.42, 0.21, 1.88, 0.13, 0.63, 0.33, 0.51, 0.91, 1.23, 4.62, 0.31, 1.01, 1.01, 1.59, 0.94, 1.11, 3.74, 6.58, 24.34, 6.57};

    // 2017
    // float exp_events[] = {0.34, 0.20, 0.35, 0.24, 1.88, 0.09, 0.64, 0.35, 0.55, 1.07, 1.60, 5.31, 0.35, 1.10, 1.14, 1.87, 1.43, 1.26, 4.37, 9.40, 27.27, 7.78};

    // 2018
    // float exp_events[] = {0.51, 0.31, 0.61, 0.35, 3.05, 0.18, 0.97, 0.53, 0.83, 1.46, 1.92, 7.89, 0.40, 1.70, 1.54, 2.50, 1.46, 1.98, 6.66, 11.34, 40.36, 11.43};

    for ( int i = 0; i < num_of_STXS_categories; i++ )
    {
       if ( exp_events[i] < 9.99) text_label.Form("%.2f", exp_events[i]);
       else text_label.Form("%.1f", exp_events[i]);
       text->SetText(0.30,0.1668 + 0.03142 * float(i),text_label);
       text->SetNDC();
       text->SetTextColor(kBlack);
       text->SetTextSize(0.2);
       text->DrawClone();
    }
    text->SetText(0.09,0.93,"Expected");
    text->DrawClone();
    text->SetText(0.16,0.90,"events");
    text->DrawClone();

    TBox  *box = new TBox(0,0.162,1.0,1.0);
    box->SetLineColor(kBlack);
    box->SetFillStyle(0);
    box->SetLineWidth(1.5);
    box->Draw();

    //=============
    // L E G E N D
    //=============

    TLegend *leg;

    legend_pad->cd();

    leg = new TLegend(0., 0., 1.0, 1.0);

    leg->SetBorderSize(1.);
    leg->SetNColumns(4);

    leg->AddEntry( Purity_Categories[Settings::bin_ggH_0J_PTH_0_10], _s_STXS_bins.at(Settings::bin_ggH_0J_PTH_0_10), "f" );//column 1
    leg->AddEntry( Purity_Categories[Settings::bin_ggH_2J_PTH_60_120], _s_STXS_bins.at(Settings::bin_ggH_2J_PTH_60_120), "f" );//column 2
    leg->AddEntry( Purity_Categories[Settings::bin_VH_had], _s_STXS_bins.at(Settings::bin_VH_had), "f" );//column 3
    leg->AddEntry( Purity_Categories[Settings::bin_VH_Lep_0_150], _s_STXS_bins.at(Settings::bin_VH_Lep_0_150), "f" );//column 4

    leg->AddEntry( Purity_Categories[Settings::bin_ggH_0J_PTH_10_200], _s_STXS_bins.at(Settings::bin_ggH_0J_PTH_10_200), "f" );//column 1
    leg->AddEntry( Purity_Categories[Settings::bin_ggH_2J_PTH_120_200], _s_STXS_bins.at(Settings::bin_ggH_2J_PTH_120_200), "f" );//column 2
    leg->AddEntry( Purity_Categories[Settings::bin_VBF_GT200], _s_STXS_bins.at(Settings::bin_VBF_GT200), "f" );//column 3
    leg->AddEntry( Purity_Categories[Settings::bin_VH_Lep_GT150], _s_STXS_bins.at(Settings::bin_VH_Lep_GT150), "f" );//column 4

    leg->AddEntry( Purity_Categories[Settings::bin_ggH_1J_PTH_0_60], _s_STXS_bins.at(Settings::bin_ggH_1J_PTH_0_60), "f" );//column 1
    leg->AddEntry( Purity_Categories[Settings::bin_ggH_PTH_200], _s_STXS_bins.at(Settings::bin_ggH_PTH_200), "f" );//column 2
    leg->AddEntry( Purity_Categories[Settings::bin_VBF_2j_mjj_350_700_2j], _s_STXS_bins.at(Settings::bin_VBF_2j_mjj_350_700_2j), "f" );//column 3
    leg->AddEntry( Purity_Categories[Settings::bin_bbH], _s_STXS_bins.at(Settings::bin_bbH), "f" );//column 4

    leg->AddEntry( Purity_Categories[Settings::bin_ggH_1J_PTH_60_120], _s_STXS_bins.at(Settings::bin_ggH_1J_PTH_60_120), "f" );//column 1
    leg->AddEntry( Purity_Categories[Settings::bin_ggH_VBF], _s_STXS_bins.at(Settings::bin_ggH_VBF), "f" );//column 2
    leg->AddEntry( Purity_Categories[Settings::bin_VBF_2j_mjj_GT700_2j], _s_STXS_bins.at(Settings::bin_VBF_2j_mjj_GT700_2j), "f" );//column 3
    leg->AddEntry( Purity_Categories[Settings::bin_ttH], _s_STXS_bins.at(Settings::bin_ttH), "f" );//column 4

    leg->AddEntry( Purity_Categories[Settings::bin_ggH_1J_PTH_120_200], _s_STXS_bins.at(Settings::bin_ggH_1J_PTH_120_200), "f" );//column 1
    leg->AddEntry((TObject*)0, "",                       "");//column 2
    leg->AddEntry( Purity_Categories[Settings::bin_VBF_2j_mjj_GT350_3j], _s_STXS_bins.at(Settings::bin_VBF_2j_mjj_GT350_3j), "f" );//column 3
    leg->AddEntry( Purity_Categories[Settings::bin_tH], _s_STXS_bins.at(Settings::bin_tH), "f" );//column 4

    leg->AddEntry( Purity_Categories[Settings::bin_ggH_2J_PTH_0_60], _s_STXS_bins.at(Settings::bin_ggH_2J_PTH_0_60), "f" );//column 1
    leg->AddEntry((TObject*)0, "",                       "");//column 2
    leg->AddEntry( Purity_Categories[Settings::bin_VBF_Rest], _s_STXS_bins.at(Settings::bin_VBF_Rest), "f" );//column 3

    leg->Draw();


    //=================
    // CMS TEXT & LUMI
    //=================

    // Draw lumi
    CMS_lumi *lumi = new CMS_lumi;
    lumi->set_simulation(c);

    _out_file_name = folder + "/STXS_Purity";
    SavePlots(c, _out_file_name, folder);
}
//=========================================================================================================




//==========================================================================================
void Histograms::plot_1D_all_cat( TString filename, TString variable_name , TString folder )
{
   int plot_index = SetPlotName( variable_name );

   TCanvas *c;

   if ( variable_name == "M4lMain" )
   {
      c = new TCanvas(variable_name, variable_name, 650, 500);
   }
   else
   {
      c = new TCanvas(variable_name, variable_name, 600, 600);
   }

   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();


   for ( int i_cat = Settings::inclusive; i_cat >= 0; i_cat-- )
   {
      is_VBF_tagged_ = (i_cat == Settings::VBF_1j_tagged || i_cat == Settings::VBF_2j_tagged);
      is_VH_tagged_  = (i_cat == Settings::VH_lepton_tagged || i_cat == Settings::VH_hadron_tagged || i_cat == Settings::VH_MET_tagged);
      is_ttH_tagged_ = (i_cat == Settings::ttH_hadron_tagged || i_cat == Settings::ttH_lepton_tagged);

      if ( variable_name == "M4lMainZoomed" && is_VBF_tagged_ )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]->SetFillColor(Cosmetics::VBF().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      }
      else if ( variable_name == "M4lMainZoomed" && is_VH_tagged_ )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]->SetFillColor(Cosmetics::VH().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]->SetFillColor(Cosmetics::Higgs_other().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      }
      else if ( variable_name == "M4lMainZoomed" && is_ttH_tagged_ )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]->SetFillColor(Cosmetics::ttH().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]->SetFillColor(Cosmetics::Higgs_other().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]->SetFillColor(Cosmetics::Higgs_other().fill_color);
      }
      else
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125]->SetFillColor(Cosmetics::Higgs_all().fill_color);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125]->SetLineColor(Cosmetics::Higgs_all().line_color);
      }

      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ]->SetFillColor(Cosmetics::qqZZ().fill_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ]->SetFillColor(Cosmetics::ggZZ().fill_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV]->SetFillColor(Cosmetics::VVV().fill_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ]->SetLineColor(Cosmetics::qqZZ().line_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ]->SetLineColor(Cosmetics::ggZZ().line_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV]->SetFillColor(Cosmetics::VVV().fill_color);


      if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
      {
         histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat]->SetFillColor(Cosmetics::ZX().fill_color);
         histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat]->SetLineColor(Cosmetics::ZX().line_color);
      }
      else
      {
         histos_1D_ZX[plot_index][Settings::fs4l][i_cat]->SetFillColor(Cosmetics::ZX().fill_color);
         histos_1D_ZX[plot_index][Settings::fs4l][i_cat]->SetLineColor(Cosmetics::ZX().line_color);
      }


      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125]->SetLineColor(Cosmetics::Higgs_all().line_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->SetLineColor(Cosmetics::Higgs_all().line_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]->SetLineColor(Cosmetics::Higgs_all().line_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]->SetLineColor(Cosmetics::Higgs_all().line_color);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]->SetLineColor(Cosmetics::Higgs_all().line_color);

      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->SetLineColor(kBlack);

      THStack *stack = new THStack();

      if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
      {
         stack->Add(histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat]);
      }
      else
      {
         stack->Add(histos_1D_ZX[plot_index][Settings::fs4l][i_cat]);
      }

      stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ]);
      stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ]);
      stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV]);

      if ( variable_name == "M4lMainZoomed" && is_VBF_tagged_ )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125bbH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125tqH]);

         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]);
         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]);

         Rebin(stack);
      }
      else if ( variable_name == "M4lMainZoomed" && is_VH_tagged_ )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125bbH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125tqH]);

         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]);
         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]);

         Rebin(stack);
      }
      else if ( variable_name == "M4lMainZoomed" && is_ttH_tagged_ )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125bbH]);
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125tqH]);

         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH]);
         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH]);

         Rebin(stack);
      }
      else
      {
         stack->Add(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125]);
      }

      stack->Draw("HIST");

      // Axis title
      stack->GetXaxis()->SetTitle(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->GetXaxis()->GetTitle());
      stack->GetYaxis()->SetTitle(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->GetYaxis()->GetTitle());


      if ( plot_index == Settings::M4lMainZoomed || plot_index == Settings::M4lMainHighMass ) stack->GetXaxis()->SetNdivisions(1005);

      if ( variable_name == "M4lMainZoomed" && (is_VBF_tagged_ || is_VH_tagged_ || is_ttH_tagged_) )
      {
         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->Rebin(2);
         ChangeYaxisTitle(stack);
      }

    float data_max = histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->GetBinContent(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->GetMaximumBin());
      float data_max_error = histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->GetBinErrorUp(histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->GetMaximumBin());

      if ( GetVarLogY( variable_name) )
      {
         stack->SetMinimum(0.5);
         stack->SetMaximum((data_max + data_max_error)*200);
      }

      else
      {
         stack->SetMinimum(1e-5);
         stack->SetMaximum((data_max + data_max_error)*1.6);
      }

      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->SetMarkerSize(0.9);
      histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data]->Draw("SAME p E1 X0");


//=============
// L E G E N D
//=============

      TLegend *legend;

      if ( variable_name == "M4lMain" || variable_name == "M4lMainHighMass" || (variable_name == "M4lMainZoomed" && (i_cat == Settings::untagged || i_cat == Settings::inclusive) ))
      {
         legend  = CreateLegend("right", histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data],
                                         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125],
                                         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ],
                                         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ],
                                         histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat],
                                         histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV]);
      }
      else if ( variable_name == "M4lMainZoomed" && (i_cat == Settings::VBF_1j_tagged || i_cat == Settings::VBF_2j_tagged) )
      {
         legend  = CreateLegendVBF("right", histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VBF],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ],
                                            histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat], 
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV], false);
      }
      else if ( variable_name == "M4lMainZoomed" && (i_cat == Settings::VH_lepton_tagged || i_cat == Settings::VH_hadron_tagged || i_cat == Settings::VH_MET_tagged) )
      {
         legend  = CreateLegendVH("right", histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data],
                                           histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125VH],
                                           histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH],
                                           histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ],
                                           histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ],
                                           histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat], 
                                           histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV], false);
      }
      else if ( variable_name == "M4lMainZoomed" && (i_cat == Settings::ttH_lepton_tagged || i_cat == Settings::ttH_hadron_tagged) )
      {
         legend  = CreateLegendttH("right", histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ttH],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125ggH],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ],
                                            histos_1D_ZX_shape[plot_index][Settings::fs4l][i_cat],
                                            histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV]);
      }
      else
      {
         legend = CreateLegend("right", histos_1D[plot_index][Settings::fs4l][i_cat][Settings::Data],
                                        histos_1D[plot_index][Settings::fs4l][i_cat][Settings::H125],
                                        histos_1D[plot_index][Settings::fs4l][i_cat][Settings::qqZZ],
                                        histos_1D[plot_index][Settings::fs4l][i_cat][Settings::ggZZ],
                                        histos_1D_ZX[plot_index][Settings::fs4l][i_cat],
                                        histos_1D[plot_index][Settings::fs4l][i_cat][Settings::VVV]);
      }

      legend->Draw();

//===========
// PLOT TEXT
//===========

      TPaveText *text;

      if ( (plot_index == Settings::M4lMainZoomed || plot_index == Settings::M4lMain) && i_cat != Settings::inclusive )
      {
         text = CreateCatText("top left", _s_category_label.at(i_cat));
         text->Draw();
      }

    if ( plot_index == Settings::M4l_110150_HighKD)
    {
      if (i_cat == Settings::VBF_2j_tagged ) text = CreateCutText("left top", "D_{bkg}^{VBF+dec} > 0.5");
        else if (i_cat == Settings::VH_hadron_tagged ) text = CreateCutText("left top", "D_{bkg}^{VH+dec} > 0.5");
        else text = CreateCutText("left top", "D_{bkg}^{kin} > 0.5");
      text->Draw();
    }

//=================
// CMS TEXT & LUMI
//=================

      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi_combination(c);

      // Draw X-axis log scale
      if ( plot_index == Settings::M4lMain )
      {
         stack->GetXaxis()->SetNdivisions(10);
         stack->GetXaxis()->SetLabelSize(0);
         DrawLogX(c, i_cat, Settings::fs4l);
      }

      _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + _s_final_state.at(Settings::fs4l) + "_" + _s_category.at(i_cat);
      SavePlots(c, _out_file_name, folder);
   }
   delete c;
}
//==========================================================================================



//=====================================================================================
void Histograms::plot_1D_all_fs( TString filename, TString variable_name , TString folder )
{
   int plot_index = SetPlotName( variable_name);

   TCanvas *c;

   if ( variable_name == "M4lMain")
      c = new TCanvas(variable_name, variable_name, 650, 500);
   else
      c = new TCanvas(variable_name, variable_name, 600, 600);

   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();

   for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
   {
      if( i_fs == Settings::fs2mu2e && _merge_2e2mu) continue;
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::H125]->SetFillColor(Cosmetics::Higgs_all().fill_color);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::qqZZ]->SetFillColor(Cosmetics::qqZZ().fill_color);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::ggZZ]->SetFillColor(Cosmetics::ggZZ().fill_color);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::VVV]->SetFillColor(Cosmetics::VVV().fill_color);

      if ( (variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass") && _merge_2e2mu )
         histos_1D_ZX_shape[plot_index][i_fs][Settings::inclusive]->SetFillColor(Cosmetics::ZX().fill_color);
      else
         histos_1D_ZX[plot_index][i_fs][Settings::inclusive]->SetFillColor(Cosmetics::ZX().fill_color);

      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::H125]->SetLineColor(Cosmetics::Higgs_all().line_color);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::qqZZ]->SetLineColor(Cosmetics::qqZZ().line_color);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::ggZZ]->SetLineColor(Cosmetics::ggZZ().line_color);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::VVV]->SetLineColor(Cosmetics::VVV().line_color);

      if ( (variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass") && _merge_2e2mu )
         histos_1D_ZX_shape[plot_index][i_fs][Settings::inclusive]->SetLineColor(Cosmetics::ZX().line_color);
      else
         histos_1D_ZX[plot_index][i_fs][Settings::inclusive]->SetLineColor(Cosmetics::ZX().line_color);

      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->SetLineColor(kBlack);


      // THStack
      THStack *stack = new THStack( "stack", "stack" );

      if ( (variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass") && _merge_2e2mu )
      {
         stack->Add(histos_1D_ZX_shape[plot_index][i_fs][Settings::inclusive]);
      }
      else
      {
         stack->Add(histos_1D_ZX[plot_index][i_fs][Settings::inclusive]);
      }

      stack->Add(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::VVV]);
      stack->Add(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::ggZZ]);
      stack->Add(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::qqZZ]);
      stack->Add(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::H125]);

      stack->Draw("HIST");

      float data_max = histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->GetBinContent(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->GetMaximumBin());
      float data_max_error = histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->GetBinErrorUp(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->GetMaximumBin());

      if ( GetVarLogY(variable_name) )
      {
         stack->SetMinimum(0.5);
         stack->SetMaximum((data_max + data_max_error)*100);
      }

      else
      {
         stack->SetMinimum(1e-5);
         stack->SetMaximum((data_max + data_max_error)*1.1);
      }


      if ( plot_index == Settings::M4lMain || plot_index == Settings::M4lMainZoomed )
      {
         if      ( i_fs == Settings::fs4e )    stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_4e);
         else if ( i_fs == Settings::fs4mu )   stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_4mu);
         else if ( i_fs == Settings::fs2e2mu ) stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_2e2mu);
         else if ( i_fs == Settings::fs2mu2e ) stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label_2mu2e);
         else                                  stack->GetXaxis()->SetTitle(Variables::M4lMain().var_X_label);;
      }
      else
      {
         stack->GetXaxis()->SetTitle(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->GetXaxis()->GetTitle());
      }

      stack->GetYaxis()->SetTitle(histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->GetYaxis()->GetTitle());


      if ( (plot_index == Settings::M4lMainZoomed) ||  (plot_index == Settings::M4lMainHighMass)) stack->GetXaxis()->SetNdivisions(1005);

    histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->SetMarkerSize(0.9);
      histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data]->Draw("SAME p E1 X0");

      TLegend *legend;

      if ( variable_name == "M4lMain" || variable_name == "M4lMainZoomed" || variable_name == "M4lMainHighMass" )
      {
         legend  = CreateLegend("right", histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data],
                                         histos_1D[plot_index][i_fs][Settings::inclusive][Settings::H125],
                                         histos_1D[plot_index][i_fs][Settings::inclusive][Settings::qqZZ],
                                         histos_1D[plot_index][i_fs][Settings::inclusive][Settings::ggZZ],
                                         histos_1D_ZX_shape[plot_index][i_fs][Settings::inclusive],
                                         histos_1D[plot_index][i_fs][Settings::inclusive][Settings::VVV]);
      }
      else
      {
         legend = CreateLegend("right", histos_1D[plot_index][i_fs][Settings::inclusive][Settings::Data],
                                        histos_1D[plot_index][i_fs][Settings::inclusive][Settings::H125],
                                        histos_1D[plot_index][i_fs][Settings::inclusive][Settings::qqZZ],
                                        histos_1D[plot_index][i_fs][Settings::inclusive][Settings::ggZZ],
                                        histos_1D_ZX[plot_index][i_fs][Settings::inclusive],
                                        histos_1D[plot_index][i_fs][Settings::inclusive][Settings::VVV]);
      }
      legend->Draw();

      //===========
    // PLOT TEXT
    //===========

    TPaveText *text;

    if ( plot_index == Settings::M4l_110150_HighKD)
    {
      text = CreateCutText("left top", "D_{bkg}^{kin} > 0.5");
      text->Draw();
    }

      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi_combination(c);

      // Draw X-axis log scale
      if ( plot_index == Settings::M4lMain )
      {
         stack->GetXaxis()->SetNdivisions(10);
         stack->GetXaxis()->SetLabelSize(0);
         DrawLogX(c, Settings::inclusive, i_fs);
      }

      _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(Settings::inclusive);
      SavePlots(c, _out_file_name, folder);
   }
   delete c;
}
//=====================================================================================



//=================================================================================================
void Histograms::plot_2D_single( TString filename, TString variable_name, TString folder, int cat )
{
   int plot_index = SetPlotName( variable_name );

   gStyle->SetPadLeftMargin(0.14);
   gStyle->SetPadRightMargin(0.15);

   TCanvas *c = new TCanvas(variable_name, variable_name, 600, 600);

   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();

   // Plot MC histogram
   TH2F* stack;

   stack = (TH2F*)histos_2D[plot_index][Settings::fs4l][cat][Settings::H125]->Clone();
   stack->Add(histos_2D[plot_index][Settings::fs4l][cat][Settings::qqZZ]);
   stack->Add(histos_2D[plot_index][Settings::fs4l][cat][Settings::ggZZ]);
   stack->Add(histos_2D[plot_index][Settings::fs4l][cat][Settings::VVV]);

   stack->SetTitleOffset(1.00, "X");
   stack->SetTitleOffset(1.25, "YZ");
   stack->SetLabelOffset(0.01, "XYZ");
   stack->SetTitleSize(0.05, "XYZ");
   stack->SetLabelSize(0.04, "XYZ");
   stack->GetZaxis()->SetTitle("Events / bin");

   stack->SetMinimum(+1e-10);
   MakeCOLZGrey(true);
   stack->Draw("COLZ");

   //Plot data histograms
   histos_2D[plot_index][Settings::fs2e2mu][cat][Settings::Data]->SetMarkerStyle(21);
   histos_2D[plot_index][Settings::fs2e2mu][cat][Settings::Data]->SetMarkerSize(0.9);
   histos_2D[plot_index][Settings::fs2e2mu][cat][Settings::Data]->SetMarkerColor(kBlue);
   histos_2D[plot_index][Settings::fs2e2mu][cat][Settings::Data]->SetLineColor(kBlue);

   histos_2D[plot_index][Settings::fs2e2mu][cat][Settings::Data]->Draw("SAME XP");

   histos_2D[plot_index][Settings::fs4mu][cat][Settings::Data]->SetMarkerStyle(20);
   histos_2D[plot_index][Settings::fs4mu][cat][Settings::Data]->SetMarkerSize(0.9);
   histos_2D[plot_index][Settings::fs4mu][cat][Settings::Data]->SetMarkerColor(kRed+1);
   histos_2D[plot_index][Settings::fs4mu][cat][Settings::Data]->SetLineColor(kRed+1);

   histos_2D[plot_index][Settings::fs4mu][cat][Settings::Data]->Draw("SAME XP");

   histos_2D[plot_index][Settings::fs4e][cat][Settings::Data]->SetMarkerStyle(22);
   histos_2D[plot_index][Settings::fs4e][cat][Settings::Data]->SetMarkerSize(0.9);
   histos_2D[plot_index][Settings::fs4e][cat][Settings::Data]->SetMarkerColor(kGreen+2);
   histos_2D[plot_index][Settings::fs4e][cat][Settings::Data]->SetLineColor(kGreen+2);


   histos_2D[plot_index][Settings::fs4e][cat][Settings::Data]->Draw("SAME XP");

   // Draw legend
   TLegend *legend;
   legend = Create2DLegend( "left", histos_2D[plot_index][Settings::fs4e][cat][Settings::Data],
                                    histos_2D[plot_index][Settings::fs4mu][cat][Settings::Data],
                                    histos_2D[plot_index][Settings::fs2e2mu][cat][Settings::Data] );
   legend->Draw();

   TPaveText *text;
   if ( plot_index == Settings::MZ1vsMZ2_M4L118130 )
   {
      text = CreateCutText("left under 2D legend", "118 < m_{4#font[12]{l}} < 130 GeV");
      text->Draw();
   }

   // Adjust color axis
   c->Update();
   TPaletteAxis *pal = (TPaletteAxis*)stack->GetListOfFunctions()->FindObject("palette");
   pal->SetX1NDC(0.855);
   pal->SetX2NDC(0.875);
   pal->SetY1NDC(0.130);
   pal->SetY2NDC(0.950);

   // Draw lumi
   CMS_lumi *lumi = new CMS_lumi;
   lumi->set_lumi_combination(c);

   _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + _s_category.at(cat);
   SavePlots(c, _out_file_name, folder);

   delete c;
}
//=================================================================================================



//=======================================================================================================
void Histograms::plot_2D_error_single( TString filename, TString variable_name, TString folder, int cat )
{
   int plot_index = SetPlotName( variable_name );

   gStyle->SetPadLeftMargin(0.14);
   gStyle->SetPadRightMargin(0.15);

   TCanvas *c = new TCanvas(variable_name, variable_name, 600, 600);

   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();

   // Plot MC histogram
   TH2F* stack;
   stack = (TH2F*)histos_2DError[plot_index][Settings::fs4l][cat][Settings::H125]->Clone();
   stack->Add(histos_2DError[plot_index][Settings::fs4l][cat][Settings::qqZZ]);
   stack->Add(histos_2DError[plot_index][Settings::fs4l][cat][Settings::ggZZ]);
   stack->Add(histos_2DError[plot_index][Settings::fs4l][cat][Settings::VVV]);

   stack->SetTitleOffset(1.00, "X");
   stack->SetTitleOffset(1.25, "YZ");
   stack->SetLabelOffset(0.01, "XYZ");
   stack->SetTitleSize(0.05, "XYZ");
   stack->SetLabelSize(0.04, "XYZ");
   stack->GetZaxis()->SetTitle("Events / bin");

   MakeCOLZGrey(false);
   stack->Draw("COLZ");

   // Plot data histograms
   histos_2DError_data[plot_index][Settings::fs2e2mu][cat]->SetMarkerStyle(21);
   histos_2DError_data[plot_index][Settings::fs2e2mu][cat]->SetMarkerSize(0.9);
   histos_2DError_data[plot_index][Settings::fs2e2mu][cat]->SetMarkerColor(kBlue);
   histos_2DError_data[plot_index][Settings::fs2e2mu][cat]->SetLineColor(kBlue);

   histos_2DError_data[plot_index][Settings::fs2e2mu][cat]->Draw("SAME P");

   histos_2DError_data[plot_index][Settings::fs4mu][cat]->SetMarkerStyle(20);
   histos_2DError_data[plot_index][Settings::fs4mu][cat]->SetMarkerSize(0.9);
   histos_2DError_data[plot_index][Settings::fs4mu][cat]->SetMarkerColor(kRed+1);
   histos_2DError_data[plot_index][Settings::fs4mu][cat]->SetLineColor(kRed+1);

   histos_2DError_data[plot_index][Settings::fs4mu][cat]->Draw("SAME P");

   histos_2DError_data[plot_index][Settings::fs4e][cat]->SetMarkerStyle(22);
   histos_2DError_data[plot_index][Settings::fs4e][cat]->SetMarkerSize(0.9);
   histos_2DError_data[plot_index][Settings::fs4e][cat]->SetMarkerColor(kGreen+2);
   histos_2DError_data[plot_index][Settings::fs4e][cat]->SetLineColor(kGreen+2);

   histos_2DError_data[plot_index][Settings::fs4e][cat]->Draw("SAME P");

   //Draw legend
   TLegend *legend;
   legend = Create2DErrorLegend( "right", histos_2DError_data[plot_index][Settings::fs4e][cat],
                                          histos_2DError_data[plot_index][Settings::fs4mu][cat],
                                          histos_2DError_data[plot_index][Settings::fs2e2mu][cat] );
   legend->Draw();

   // Adjust color axis
   c->Update();

   TPaletteAxis *pal = (TPaletteAxis*)stack->GetListOfFunctions()->FindObject("palette");
   pal->SetX1NDC(0.855);
   pal->SetX2NDC(0.875);
   pal->SetY1NDC(0.130);
   pal->SetY2NDC(0.840);

   // Draw lumi
   CMS_lumi *lumi = new CMS_lumi;
   lumi->set_lumi_combination(c);

   _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + _s_category.at(cat);
   SavePlots(c, _out_file_name, folder);

   delete c;
}
//=======================================================================================================



//===============================================================================================
void Histograms::plot_2D_error_all_cat( TString filename, TString variable_name, TString folder )
{
   int plot_index = SetPlotName( variable_name );

   gStyle->SetPadLeftMargin(0.14);
   gStyle->SetPadRightMargin(0.15);

   TCanvas *c = new TCanvas(variable_name, variable_name, 600, 600);

   TPad* plot_pad = new TPad("plot_pad", "plot_pad", 0.0, 0.05, 1.0, 0.85);
   TPad* legend_pad = new TPad("legend_pad", "legend_pad", gStyle->GetPadLeftMargin(), 0.81, 1 - gStyle->GetPadRightMargin(), 0.95);

   if ( GetVarLogX( variable_name) ) plot_pad->SetLogx();
   if ( GetVarLogY( variable_name) ) plot_pad->SetLogy();


   plot_pad->Draw();
   legend_pad->Draw();

   plot_pad->cd();

   // Plot MC histogram
   TH2F* stack;
   stack = (TH2F*)histos_2DError[plot_index][Settings::fs4l][Settings::inclusive][Settings::H125]->Clone();
   stack->Add(histos_2DError[plot_index][Settings::fs4l][Settings::inclusive][Settings::qqZZ]);
   stack->Add(histos_2DError[plot_index][Settings::fs4l][Settings::inclusive][Settings::ggZZ]);
   stack->Add(histos_2DError[plot_index][Settings::fs4l][Settings::inclusive][Settings::VVV]);
   histos_2DError_ZX[plot_index][Settings::fs4l][Settings::inclusive]->Smooth();
   stack->Add(histos_2DError_ZX[plot_index][Settings::fs4l][Settings::inclusive]);

   stack->SetTitleOffset(1.00, "X");
   stack->SetTitleOffset(1.25, "YZ");
   stack->SetLabelOffset(0.01, "XYZ");
   stack->SetTitleSize(0.05, "XYZ");
   stack->SetLabelSize(0.04, "XYZ");
   stack->GetZaxis()->SetTitle("Events / bin");

//   stack->GetYaxis()->SetRangeUser(0.0, 1.0);

   MakeCOLZGrey(false);
   stack->Draw("COLZ");

   // Plot data histograms
   for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
   {
      histos_2DError_data[plot_index][i_fs][Settings::untagged]->SetMarkerStyle(20);
      histos_2DError_data[plot_index][i_fs][Settings::untagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::VBF_1j_tagged]->SetMarkerStyle(26);
      histos_2DError_data[plot_index][i_fs][Settings::VBF_1j_tagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::VBF_2j_tagged]->SetMarkerStyle(32);
      histos_2DError_data[plot_index][i_fs][Settings::VBF_2j_tagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::VH_lepton_tagged]->SetMarkerStyle(28);
      histos_2DError_data[plot_index][i_fs][Settings::VH_lepton_tagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::VH_hadron_tagged]->SetMarkerStyle(27);
      histos_2DError_data[plot_index][i_fs][Settings::VH_hadron_tagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::ttH_hadron_tagged]->SetMarkerStyle(30);
      histos_2DError_data[plot_index][i_fs][Settings::ttH_hadron_tagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::ttH_lepton_tagged]->SetMarkerStyle(25);
      histos_2DError_data[plot_index][i_fs][Settings::ttH_lepton_tagged]->SetMarkerSize(0.9);

      histos_2DError_data[plot_index][i_fs][Settings::VH_MET_tagged]->SetMarkerStyle(25);
      histos_2DError_data[plot_index][i_fs][Settings::VH_MET_tagged]->SetMarkerSize(0.9);
   }

   for (int i_cat = 0; i_cat < Settings::inclusive; i_cat++)
   {

      histos_2DError_data[plot_index][Settings::fs4e][i_cat]->SetMarkerColor(kGreen+2);
      histos_2DError_data[plot_index][Settings::fs4e][i_cat]->SetLineColor(kGreen+2);

      histos_2DError_data[plot_index][Settings::fs4mu][i_cat]->SetMarkerColor(kRed+1);
      histos_2DError_data[plot_index][Settings::fs4mu][i_cat]->SetLineColor(kRed+1);

      histos_2DError_data[plot_index][Settings::fs2e2mu][i_cat]->SetMarkerColor(kBlue);
      histos_2DError_data[plot_index][Settings::fs2e2mu][i_cat]->SetLineColor(kBlue);

      histos_2DError_data[plot_index][Settings::fs4l][i_cat]->SetMarkerColor(kBlack);
      histos_2DError_data[plot_index][Settings::fs4l][i_cat]->SetLineColor(kBlack);
   }

   for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat < Settings::inclusive; i_cat++)
      {
         if(plot_index == Settings::DVBFDECvsM4lZoomed && i_cat == Settings::VBF_2j_tagged) histos_2DError_data[plot_index][i_fs][i_cat]->Draw("SAME P");
      else if(plot_index == Settings::DVHDECvsM4lZoomed && i_cat == Settings::VH_hadron_tagged) histos_2DError_data[plot_index][i_fs][i_cat]->Draw("SAME P");
      else histos_2DError_data[plot_index][i_fs][i_cat]->Draw("SAME P");
      }
   }

   // Adjust axis color
   c->Update();
   TPaletteAxis* pal = (TPaletteAxis*)stack->GetListOfFunctions()->FindObject("palette");
   pal->SetX1NDC(0.855);
   pal->SetX2NDC(0.875);
   pal->SetY1NDC(0.130);
   pal->SetY2NDC(0.950);

   if ( plot_index == Settings::DVHvsM4lZoomed || plot_index == Settings::DWHvsM4lZoomed || plot_index == Settings::DZHvsM4lZoomed ||
        plot_index == Settings::D1jetvsM4lZoomed || plot_index == Settings::D2jetvsM4lZoomed )
   {
      TLine *wp_line;
      wp_line = CreateDashedLine(100., 0.5, 170., 0.5);
      wp_line->Draw();
   }

   //Draw legend
   legend_pad->cd();
   TLegend *legend = new TLegend(0.,0.,0.,0.);
   if ( plot_index == Settings::KDvsM4lZoomed) legend = Create2DLegendAllCat_KD( "top", histos_2DError_data[plot_index][Settings::fs4e][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4mu][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs2e2mu][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::VBF_1j_tagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::VH_lepton_tagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::ttH_hadron_tagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::ttH_lepton_tagged]);

  if ( plot_index == Settings::DVBFDECvsM4lZoomed) legend = Create2DLegendAllCat_DVBFDEC( "top", histos_2DError_data[plot_index][Settings::fs4e][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4mu][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs2e2mu][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::VBF_2j_tagged]);
   if ( plot_index == Settings::DVHDECvsM4lZoomed) legend = Create2DLegendAllCat_DVHDEC( "top", histos_2DError_data[plot_index][Settings::fs4e][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4mu][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs2e2mu][Settings::untagged],
                                         histos_2DError_data[plot_index][Settings::fs4l][Settings::VH_hadron_tagged]);
   legend->Draw();

   // Draw lumi
   CMS_lumi *lumi = new CMS_lumi;
   lumi->set_lumi_combination(c);

   _out_file_name = folder + "/" + variable_name + "_" + filename + "_" + "all_categories";
   SavePlots(c, _out_file_name, folder);

   delete c;
}
//===============================================================================================



//==============================================================
void Histograms::FillYieldGraphs( float M4l_down, float M4l_up , TString fit_option = "")
{
   system("mkdir -p ROOT_files");
   TFile *f_signal_fits = new TFile("ROOT_files/Signal_fits.root", "recreate");

   gStyle->SetOptFit(0);

  TCanvas *c = new TCanvas("c1", "c1", 1500, 1500);

   TPaveText* pav = new TPaveText(0.13,0.93,0.8,1.,"brNDC");
   pav->SetFillStyle(0);
   pav->SetBorderSize(0);
   pav->SetTextAlign(11);
   pav->SetTextSize(0.04);

   TLegend* lgd;
   lgd = new TLegend(0.2,0.73,0.4,0.88);
   lgd->SetFillStyle(0);



   int process = 0;
   double temp_yield;
   double temp_error;
   vector<double> yield;
   vector<double> error;

   int fs_marker[num_of_final_states] = {20, 22, 21, 33, 29};
   Color_t catColor[num_of_categories] = {kBlue-9, kCyan-6, kGreen-6, kRed-7, kOrange+6, kMagenta-6, kYellow - 3 ,kBlack};

   for ( int i_prod_mode = 0; i_prod_mode < num_of_production_modes - 2; i_prod_mode++ ) // -2 because we want only signal
   {
      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
         for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
         {
            if( i_fs == Settings::fs2mu2e ) continue;

            for( int i_mass_point = 0; i_mass_point < 5; i_mass_point++)
            {
               process = SetProcess(i_mass_point, i_prod_mode);

               temp_yield = histos_1D[Settings::M4lYields][i_fs][i_cat][process]->IntegralAndError(
                            histos_1D[Settings::M4lYields][i_fs][i_cat][process]->FindBin(M4l_down),
                            histos_1D[Settings::M4lYields][i_fs][i_cat][process]->FindBin(M4l_up) - 1, temp_error);

               //cout << " " << temp_yield;

               yield.push_back(temp_yield);
               error.push_back(temp_error);
               mass_points.push_back(SetMassPoint(i_mass_point));

            } // i_mass_point

            yields_graph[i_fs][i_cat][i_prod_mode] = new TGraphErrors(yield.size(), &(mass_points[0]), &(yield[0]), 0, &(error[0]));

            yield.clear();
            error.clear();
            mass_points.clear();

            // Cosmetics
            yields_graph[i_fs][i_cat][i_prod_mode]->GetXaxis()->SetLimits( 115, 135);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetMinimum(0.);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetMarkerStyle(22);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetMarkerSize(3);
            yields_graph[i_fs][i_cat][i_prod_mode]->GetXaxis()->SetTitle("generated m_{H}");
            yields_graph[i_fs][i_cat][i_prod_mode]->GetYaxis()->SetTitle(Form("expected yield in %.1f fb^{-1}", _lumi));
            yields_graph[i_fs][i_cat][i_prod_mode]->SetMarkerStyle(fs_marker[i_fs]);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetMarkerColor(catColor[i_cat]);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetLineColor(kBlack);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetLineWidth(2);


            // Set fit model
            TString fit_model;
            fit_model = "pol2";

            if( i_prod_mode == Settings::bbH
            || (i_prod_mode == Settings::WH_had && i_cat != Settings::untagged)
            || (i_prod_mode == Settings::WH_lep && i_cat != Settings::untagged)
            || (i_prod_mode == Settings::ZH_had && i_cat != Settings::untagged)
            || (i_prod_mode == Settings::ZH_lep)
            || (i_prod_mode == Settings::ttH_had && i_cat != Settings::untagged)
            || (i_prod_mode == Settings::ttH_lep && i_cat != Settings::untagged)
            || (i_prod_mode == Settings::ggH && (i_cat != Settings::untagged && i_cat != Settings::VBF_1j_tagged && i_cat != Settings::VBF_2j_tagged))
            || (i_prod_mode == Settings::VBF && (i_cat != Settings::untagged && i_cat != Settings::VBF_1j_tagged && i_cat != Settings::VBF_2j_tagged))) fit_model = "pol1";

            if( i_prod_mode == Settings::tqH
            || (i_prod_mode == Settings::WH_had && i_cat == Settings::ttH_lepton_tagged)
//            || (i_prod_mode == Settings::WH_had && i_cat == Settings::ttH_hadron_tagged)
//            || (i_prod_mode == Settings::WH_had && i_cat == Settings::VH_lepton_tagged)
            || (i_prod_mode == Settings::ZH_had && i_cat == Settings::ttH_lepton_tagged) ) fit_model = "pol0";
//            || (i_prod_mode == Settings::ZH_had && i_cat == Settings::ttH_hadron_tagged)
//            || (i_prod_mode == Settings::ZH_had && i_cat == Settings::VH_lepton_tagged)
//        || (i_prod_mode == Settings::WH_lep && i_cat == Settings::ttH_lepton_tagged)
//            || (i_prod_mode == Settings::WH_lep && i_cat == Settings::ttH_hadron_tagged)
//            || (i_prod_mode == Settings::WH_lep && i_cat == Settings::VH_hadron_tagged)
//            || (i_prod_mode == Settings::ZH_lep && i_cat == Settings::ttH_lepton_tagged)
//            || (i_prod_mode == Settings::ZH_lep && i_cat == Settings::ttH_hadron_tagged)
//            || (i_prod_mode == Settings::ZH_lep && i_cat == Settings::VH_hadron_tagged)

            TF1 *fit_function = new TF1("fit_function", fit_model, 119, 131);
            fit_function->SetLineColor(catColor[i_cat]);
            fit_function->SetLineWidth(2);

            yields_graph[i_fs][i_cat][i_prod_mode]->Fit(fit_function, "R" + fit_option);

            _fit_funct_name = "f_" + _s_production_mode.at(i_prod_mode) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
            fit_function->SetName(_fit_funct_name);
            fit_function->Write();

             _graph_name = "g_" + _s_production_mode.at(i_prod_mode) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
            yields_graph[i_fs][i_cat][i_prod_mode]->SetName(_graph_name);
            yields_graph[i_fs][i_cat][i_prod_mode]->Write();

            lgd->AddEntry(yields_graph[i_fs][i_cat][i_prod_mode],_s_final_state[i_fs],"P");

         } // i_fs

      c->cd();
      pav->AddText(_s_production_mode[i_prod_mode] + ", " + _s_category[i_cat]);

      yields_graph[Settings::fs2e2mu][i_cat][i_prod_mode]->Draw("AP");
      yields_graph[Settings::fs4e][i_cat][i_prod_mode]->Draw("P");
      yields_graph[Settings::fs4mu][i_cat][i_prod_mode]->Draw("P");
      pav->Draw();
      lgd->Draw();

      TString fit_name = _s_production_mode.at(i_prod_mode) + "_" + _s_category.at(i_cat);

       _out_file_name = "Fits/" + _s_production_mode.at(i_prod_mode) + "_" + _s_category.at(i_cat);
       SavePlots(c, _out_file_name, "Fits");

      pav->Clear();
      lgd->Clear();
      c->Clear();

      } // i_cat
   } // i_prod_mode

   f_signal_fits->Close();
   delete f_signal_fits;
}
//==============================================================



//========================================================================================================
void Histograms::PrepareYamlFiles( TString sqrt, float M4l_down, float M4l_up, vector< vector <float> > _expected_yield_SR)
{

   TString out_file_name[Settings::fs4l];
   ofstream out_file[Settings::fs4l];
   TString folder_name = "YAML_files";
   system("mkdir -p " + folder_name);

   TFile *f_signal_fits;
   f_signal_fits = TFile::Open("ROOT_files/Signal_fits.root");

   TF1* fit_function;
   M4lZX *ZXYields = new M4lZX();

   int num_of_parameters;

   for ( int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++ )
   {
      if( i_fs == Settings::fs2mu2e ) continue;

      out_file_name[i_fs] = folder_name + "/yields_per_tag_category_" + sqrt + "TeV_" + _s_final_state.at(i_fs) + ".yaml";
      out_file[i_fs].open(out_file_name[i_fs]);

      out_file[i_fs] << "---" << endl;
      out_file[i_fs] << "# sqrt(s) = " << sqrt << " TeV" << endl;
      out_file[i_fs] << "# integrated luminosity = " << _lumi << " fb-1" << endl;
      out_file[i_fs] << endl;
      out_file[i_fs] << "mass_range: '" << M4l_down << ", " << M4l_up << "'" << endl;
      out_file[i_fs] << "kd_range: '0, 1'"<<endl;
      out_file[i_fs] << endl;

      for ( int i_cat = 0; i_cat < num_of_categories - 1; i_cat++ )
      {
         out_file[i_fs] << _s_category.at(i_cat) << ":" << endl;

         for ( int i_prod_mode = 0; i_prod_mode < num_of_production_modes - 2; i_prod_mode++ )
         {
            out_file[i_fs] <<  "    " << _s_production_mode.at(i_prod_mode) << "_hzz: ";

            _fit_funct_name = "f_" + _s_production_mode.at(i_prod_mode) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
            fit_function = (TF1*)f_signal_fits->Get(_fit_funct_name);
            num_of_parameters = fit_function->GetNpar();
            out_file[i_fs] << "TMath::Max(0,";
            for ( int i_par = 0; i_par <= 2; i_par++ )
            {
               double parameter = 0;
               if(i_par < num_of_parameters) parameter = fit_function->GetParameter(i_par);

               out_file[i_fs] << "(" << (parameter);

               for ( int i_par_2 = 0; i_par_2 <= i_par-1; i_par_2++)
                  out_file[i_fs] << "*@0";

               out_file[i_fs]<<")";

               if ( i_par < 2)
               {
                  out_file[i_fs] << "+";
               }
            }
            out_file[i_fs] << ")";
            out_file[i_fs] << endl;
         }
         out_file[i_fs] <<  "    "  << _s_process.at(Settings::yqqZZ) << "_hzz: '" << histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->Integral(
                                                                      histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->FindBin(M4l_down),
                                                                      histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->FindBin(M4l_up) - 1) << "'" << endl;

         out_file[i_fs] <<  "    "  << _s_process.at(Settings::yggZZ) << "_hzz: '" << histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->Integral(
                                                                     histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->FindBin(M4l_down),
                                                                     histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->FindBin(M4l_up) - 1) << "'" << endl;

         out_file[i_fs] <<  "    "  << "zjets_hzz: '" << ZXYields->GetM4lZX_Yields(_expected_yield_SR, M4l_down, M4l_up, i_fs, i_cat) << "'" << endl;
         out_file[i_fs] << endl;

      } // end i_cat

      out_file[i_fs].close();

   } // end i_fs
}
//========================================================================================================




//========================================================================================================
void Histograms::PrintYields( vector< vector <float> > _expected_yield_SR)
{

   cout << std::setprecision(2) << fixed;

   cout << endl;
   cout << "=========================" << endl;
   cout << "Yields in full mass range" << endl;
   cout << "=========================" << endl;
   cout << _s_final_state.at(0) << "   " << _s_final_state.at(1) << "   " << _s_final_state.at(2) << "   " << _s_final_state.at(4) << endl;

   for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
   {

      if ((i_proc > Settings::yH130 && i_proc < Settings::yqqZZ) || i_proc == Settings::yH120 || i_proc == Settings::yH124
          || i_proc == Settings::yH126 || i_proc == Settings::yH130 || i_proc == Settings::yDY || i_proc == Settings::yttbar) continue;

      cout << endl;

      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         cout << _s_process.at(i_proc) << "   " << _s_category.at(i_cat) << "   ";

         for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
         {
            if( _merge_2e2mu && i_fs == Settings::fs2mu2e ) continue;

            cout << histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]->Integral() << "   ";
         }
         cout << endl;
      }
   }

   // Z+X

   M4lZX *ZXYields = new M4lZX();

   cout << endl;

   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      cout << "Z+X" << "   " << _s_category.at(i_cat) << "   ";

      for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
      {
         if( _merge_2e2mu && i_fs == Settings::fs2mu2e ) continue;

          cout << ZXYields->GetM4lZX_Yields(_expected_yield_SR, 0, 3000, i_fs, i_cat) << "   ";
      }
      cout << endl;
   }

  // Total MC
   cout << endl;

   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      cout << "Total" << "   " << _s_category.at(i_cat) << "   ";

      for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
      {
         if( _merge_2e2mu && i_fs == Settings::fs2mu2e ) continue;

          cout << ZXYields->GetM4lZX_Yields(_expected_yield_SR, 0, 3000, i_fs, i_cat) + histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Integral() + histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->Integral() + histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->Integral() << "   ";
      }
      cout << endl;
   }

   delete ZXYields;
}
//========================================================================================================



//========================================================
void Histograms::PrintYields(float M4l_down, float M4l_up, vector< vector <float> > _expected_yield_SR)
{
   //cout << std::setprecision(2) << fixed;

   cout << endl;
   cout << "===============================" << endl;
   cout << "Yields in [" << M4l_down << ", " << M4l_up << "] mass range" << endl;
   cout << "===============================" << endl;
   cout << _s_final_state.at(0) << "   " << _s_final_state.at(1) << "   " << _s_final_state.at(2) << "   " << _s_final_state.at(4) << endl;

   for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
   {

      if ((i_proc > Settings::yH130 && i_proc < Settings::yqqZZ) || i_proc == Settings::yH120 || i_proc == Settings::yH124
          || i_proc == Settings::yH126 || i_proc == Settings::yH130 || i_proc == Settings::yDY || i_proc == Settings::yttbar) continue;

      cout << endl;

      for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
      {
         cout << _s_process.at(i_proc) << "   " << _s_category.at(i_cat) << "   ";

         for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
         {
            if( _merge_2e2mu && i_fs == Settings::fs2mu2e ) continue;

            cout << histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]->Integral(
                    histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]->FindBin(M4l_down),
                    histos_1D[Settings::M4lYields][i_fs][i_cat][i_proc]->FindBin(M4l_up) - 1)
            << "   ";
         }
         cout << endl;
      }
   }

   // Z+X
   M4lZX *ZXYields = new M4lZX();

   cout << endl;

   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      cout << "Z+X" << "   " << _s_category.at(i_cat) << "   ";

      for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
      {
         if( _merge_2e2mu && i_fs == Settings::fs2mu2e ) continue;

         cout << ZXYields->GetM4lZX_Yields(_expected_yield_SR, M4l_down, M4l_up, i_fs, i_cat) << "   ";
      }
      cout << endl;
   }

    // Total MC
   cout << endl;

   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      cout << "Total" << "   " << _s_category.at(i_cat) << "   ";

      for ( int i_fs = 0; i_fs < num_of_final_states; i_fs++ )
      {
         if( _merge_2e2mu && i_fs == Settings::fs2mu2e ) continue;

          cout << ZXYields->GetM4lZX_Yields(
          _expected_yield_SR, M4l_down, M4l_up, i_fs, i_cat) +
          histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->Integral(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->FindBin(M4l_down),histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yH125]->FindBin(M4l_up) - 1) +
          histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->Integral(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->FindBin(M4l_down),histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yqqZZ]->FindBin(M4l_up) - 1) +
          histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->Integral(histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->FindBin(M4l_down),histos_1D[Settings::M4lYields][i_fs][i_cat][Settings::yggZZ]->FindBin(M4l_up) - 1) << "   ";
      }
      cout << endl;
   }

   delete ZXYields;
}
//========================================================



//========================================================
void Histograms::PrintLatexTables(float M4l_down, float M4l_up, vector< vector <float> > _expected_yield_SR)
{
   cout << std::setprecision(2) << fixed;

   map<int, vector<float>> yields_map;
   vector<float> yields_ZX;
   float temp_yield;

   // Z+X
   M4lZX *ZXYields = new M4lZX();

   cout << endl << endl;
   cout << "//==============" << endl;
   cout << "// T a b l e  1 " << endl;
   cout << "//==============" << endl;
   cout << endl;

//==============
// T a b l e  1
//==============

  for (int i_proc = 0; i_proc < num_of_processes_yields; i_proc++)
   {
      if (i_proc == Settings::yH125 || i_proc == Settings::yqqZZ || i_proc == Settings::yggZZ)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            if ( _merge_2e2mu && i_fs == Settings::fs2mu2e) continue;
            temp_yield = histos_1D[Settings::M4lYields][i_fs][Settings::inclusive][i_proc]->Integral();
            yields_map[i_proc].push_back(temp_yield);
         }
      }
   }

   // Z+X
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      if ( _merge_2e2mu && i_fs == Settings::fs2mu2e) continue;
      temp_yield = ZXYields->GetM4lZX_Yields(_expected_yield_SR, 30., 3000., i_fs, Settings::inclusive);
      yields_ZX.push_back(temp_yield);
   }


   cout << "\\begin{tabular}{l|c|c|c|c}" << endl;
   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "\\textbf{Channel} & $4\\Pe$ & $4\\Pgm$ & $2\\Pe2\\Pgm$ & $4\\ell$ \\\\" << endl;
   cout << "\\hline" << endl;


   cout << "\\qqZZ ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yqqZZ].at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\ggZZ ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yggZZ].at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\cPZ\\ + X ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_ZX.at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;

   cout << "Sum of backgrounds ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yqqZZ].at(i_fs) + yields_map[Settings::yggZZ].at(i_fs)+ yields_ZX.at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;

   cout << "Signal ($\\mH=125~\\GeV$) ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yH125].at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;

   cout << "Total expected ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yH125].at(i_fs) + yields_map[Settings::yqqZZ].at(i_fs) + yields_map[Settings::yggZZ].at(i_fs)+ yields_ZX.at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;
   cout << "Observed & XXX & XXX & XXX & XXX \\\\" << endl;

   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "\\end{tabular}" << endl;


   cout << endl << endl;
   cout << "//==============" << endl;
   cout << "// T a b l e  2 " << endl;
   cout << "//==============" << endl;
   cout << endl;

//==============
// T a b l e  2
//==============

   yields_map.clear();
   yields_ZX.clear();


   for (int i_proc = 0; i_proc < num_of_processes_yields; i_proc++)
   {
      if (i_proc == Settings::yH125 || i_proc == Settings::yqqZZ || i_proc == Settings::yggZZ)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            if ( _merge_2e2mu && i_fs == Settings::fs2mu2e) continue;
            temp_yield = histos_1D[Settings::M4lYields][i_fs][Settings::inclusive][i_proc]->Integral(
                         histos_1D[Settings::M4lYields][i_fs][Settings::inclusive][i_proc]->FindBin(M4l_down),
                         histos_1D[Settings::M4lYields][i_fs][Settings::inclusive][i_proc]->FindBin(M4l_up)-1);
            yields_map[i_proc].push_back(temp_yield);
         }
      }
   }

   // Z+X
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      if ( _merge_2e2mu && i_fs == Settings::fs2mu2e) continue;
      temp_yield = ZXYields->GetM4lZX_Yields(_expected_yield_SR, M4l_down, M4l_up, i_fs, Settings::inclusive);
      yields_ZX.push_back(temp_yield);
   }


   cout << "\\begin{tabular}{l|c|c|c|c}" << endl;
   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "\\textbf{Channel} & $4\\Pe$ & $4\\Pgm$ & $2\\Pe2\\Pgm$ & $4\\ell$ \\\\" << endl;
   cout << "\\hline" << endl;


   cout << "\\qqZZ ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yqqZZ].at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\ggZZ ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yggZZ].at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\cPZ\\ + X ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_ZX.at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;

   cout << "Sum of backgrounds ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yqqZZ].at(i_fs) + yields_map[Settings::yggZZ].at(i_fs)+ yields_ZX.at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;

   cout << "Signal ($\\mH=125~\\GeV$) ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yH125].at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;

   cout << "Total expected ";
   for (int i_fs = 0; i_fs < num_of_final_states - 1; i_fs++)
   {
      cout << "& $" << yields_map[Settings::yH125].at(i_fs) + yields_map[Settings::yqqZZ].at(i_fs) + yields_map[Settings::yggZZ].at(i_fs)+ yields_ZX.at(i_fs) << "^{+ y.y}_{- z.z}$ ";
   }
   cout << "\\\\" << endl;

   cout << "\\hline" << endl;
   cout << "Observed & XXX & XXX & XXX & XXX \\\\" << endl;

   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "\\end{tabular}" << endl;


   cout << endl << endl;
   cout << "//==============" << endl;
   cout << "// T a b l e  3 " << endl;
   cout << "//==============" << endl;
   cout << endl;

//==============
// T a b l e  3
//==============

   yields_map.clear();
   yields_ZX.clear();
   double total;

   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      for ( int i_proc = 0; i_proc < num_of_processes_yields; i_proc++ )
      {
         if (i_proc == Settings::yH125ggH || i_proc == Settings::yH125VBF || i_proc == Settings::yH125WHlep || i_proc == Settings::yH125WHhad || i_proc == Settings::yH125ZHlep || i_proc == Settings::yH125ZHhad
         || i_proc == Settings::yH125ttHlep || i_proc == Settings::yH125ttHhad || i_proc == Settings::yH125tqH || i_proc == Settings::yH125bbH || i_proc == Settings::yqqZZ || i_proc == Settings::yggZZ)
         {
            temp_yield = histos_1D[Settings::M4lYields][Settings::fs4l][i_cat][i_proc]->Integral(
                         histos_1D[Settings::M4lYields][Settings::fs4l][i_cat][i_proc]->FindBin(M4l_down),
                         histos_1D[Settings::M4lYields][Settings::fs4l][i_cat][i_proc]->FindBin(M4l_up) - 1);

            yields_map[i_cat].push_back(temp_yield);
         }
      }
   }

   // Z+X
   for ( int i_cat = 0; i_cat < num_of_categories; i_cat++ )
   {
      temp_yield = ZXYields->GetM4lZX_Yields(_expected_yield_SR, M4l_down, M4l_up, Settings::fs4l, i_cat);
      yields_ZX.push_back(temp_yield);
   }


   cout << "\\begin{tabular}{c|cccccccccc|ccc|c|c}" << endl;
   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "Event & \\multicolumn{10}{c|}{Signal} & \\multicolumn{3}{c|}{Background} & Total & Observed \\\\" << endl;
   cout << "category & $\\ggH$ & VBF & WH-lep & WH-had & ZH-lep & ZH-had & ttH-lep & ttH-had & bbH & tqH & \\qqZZ & \\ggZZ & \\cPZ\\ + X & expected & \\\\" << endl;
   cout << "\\hline" << endl;


   cout << "Untagged ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::untagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::untagged].at(i_proc) << "$ ";
      total += yields_map[Settings::untagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::untagged) << "$ & $" << total + yields_ZX.at(Settings::untagged) << "$ & XXX \\\\" << endl;


   cout << "VBF-1j ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::VBF_1j_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::VBF_1j_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::VBF_1j_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::VBF_1j_tagged) << "$ & $" << total + yields_ZX.at(Settings::VBF_1j_tagged) << "$ & XXX \\\\" << endl;


   cout << "VBF-2j ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::VBF_2j_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::VBF_2j_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::VBF_2j_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::VBF_2j_tagged) << "$ & $" << total + yields_ZX.at(Settings::VBF_1j_tagged) << "$ & XXX \\\\" << endl;


   cout << "VH-lept. ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::VH_lepton_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::VH_lepton_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::VH_lepton_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::VH_lepton_tagged) << "$ & $" << total + yields_ZX.at(Settings::VH_lepton_tagged) << "$ & XXX \\\\" << endl;


   cout << "VH-hadr. ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::VH_hadron_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::VH_hadron_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::VH_hadron_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::VH_hadron_tagged) << "$ & $" << total + yields_ZX.at(Settings::VH_hadron_tagged) << "$ & XXX \\\\" << endl;


   cout << "ttH-lept. ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::ttH_lepton_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::ttH_lepton_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::ttH_lepton_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::ttH_lepton_tagged) << "$ & $" << total + yields_ZX.at(Settings::ttH_lepton_tagged) << "$ & XXX \\\\" << endl;

  cout << "ttH-hadr. ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::ttH_hadron_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::ttH_hadron_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::ttH_hadron_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::ttH_hadron_tagged) << "$ & $" << total + yields_ZX.at(Settings::ttH_hadron_tagged) << "$ & XXX \\\\" << endl;


   cout << "VH-MET ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::VH_MET_tagged].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::VH_MET_tagged].at(i_proc) << "$ ";
      total += yields_map[Settings::VH_MET_tagged].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::VH_MET_tagged) << "$ & $" << total + yields_ZX.at(Settings::VH_MET_tagged) << "$ & XXX \\\\" << endl;

   cout << "\\hline" << endl;


   cout << "Total ";
   total = 0;
   for (unsigned int i_proc = 0; i_proc < yields_map[Settings::inclusive].size(); i_proc++ )
   {
      cout << "& $" << yields_map[Settings::inclusive].at(i_proc) << "$ ";
      total += yields_map[Settings::inclusive].at(i_proc);
   }
   cout << "& $" << yields_ZX.at(Settings::inclusive) << "$ & $" << total + yields_ZX.at(Settings::inclusive) << "$ & XXX \\\\" << endl;

   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "\\end{tabular}" << endl;

   delete ZXYields;
}
//========================================================



//=======================================
void Histograms::SavePlots( TCanvas *c, TString name, TString folder)
{
   system("mkdir -p " + folder);
   c->SaveAs(name + ".pdf");
   c->SaveAs(name + ".root");
   c->SaveAs(name + ".eps");
   c->SaveAs(name + ".png");
   //gSystem->Exec("convert -density 150 -quality 100 " + name + ".eps " + name + ".png");
}
//=======================================



//==================================================
int Histograms::SetPlotName( TString variable_name )
{
   //=============
   // M4l
   //=============
   if ( variable_name == "M4lMain" )              return Settings::M4lMain;
   else if ( variable_name == "M4lMainZoomed" )   return Settings::M4lMainZoomed;
   else if ( variable_name == "M4lMainHighMass" )   return Settings::M4lMainHighMass;

   //=============
   // MZ1
   //=============
   else if ( variable_name == "MZ1" )             return Settings::MZ1;
   else if ( variable_name == "MZ1_M4L118130" )   return Settings::MZ1_M4L118130;

   //=============
   // MZ2
   //=============
   else if ( variable_name == "MZ2" )             return Settings::MZ2;
   else if ( variable_name == "MZ2_M4L118130" )   return Settings::MZ2_M4L118130;

   //=============
   // KD
   //=============
   else if ( variable_name == "KD" )                return Settings::KD;
   else if ( variable_name == "KD_M4L118130" )      return Settings::KD_M4L118130;
   else if ( variable_name == "DVBFDEC" )           return Settings::DVBFDEC;
   else if ( variable_name == "DVBFDEC_M4L118130" ) return Settings::DVBFDEC_M4L118130;
   else if ( variable_name == "DVHDEC" )            return Settings::DVHDEC;
   else if ( variable_name == "DVHDEC_M4L118130" )  return Settings::DVHDEC_M4L118130;
   else if ( variable_name == "D1jet" )             return Settings::D1jet;
   else if ( variable_name == "D1jet_M4L118130" )   return Settings::D1jet_M4L118130;
   else if ( variable_name == "D2jet" )             return Settings::D2jet;
   else if ( variable_name == "D2jet_M4L118130" )   return Settings::D2jet_M4L118130;
   else if ( variable_name == "DWH" )               return Settings::DWH;
   else if ( variable_name == "DWH_M4L118130" )     return Settings::DWH_M4L118130;
   else if ( variable_name == "DZH" )               return Settings::DZH;
   else if ( variable_name == "DZH_M4L118130" )     return Settings::DZH_M4L118130;
   else if ( variable_name == "DVH" )               return Settings::DVH;
   else if ( variable_name == "DVH_M4L118130" )     return Settings::DVH_M4L118130;

   //=============
   // MZ1vsMZ2
   //=============
   else if ( variable_name == "MZ1vsMZ2" )          return Settings::MZ1vsMZ2;
   else if ( variable_name == "MZ1vsMZ2_M4L118130" )return Settings::MZ1vsMZ2_M4L118130;

   //=============
   // KDvsM4l
   //=============
   else if ( variable_name == "KDvsM4l" )           return Settings::KDvsM4l;
   else if ( variable_name == "KDvsM4lZoomed" )     return Settings::KDvsM4lZoomed;
   else if ( variable_name == "KDvsM4lHighMass" )   return Settings::KDvsM4lHighMass;
   else if ( variable_name == "DVBFDECvsM4lZoomed" )return Settings::DVBFDECvsM4lZoomed;
   else if ( variable_name == "DVHDECvsM4lZoomed" ) return Settings::DVHDECvsM4lZoomed;
   else if ( variable_name == "D1jetvsM4lZoomed" )  return Settings::D1jetvsM4lZoomed;
   else if ( variable_name == "D2jetvsM4lZoomed" )  return Settings::D2jetvsM4lZoomed;
   else if ( variable_name == "DWHvsM4lZoomed" )    return Settings::DWHvsM4lZoomed;
   else if ( variable_name == "DZHvsM4lZoomed" )    return Settings::DZHvsM4lZoomed;
   else if ( variable_name == "DVHvsM4lZoomed" )    return Settings::DVHvsM4lZoomed;


  //=============
   // Others
   //=============
   else if ( variable_name == "PFMET" )              return Settings::PFMET;
   else if ( variable_name == "Pt4l" )               return Settings::Pt4l;
   else if ( variable_name == "Eta4l" )              return Settings::Eta4l;
   else if ( variable_name == "Pt_leading" )         return Settings::Pt_leading;
   else if ( variable_name == "Pt_trailing" )        return Settings::Pt_trailing;
   else if ( variable_name == "Eta_leading" )         return Settings::Eta_leading;
   else if ( variable_name == "Eta_trailing" )        return Settings::Eta_trailing;
   else if ( variable_name == "SIP_leading" )        return Settings::SIP_leading;
   else if ( variable_name == "SIP_trailing" )       return Settings::SIP_trailing;
   else if ( variable_name == "ISO_leading" )        return Settings::ISO_leading;
   else if ( variable_name == "ISO_trailing" )       return Settings::ISO_trailing;
   else if ( variable_name == "NExtraLep" )          return Settings::NExtraLep;
  else if ( variable_name == "NJets" )              return Settings::NJets;
   else if ( variable_name == "NJetsBTagged" )       return Settings::NJetsBTagged;
   else if ( variable_name == "M4l_110150_HighKD" )  return Settings::M4l_110150_HighKD;

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      //abort;
      return Settings::M4lMain;
   }
}
//==================================================



//===================================================
bool Histograms::GetVarLogX ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4lMain")                return bool(Variables::M4lMain().var_log_x);
   else if (variable_name == "M4lMainZoomed")    return bool(Variables::M4lMainZoomed().var_log_x);
   else if (variable_name == "M4lMainHighMass")  return bool(Variables::M4lMainHighMass().var_log_x);

   //=============
   // MZ1
   //=============
   else if (variable_name == "MZ1")              return bool(Variables::MZ1().var_log_x);
   else if (variable_name == "MZ1_M4L118130")    return bool(Variables::MZ1_M4L118130().var_log_x);

   //=============
   // MZ2
   //=============
   else if (variable_name == "MZ2")              return bool(Variables::MZ2().var_log_x);
   else if (variable_name == "MZ2_M4L118130")    return bool(Variables::MZ2_M4L118130().var_log_x);

   //=============
   // KD
   //=============
   else if (variable_name == "KD")               return bool(Variables::KD().var_log_x);
   else if (variable_name == "KD_M4L118130")     return bool(Variables::KD_M4L118130().var_log_x);
   else if (variable_name == "DVBFDEC")          return bool(Variables::DVBFDEC().var_log_x);
   else if (variable_name == "DVBFDEC_M4L118130")return bool(Variables::DVBFDEC_M4L118130().var_log_x);
   else if (variable_name == "DVHDEC")           return bool(Variables::DVHDEC().var_log_x);
   else if (variable_name == "DVHDEC_M4L118130") return bool(Variables::DVHDEC_M4L118130().var_log_x);
   else if (variable_name == "D1jet")            return bool(Variables::D1jet().var_log_x);
   else if (variable_name == "D1jet_M4L118130")  return bool(Variables::D1jet_M4L118130().var_log_x);
   else if (variable_name == "D2jet")            return bool(Variables::D2jet().var_log_x);
   else if (variable_name == "D2jet_M4L118130")  return bool(Variables::D2jet_M4L118130().var_log_x);
   else if (variable_name == "DWH")              return bool(Variables::DWH().var_log_x);
   else if (variable_name == "DWH_M4L118130")    return bool(Variables::DWH_M4L118130().var_log_x);
   else if (variable_name == "DZH")              return bool(Variables::DZH().var_log_x);
   else if (variable_name == "DZH_M4L118130")    return bool(Variables::DZH_M4L118130().var_log_x);
   else if (variable_name == "DVH")              return bool(Variables::DVH().var_log_x);
   else if (variable_name == "DVH_M4L118130")    return bool(Variables::DVH_M4L118130().var_log_x);

   //=============
   // MZ1vsMZ2
   //=============
   else if (variable_name == "MZ1vsMZ2")         return bool(Variables::MZ1vsMZ2().var_log_x);
   else if (variable_name == "MZ1vsMZ2_M4L118130")return bool(Variables::MZ1vsMZ2_M4L118130().var_log_x);

   //=============
   // KDvsM4l
   //=============
   else if (variable_name == "KDvsM4l")           return bool(Variables::KDvsM4l().var_log_x);
   else if (variable_name == "KDvsM4lZoomed")     return bool(Variables::KDvsM4lZoomed().var_log_x);
   else if (variable_name == "KDvsM4lHighMass")   return bool(Variables::KDvsM4lHighMass().var_log_x);
   else if (variable_name == "DVBFDECvsM4lZoomed")return bool(Variables::DVBFDECvsM4lZoomed().var_log_x);
   else if (variable_name == "DVHDECvsM4lZoomed") return bool(Variables::DVHDECvsM4lZoomed().var_log_x);
   else if (variable_name == "D1jetvsM4lZoomed")  return bool(Variables::D1jetvsM4lZoomed().var_log_x);
   else if (variable_name == "D2jetvsM4lZoomed")  return bool(Variables::D2jetvsM4lZoomed().var_log_x);
   else if (variable_name == "DWHvsM4lZoomed")    return bool(Variables::DWHvsM4lZoomed().var_log_x);
   else if (variable_name == "DZHvsM4lZoomed")    return bool(Variables::DZHvsM4lZoomed().var_log_x);
   else if (variable_name == "DVHvsM4lZoomed")    return bool(Variables::DVHvsM4lZoomed().var_log_x);

  //=============
   // Others
   //=============
   else if ( variable_name == "PFMET" )              return bool(Variables::PFMET().var_log_x);
   else if ( variable_name == "Pt4l" )               return bool(Variables::Pt4l().var_log_x);
   else if ( variable_name == "Eta4l" )              return bool(Variables::Eta4l().var_log_x);
   else if ( variable_name == "Pt_leading" )         return bool(Variables::Pt_leading().var_log_x);
   else if ( variable_name == "Pt_trailing" )        return bool(Variables::Pt_trailing().var_log_x);
   else if ( variable_name == "Eta_leading" )         return bool(Variables::Eta_leading().var_log_x);
   else if ( variable_name == "Eta_trailing" )        return bool(Variables::Eta_trailing().var_log_x);
   else if ( variable_name == "SIP_leading" )        return bool(Variables::SIP_leading().var_log_x);
   else if ( variable_name == "SIP_trailing" )       return bool(Variables::SIP_trailing().var_log_x);
   else if ( variable_name == "ISO_leading" )        return bool(Variables::ISO_leading().var_log_x);
   else if ( variable_name == "ISO_trailing" )       return bool(Variables::ISO_trailing().var_log_x);
   else if ( variable_name == "NExtraLep" )          return bool(Variables::NExtraLep().var_log_x);
  else if ( variable_name == "NJets" )              return bool(Variables::NJets().var_log_x);
   else if ( variable_name == "NJetsBTagged" )       return bool(Variables::NJetsBTagged().var_log_x);
   else if ( variable_name == "M4l_110150_HighKD" )  return bool(Variables::M4l_110150_HighKD().var_log_x);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
//      abort;
      return bool(Variables::M4lMain().var_log_x);
   }
}
//===================================================



//===================================================
bool Histograms::GetVarLogY ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4lMain")                return bool(Variables::M4lMain().var_log_y);
   else if (variable_name == "M4lMainZoomed")    return bool(Variables::M4lMainZoomed().var_log_y);
   else if (variable_name == "M4lMainHighMass")  return bool(Variables::M4lMainHighMass().var_log_y);

   //=============
   // MZ1
   //=============
   else if (variable_name == "MZ1")              return bool(Variables::MZ1().var_log_y);
   else if (variable_name == "MZ1_M4L118130")    return bool(Variables::MZ1_M4L118130().var_log_y);

   //=============
   // MZ2
   //=============
   else if (variable_name == "MZ2")              return bool(Variables::MZ2().var_log_y);
   else if (variable_name == "MZ2_M4L118130")    return bool(Variables::MZ2_M4L118130().var_log_y);

   //=============
   // KD
   //=============
   else if (variable_name == "KD")               return bool(Variables::KD().var_log_y);
   else if (variable_name == "KD_M4L118130")     return bool(Variables::KD_M4L118130().var_log_y);
   else if (variable_name == "DVBFDEC")          return bool(Variables::DVBFDEC().var_log_y);
   else if (variable_name == "DVBFDEC_M4L118130")return bool(Variables::DVBFDEC_M4L118130().var_log_y);
   else if (variable_name == "DVHDEC")           return bool(Variables::DVHDEC().var_log_y);
   else if (variable_name == "DVHDEC_M4L118130") return bool(Variables::DVHDEC_M4L118130().var_log_y);
   else if (variable_name == "D1jet")            return bool(Variables::D1jet().var_log_y);
   else if (variable_name == "D1jet_M4L118130")  return bool(Variables::D1jet_M4L118130().var_log_y);
   else if (variable_name == "D2jet")            return bool(Variables::D2jet().var_log_y);
   else if (variable_name == "D2jet_M4L118130")  return bool(Variables::D2jet_M4L118130().var_log_y);
   else if (variable_name == "DWH")              return bool(Variables::DWH().var_log_y);
   else if (variable_name == "DWH_M4L118130")    return bool(Variables::DWH_M4L118130().var_log_y);
   else if (variable_name == "DZH")              return bool(Variables::DZH().var_log_y);
   else if (variable_name == "DZH_M4L118130")    return bool(Variables::DZH_M4L118130().var_log_y);
   else if (variable_name == "DVH")              return bool(Variables::DVH().var_log_y);
   else if (variable_name == "DVH_M4L118130")    return bool(Variables::DVH_M4L118130().var_log_y);

   //=============
   // MZ1vsMZ2
   //=============
   else if (variable_name == "MZ1vsMZ2")         return bool(Variables::MZ1vsMZ2().var_log_y);
   else if (variable_name == "MZ1vsMZ2_M4L118130")return bool(Variables::MZ1vsMZ2_M4L118130().var_log_y);

   //=============
   // KDvsM4l
   //=============
   else if (variable_name == "KDvsM4l")           return bool(Variables::KDvsM4l().var_log_y);
   else if (variable_name == "KDvsM4lZoomed")     return bool(Variables::KDvsM4lZoomed().var_log_y);
   else if (variable_name == "KDvsM4lHighMass")   return bool(Variables::KDvsM4lHighMass().var_log_y);
   else if (variable_name == "DVBFDECvsM4lZoomed")return bool(Variables::DVBFDECvsM4lZoomed().var_log_y);
   else if (variable_name == "DVHDECvsM4lZoomed") return bool(Variables::DVHDECvsM4lZoomed().var_log_y);
   else if (variable_name == "D1jetvsM4lZoomed")  return bool(Variables::D1jetvsM4lZoomed().var_log_y);
   else if (variable_name == "D2jetvsM4lZoomed")  return bool(Variables::D2jetvsM4lZoomed().var_log_y);
   else if (variable_name == "DWHvsM4lZoomed")    return bool(Variables::DWHvsM4lZoomed().var_log_y);
   else if (variable_name == "DZHvsM4lZoomed")    return bool(Variables::DZHvsM4lZoomed().var_log_y);
   else if (variable_name == "DVHvsM4lZoomed")    return bool(Variables::DVHvsM4lZoomed().var_log_y);

  //=============
   // Others
   //=============
   else if ( variable_name == "PFMET" )              return bool(Variables::PFMET().var_log_y);
   else if ( variable_name == "Pt4l" )               return bool(Variables::Pt4l().var_log_y);
   else if ( variable_name == "Eta4l" )              return bool(Variables::Eta4l().var_log_y);
   else if ( variable_name == "Pt_leading" )         return bool(Variables::Pt_leading().var_log_y);
   else if ( variable_name == "Pt_trailing" )        return bool(Variables::Pt_trailing().var_log_y);
   else if ( variable_name == "Eta_leading" )         return bool(Variables::Eta_leading().var_log_y);
   else if ( variable_name == "Eta_trailing" )        return bool(Variables::Eta_trailing().var_log_y);
   else if ( variable_name == "SIP_leading" )        return bool(Variables::SIP_leading().var_log_y);
   else if ( variable_name == "SIP_trailing" )       return bool(Variables::SIP_trailing().var_log_y);
   else if ( variable_name == "ISO_leading" )        return bool(Variables::ISO_leading().var_log_y);
   else if ( variable_name == "ISO_trailing" )       return bool(Variables::ISO_trailing().var_log_y);
   else if ( variable_name == "NExtraLep" )          return bool(Variables::NExtraLep().var_log_y);
   else if ( variable_name == "NJets" )              return bool(Variables::NJets().var_log_y);
   else if ( variable_name == "NJetsBTagged" )       return bool(Variables::NJetsBTagged().var_log_y);
   else if ( variable_name == "M4l_110150_HighKD" )  return bool(Variables::M4l_110150_HighKD().var_log_y);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
//      abort;
      return bool(Variables::M4lMain().var_log_y);
   }
}
//===================================================



//============================================================================================================
TLegend* Histograms::CreateLegend( string position, TH1F *data, TH1F *h125, TH1F *qqZZ, TH1F *ggZZ, TH1F *ZX, TH1F *VVV )
{
   TLegend *leg;

   if ( position == "right" )
   {
      leg = new TLegend(.71, .71, .91, .91);
   }
   else
   {
      leg = new TLegend(.21, .71, .41, .91);
   }

   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);

   leg->AddEntry( data, "Data", "p E" );
   leg->AddEntry( h125, "H(125)","f");
   leg->AddEntry( qqZZ, "q#bar{q}#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ggZZ, "gg#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ZX,   "Z+X", "f" );
   leg->AddEntry( VVV,   "VVV, tt + (V)V", "f" );

   return leg;
}
//============================================================================================================



//==============================================================================================================================
TLegend* Histograms::CreateLegendVBF( string position, TH1F *data, TH1F *h125VBF, TH1F *h125_other, TH1F *qqZZ, TH1F *ggZZ, TH1F *ZX , TH1F* VVV, bool mask)
{
   TLegend *leg;

   if ( position == "right" )
   {
      leg = new TLegend(.71, .67, .91, .91);
   }
   else
   {
      leg = new TLegend(.21, .67, .41, .91);
   }

   if(!mask)
   {
      leg->SetFillColor(0);
    leg->SetFillStyle(0);
  }
  leg->SetBorderSize(0);

   leg->AddEntry( data, "Data", "p E" );
   leg->AddEntry( h125VBF,"H(125), VBF","f");
   leg->AddEntry( h125_other,"H(125), other","f");
   leg->AddEntry( qqZZ, "q#bar{q}#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ggZZ, "gg#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ZX, "Z+X", "f" );
   leg->AddEntry( VVV,   "VVV, tt + (V)V", "f" );

   return leg;
}
//==============================================================================================================================



//==============================================================================================================================
TLegend* Histograms::CreateLegendVH( string position, TH1F *data, TH1F *h125VH, TH1F *h125_other, TH1F *qqZZ, TH1F *ggZZ, TH1F *ZX, TH1F* VVV, bool mask )
{
   TLegend *leg;

   if ( position == "right" )
   {
      leg = new TLegend(.71, .67, .91, .91);
   }
   else
   {
      leg = new TLegend(.21, .67, .41, .91);
   }

   if(!mask)
   {
      leg->SetFillColor(0);
    leg->SetFillStyle(0);
  }
  leg->SetBorderSize(0);

   leg->AddEntry( data, "Data", "p E" );
   leg->AddEntry( h125VH,"H(125), VH","f");
   leg->AddEntry( h125_other,"H(125), other","f");
   leg->AddEntry( qqZZ, "q#bar{q}#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ggZZ, "gg#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ZX, "Z+X", "f" );
   leg->AddEntry( VVV,   "VVV, tt + (V)V", "f" );

   return leg;
}
//==============================================================================================================================



//==============================================================================================================================
TLegend* Histograms::CreateLegendttH( string position, TH1F *data, TH1F *h125ttH, TH1F *h125_other, TH1F *qqZZ, TH1F *ggZZ, TH1F *ZX, TH1F* VVV )
{
   TLegend *leg;

   if ( position == "right" )
   {
      leg = new TLegend(.71, .67, .91, .91);
   }
   else
   {
      leg = new TLegend(.21, .67, .41, .91);
   }

   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);

   leg->AddEntry( data, "Data", "p E" );
   leg->AddEntry( h125ttH,"H(125), ttH","f");
   leg->AddEntry( h125_other,"H(125), other","f");
   leg->AddEntry( qqZZ, "q#bar{q}#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ggZZ, "gg#rightarrowZZ, Z#gamma*", "f" );
   leg->AddEntry( ZX, "Z+X", "f" );
   leg->AddEntry( VVV,   "VVV, tt + (V)V", "f" );
   return leg;
}
//==============================================================================================================================


//==============================================================================================================================
TLegend* Histograms::CreateLegendSTXS( string position, TH1F *data, TH1F *h125ggH, TH1F *h125bbH, TH1F *h125VBF, TH1F *h125VH, TH1F *h125ttH, TH1F *h125tqH, TH1F *qqZZ, TH1F *ggZZ, TH1F *ZX , TH1F *VVV , bool mask)
{
    TLegend *leg;

    if ( position == "right" )
    {
        leg = new TLegend(.61, .77, .97, .91);
    }
    else
    {
        leg = new TLegend(.21, .67, .41, .91);
    }

    if(!mask)
    {
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
    }
    leg->SetBorderSize(0);
    leg->SetNColumns(2);

    leg->AddEntry( data, "Data", "p E" );
    leg->AddEntry( h125ttH,"H(125), ttH","f");
    leg->AddEntry( h125ggH,"H(125), ggH","f");
    leg->AddEntry( h125tqH,"H(125), tH","f");
    leg->AddEntry( h125bbH,"H(125), bbH","f");
    leg->AddEntry( qqZZ, "q#bar{q}#rightarrowZZ, Z#gamma*", "f" );
    leg->AddEntry( h125VBF,"H(125), VBF","f");
    leg->AddEntry( ggZZ, "gg#rightarrowZZ, Z#gamma*", "f" );
    leg->AddEntry( h125VH,"H(125), VH","f");
    leg->AddEntry( ZX, "Z+X", "f" );
    leg->AddEntry( VVV, "EW bkg", "f" );

    return leg;
}
//==============================================================================================================================


//===========================================================================================
TLegend* Histograms::Create2DLegend( string position, TH2F *fs4e, TH2F *fs4mu, TH2F *fs2e2mu)
{
   TLegend *leg;
   leg = new TLegend(0.18, 0.80, 0.27, 0.92); // Default initalisiation to avoid scram compilation error
   if ( position == "left" ) leg = new TLegend(0.18, 0.80, 0.27, 0.92);

   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.032);

   leg->AddEntry(fs4e,   "4e",    "p");
   leg->AddEntry(fs4mu,  "4#mu",   "p");
   leg->AddEntry(fs2e2mu,"2e2#mu", "p");

//   leg->Draw();

   return leg;
}
//===========================================================================================



//=======================================================================================================================================================================
TLegend* Histograms::Create2DLegendAllCat_KD( string position, TGraphErrors *fs4e, TGraphErrors *fs4mu, TGraphErrors *fs2e2mu, TGraphErrors *untagged, TGraphErrors *VBF1jet, TGraphErrors *VHlep, TGraphErrors *ttHhad, TGraphErrors *ttHlep)
{
   TLegend *leg;
   leg = new TLegend(0.00, 0.00, 1.00, 1.00);

   leg->SetFillStyle(0);
   leg->SetBorderSize(1);
   leg->SetTextFont(42);
   leg->SetTextSize(0.18);
   leg->SetNColumns(3);

   leg->AddEntry(fs4e,    "4e",              "lp");
   leg->AddEntry(fs4mu,   "4#mu",            "lp");
   leg->AddEntry(fs2e2mu, "2e2#mu",          "lp");
   leg->AddEntry(untagged, "untagged",        "lp");
   leg->AddEntry(VHlep,    "VH-lept. tagged", "lp");
   leg->AddEntry(ttHhad,   "t#bar{t}H-hadr. tagged","lp");
   leg->AddEntry(VBF1jet, "VBF-1jet tagged",   "lp");
   leg->AddEntry(ttHlep,  "t#bar{t}H-lept. tagged", "lp");
   leg->AddEntry((TObject*)0, "",                       "");


   return leg;
}
//=======================================================================================================================================================================

//=======================================================================================================================================================================
TLegend* Histograms::Create2DLegendAllCat_DVBFDEC( string position, TGraphErrors *fs4e, TGraphErrors *fs4mu, TGraphErrors *fs2e2mu, TGraphErrors *VBF2jet)
{
   TLegend *leg;
   leg = new TLegend(0.00, 0.00, 1.00, 1.00);

   leg->SetFillStyle(0);
   leg->SetBorderSize(1);
   leg->SetTextFont(42);
   leg->SetTextSize(0.20);
   leg->SetNColumns(3);

   leg->AddEntry(fs4e,    "4e                   ",  "lp");
   leg->AddEntry(fs4mu,   "4#mu                 ",  "lp");
   leg->AddEntry(fs2e2mu, "2e2#mu        ",  "lp");
   leg->AddEntry((TObject*)0, "",    "");
   leg->AddEntry(VBF2jet,   "VBF-2jet tagged", "lp");
   leg->AddEntry((TObject*)0, "",    "");



   return leg;
}
//=======================================================================================================================================================================

//=======================================================================================================================================================================
TLegend* Histograms::Create2DLegendAllCat_DVHDEC( string position, TGraphErrors *fs4e, TGraphErrors *fs4mu, TGraphErrors *fs2e2mu, TGraphErrors *VHhad)
{
   TLegend *leg;
   leg = new TLegend(0.00, 0.00, 1.00, 1.00);

   leg->SetFillStyle(0);
   leg->SetBorderSize(1);
   leg->SetTextFont(42);
   leg->SetTextSize(0.20);
   leg->SetNColumns(3);

   leg->AddEntry(fs4e,    "4e                   ",  "lp");
   leg->AddEntry(fs4mu,   "4#mu                 ",  "lp");
   leg->AddEntry(fs2e2mu, "2e2#mu     ",  "lp");
   leg->AddEntry((TObject*)0, "",    "");
   leg->AddEntry(VHhad,   "VH-hadr. tagged", "lp");
   leg->AddEntry((TObject*)0, "",     "");


   return leg;
}
//=======================================================================================================================================================================



//=======================================================================================================================
TLegend* Histograms::Create2DErrorLegend( string position, TGraphErrors *fs4e, TGraphErrors *fs4mu,TGraphErrors *fs2e2mu)
{
   TLegend *leg;
   leg = new TLegend(0.855, 0.855, 0.955, 0.955);// Default initalisiation to avoid scram compilation error
   if(position == "right") leg = new TLegend(0.855, 0.855, 0.955, 0.955);

   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextFont(42);
   leg->SetTextSize(0.032);

   leg->AddEntry(fs4e,    "4e",    "lp");
   leg->AddEntry(fs4mu,   "4#mu",   "lp");
   leg->AddEntry(fs2e2mu, "2e2#mu", "lp");

   return leg;
}
//=======================================================================================================================



//========================================================================
TPaveText* Histograms::CreateCutText( string position, TString cut_label )
{
   TPaveText *pav;
   pav = new TPaveText(.21, .81, .51, .91 ,"brNDC");// Default initalisiation to avoid scram compilation error

   if ( position == "left top" )    pav = new TPaveText(.21, .81, .51, .91 ,"brNDC");
   if ( position == "right top" )   pav = new TPaveText(.61, .81, .91, .91 ,"brNDC");
   if ( position == "middle top" )  pav = new TPaveText(.41, .81, .71, .91 ,"brNDC");
   if ( position == "left under 2D legend" ) pav = new TPaveText(.18, .72, .31, .78 ,"brNDC");

   pav->SetFillStyle(0);
   pav->SetBorderSize(0);
   pav->SetTextAlign(11);
   pav->SetTextSize(0.032);
   pav->SetTextFont(42);
   pav->AddText(cut_label);

   return pav;
}
//========================================================================



//=======================================================================================================================
TPaveText* Histograms::CreateCatText( string position, TString cat_label)
{
   TPaveText *pav;
   pav = new TPaveText(.21, .81, .51, .91,"brNDC");// Default initalisiation to avoid scram compilation error
   if (position == "top left")  pav = new TPaveText(.21, .81, .51, .91,"brNDC");
   if (position == "left under cut text" )    pav = new TPaveText(.21, .76, .51, .86 ,"brNDC");
   pav->SetFillStyle(0);
   pav->SetBorderSize(0);
   pav->SetTextAlign(11);
   pav->SetTextSize(0.032);
   pav->SetTextFont(42);
   pav->AddText(cat_label);

   return pav;
}
//=======================================================================================================================



//=======================================================================================================================
TLine* Histograms::CreateDashedLine(float x1, float y1, float x2, float y2)
{
   TLine *line;
   line = new TLine(x1, y1, x2, y2);

   line->SetLineStyle(9);
   line->SetLineWidth(2);
   line->SetLineColor(kBlack);

   return line;
}
//=======================================================================================================================



//======================================================
void Histograms::DrawLogX( TCanvas *c, int cat, int fs )
{
   int x_low = 100;
   int x_up  = 1100;
   int step  = 100;

   float label_margin = -1.8; //Change this if you want to move X-axis numbers up/dn

   TLatex *latex_80 = new TLatex(80, label_margin , "80");
   latex_80->SetTextAlign(23);
   latex_80->SetTextFont (42);
   latex_80->SetTextSize (0.04);
   latex_80->Draw();

   for ( int i = x_low; i < x_up; i += step )
   {
      if (i == 600 || i == 800 || i == 900) continue;
      float i_x = i;

      TLatex *latex = new TLatex(i, label_margin , Form("%.0f", i_x));

      latex->SetTextAlign(23);
      latex->SetTextFont (42);
      latex->SetTextSize (0.04);
      latex->Draw();
   }
}
//======================================================



//=======================================
void Histograms::MakeCOLZGrey(bool shift)
{
//   int col = 1;
   const Int_t NRGBs = 2;
   const Int_t NCont = 255;
   Double_t stops[NRGBs] = { 0.00, 1.00 };
   Double_t red[NRGBs]   = { 1., 0.5 };
   Double_t green[NRGBs] = { 1., 0.5 };
   Double_t blue[NRGBs]  = { 1., 0.5 };

   if(shift)
   {
      for(int i=0; i<NRGBs; i++){
         red  [i] -= 0.10;
         green[i] -= 0.10;
         blue [i] -= 0.10;
      }
   }
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}
//=======================================



//=======================================
float Histograms::SetMassPoint(int point)
{
   switch (point) {
      case 0:
         return 120.;
         break;

      case 1:
         return 124.;
         break;

      case 2:
         return 125;
         break;

      case 3:
         return 126;
         break;

      case 4:
         return 130;
         break;

      default:
         cout << "[ERROR] Mass point index out of range: " << point << endl;
         abort();
         break;
   }
}
//=======================================



//=======================================
int Histograms::SetProcess(int point, int production_mode)
{
   if (point == 0)
   {
      switch (production_mode) {
         case 0:
            return Settings::yH120ggH;
            break;

         case 1:
            return Settings::yH120VBF;
            break;

         case 2:
            return Settings::yH120WHlep;
            break;

         case 3:
            return Settings::yH120WHhad;
            break;

      case 4:
        return Settings::yH120ZHlep;
        break;

      case 5:
        return Settings::yH120ZHhad;
        break;

         case 6:
            return Settings::yH120ttHlep;
            break;

         case 7:
            return Settings::yH120ttHhad;
            break;

         case 8:
            return Settings::yH120bbH;
            break;

         case 9:
            return Settings::yH120tqH;
            break;

         default:
            cout << "[ERROR] Mass point index out of range: " << production_mode << endl;
            abort();
            break;
      }
   }

   else if (point == 1)
   {
      switch (production_mode) {
      case 0:
        return Settings::yH124ggH;
        break;

      case 1:
        return Settings::yH124VBF;
        break;

      case 2:
        return Settings::yH124WHlep;
        break;

      case 3:
        return Settings::yH124WHhad;
        break;

      case 4:
        return Settings::yH124ZHlep;
        break;

      case 5:
        return Settings::yH124ZHhad;
        break;

      case 6:
        return Settings::yH124ttHlep;
        break;

      case 7:
        return Settings::yH124ttHhad;
        break;

         case 8:
            return Settings::yH124bbH;
            break;

         case 9:
            return Settings::yH124tqH;
            break;

         default:
            cout << "[ERROR] Mass point index out of range: " << production_mode << endl;
            abort();
            break;

      }
   }

   else if (point == 2)
   {
      switch (production_mode) {
      case 0:
        return Settings::yH125ggH;
        break;

      case 1:
        return Settings::yH125VBF;
        break;

      case 2:
        return Settings::yH125WHlep;
        break;

      case 3:
        return Settings::yH125WHhad;
        break;

      case 4:
        return Settings::yH125ZHlep;
        break;

      case 5:
        return Settings::yH125ZHhad;
        break;

      case 6:
        return Settings::yH125ttHlep;
        break;

      case 7:
        return Settings::yH125ttHhad;
        break;

         case 8:
            return Settings::yH125bbH;
            break;

         case 9:
            return Settings::yH125tqH;
            break;

         default:
            cout << "[ERROR] Mass point index out of range: " << production_mode << endl;
            abort();
            break;

      }
   }

   else if (point == 3)
   {
      switch (production_mode) {
      case 0:
        return Settings::yH126ggH;
        break;

      case 1:
        return Settings::yH126VBF;
        break;

      case 2:
        return Settings::yH126WHlep;
        break;

      case 3:
        return Settings::yH126WHhad;
        break;

      case 4:
        return Settings::yH126ZHlep;
        break;

      case 5:
        return Settings::yH126ZHhad;
        break;

      case 6:
        return Settings::yH126ttHlep;
        break;

      case 7:
        return Settings::yH126ttHhad;
        break;

         case 8:
            return Settings::yH126bbH;
            break;

         case 9:
            return Settings::yH126tqH;
            break;

         default:
            cout << "[ERROR] Mass point index out of range: " << production_mode << endl;
            abort();
            break;

      }
   }

   else if (point == 4)
   {
      switch (production_mode) {
      case 0:
        return Settings::yH130ggH;
        break;

      case 1:
        return Settings::yH130VBF;
        break;

      case 2:
        return Settings::yH130WHlep;
        break;

      case 3:
        return Settings::yH130WHhad;
        break;

      case 4:
        return Settings::yH130ZHlep;
        break;

      case 5:
        return Settings::yH130ZHhad;
        break;

      case 6:
        return Settings::yH130ttHlep;
        break;

      case 7:
        return Settings::yH130ttHhad;
        break;

         case 8:
            return Settings::yH130bbH;
            break;

         case 9:
            return Settings::yH130tqH;
            break;

         default:
            cout << "[ERROR] Mass point index out of range: " << production_mode << endl;
            abort();
            break;
      }
   }
   else
   {
      cout << "[ERROR] Mass point index out of range: " << point << endl;
      abort();
   }
}
//=======================================



//=====================================
void Histograms::Rebin(THStack *stack)
{
   for ( int i = 0; i < stack->GetHists()->GetSize(); i++ )
   {
      TH1F *h = (TH1F*)stack->GetHists()->At(i);
      h->Rebin(2);
   }
}
//=====================================



//========================================
void Histograms::ChangeYaxisTitle(THStack *stack)
{
   stack->GetYaxis()->SetTitle("Events / 4 GeV");
}
//========================================

//========================================
void Histograms::Split_2e2mu()
{
   _merge_2e2mu = false;
}
//========================================
