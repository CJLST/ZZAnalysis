#ifndef Settings_h
#define Settings_h

using namespace std;

class Settings
{

public:
	
	Settings();
	~Settings();
   
   enum _blindings { fullyblind = 0, blindabove_110 = 1, blindbelow_150 = 2, blind_110_150 = 3, unblinded = 4 };
   enum _final_state { fs4e = 0, fs4mu = 1, fs2e2mu = 2, fs2mu2e = 3, fs4l = 4, MAX_NUM_OF_FINAL_STATES };
   
   enum _category
   {
      untagged          = 0,
      VBF_1j_tagged     = 1,
      VBF_2j_tagged     = 2,
      VH_lepton_tagged  = 3,
      VH_hadron_tagged  = 4,
      ttH_lepton_tagged = 5,
      ttH_hadron_tagged = 6,
      VH_MET_tagged     = 7,
      inclusive         = 8,
      MAX_NUM_OF_CATEGORIES
   };
   
   enum _1D_plot_name
   {
      M4lMain         = 0,
      M4lMainZoomed   = 1,
      M4lMainHighMass = 2,
      MZ1             = 3,
      MZ1_M4L118130   = 4,
      MZ2             = 5,
      MZ2_M4L118130   = 6,
      KD              = 7,
      KD_M4L118130    = 8,
      D1jet           = 9,
      D1jet_M4L118130 = 10,
      D2jet           = 11,
      D2jet_M4L118130 = 12,
      DWH             = 13,
      DWH_M4L118130   = 14,
      DZH             = 15,
      DZH_M4L118130   = 16,
      DVH             = 17,
      DVH_M4L118130   = 18,
      M4lYields       = 19,
      MAX_NUM_OF_1D_PLOT_NAMES
   };
   
   enum _2D_plot_name
   {
      MZ1vsMZ2           = 0,
      MZ1vsMZ2_M4L118130 = 1,
      MAX_NUM_OF_2D_PLOT_NAMES
   };
   
   enum _2D_error_plot_name 
   {
      KDvsM4l          = 0,
      KDvsM4lZoomed    = 1,
      KDvsM4lHighMass  = 2,
      D1jetvsM4lZoomed = 3,
      D2jetvsM4lZoomed = 4,
      DWHvsM4lZoomed   = 5,
      DZHvsM4lZoomed   = 6,
      DVHvsM4lZoomed   = 7,
      MAX_NUM_OF_2D_ERROR_PLOT_NAMES
   };
   
   enum _production_modes
	{
		ggH     = 0,
		VBF     = 1,
		WH_lep  = 2,
		WH_had  = 3,
		ZH_lep  = 4,
		ZH_had  = 5,
		ttH_lep = 6,
		ttH_had = 7,
		bbH     = 8,
		tqH     = 9,
		MAX_NUM_OF_PRODUCTION_MODES
		
	};
   
   enum _process
   { 
      Data    = 0,
      H125    = 1,
      H125ggH = 2,
      H125VBF = 3,
      H125VH  = 4,
      H125ttH = 5,
      H125bbH = 6,
      H125tqH = 7,
      qqZZ    = 8,
      ggZZ    = 9,
      DY      = 10,
      ttbar   = 11,
      MAX_NUM_OF_PROCESSES
   };
   
   enum _process_yields
   {
      yData      = 0,
      yH120      = 1,
      yH124      = 2,
      yH125      = 3,
      yH126      = 4,
      yH130      = 5,
      yH120ggH   = 6,
      yH124ggH   = 7,
      yH125ggH   = 8,
      yH126ggH   = 9,
      yH130ggH   = 10,
      yH120VBF   = 11,
      yH124VBF   = 12,
      yH125VBF   = 13,
      yH126VBF   = 14,
      yH130VBF   = 15,
      yH120WHlep = 16,
      yH124WHlep = 17,
      yH125WHlep = 18,
      yH126WHlep = 19,
      yH130WHlep = 20,
		yH120WHhad = 21,
		yH124WHhad = 22,
		yH125WHhad = 23,
		yH126WHhad = 24,
		yH130WHhad = 25,
      yH120ZHlep = 26,
      yH124ZHlep = 27,
      yH125ZHlep = 28,
      yH126ZHlep = 29,
      yH130ZHlep = 30,
		yH120ZHhad = 31,
		yH124ZHhad = 32,
		yH125ZHhad = 33,
		yH126ZHhad = 34,
		yH130ZHhad = 35,
      yH120ttHlep= 36,
      yH124ttHlep= 37,
      yH125ttHlep= 38,
      yH126ttHlep= 39,
      yH130ttHlep= 40,
		yH120ttHhad= 41,
      yH124ttHhad= 42,
      yH125ttHhad= 43,
      yH126ttHhad= 44,
      yH130ttHhad= 45,
		yH120bbH   = 46,
      yH124bbH   = 47,
      yH125bbH   = 48,
      yH126bbH   = 49,
      yH130bbH   = 50,
		yH120tqH   = 51,
      yH124tqH   = 52,
      yH125tqH   = 53,
      yH126tqH   = 54,
      yH130tqH   = 55,
      yqqZZ      = 56,
      yggZZ      = 57,
      yDY        = 58,
      yttbar     = 59,
      MAX_NUM_OF_PROCESSES_YIELDS
   };
   
   static const int num_of_production_modes    = MAX_NUM_OF_PRODUCTION_MODES;
   static const int num_of_processes           = MAX_NUM_OF_PROCESSES;
   static const int num_of_processes_yields    = MAX_NUM_OF_PROCESSES_YIELDS;
   static const int num_of_final_states        = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_categories          = MAX_NUM_OF_CATEGORIES;
   static const int num_of_1D_plot_names       = MAX_NUM_OF_1D_PLOT_NAMES;
   static const int num_of_2D_plot_names       = MAX_NUM_OF_2D_PLOT_NAMES;
   static const int num_of_2D_error_plot_names = MAX_NUM_OF_2D_ERROR_PLOT_NAMES;
   
   private:
      
};
#endif
