#ifndef Settings_h
#define Settings_h

using namespace std;

class Settings
{

public:
	
	Settings();
	~Settings();
   
   enum _process
   {
      Data = 0,
      WZ = 1,
      qqZZ = 2,
      DY = 3,
      ttbar = 4,
      Total = 5,
      MAX_NUM_OF_PROCESSES
   };
   
   enum _flavour
	{
		ele = 0,
		mu = 1,
		MAX_NUM_OF_FLAVOURS
		
	};
	
   enum _final_state
	{
		fs4mu = 0,
		fs4e = 1,
		fs2e2mu = 2,
		fs2mu2e = 3,
		fs4l = 4,
		MAX_NUM_OF_FINAL_STATES
	};

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
   
   enum _category_STXS
   {
      ggH_0J_PTH_0_10     = 0,
      ggH_0J_PTH_10_200   = 1,
      ggH_1J_PTH_0_60     = 2,
      ggH_1J_PTH_60_120   = 3,
      ggH_1J_PTH_120_200  = 4,
      ggH_2J_PTH_0_60     = 5,
      ggH_2J_PTH_60_120   = 6,
      ggH_2J_PTH_120_200  = 7,
      ggH_PTH_200         = 8,
      ggH_VBF             = 9,
      VBF_1j              = 10,
      VBF_2j              = 11,
      VBF_2j_mjj_350_700_2j = 12,
      VBF_2j_mjj_GT700_2j   = 13,
      VBF_2j_mjj_GT350_3j   = 14,
      VBF_GT200_2J          = 15,
      VH_Had                = 16,
      VBF_rest_VH           = 17,
      VH_lep_0_150          = 18,
      VH_Lep_GT150          = 19,
      ttH_Lep               = 20,
      ttH_Had               = 21,
      inclusive_stxs        = 22,
      MAX_NUM_OF_CATEGORIES_STXS
   };

   enum _eta_bins
	{
		EB = 0,
		EE = 1,
		MAX_NUM_OF_ETA_BINS
	};
   
   enum _regions_OS
	{
		reg2P2F = 0,
		reg3P1F = 1,
		regOS   = 2,
		MAX_NUM_OF_REGIONS_OS
	};
   
   enum _regions_SS
	{
		regZLL = 0,
		MAX_NUM_OF_REGIONS_SS
	};
   
   enum _fake_rates
	{
		corrected = 0 ,
		uncorrected = 1,
		MAX_NUM_OF_FAKE_RATES
		
	};
	
	enum _FR_variations
	{
		nominal = 0,
		Up      = 1,
		Dn      = 2,
		MAX_NUM_OF_FR_VARIATIONS
	};
	
	enum _Z_mass_windows
	{
		_40_MZ1_120 = 0 ,
		_MZ1mMZtrue_7 = 1,
		_60_MZ1_120 = 2,
		_MZ1EmMZtrue_5 = 3,
		MAX_NUM_OF_Z_MASS_WINDOWS
	};
   static const int num_of_processes         = MAX_NUM_OF_PROCESSES;
   static const int num_of_flavours          = MAX_NUM_OF_FLAVOURS;
   static const int num_of_final_states      = MAX_NUM_OF_FINAL_STATES;
   static const int num_of_categories        = MAX_NUM_OF_CATEGORIES;
   static const int num_of_categories_stxs   = MAX_NUM_OF_CATEGORIES_STXS;
   static const int num_of_eta_bins          = MAX_NUM_OF_ETA_BINS;
   static const int num_of_regions_os        = MAX_NUM_OF_REGIONS_OS;
   static const int num_of_regions_ss        = MAX_NUM_OF_REGIONS_SS;
   static const int num_of_fake_rates        = MAX_NUM_OF_FAKE_RATES;
   static const int num_of_fr_variations     = MAX_NUM_OF_FR_VARIATIONS;
	static const int num_of_z_mass_windows    = MAX_NUM_OF_Z_MASS_WINDOWS;

   private:
      
};
#endif
