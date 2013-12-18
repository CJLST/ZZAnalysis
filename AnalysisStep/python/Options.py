import FWCore.ParameterSet.VarParsing as VarParsing
import sys

# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#Passing_Command_Line_Arguments_T



############################################################################################
#  __________   _                _                       ___        _   _                  #
# |__  /__  /  / \   _ __   __ _| |_   _ _______ _ __   / _ \ _ __ | |_(_) ___  _ __  ___  #
#   / /  / /  / _ \ | '_ \ / _` | | | | |_  / _ \ '__| | | | | '_ \| __| |/ _ \| '_ \/ __| #
#  / /_ / /_ / ___ \| | | | (_| | | |_| |/ /  __/ |    | |_| | |_) | |_| | (_) | | | \__ \ #
# /____/____/_/   \_\_| |_|\__,_|_|\__, /___\___|_|     \___/| .__/ \__|_|\___/|_| |_|___/ #
#                                  |___/                     |_|                           #
############################################################################################

def getAnalysisOptions(options=None):   
    if options == None:
        #options = VarParsing.VarParsing()
        #options = VarParsing.VarParsing ('analysis')
        options = VarParsing.VarParsing ('standard')

    #  default, -1 means all events
    options.maxEvents = -1

    #  Setting isMC as option
    options.register( 'isMC',
		  0, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,            # string, int, or float
		  "Tells if the sample is MC or not (default = 0)." )

    options.register( 'tune',
		  '', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,         # string, int, or float
		  "Can be DoubleEle, DoubleMu, MuEG or MCAllEvents." )
		  

    options.register( 'elecorrtype',
		  '', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,         # string, int, or float
		  "For electron corrections. Possible for DATA ( 13JulReReco,06AugReReco, 24AugReReco, Prompt2012C, Jan16ReReco ), for MC (Summer12_53X, Fall11) " )

    options.register ('setup',
		  2012, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,            # string, int, or float
		  "Tells about the tags used for rho corrections and possibly other things." )                 

    options.register ('filePrepend',
		  '', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,         # string, int, or float
		  "Prefix added to output file names. Could be directory path." )
    options.register ('tag',
		  '', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,         # string, int, or float
		  "Tag to add to file names." )

    options.register ('goodlumi',
		  'NO_JSON_FILE', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "JSON  good lumisection file to be passed for data. No default." )
		  
    #Setting input tag as option
    options.register ('globaltag',
		  'NOT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Global tag from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions")

    # Protection in case sys.argv is missing due to various edm tools
    if not hasattr(sys, "argv"):
        return options

    # Hack to be able to pass multiple arguments from multicrab.cfg
    print sys.argv
    if len(sys.argv) > 0:
        last = sys.argv.pop()
        sys.argv.extend(last.split(":"))
        print sys.argv

    #options.parseArguments()

    return options



############################################################
#  _____                 ___        _   _                  #
# |_   _| __ ___  ___   / _ \ _ __ | |_(_) ___  _ __  ___  #
#   | || '__/ _ \/ _ \ | | | | '_ \| __| |/ _ \| '_ \/ __| #
#   | || | |  __/  __/ | |_| | |_) | |_| | (_) | | | \__ \ #
#   |_||_|  \___|\___|  \___/| .__/ \__|_|\___/|_| |_|___/ #
#                            |_|                           #
############################################################

def getTreeOptions(options=None):
    if options == None:
        #options = VarParsing.VarParsing()
        #options = VarParsing.VarParsing ('analysis')
        options = VarParsing.VarParsing ('standard')
        
    options.maxEvents = -1 # default, -1 means all events

    #Setting input tag as option
    options.register ('globaltag',
		  'NOT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Global tag from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions")
    #Setting applyskim as option
    options.register ('applySkim',
		  1, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Filters skimed events")
    #Setting isMC as option
    options.register ('isMC',
		  0, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Tells if the sample is MC or not (default = 0)")


    options.register ('sample',
		  '', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Sample type for MC  (isMC=True must be provided ). Posible options: [DYjets,ZZ,ZZtau]. Other options are not used in the analysis."
		  )
    options.register ('allowB',
		  'NOT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Set to 0 to filter out heavy flavors and 1 to filter out light flavors."
		  )		  
    options.register ('finalstate',
		  'LLLL', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Final state to analyze. Possible options are [default=LLLL, EEEE, MMMM, EEMM]."
		  )
    options.register ('dataset',
		  'NOT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Name of data set to be used for electron calibration (Claude provides it)")
		  
    options.register ('leptonSetup',
		  2011, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Tells about the tags used for rho corrctions and possibly other things.")		  
		  
    options.register ('goodlumi',
		  'NO_JSON_FILE', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "JSON  good lumisection file to be passed for data. No default."
		  )
		  
    options.register ('useHLTdata',
		  '',
		  #'NO_HLT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.list, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "HLT paths that can be provided to configuration (example: useHLT=first_path,second_path...). Default defined in configuration."
		  )
		  
    options.register ('vetoHLTdata',
		  '',
		  #'NO_VETO_HLT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.list, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Veto HLT paths that can be provided to configuration as list [example:vetoHLT=firs_path,second_path...]. Default defined in configuration."
		  )
	  
    options.register ('eleHLTmc',
		  '',
		  #'NO_HLT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.list, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "HLT paths that can be provided to configuration (example: eleHLTmc=first_path,second_path...). Used in case of MC. Default defined in configuration."
		  )
		  
    options.register ('muHLTmc',
		  '',
		  #'NO_VETO_HLT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.list, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Veto HLT paths that can be provided to configuration as list [example:muHLTmc=firs_path,second_path...]. Used in case of MC. Default defined in configuration."
		  )
			  
    options.register ('llrTree',
		  0, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Tells if the LLR simple tree should be produced(default = 0)"
		  )
		  
    options.register ('skim',
		  0, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Tells if the skimming should be applyed (default = 0)")		  
		  
    options.register ('bdtEleId',
		  0, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Switch on to do BDT electron ID (default = 0)")		  		  

    options.register ('pfIso',
		  0, # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.int,          # string, int, or float
		  "Switch on to do use ParticleFlow isolation (default = 0)")		  		  
		  
    options.register ('tag',
		  '',
		  #'NO_VETO_HLT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Tag to add to file names."
		  )
		  
    options.register ('filePrepend',
		  '',
		  #'NO_VETO_HLT_PROVIDED', # default value
		  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
		  VarParsing.VarParsing.varType.string,          # string, int, or float
		  "Prefix added to output file names. Could be directory path."
		  )

    # Protection in case sys.argv is missing due to various edm tools
    if not hasattr(sys, "argv"):
        return options

    # Hack to be able to pass multiple arguments from multicrab.cfg
    print sys.argv
    if len(sys.argv) > 0:
        last = sys.argv.pop()
        sys.argv.extend(last.split(":"))
        print sys.argv

    #options.parseArguments()

    return options

