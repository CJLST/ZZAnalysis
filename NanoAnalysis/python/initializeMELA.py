import Mela



### initializes a pointer to a MELA object given parameters set by the user. 
def initializeMELA(runMELA, year, matrixelement = None, process = None, couplings = None, production = None, prod = None, dec = None, computeprop = None, propscheme = Mela.ResonancePropagatorScheme.FixedWidth, dividep = None, separatewwzz = False, useconstant = False, match_MX = False, lepton_interference = Mela.LeptonInterference.DefaultLeptonInterf): 

   
    if runMELA == True:
        sqrts = 13. 
        if year >= 2022:
            sqrts = 13.6
        import sys, os
        devnull = os.open(os.devnull, os.O_WRONLY)
        stdout_fd = sys.stdout.fileno()
        saved_stdout = os.dup(stdout_fd)
        os.dup2(devnull, stdout_fd)
        m = Mela.Mela(sqrts, 125, Mela.VerbosityLevel.ERROR) 
        m.setCandidateDecayMode(Mela.CandidateDecayMode.CandidateDecay_ZZ)  
        os.dup2(saved_stdout, stdout_fd)
        os.close(devnull)
        os.close(saved_stdout)
        print(f"initializeMELA: created Mela({sqrts:.1f},125,TVar.CandidateDecay_ZZ)", flush=True)
        print("initializeMELA: ", "MATRIXELEMENT: ", matrixelement, " PROCESS: ", process) 

        requiredSettings = [matrixelement, process, couplings, production, prod, dec, computeprop]
    
        ### If at least one of the essential settings to calculate probabilities are set, make clear that only angles will be computed. 
        if any(elem is None for elem in requiredSettings): 
            print("initializeMELA: The required settings to compute probabilites have not been given. Only angles will be computed.")
            probSettingsDict = None
        else: 
            probSettingsDict = settingsParser(matrixelement, process, couplings, production, prod, dec, computeprop, propscheme, dividep, separatewwzz, useconstant, match_MX, lepton_interference)
            m.setProcess(probSettingsDict["process"], probSettingsDict["matrixelement"], probSettingsDict["production"])

        

    else: 
        m = None
        probSettingsDict = None
    return m, probSettingsDict



def check_enum(entry, enum):
    found = False
    i = 0
    possible_value = tuple(enum.__members__.keys())
    mapping = enum.__members__
    while (not found) and (i < len(possible_value)):
        if entry.lower() == possible_value[i].lower():
            found = True
            entry = possible_value[i]
        i += 1
    if not found:
        possible_value = tuple(map(str.lower, possible_value))
        errortext = "Unknown matrix element given!"
        errortext += "\nThe following are valid matrix elements"
        errortext += "\n" + "\n".join(possible_value)
        errortext = help.print_msg_box(errortext, title="ERROR")
        raise ValueError("\n" + errortext)
    return mapping[entry]

def couplingsParser(couplings):
    couplings = couplings.replace("|", ",")
    couplings = couplings.replace("+", ",")

    coupling_list = eval(couplings)
    
    return coupling_list


def settingsParser(matrixelement, process, couplings, production, prod, dec, computeprop, propscheme, dividep, separatewwzz, useconstant, match_MX, lepton_interference):
    settingsDict = {}
    settingsDict["matrixelement"] = check_enum(matrixelement, Mela.MatrixElement)
    settingsDict["process"] = check_enum(process, Mela.Process)
    settingsDict["couplings"] = couplingsParser(couplings)
    settingsDict["production"] = check_enum(production, Mela.Production)
    settingsDict["prod"] = prod
    settingsDict["dec"] = dec
    settingsDict["computeprop"] = computeprop 
    if type(propscheme) == str: 
        settingsDict["propscheme"] = check_enum(propscheme, Mela.ResonancePropagatorScheme)
    settingsDict["dividep"] = dividep ###TODO is this really useful if we're only just calculating the native probability? 
    settingsDict["separatewwzz"] = separatewwzz 
    settingsDict["useconstant"] = useconstant
    settingsDict["match_MX"] = match_MX
    if type (lepton_interference) == str:
        settingsDict["lepton_interference"] = check_enum(lepton_interference, Mela.LeptonInterference)

    return settingsDict
