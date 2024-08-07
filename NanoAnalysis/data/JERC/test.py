import os
import correctionlib._core as core

# path to directory of this script
__this_dir__ = os.path.dirname(__file__)


example_value_dict = {
    # jet transverse momentum
    "JetPt": 100.0,
    # jet pseudorapidity
    "JetEta": 0.0,
    # jet azimuthal angle
    "JetPhi": 0.2,
    # jet area
    "JetA": 0.5,
    # median energy density (pileup)
    "Rho": 15.0,
    # systematic variation (only for JER SF)
    "systematic": "nom",
    # pT of matched gen-level jet (only for JER smearing)
    "GenPt": 80.0,  # or -1 if no match
    # unique event ID used for deterministic
    # pseudorandom number generation (only for JER smearing)
    "EventID": 12345,
}



#
# JEC-related examples
#

# JEC base tag
jec = "Summer19UL16_V7_MC"

# jet algorithms
algo = "AK4PFchs"
algo_ak8 = "AK8PFPuppi"

# jet energy correction level
lvl = "L2Relative"

# jet energy correction level
lvl_compound = "L1L2L3Res"

# jet energy uncertainty
unc = "Total"

# print input information
print("\n\nJEC parameters")
print("##############\n")

print("jec = {}".format(jec))
print("algo = {}".format(algo))
print("algo_ak8 = {}".format(algo_ak8))
for v in ("JetPt", "JetEta", "JetA", "JetPhi", "JetA", "Rho"):
    print("{} = {}".format(v, example_value_dict[v]))

# tool for JER smearing
fname_jersmear = os.path.join(__this_dir__, "jer_smear.json")
print("\nLoading JSON file: {}".format(fname_jersmear))
cset_jersmear = core.CorrectionSet.from_file(os.path.join(fname_jersmear))


