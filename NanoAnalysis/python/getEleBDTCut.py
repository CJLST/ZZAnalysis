##
# Return the function to be used to check if an ele passes the BDT cut, pre-configured for:
# - era (2018, 2022...)
# - dataTag ('UL', etc)
# - nanoAOD version (9, ecc)
##

def getEleBDTCut(era, dataTag, nanoVersion) :
    # nanoAODv9 includes mvaFall17V2Iso = 2017 WP and training (ElectronMVAEstimatorRun2Fall17IsoV2Values)

    def eleBDTCut_RunIIpreUL_v9(ele) :
        # pre-UL WP for Run II (miniAOD branch: Run2_CutBased_BTag16)
        fSCeta = abs(ele.eta + ele.deltaEtaSC)
        BDT = ele.mvaFall17V2Iso
        return (ele.pt<=10. and     ((fSCeta<0.8                   and BDT > 0.85216885148) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.82684550976) or \
                                     (fSCeta>=1.479                and BDT > 0.86937630022))) \
                or (ele.pt>10. and  ((fSCeta<0.8                   and BDT > 0.98248928759) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.96919224579) or \
                                     (fSCeta>=1.479                and BDT > 0.79349796445)))
  
    def eleBDTCut_RunIIUL_v9(ele) :
        # UL WP (miniAOD branch Run2_CutBased_UL)
        fSCeta = abs(ele.eta + ele.deltaEtaSC)
        BDT = ele.mvaFall17V2Iso
        return (ele.pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9128577458) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9056792368) or \
                                     (fSCeta>=1.479                and BDT > 0.9439440575))) \
                or (ele.pt>10. and  ((fSCeta<0.8                   and BDT > 0.1559788054) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0273863727) or \
                                     (fSCeta>=1.479                and BDT > -0.5532483665)))


    def eleBDTCut_RunIII_ULTraining(ele) :
        fSCeta = abs(ele.eta + ele.deltaEtaSC)
        BDT = ele.mvaHZZIso #In Run3 samples: 2018 UL tuning (ElectronMVAEstimatorRun2Summer18ULIdIsoValues)
        return (ele.pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9044286167) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9094166886) or \
                                     (fSCeta>=1.479                and BDT > 0.9443653660))) \
                or (ele.pt>10. and  ((fSCeta<0.8                   and BDT > 0.1968600840) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0759172100) or \
                                     (fSCeta>=1.479                and BDT > -0.5169136775)))

    if era == 2017 or era == 2018 :
        if "UL" in dataTag :
            if nanoVersion <10 :
                return eleBDTCut_RunIIUL_v9
        else :
            if nanoVersion <10 :
                return eleBDTCut_RunIIpreUL_v9

    elif era >=2022 :
        return eleBDTCut_RunIII_ULTraining

    # Fallback: combination not supported
    raise ValueError('getEleBDTCut: era '+ str(era)+', dataTag ' + dataTag + ', nanoVersion ' + nanoVersion + ' not supported')


