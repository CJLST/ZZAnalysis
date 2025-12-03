##
# Return the function to be used to check if an ele passes the BDT cut, pre-configured for:
# - era (2018, 2022...)
# - dataTag ('UL', etc)
# - nanoAOD version (9, ecc)
##

def getEleBDTCut(era, dataTag, nanoVersion, useUncorrPt=False) :

    # nanoAODv9 and older include only mvaFall17V2Iso = 2017 WP and training (ElectronMVAEstimatorRun2Fall17IsoV2Values), for all years
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

    def eleBDTCut_RunII2016UL_v15(ele) :
        fSCeta = ele.superclusterEta
        BDT = ele.mvaHZZIso
        return (ele.pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9557993256) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9475406570) or \
                                     (fSCeta>=1.479                and BDT > 0.9285158721))) \
                or (ele.pt>10. and  ((fSCeta<0.8                   and BDT > 0.3272075608) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.2468345995) or \
                                     (fSCeta>=1.479                and BDT > -0.5955762814)))

    def eleBDTCut_RunII2017UL_v15(ele) :
        fSCeta = ele.superclusterEta
        BDT = ele.mvaHZZIso
        return (ele.pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9128577458 ) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9056792368 ) or \
                                     (fSCeta>=1.479                and BDT > 0.9439440575 ))) \
                or (ele.pt>10. and  ((fSCeta<0.8                   and BDT > 0.1559788054 ) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0273863727 ) or \
                                     (fSCeta>=1.479                and BDT > -0.5532483665)))

    def eleBDTCut_RunII2018UL_v15(ele) :
        fSCeta = ele.superclusterEta
        BDT = ele.mvaHZZIso
        return (ele.pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9044286167) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9094166886) or \
                                     (fSCeta>=1.479                and BDT > 0.9443653660))) \
                or (ele.pt>10. and  ((fSCeta<0.8                   and BDT > 0.1968600840) or \
                                     (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0759172100) or \
                                     (fSCeta>=1.479                and BDT > -0.5169136775)))

    # Run3 nanoAODv12 samples have the 2018 UL tuning (ElectronMVAEstimatorRun2Summer18ULIdIsoValues)
    # The WP was derived before scale corrections, so the uncorrected pt should be used when available.
    def eleBDTCut_RunIII_ULTraining_def(ele) :
        return(eleBDTCut_RunIII_ULTraining(ele.pt, abs(ele.eta + ele.deltaEtaSC), ele.mvaHZZIso))
        
    def eleBDTCut_RunIII_ULTraining_uncorr(ele) :
        return(eleBDTCut_RunIII_ULTraining(ele.uncorrected_pt, abs(ele.eta + ele.deltaEtaSC), ele.mvaHZZIso))
               
    def eleBDTCut_RunIII_ULTraining(pt, fSCeta, BDT) :
        return (pt<=10. and     ((fSCeta<0.8                   and BDT > 0.9044286167) or \
                                 (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.9094166886) or \
                                 (fSCeta>=1.479                and BDT > 0.9443653660))) \
               or (pt>10. and  ((fSCeta<0.8                   and BDT > 0.1968600840) or \
                                (fSCeta>=0.8 and fSCeta<1.479 and BDT > 0.0759172100) or \
                                (fSCeta>=1.479                and BDT > -0.5169136775)))

    def eleBDTCut_RunIII_2022Training_WP(ele):
        return ele.mvaIso_WPHZZ

    
    if era >= 2016 and era <= 2018 :
        if nanoVersion < 10 :
            if "UL" in dataTag :
                return eleBDTCut_RunIIUL_v9
            else :
                return eleBDTCut_RunIIpreUL_v9
        elif nanoVersion == 15 :
            cutsv15 = {2016: eleBDTCut_RunII2016UL_v15, 2017: eleBDTCut_RunII2017UL_v15, 2018: eleBDTCut_RunII2018UL_v15}
            return cutsv15[era]
        else :
            raise ValueError('getEleBDTCut: era '+ str(era)+', dataTag ' + dataTag + ', nanoVersion ' + str(nanoVersion) + ' not supported')
          

    elif era >=2022 :
        if nanoVersion <14:
            if useUncorrPt:
                return eleBDTCut_RunIII_ULTraining_uncorr
            else :
                return eleBDTCut_RunIII_ULTraining_def
        else:
            return eleBDTCut_RunIII_2022Training_WP

    # Fallback: combination not supported
    raise ValueError('getEleBDTCut: era '+ str(era)+', dataTag ' + dataTag + ', nanoVersion ' + str(nanoVersion) + ' not supported')


