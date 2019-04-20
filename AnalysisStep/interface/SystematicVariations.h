#ifndef SYSTEMATICVARIATIONS_H
#define SYSTEMATICVARIATIONS_H


namespace SystematicsHelpers{

  enum SystematicVariationTypes{
    sNominal,
    tPDFScaleDn, tPDFScaleUp,
    tQCDScaleDn, tQCDScaleUp,
    tAsMZDn, tAsMZUp,
    tPDFReplicaDn, tPDFReplicaUp,
    tPythiaScaleDn, tPythiaScaleUp,
    tPythiaTuneDn, tPythiaTuneUp,
    eLepSFEleDn, eLepSFEleUp,
    eLepSFMuDn, eLepSFMuUp,
    eLepScaleEleDn, eLepScaleEleUp,
    eLepScaleMuDn, eLepScaleMuUp,
    eLepResEleDn, eLepResEleUp,
    eLepResMuDn, eLepResMuUp,
    eMETDn, eMETUp,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    eBTagSFDn, eBTagSFUp,
    nSystematicVariations
  };

}

#endif
