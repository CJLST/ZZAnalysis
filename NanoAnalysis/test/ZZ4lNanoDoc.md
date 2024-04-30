# H4l nanoAOD file format documentation
This document is generated automatically with:  
`./inspectNanoFile.py ZZ4lAnalysis.root --docmd ZZ4lNanoDoc.md`

Jump to:

* [Events tree](#events-tree-content)
* [AllEvents tree](#allevents-tree-content)
* [Runs tree](#runs-tree-content)
* [LuminosityBlocks tree](#luminosityblocks-tree-content)
## Events tree content

| Collection | Description |
| - | - |
| [**Electron**](#electron) | slimmedElectrons after basic selection (pt > 5 ) |
| [**FidDressedLeps**](#fiddressedleps) | gen dressed leps for fiducial analysis |
| [**FidZ**](#fidz) | Zs made of FidDressedLeps |
| [**FidZ1**](#fidz1) | FidZ1_mass/F |
| [**FidZ2**](#fidz2) | FidZ2_mass/F |
| [**FidZZ**](#fidzz) | Index of 1st Z1 daughter in FidDressedLeps collection |
| [**Flag**](#flag) | Trigger/flag bit (process: PAT) |
| [**FsrPhoton**](#fsrphoton) | Final state radiation photons emitted by muons or electrons |
| [**GenDressedLepton**](#gendressedlepton) | Dressed leptons from Rivet-based ParticleLevelProducer |
| [**GenPart**](#genpart) | interesting gen particles  |
| [**GenZZ**](#genzz) | product of pdgId of the four gen leptons from ZZ decay |
| [**Generator**](#generator) | MC generator weight |
| [**HLT**](#hlt) | Trigger/flag bit (process: HLT) |
| [**HTXS**](#htxs) | number of jets with pt>30 GeV as identified in HTXS |
| [**Jet**](#jet) | slimmedJetsPuppi, i.e. ak4 PFJets Puppi with JECs applied, after basic selection (pt > 15) |
| [**JetLeadingIdx**](#jetleadingidx) | JetLeadingIdx/S |
| [**JetSubleadingIdx**](#jetsubleadingidx) | JetSubleadingIdx/S |
| [**LHEPdfWeight**](#lhepdfweight) | LHE pdf variation weights (w_var / w_nominal) for LHA IDs 325300 - 325402 |
| [**LHEReweightingWeight**](#lhereweightingweight) |  |
| [**LHEScaleWeight**](#lhescaleweight) | LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0  |
| [**MET**](#met) | pt |
| [**Muon**](#muon) | slimmedMuons after basic selection (pt > 15 \|\| (pt > 3 && (passed("CutBasedIdLoose") \|\| passed("SoftCutBasedId") \|\| passed("SoftMvaId") \|\| passed("CutBasedIdGlobalHighPt") \|\| passed("CutBasedIdTrkHighPt")))) |
| [**PSWeight**](#psweight) | PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5; |
| [**Pileup**](#pileup) | the number of pileup interactions that have been added to the event in the current bunch crossing |
| [**ZCand**](#zcand) | Z candidates passing the full H4l selection |
| [**ZLCand**](#zlcand) | Index of extra lep for the ZL CR |
| [**ZLLCand**](#zllcand) | Z+LL control region candidates |
| [**ZLLbest2P2FIdx**](#zllbest2p2fidx) | best candidate for the 2P2F CR |
| [**ZLLbest3P1FIdx**](#zllbest3p1fidx) | best candidate for the 3P1F CR |
| [**ZLLbestSIPCRIdx**](#zllbestsipcridx) | best candidate for the SIP CR |
| [**ZLLbestSSIdx**](#zllbestssidx) | best candidate for the SS CR |
| [**ZZCand**](#zzcand) | ZZ candidates passing the full H4l selection |
| [**bestCandIdx**](#bestcandidx) | Seleced ZZ candidate in the event |
| [**bestZIdx**](#bestzidx) | Best Z in the event (mass closest to mZ) |
| [**event**](#event) | event/l |
| [**genWeight**](#genweight) | generator weight |
| [**ggH**](#ggh) | Reweighting for ggH as a function of njets and pT |
| [**luminosityBlock**](#luminosityblock) | luminosityBlock/i |
| [**nCleanedJetsPt30**](#ncleanedjetspt30) | nCleanedJetsPt30/B |
| [**overallEventWeight**](#overalleventweight) | overallEventWeight/F |
| [**passedFiducial**](#passedfiducial) | event passes fiducial selection at gen level |
| [**puWeight**](#puweight) | puWeight/F |
| [**run**](#run) | run/i |

## Events tree detail

### <a id='electron'></a>Electron [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Electron_ZZFullId** | Bool_t| pass H4l full ID selection |
| **Electron_ZZFullSel** | Bool_t| pass H4l full SR selection (for electrons,=FullID) |
| **Electron_ZZFullSelNoSIP** | Bool_t| pass H4l full ID without SIP (base for CR SIP  method) |
| **Electron_ZZRelaxedId** | Bool_t| pass H4l relaxed ID including SIP (base for SS and OS CRs) |
| **Electron_ZZRelaxedIdNoSIP** | Bool_t| pass H4l relaxed ID without SIP (widest subset of cuts commont to all CRs) |
| **Electron_charge** | Int_t| electric charge |
| **Electron_convVeto** | Bool_t| pass conversion veto |
| **Electron_cutBased** | UChar_t| cut-based ID RunIII Winter22 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight) |
| **Electron_cutBased_HEEP** | Bool_t| cut-based HEEP ID |
| **Electron_deltaEtaSC** | Float_t| delta eta (SC,ele) with sign |
| **Electron_dr03EcalRecHitSumEt** | Float_t| Non-PF Ecal isolation within a delta R cone of 0.3 with electron pt > 35 GeV |
| **Electron_dr03HcalDepth1TowerSumEt** | Float_t| Non-PF Hcal isolation within a delta R cone of 0.3 with electron pt > 35 GeV |
| **Electron_dr03TkSumPt** | Float_t| Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV |
| **Electron_dr03TkSumPtHEEP** | Float_t| Non-PF track isolation within a delta R cone of 0.3 with electron pt > 35 GeV used in HEEP ID |
| **Electron_dxy** | Float_t| dxy (with sign) wrt first PV, in cm |
| **Electron_dxyErr** | Float_t| dxy uncertainty, in cm |
| **Electron_dz** | Float_t| dz (with sign) wrt first PV, in cm |
| **Electron_dzErr** | Float_t| dz uncertainty, in cm |
| **Electron_eInvMinusPInv** | Float_t| 1/E_SC - 1/p_trk |
| **Electron_energyErr** | Float_t| energy error of the cluster-track combination |
| **Electron_eta** | Float_t| eta |
| **Electron_fsrPhotonIdx** | Short_t(index to Fsrphoton)| Index of the lowest-dR/ET2 among associated FSR photons |
| **Electron_genPartFlav** | UChar_t| Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 electrons or photons: 1 = prompt electron (including gamma*->mu mu), 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched |
| **Electron_genPartIdx** | Short_t(index to Genpart)| Index into genParticle list for MC matching to status==1 electrons or photons |
| **Electron_hoe** | Float_t| H over E |
| **Electron_ip3d** | Float_t| 3D impact parameter wrt first PV, in cm |
| **Electron_isPFcand** | Bool_t| electron is PF candidate |
| **Electron_jetIdx** | Short_t(index to Jet)| index of the associated jet (-1 if none) |
| **Electron_jetNDauCharged** | UChar_t| number of charged daughters of the closest jet |
| **Electron_jetPtRelv2** | Float_t| Relative momentum of the lepton with respect to the closest jet after subtracting the lepton |
| **Electron_jetRelIso** | Float_t| Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet) |
| **Electron_lostHits** | UChar_t| number of missing inner hits |
| **Electron_mass** | Float_t| mass |
| **Electron_miniPFRelIso_all** | Float_t| mini PF relative isolation, total (with scaled rho*EA PU Winter22V1 corrections) |
| **Electron_miniPFRelIso_chg** | Float_t| mini PF relative isolation, charged component |
| **Electron_mvaHZZIso** | Float_t| HZZ MVA Iso ID score |
| **Electron_mvaIso** | Float_t| MVA Iso ID score, Winter22V1 |
| **Electron_mvaIso_WP80** | Bool_t| MVA Iso ID WP80, Winter22V1 |
| **Electron_mvaIso_WP90** | Bool_t| MVA Iso ID WP90, Winter22V1 |
| **Electron_mvaNoIso** | Float_t| MVA noIso ID score, Winter22V1 |
| **Electron_mvaNoIso_WP80** | Bool_t| MVA noIso ID WP80, Winter22V1 |
| **Electron_mvaNoIso_WP90** | Bool_t| MVA noIso ID WP90, Winter22V1 |
| **Electron_mvaTTH** | Float_t| TTH MVA lepton ID score |
| **Electron_passBDT** | Bool_t| pass H4l BDT cut |
| **Electron_passIso** | Bool_t| always True for electrons |
| **Electron_pdgId** | Int_t| PDG code assigned by the event reconstruction (not by MC truth) |
| **Electron_pfRelIso03FsrCorr** | Float_t| FSR-subtracted pfRelIso03 |
| **Electron_pfRelIso03_all** | Float_t| PF relative isolation dR=0.3, total (with rho*EA PU Winter22V1 corrections) |
| **Electron_pfRelIso03_chg** | Float_t| PF relative isolation dR=0.3, charged component |
| **Electron_phi** | Float_t| phi |
| **Electron_photonIdx** | Short_t(index to Photon)| index of the first associated photon (-1 if none) |
| **Electron_pt** | Float_t| pt |
| **Electron_r9** | Float_t| R9 of the supercluster, calculated with full 5x5 region |
| **Electron_scEtOverPt** | Float_t| (supercluster transverse energy)/pt-1 |
| **Electron_scaleDn_pt** | Float_t| Electron_scaleDn_pt[nElectron]/F |
| **Electron_scaleUp_pt** | Float_t| Electron_scaleUp_pt[nElectron]/F |
| **Electron_seedGain** | UChar_t| Gain of the seed crystal |
| **Electron_seediEtaOriX** | Char_t| iEta or iX of seed crystal. iEta is barrel-only, iX is endcap-only. iEta runs from -85 to +85, with no crystal at iEta=0. iX runs from 1 to 100. |
| **Electron_seediPhiOriY** | Int_t| iPhi or iY of seed crystal. iPhi is barrel-only, iY is endcap-only. iPhi runs from 1 to 360. iY runs from 1 to 100. |
| **Electron_sieie** | Float_t| sigma_IetaIeta of the supercluster, calculated with full 5x5 region |
| **Electron_sip3d** | Float_t| 3D impact parameter significance wrt first PV, in cm |
| **Electron_smearDn_pt** | Float_t| Electron_smearDn_pt[nElectron]/F |
| **Electron_smearUp_pt** | Float_t| Electron_smearUp_pt[nElectron]/F |
| **Electron_svIdx** | Short_t(index to Sv)| index of matching secondary vertex |
| **Electron_tightCharge** | UChar_t| Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent) |
| **Electron_uncorrected_pt** | Float_t| Electron_uncorrected_pt[nElectron]/F |
| **Electron_vidNestedWPBitmap** | Int_t| VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleEBEECut,GsfEleHadronicOverEMEnergyScaledCut,GsfEleEBEECut,GsfEleRelPFIsoScaledCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 3 bits per cut |
| **Electron_vidNestedWPBitmapHEEP** | Int_t| VID compressed bitmap (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleEBEECut,GsfEleEBEECut,GsfEleFull5x5SigmaIEtaIEtaWithSatCut,GsfEleFull5x5E2x5OverE5x5WithSatCut,GsfEleHadronicOverEMLinearCut,GsfEleTrkPtIsoCut,GsfEleEmHadD1IsoRhoCut,GsfEleDxyCut,GsfEleMissingHitsCut,GsfEleEcalDrivenCut), 1 bits per cut |
| **nElectron** | Int_t| slimmedElectrons after basic selection (pt > 5 ) |

### <a id='fiddressedleps'></a>FidDressedLeps [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidDressedLeps_RelIso** | Float_t| FidDressedLeps_RelIso[nFidDressedLeps]/F |
| **FidDressedLeps_eta** | Float_t| FidDressedLeps_eta[nFidDressedLeps]/F |
| **FidDressedLeps_id** | Float_t| FidDressedLeps_id[nFidDressedLeps]/F |
| **FidDressedLeps_mass** | Float_t| FidDressedLeps_mass[nFidDressedLeps]/F |
| **FidDressedLeps_momid** | Float_t| FidDressedLeps_momid[nFidDressedLeps]/F |
| **FidDressedLeps_mommomid** | Float_t| FidDressedLeps_mommomid[nFidDressedLeps]/F |
| **FidDressedLeps_phi** | Float_t| FidDressedLeps_phi[nFidDressedLeps]/F |
| **FidDressedLeps_pt** | Float_t| FidDressedLeps_pt[nFidDressedLeps]/F |
| **nFidDressedLeps** | Int_t| gen dressed leps for fiducial analysis |

### <a id='fidz'></a>FidZ [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZ_DauPdgId** | Int_t| FidZ decay flavour |
| **FidZ_MomPdgId** | Int_t| FidZ mother Id |
| **nFidZ** | Int_t| Zs made of FidDressedLeps |

### <a id='fidz1'></a>FidZ1 [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZ1_mass** | Float_t| FidZ1_mass/F |

### <a id='fidz2'></a>FidZ2 [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZ2_mass** | Float_t| FidZ2_mass/F |

### <a id='fidzz'></a>FidZZ [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZZ_Z1l1Idx** | Short_t(index to Z1L1)| Index of 1st Z1 daughter in FidDressedLeps collection |
| **FidZZ_Z1l2Idx** | Short_t(index to Z1L2)| Index of 2nd Z1 daughter in FidDressedLeps collection |
| **FidZZ_Z2l1Idx** | Short_t(index to Z2L1)| Index of 1st Z2 daughter in FidDressedLeps collection |
| **FidZZ_Z2l2Idx** | Short_t(index to Z2L2)| Index of 2nd Z1 daughter in FidDressedLeps collection |
| **FidZZ_eta** | Float_t| FidZZ_eta/F |
| **FidZZ_mass** | Float_t| mass of gen ZZ made with FidDressedLeps |
| **FidZZ_phi** | Float_t| FidZZ_phi/F |
| **FidZZ_pt** | Float_t| FidZZ_pt/F |
| **FidZZ_rapidity** | Float_t| FidZZ_rapidity/F |

### <a id='flag'></a>Flag [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Flag_BadChargedCandidateFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadChargedCandidateSummer16Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadPFMuonDzFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadPFMuonFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_BadPFMuonSummer16Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_CSCTightHalo2015Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_CSCTightHaloFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_CSCTightHaloTrkMuUnvetoFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_EcalDeadCellBoundaryEnergyFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_EcalDeadCellTriggerPrimitiveFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_HBHENoiseFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_HBHENoiseIsoFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_HcalStripHaloFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_METFilters** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_chargedHadronTrackResolutionFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_ecalBadCalibFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_ecalLaserCorrFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_eeBadScFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_globalSuperTightHalo2016Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_globalTightHalo2016Filter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_goodVertices** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_hcalLaserEventFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_hfNoisyHitsFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_muonBadTrackFilter** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOGFilters** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOG_logErrorTooManyClusters** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOG_manystripclus53X** | Bool_t| Trigger/flag bit (process: PAT) |
| **Flag_trkPOG_toomanystripclus53X** | Bool_t| Trigger/flag bit (process: PAT) |

### <a id='fsrphoton'></a>FsrPhoton [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FsrPhoton_dROverEt2** | Float_t| deltaR to associated muon divided by photon et2 |
| **FsrPhoton_electronIdx** | Short_t(index to Electron)| index of associated electron |
| **FsrPhoton_eta** | Float_t| eta |
| **FsrPhoton_genFsrIdx** | Short_t(index to Genfsr)| index of GenPart matched to FSRPhoton (ie closest gen FSR of a lepton from Z decay) |
| **FsrPhoton_mass** | Float_t| FsrPhoton_mass[nFsrPhoton]/F |
| **FsrPhoton_muonIdx** | Short_t(index to Muon)| index of associated muon |
| **FsrPhoton_phi** | Float_t| phi |
| **FsrPhoton_pt** | Float_t| pt |
| **FsrPhoton_relIso03** | Float_t| relative isolation in a 0.3 cone without CHS |
| **nFsrPhoton** | Int_t| Final state radiation photons emitted by muons or electrons |

### <a id='gendressedlepton'></a>GenDressedLepton [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **GenDressedLepton_eta** | Float_t| eta |
| **GenDressedLepton_hasTauAnc** | Bool_t| true if Dressed lepton has a tau as ancestor |
| **GenDressedLepton_mass** | Float_t| mass |
| **GenDressedLepton_pdgId** | Int_t| PDG id |
| **GenDressedLepton_phi** | Float_t| phi |
| **GenDressedLepton_pt** | Float_t| pt |
| **nGenDressedLepton** | Int_t| Dressed leptons from Rivet-based ParticleLevelProducer |

### <a id='genpart'></a>GenPart [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **GenPart_eta** | Float_t| eta |
| **GenPart_genPartIdxMother** | Short_t(index to Genpart)| index of the mother particle |
| **GenPart_mass** | Float_t| Mass stored for all particles with the exception of quarks (except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411). For these particles, you can lookup the value from PDG. |
| **GenPart_pdgId** | Int_t| PDG id |
| **GenPart_phi** | Float_t| phi |
| **GenPart_pt** | Float_t| pt |
| **GenPart_status** | Int_t| Particle status. 1=stable |
| **GenPart_statusFlags** | UShort_t| gen status flags stored bitwise, bits are: 0 : isPrompt, 1 : isDecayedLeptonHadron, 2 : isTauDecayProduct, 3 : isPromptTauDecayProduct, 4 : isDirectTauDecayProduct, 5 : isDirectPromptTauDecayProduct, 6 : isDirectHadronDecayProduct, 7 : isHardProcess, 8 : fromHardProcess, 9 : isHardProcessTauDecayProduct, 10 : isDirectHardProcessTauDecayProduct, 11 : fromHardProcessBeforeFSR, 12 : isFirstCopy, 13 : isLastCopy, 14 : isLastCopyBeforeFSR,  |
| **nGenPart** | Int_t| interesting gen particles  |

### <a id='genzz'></a>GenZZ [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **GenZZ_FinalState** | Int_t| product of pdgId of the four gen leptons from ZZ decay |
| **GenZZ_Z1l1Idx** | Short_t(index to Z1L1)| index of the 1st gen lepton from ZZ decay in the GenPart collection (pre-FSR) |
| **GenZZ_Z1l2Idx** | Short_t(index to Z1L2)| index of the 2nd gen lepton from ZZ decay in the GenPart collection (pre-FSR) |
| **GenZZ_Z2l1Idx** | Short_t(index to Z2L1)| index of the 3rd gen lepton from ZZ decay in the GenPart collection (pre-FSR) |
| **GenZZ_Z2l2Idx** | Short_t(index to Z2L2)| index of the 4th gen lepton from ZZ decay in the GenPart collection (pre-FSR) |
| **GenZZ_mass** | Float_t| mass of the 4l system (pre-FSR) |

### <a id='generator'></a>Generator [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Generator_weight** | Float_t| MC generator weight |

### <a id='hlt'></a>HLT [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiMu9_Ele9_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle10_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle24_eta2p1_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle25_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle27_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle33_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle4_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle4p5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle5p5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle6_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle6p5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle7_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle7p5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle8_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle8p5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle9_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_DoubleEle9p5_eta1p22_mMax6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele115_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele135_CaloIdVT_GsfTrkIdT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT450_PFMET50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_IsoVVVL_PFHT600** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele15_WPLoose_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele17_CaloIdM_TrackIdM_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele20_WPLoose_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_CaloIdM_TrackIdM_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele27_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele28_HighEta_SC20_Mass55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele28_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele28_eta2p1_WPTight_Gsf_HT150** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele30_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele32_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele32_WPTight_Gsf_L1DoubleEG** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele35_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele35_WPTight_Gsf_L1EGMT** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele38_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele40_WPTight_Gsf** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele50_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Ele8_CaloIdM_TrackIdM_PFJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_TwoProngs35** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS180_eta2p1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_LooseDeepTauPFTauHPS30_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS20_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS35_L2NN_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu24_eta2p1_MediumDeepTauPFTauHPS45_L2NN_eta2p1_CrossL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu27_MediumDeepTauPFTauHPS20_eta2p1_SingleL1** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu50_AK8PFJet230_SoftDropMass40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_IsoMu50_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetBB0p35** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu0_L1DoubleMu** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets100_PFBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets100_PFBTagDeepJet_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets200_PFBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets200_PFBTagDeepJet_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets350_PFBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets350_PFBTagDeepJet_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets40MaxDeta1p6_DoublePFBTagDeepJet_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets40_PFBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets40_PFBTagDeepJet_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepCSV_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_DoublePFJets54MaxDeta1p6_DoublePFBTagDeepJet_p71** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_IP6** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12eta2p3** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu12eta2p3_PFJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT450_PFMET50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu15_IsoVVVL_PFHT600** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_Photon30_IsoCaloId** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu18_Mu9_SameSign** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20NoFiltersNoVtxDisplaced_Photon20_CaloCustomId** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu20_TkMu0_Phi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu25_TkMu0_Onia** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu25_TkMu0_Phi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu27** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu27_Ele37_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu30_TkMu0_Psi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu30_TkMu0_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu37_Ele27_CaloIdL_MW** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu37_TkMu27** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3_L1SingleMu5orSingleMu7** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3_PFJet40** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu4_L1DoubleMu** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu50** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu50_IsoVVVL_PFHT450** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu50_L1SingleMuShower** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu55** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu6HT240_DisplacedDijet30_Inclusive0PtrkShortSig5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu6HT240_DisplacedDijet30_Inclusive1PtrkShortSig5_DisplacedLoose** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu6HT240_DisplacedDijet35_Inclusive0PtrkShortSig5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu6HT240_DisplacedDijet35_Inclusive1PtrkShortSig5_DisplacedLoose** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu6HT240_DisplacedDijet40_Inclusive0PtrkShortSig5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu6HT240_DisplacedDijet40_Inclusive1PtrkShortSig5_DisplacedLoose** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_L2Mu2_Jpsi** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu7p5_L2Mu2_Upsilon** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_DiEle12_CaloIdL_TrackIdL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBTagParticleNet_2BTagSum0p65** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepJet_1p5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_10_5_5_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_12_10_5** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_5_3_3_Mass3p8_DCA** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_TripleMu_5_3_3_Mass3p8_DZ** | Bool_t| Trigger/flag bit (process: HLT) |
| **HLT_passZZ4l** | Bool_t| HLT_passZZ4l/O |
| **HLT_passZZ4lEle** | Bool_t| HLT_passZZ4lEle/O |
| **HLT_passZZ4lMu** | Bool_t| HLT_passZZ4lMu/O |
| **HLT_passZZ4lMuEle** | Bool_t| HLT_passZZ4lMuEle/O |

### <a id='htxs'></a>HTXS [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **HTXS_Higgs_pt** | Float_t| pt of the Higgs boson as identified in HTXS |
| **HTXS_Higgs_y** | Float_t| rapidity of the Higgs boson as identified in HTXS |
| **HTXS_njets30** | UChar_t| number of jets with pt>30 GeV as identified in HTXS |

### <a id='jet'></a>Jet [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Jet_PNetRegPtRawCorr** | Float_t| ParticleNet universal flavor-aware visible pT regression (no neutrinos), correction relative to raw jet pT |
| **Jet_PNetRegPtRawCorrNeutrino** | Float_t| ParticleNet universal flavor-aware pT regression neutrino correction, relative to visible. To apply full regression, multiply raw jet pT by both PNetRegPtRawCorr and PNetRegPtRawCorrNeutrino. |
| **Jet_PNetRegPtRawRes** | Float_t| ParticleNet universal flavor-aware jet pT resolution estimator, (q84 - q16)/2 |
| **Jet_ZZMask** | Bool_t| Jet_ZZMask[nJet]/O |
| **Jet_area** | Float_t| jet catchment area, for JECs |
| **Jet_btagDeepFlavB** | Float_t| DeepJet b+bb+lepb tag discriminator |
| **Jet_btagDeepFlavCvB** | Float_t| DeepJet c vs b+bb+lepb discriminator |
| **Jet_btagDeepFlavCvL** | Float_t| DeepJet c vs uds+g discriminator |
| **Jet_btagDeepFlavQG** | Float_t| DeepJet g vs uds discriminator |
| **Jet_btagPNetB** | Float_t| ParticleNet b vs. udscg |
| **Jet_btagPNetCvB** | Float_t| ParticleNet c vs. b |
| **Jet_btagPNetCvL** | Float_t| ParticleNet c vs. udsg |
| **Jet_btagPNetQvG** | Float_t| ParticleNet q (udsbc) vs. g |
| **Jet_btagPNetTauVJet** | Float_t| ParticleNet tau vs. jet |
| **Jet_btagRobustParTAK4B** | Float_t| RobustParTAK4 b+bb+lepb tag discriminator |
| **Jet_btagRobustParTAK4CvB** | Float_t| RobustParTAK4 c vs b+bb+lepb discriminator |
| **Jet_btagRobustParTAK4CvL** | Float_t| RobustParTAK4 c vs uds+g discriminator |
| **Jet_btagRobustParTAK4QG** | Float_t| RobustParTAK4 g vs uds discriminator |
| **Jet_chEmEF** | Float_t| charged Electromagnetic Energy Fraction |
| **Jet_chHEF** | Float_t| charged Hadron Energy Fraction |
| **Jet_electronIdx1** | Short_t(index to Electron)| index of first matching electron |
| **Jet_electronIdx2** | Short_t(index to Electron)| index of second matching electron |
| **Jet_eta** | Float_t| eta |
| **Jet_genJetIdx** | Short_t(index to Genjet)| index of matched gen jet |
| **Jet_hadronFlavour** | UChar_t| flavour from hadron ghost clustering |
| **Jet_hfadjacentEtaStripsSize** | Int_t| eta size of the strips next to the central tower strip in HF (noise discriminating variable) |
| **Jet_hfcentralEtaStripSize** | Int_t| eta size of the central tower strip in HF (noise discriminating variable) |
| **Jet_hfsigmaEtaEta** | Float_t| sigmaEtaEta for HF jets (noise discriminating variable) |
| **Jet_hfsigmaPhiPhi** | Float_t| sigmaPhiPhi for HF jets (noise discriminating variable) |
| **Jet_jetId** | UChar_t| Jet ID flag: bit2 is tight, bit3 is tightLepVeto |
| **Jet_mass** | Float_t| mass |
| **Jet_muEF** | Float_t| muon Energy Fraction |
| **Jet_muonIdx1** | Short_t(index to Muon)| index of first matching muon |
| **Jet_muonIdx2** | Short_t(index to Muon)| index of second matching muon |
| **Jet_muonSubtrFactor** | Float_t| 1-(muon-subtracted raw pt)/(raw pt) |
| **Jet_nConstituents** | UChar_t| Number of particles in the jet |
| **Jet_nElectrons** | UChar_t| number of electrons in the jet |
| **Jet_nMuons** | UChar_t| number of muons in the jet |
| **Jet_nSVs** | UChar_t| number of secondary vertices in the jet |
| **Jet_neEmEF** | Float_t| neutral Electromagnetic Energy Fraction |
| **Jet_neHEF** | Float_t| neutral Hadron Energy Fraction |
| **Jet_partonFlavour** | Short_t| flavour from parton matching |
| **Jet_phi** | Float_t| phi |
| **Jet_pt** | Float_t| pt |
| **Jet_rawFactor** | Float_t| 1 - Factor to get back to raw pT |
| **Jet_svIdx1** | Short_t(index to Sv)| index of first matching secondary vertex |
| **Jet_svIdx2** | Short_t(index to Sv)| index of second matching secondary vertex |
| **nJet** | Int_t| slimmedJetsPuppi, i.e. ak4 PFJets Puppi with JECs applied, after basic selection (pt > 15) |

### <a id='jetleadingidx'></a>JetLeadingIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **JetLeadingIdx** | Short_t(index to Jetleading)| JetLeadingIdx/S |

### <a id='jetsubleadingidx'></a>JetSubleadingIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **JetSubleadingIdx** | Short_t(index to Jetsubleading)| JetSubleadingIdx/S |

### <a id='lhepdfweight'></a>LHEPdfWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **LHEPdfWeight** | Float_t| LHE pdf variation weights (w_var / w_nominal) for LHA IDs 325300 - 325402 |
| **nLHEPdfWeight** | Int_t|  |

### <a id='lhereweightingweight'></a>LHEReweightingWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **LHEReweightingWeight** | Float_t|  |
| **nLHEReweightingWeight** | Int_t|  |

### <a id='lhescaleweight'></a>LHEScaleWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **LHEScaleWeight** | Float_t| LHE scale variation weights (w_var / w_nominal); [0] is renscfact=0.5d0 facscfact=0.5d0 ; [1] is renscfact=0.5d0 facscfact=1d0 ; [2] is renscfact=0.5d0 facscfact=2d0 ; [3] is renscfact=1d0 facscfact=0.5d0 ; [4] is renscfact=1d0 facscfact=1d0 ; [5] is renscfact=1d0 facscfact=2d0 ; [6] is renscfact=2d0 facscfact=0.5d0 ; [7] is renscfact=2d0 facscfact=1d0 ; [8] is renscfact=2d0 facscfact=2d0  |
| **nLHEScaleWeight** | Int_t|  |

### <a id='met'></a>MET [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **MET_pt** | Float_t| pt |

### <a id='muon'></a>Muon [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Muon_ZZFullId** | Bool_t| pass H4l full ID selection |
| **Muon_ZZFullSel** | Bool_t| pass H4l full SR selection (FullID + isolation) |
| **Muon_ZZFullSelNoSIP** | Bool_t| pass H4l full ID without SIP (base for CR SIP  method) |
| **Muon_ZZRelaxedId** | Bool_t| pass H4l relaxed ID including SIP (base for SS and OS CRs) |
| **Muon_ZZRelaxedIdNoSIP** | Bool_t| pass H4l relaxed ID without SIP (widest subset of cuts commont to all CRs) |
| **Muon_bsConstrainedChi2** | Float_t| chi2 of beamspot constraint |
| **Muon_bsConstrainedPt** | Float_t| pT with beamspot constraint |
| **Muon_bsConstrainedPtErr** | Float_t| pT error with beamspot constraint  |
| **Muon_charge** | Int_t| electric charge |
| **Muon_dxy** | Float_t| dxy (with sign) wrt first PV, in cm |
| **Muon_dxyErr** | Float_t| dxy uncertainty, in cm |
| **Muon_dxybs** | Float_t| dxy (with sign) wrt the beam spot, in cm |
| **Muon_dz** | Float_t| dz (with sign) wrt first PV, in cm |
| **Muon_dzErr** | Float_t| dz uncertainty, in cm |
| **Muon_eta** | Float_t| eta |
| **Muon_fsrPhotonIdx** | Short_t(index to Fsrphoton)| Index of the lowest-dR/ET2 among associated FSR photons |
| **Muon_genPartFlav** | UChar_t| Flavour of genParticle (DressedLeptons for electrons) for MC matching to status==1 muons: 1 = prompt muon (including gamma*->mu mu), 15 = muon from prompt tau, 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched |
| **Muon_genPartIdx** | Short_t(index to Genpart)| Index into genParticle list for MC matching to status==1 muons |
| **Muon_highPtId** | UChar_t| high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT) |
| **Muon_highPurity** | Bool_t| inner track is high purity |
| **Muon_inTimeMuon** | Bool_t| inTimeMuon ID |
| **Muon_ip3d** | Float_t| 3D impact parameter wrt first PV, in cm |
| **Muon_isGlobal** | Bool_t| muon is global muon |
| **Muon_isPFcand** | Bool_t| muon is PF candidate |
| **Muon_isStandalone** | Bool_t| muon is a standalone muon |
| **Muon_isTracker** | Bool_t| muon is tracker muon |
| **Muon_jetIdx** | Short_t(index to Jet)| index of the associated jet (-1 if none) |
| **Muon_jetNDauCharged** | UChar_t| number of charged daughters of the closest jet |
| **Muon_jetPtRelv2** | Float_t| Relative momentum of the lepton with respect to the closest jet after subtracting the lepton |
| **Muon_jetRelIso** | Float_t| Relative isolation in matched jet (1/ptRatio-1, pfRelIso04_all if no matched jet) |
| **Muon_looseId** | Bool_t| muon is loose muon |
| **Muon_mass** | Float_t| mass |
| **Muon_mediumId** | Bool_t| cut-based ID, medium WP |
| **Muon_mediumPromptId** | Bool_t| cut-based ID, medium prompt WP |
| **Muon_miniIsoId** | UChar_t| MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight) |
| **Muon_miniPFRelIso_all** | Float_t| mini PF relative isolation, total (with scaled rho*EA PU corrections) |
| **Muon_miniPFRelIso_chg** | Float_t| mini PF relative isolation, charged component |
| **Muon_multiIsoId** | UChar_t| MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium) |
| **Muon_mvaLowPt** | Float_t| Low pt muon ID score |
| **Muon_mvaMuID** | Float_t| MVA-based ID score  |
| **Muon_mvaMuID_WP** | UChar_t| MVA-based ID selector WPs (1=MVAIDwpMedium,2=MVAIDwpTight) |
| **Muon_mvaTTH** | Float_t| TTH MVA lepton ID score |
| **Muon_nStations** | UChar_t| number of matched stations with default arbitration (segment & track) |
| **Muon_nTrackerLayers** | UChar_t| number of layers in the tracker |
| **Muon_passIso** | Bool_t| Pass ZZ isolation cut |
| **Muon_pdgId** | Int_t| PDG code assigned by the event reconstruction (not by MC truth) |
| **Muon_pfIsoId** | UChar_t| PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight) |
| **Muon_pfRelIso03FsrCorr** | Float_t| FSR-subtracted pfRelIso03 |
| **Muon_pfRelIso03_all** | Float_t| PF relative isolation dR=0.3, total (deltaBeta corrections) |
| **Muon_pfRelIso03_chg** | Float_t| PF relative isolation dR=0.3, charged component |
| **Muon_pfRelIso04_all** | Float_t| PF relative isolation dR=0.4, total (deltaBeta corrections) |
| **Muon_phi** | Float_t| phi |
| **Muon_pt** | Float_t| pt |
| **Muon_ptErr** | Float_t| ptError of the muon track |
| **Muon_puppiIsoId** | UChar_t| PuppiIsoId from miniAOD selector (1=Loose, 2=Medium, 3=Tight) |
| **Muon_segmentComp** | Float_t| muon segment compatibility |
| **Muon_sip3d** | Float_t| 3D impact parameter significance wrt first PV |
| **Muon_softId** | Bool_t| soft cut-based ID |
| **Muon_softMva** | Float_t| soft MVA ID score |
| **Muon_softMvaId** | Bool_t| soft MVA ID |
| **Muon_svIdx** | Short_t(index to Sv)| index of matching secondary vertex |
| **Muon_tightCharge** | UChar_t| Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass) |
| **Muon_tightId** | Bool_t| cut-based ID, tight WP |
| **Muon_tkIsoId** | UChar_t| TkIso ID (1=TkIsoLoose, 2=TkIsoTight) |
| **Muon_tkRelIso** | Float_t| Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt |
| **Muon_triggerIdLoose** | Bool_t| TriggerIdLoose ID |
| **Muon_tunepRelPt** | Float_t| TuneP relative pt, tunePpt/pt |
| **nMuon** | Int_t| slimmedMuons after basic selection (pt > 15 \|\| (pt > 3 && (passed("CutBasedIdLoose") \|\| passed("SoftCutBasedId") \|\| passed("SoftMvaId") \|\| passed("CutBasedIdGlobalHighPt") \|\| passed("CutBasedIdTrkHighPt")))) |

### <a id='psweight'></a>PSWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **PSWeight** | Float_t| PS weights (w_var / w_nominal);   [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5; |
| **nPSWeight** | Int_t|  |

### <a id='pileup'></a>Pileup [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Pileup_gpudensity** | Float_t| Generator-level PU vertices / mm |
| **Pileup_nPU** | Int_t| the number of pileup interactions that have been added to the event in the current bunch crossing |
| **Pileup_nTrueInt** | Float_t| the true mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled |
| **Pileup_pudensity** | Float_t| PU vertices / mm |
| **Pileup_sumEOOT** | Int_t| number of early out of time pileup |
| **Pileup_sumLOOT** | Int_t| number of late out of time pileup |

### <a id='zcand'></a>ZCand [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZCand_eta** | Float_t| ZCand_eta[nZCand]/F |
| **ZCand_flav** | Float_t| Product of the pdgIds of the 2 daughters |
| **ZCand_fsr1Idx** | Short_t(index to Fsr1)| index of FSR associated to l1 (-1 if none) |
| **ZCand_fsr2Idx** | Short_t(index to Fsr2)| index of FSR associated to l2 (-1 if none) |
| **ZCand_l1Idx** | Short_t(index to L1)| index of 1st daughter in Electron+Muon merged collection |
| **ZCand_l2Idx** | Short_t(index to L2)| index of 2nd daughter in Electron+Muon merged collection |
| **ZCand_mass** | Float_t| mass |
| **ZCand_phi** | Float_t| ZCand_phi[nZCand]/F |
| **ZCand_pt** | Float_t| ZCand_pt[nZCand]/F |
| **ZCand_rapidity** | Float_t| ZCand_rapidity[nZCand]/F |
| **nZCand** | Int_t| Z candidates passing the full H4l selection |

### <a id='zlcand'></a>ZLCand [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZLCand_lepIdx** | Short_t(index to Lep)| Index of extra lep for the ZL CR |

### <a id='zllcand'></a>ZLLCand [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZLLCand_KD** | Float_t| ZLLCand_KD[nZLLCand]/F |
| **ZLLCand_Z1flav** | Int_t| ZLLCand_Z1flav[nZLLCand]/I |
| **ZLLCand_Z1l1Idx** | Short_t(index to Z1L1)| ZLLCand_Z1l1Idx[nZLLCand]/S |
| **ZLLCand_Z1l2Idx** | Short_t(index to Z1L2)| ZLLCand_Z1l2Idx[nZLLCand]/S |
| **ZLLCand_Z1mass** | Float_t| ZLLCand_Z1mass[nZLLCand]/F |
| **ZLLCand_Z2flav** | Short_t| ZLLCand_Z2flav[nZLLCand]/S |
| **ZLLCand_Z2l1Idx** | Short_t(index to Z2L1)| ZLLCand_Z2l1Idx[nZLLCand]/S |
| **ZLLCand_Z2l2Idx** | Short_t(index to Z2L2)| ZLLCand_Z2l2Idx[nZLLCand]/S |
| **ZLLCand_Z2mass** | Float_t| ZLLCand_Z2mass[nZLLCand]/F |
| **ZLLCand_eta** | Float_t| ZLLCand_eta[nZLLCand]/F |
| **ZLLCand_mass** | Float_t| ZLLCand_mass[nZLLCand]/F |
| **ZLLCand_massPreFSR** | Float_t| ZLLCand_massPreFSR[nZLLCand]/F |
| **ZLLCand_phi** | Float_t| ZLLCand_phi[nZLLCand]/F |
| **ZLLCand_pt** | Float_t| ZLLCand_pt[nZLLCand]/F |
| **ZLLCand_rapidity** | Float_t| ZLLCand_rapidity[nZLLCand]/F |
| **nZLLCand** | Int_t| Z+LL control region candidates |

### <a id='zllbest2p2fidx'></a>ZLLbest2P2FIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZLLbest2P2FIdx** | Short_t(index to Zllbest2P2F)| best candidate for the 2P2F CR |

### <a id='zllbest3p1fidx'></a>ZLLbest3P1FIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZLLbest3P1FIdx** | Short_t(index to Zllbest3P1F)| best candidate for the 3P1F CR |

### <a id='zllbestsipcridx'></a>ZLLbestSIPCRIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZLLbestSIPCRIdx** | Short_t(index to Zllbestsipcr)| best candidate for the SIP CR |

### <a id='zllbestssidx'></a>ZLLbestSSIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZLLbestSSIdx** | Short_t(index to Zllbestss)| best candidate for the SS CR |

### <a id='zzcand'></a>ZZCand [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ZZCand_KD** | Float_t| Kinematic discriminant for the choice of best candidate |
| **ZZCand_Z1flav** | Int_t| Product of the pdgIds of the 2 Z1 daughters |
| **ZZCand_Z1l1Idx** | Short_t(index to Z1L1)| Index of 1st Z1 daughter in the Electron+Muon merged collection |
| **ZZCand_Z1l2Idx** | Short_t(index to Z1L2)| Index of 2nd Z1 daughter in the Electron+Muon merged collection |
| **ZZCand_Z1mass** | Float_t| Z1 mass |
| **ZZCand_Z2flav** | Int_t| Product of the pdgIds of the 2 Z2 daughters |
| **ZZCand_Z2l1Idx** | Short_t(index to Z2L1)| Index of 1st Z2 daughter in the Electron+Muon merged collection |
| **ZZCand_Z2l2Idx** | Short_t(index to Z2L2)| Index of 2nd Z2 daughter in the Electron+Muon merged collection |
| **ZZCand_Z2mass** | Float_t| Z2 mass |
| **ZZCand_Z2sumpt** | Float_t| sum of Z2 daughter pts (used in the choice of best candidate) |
| **ZZCand_dataMCWeight** | Float_t| data/MC efficiency correction weight |
| **ZZCand_eta** | Float_t| ZZCand_eta[nZZCand]/F |
| **ZZCand_mass** | Float_t| mass |
| **ZZCand_massPreFSR** | Float_t| mass without FSR photons |
| **ZZCand_nExtraLep** | Int_t| number of extra leptons passing H4l full sel |
| **ZZCand_nExtraZ** | Int_t| number of extra Zs passing H4l full sel |
| **ZZCand_phi** | Float_t| ZZCand_phi[nZZCand]/F |
| **ZZCand_pt** | Float_t| ZZCand_pt[nZZCand]/F |
| **ZZCand_rapidity** | Float_t| ZZCand_rapidity[nZZCand]/F |
| **nZZCand** | Int_t| ZZ candidates passing the full H4l selection |

### <a id='bestcandidx'></a>bestCandIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **bestCandIdx** | Short_t(index to Bestcand)| Seleced ZZ candidate in the event |

### <a id='bestzidx'></a>bestZIdx [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **bestZIdx** | Short_t(index to Bestz)| Best Z in the event (mass closest to mZ) |

### <a id='event'></a>event [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **event** | ULong64_t| event/l |

### <a id='genweight'></a>genWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **genWeight** | Float_t| generator weight |

### <a id='ggh'></a>ggH [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ggH_NNLOPS_Weight** | Float_t| Reweighting for ggH as a function of njets and pT |

### <a id='luminosityblock'></a>luminosityBlock [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **luminosityBlock** | UInt_t| luminosityBlock/i |

### <a id='ncleanedjetspt30'></a>nCleanedJetsPt30 [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **nCleanedJetsPt30** | Char_t| nCleanedJetsPt30/B |
| **nCleanedJetsPt30_jesDn** | Char_t| nCleanedJetsPt30_jesDn/B |
| **nCleanedJetsPt30_jesUp** | Char_t| nCleanedJetsPt30_jesUp/B |

### <a id='overalleventweight'></a>overallEventWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **overallEventWeight** | Float_t| overallEventWeight/F |

### <a id='passedfiducial'></a>passedFiducial [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **passedFiducial** | Char_t| event passes fiducial selection at gen level |

### <a id='puweight'></a>puWeight [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **puWeight** | Float_t| puWeight/F |

### <a id='run'></a>run [<sup>[back to top]</sup>](#events-tree-content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

## AllEvents tree content

| Collection | Description |
| - | - |
| [**FidDressedLeps**](#fiddressedleps) | gen dressed leps for fiducial analysis |
| [**FidZ**](#fidz) | Zs made of FidDressedLeps |
| [**FidZ1**](#fidz1) | FidZ1_mass/F |
| [**FidZ2**](#fidz2) | FidZ2_mass/F |
| [**FidZZ**](#fidzz) | Index of 1st Z1 daughter in FidDressedLeps collection |
| [**GenDressedLepton**](#gendressedlepton) | Dressed leptons from Rivet-based ParticleLevelProducer |
| [**Generator**](#generator) | MC generator weight |
| [**event**](#event) | event/l |
| [**ggH**](#ggh) | Reweighting for ggH as a function of njets and pT |
| [**luminosityBlock**](#luminosityblock) | luminosityBlock/i |
| [**overallEventWeight**](#overalleventweight) | overallEventWeight/F |
| [**passedFiducial**](#passedfiducial) | event passes fiducial selection at gen level |
| [**puWeight**](#puweight) | puWeight/F |
| [**run**](#run) | run/i |

## AllEvents tree detail

### <a id='fiddressedleps'></a>FidDressedLeps [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidDressedLeps_RelIso** | Float_t| FidDressedLeps_RelIso[nFidDressedLeps]/F |
| **FidDressedLeps_eta** | Float_t| FidDressedLeps_eta[nFidDressedLeps]/F |
| **FidDressedLeps_id** | Float_t| FidDressedLeps_id[nFidDressedLeps]/F |
| **FidDressedLeps_mass** | Float_t| FidDressedLeps_mass[nFidDressedLeps]/F |
| **FidDressedLeps_momid** | Float_t| FidDressedLeps_momid[nFidDressedLeps]/F |
| **FidDressedLeps_mommomid** | Float_t| FidDressedLeps_mommomid[nFidDressedLeps]/F |
| **FidDressedLeps_phi** | Float_t| FidDressedLeps_phi[nFidDressedLeps]/F |
| **FidDressedLeps_pt** | Float_t| FidDressedLeps_pt[nFidDressedLeps]/F |
| **nFidDressedLeps** | Int_t| gen dressed leps for fiducial analysis |

### <a id='fidz'></a>FidZ [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZ_DauPdgId** | Int_t| FidZ decay flavour |
| **FidZ_MomPdgId** | Int_t| FidZ mother Id |
| **nFidZ** | Int_t| Zs made of FidDressedLeps |

### <a id='fidz1'></a>FidZ1 [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZ1_mass** | Float_t| FidZ1_mass/F |

### <a id='fidz2'></a>FidZ2 [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZ2_mass** | Float_t| FidZ2_mass/F |

### <a id='fidzz'></a>FidZZ [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **FidZZ_Z1l1Idx** | Short_t(index to Z1L1)| Index of 1st Z1 daughter in FidDressedLeps collection |
| **FidZZ_Z1l2Idx** | Short_t(index to Z1L2)| Index of 2nd Z1 daughter in FidDressedLeps collection |
| **FidZZ_Z2l1Idx** | Short_t(index to Z2L1)| Index of 1st Z2 daughter in FidDressedLeps collection |
| **FidZZ_Z2l2Idx** | Short_t(index to Z2L2)| Index of 2nd Z1 daughter in FidDressedLeps collection |
| **FidZZ_eta** | Float_t| FidZZ_eta/F |
| **FidZZ_mass** | Float_t| mass of gen ZZ made with FidDressedLeps |
| **FidZZ_phi** | Float_t| FidZZ_phi/F |
| **FidZZ_pt** | Float_t| FidZZ_pt/F |
| **FidZZ_rapidity** | Float_t| FidZZ_rapidity/F |

### <a id='gendressedlepton'></a>GenDressedLepton [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **GenDressedLepton_eta** | Float_t| eta |
| **GenDressedLepton_hasTauAnc** | Bool_t| true if Dressed lepton has a tau as ancestor |
| **GenDressedLepton_mass** | Float_t| mass |
| **GenDressedLepton_pdgId** | Int_t| PDG id |
| **GenDressedLepton_phi** | Float_t| phi |
| **GenDressedLepton_pt** | Float_t| pt |
| **nGenDressedLepton** | Int_t| Dressed leptons from Rivet-based ParticleLevelProducer |

### <a id='generator'></a>Generator [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **Generator_weight** | Float_t| MC generator weight |

### <a id='event'></a>event [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **event** | ULong64_t| event/l |

### <a id='ggh'></a>ggH [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **ggH_NNLOPS_Weight** | Float_t| Reweighting for ggH as a function of njets and pT |

### <a id='luminosityblock'></a>luminosityBlock [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **luminosityBlock** | UInt_t| luminosityBlock/i |

### <a id='overalleventweight'></a>overallEventWeight [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **overallEventWeight** | Float_t| overallEventWeight/F |

### <a id='passedfiducial'></a>passedFiducial [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **passedFiducial** | Char_t| event passes fiducial selection at gen level |

### <a id='puweight'></a>puWeight [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **puWeight** | Float_t| puWeight/F |

### <a id='run'></a>run [<sup>[back to top]</sup>](#allevents-tree-content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

## Runs tree content

| Collection | Description |
| - | - |
| [**LHEPdfSumw**](#lhepdfsumw) | Sum of genEventWeight * LHEPdfWeight[i], divided by genEventSumw |
| [**LHEScaleSumw**](#lhescalesumw) | Sum of genEventWeight * LHEScaleWeight[i], divided by genEventSumw |
| [**PSSumw**](#pssumw) | Sum of genEventWeight * PSWeight[i], divided by genEventSumw |
| [**genEventCount**](#geneventcount) | event count |
| [**genEventSumw**](#geneventsumw) | sum of gen weights |
| [**genEventSumw2**](#geneventsumw2) | sum of gen (weight^2) |
| [**run**](#run) | run/i |

## Runs tree detail

### <a id='lhepdfsumw'></a>LHEPdfSumw [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **LHEPdfSumw** | Double_t| Sum of genEventWeight * LHEPdfWeight[i], divided by genEventSumw |
| **nLHEPdfSumw** | Int_t| Number of entries in LHEPdfSumw |

### <a id='lhescalesumw'></a>LHEScaleSumw [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **LHEScaleSumw** | Double_t| Sum of genEventWeight * LHEScaleWeight[i], divided by genEventSumw |
| **nLHEScaleSumw** | Int_t| Number of entries in LHEScaleSumw |

### <a id='pssumw'></a>PSSumw [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **PSSumw** | Double_t| Sum of genEventWeight * PSWeight[i], divided by genEventSumw |
| **nPSSumw** | Int_t| Number of entries in PSSumw |

### <a id='geneventcount'></a>genEventCount [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **genEventCount** | Long64_t| event count |

### <a id='geneventsumw'></a>genEventSumw [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **genEventSumw** | Double_t| sum of gen weights |

### <a id='geneventsumw2'></a>genEventSumw2 [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **genEventSumw2** | Double_t| sum of gen (weight^2) |

### <a id='run'></a>run [<sup>[back to top]</sup>](#runs-tree-content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

## LuminosityBlocks tree content

| Collection | Description |
| - | - |
| [**GenFilter**](#genfilter) | generator filter: efficiency |
| [**luminosityBlock**](#luminosityblock) | luminosityBlock/i |
| [**run**](#run) | run/i |

## LuminosityBlocks tree detail

### <a id='genfilter'></a>GenFilter [<sup>[back to top]</sup>](#luminosityblocks-tree-content)
| Object property | Type | Description |
| - | - | - |
| **GenFilter_filterEfficiency** | Float_t| generator filter: efficiency |
| **GenFilter_filterEfficiencyError** | Float_t| generator filter: efficiency error |
| **GenFilter_numEventsPassed** | Int_t| generator filter: passed number of events |
| **GenFilter_numEventsTotal** | Int_t| generator filter: total number of events |

### <a id='luminosityblock'></a>luminosityBlock [<sup>[back to top]</sup>](#luminosityblocks-tree-content)
| Object property | Type | Description |
| - | - | - |
| **luminosityBlock** | UInt_t| luminosityBlock/i |

### <a id='run'></a>run [<sup>[back to top]</sup>](#luminosityblocks-tree-content)
| Object property | Type | Description |
| - | - | - |
| **run** | UInt_t| run/i |

