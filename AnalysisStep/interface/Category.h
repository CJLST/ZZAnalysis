#ifndef CATEGORY_H
#define CATEGORY_H

enum Category {
  Untagged     = 0,
  OneJetTagged = 1,
  VBFTagged    = 2, 
  VHLeptTagged = 3, 
  VHHadrTagged = 4, 
  ttHTagged    = 5
};

//int category(
extern "C" int category(
	     int nExtraLeptons,
	     float ZZPt,
	     float ZZMass,
	     int nJets, 
	     int nBTaggedJets,
	     float* jetpt,
	     float* jeteta,
	     float* jetphi,
	     float* jetmass,
	     float Fisher
	     );

#endif
