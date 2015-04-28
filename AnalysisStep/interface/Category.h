#ifndef CATEGORY_H
#define CATEGORY_H

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
