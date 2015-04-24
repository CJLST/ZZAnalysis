#ifndef CATEGORY_H
#define CATEGORY_H

//int category(
extern "C" int category(
	     int nExtraLeptons,
	     float ZZPt,
	     float ZZMass,
	     int nJets, 
	     int nBTaggedJets,
	     float ptj1,
	     float ptj2,
	     float etaj1,
	     float etaj2,
	     float mjj,
	     float Fisher
	     );

#endif
