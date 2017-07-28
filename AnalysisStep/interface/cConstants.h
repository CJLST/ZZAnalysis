#ifndef CCONSTANTS_H
#define CCONSTANTS_H

#include "TSpline3.h"
#include "TString.h"

class cConstantSpline {
  public:
    cConstantSpline(const TString& filename);
    void initspline();
    double eval(double ZZMass);
  private:
    const TString filename_;
    unique_ptr<TFile> f_;
    TSpline3* spline_;
};

extern "C" float getDVBF2jetsConstant(float ZZMass);
extern "C" float getDVBF1jetConstant(float ZZMass);
extern "C" float getDWHhConstant(float ZZMass);
extern "C" float getDZHhConstant(float ZZMass);

extern "C" float getDVBF2jetsWP(float ZZMass, bool useQGTagging);
extern "C" float getDVBF1jetWP(float ZZMass, bool useQGTagging);
extern "C" float getDWHhWP(float ZZMass, bool useQGTagging);
extern "C" float getDZHhWP(float ZZMass, bool useQGTagging);

extern "C" float getDVBF2jetsConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP);
extern "C" float getDVBF1jetConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP);
extern "C" float getDWHhConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP);
extern "C" float getDZHhConstant_shiftWP(float ZZMass, bool useQGTagging, float newWP);

extern "C" float getDbkgkinConstant(int ZZflav, float ZZMass);
extern "C" float getDbkgConstant(int ZZflav, float ZZMass);

#endif
