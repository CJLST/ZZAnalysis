#ifndef CCONSTANTS_H
#define CCONSTANTS_H

#include <TSpline.h>
#include <TString.h>
#include <memory>

class cConstantSpline {
  public:
    cConstantSpline(const TString& filename);
    void initspline(bool isDbkg);
    double eval(double ZZMass, bool isDbkg);
  private:
    const TString filename_;
    std::unique_ptr<TSpline3> spline_;
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

extern "C" float getDbkgVBFdecConstant(int ZZflav, float ZZMass);
extern "C" float getDbkgVHdecConstant(int ZZflav, float ZZMass);

extern "C" float getDbkgkinConstant(int ZZflav, float ZZMass);
extern "C" float getDbkgConstant(int ZZflav, float ZZMass);

#endif
