#ifndef METOBJECT_H
#define METOBJECT_H


class METVariables{
public:
  float met_raw;
  float phi_raw;
  float met_original;
  float phi_original;
  float met;
  float phi;
  float met_METup;
  float phi_METup;
  float met_METdn;
  float phi_METdn;
  float met_JECup;
  float phi_JECup;
  float met_JECdn;
  float phi_JECdn;
  float met_JERup;
  float phi_JERup;
  float met_JERdn;
  float phi_JERdn;
  float met_PUup;
  float phi_PUup;
  float met_PUdn;
  float phi_PUdn;

  METVariables();
  METVariables(METVariables const& other);
  METVariables& operator=(const METVariables& other);

  void swap(METVariables& other);

};

class METObject{
public:
  METVariables extras;

  METObject();
  METObject(const METObject& other);
  METObject& operator=(const METObject& other);
  ~METObject();

  void swap(METObject& other);

};

#endif
