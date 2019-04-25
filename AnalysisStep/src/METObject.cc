#include <algorithm>
#include <utility>
#include <ZZAnalysis/AnalysisStep/interface/METObject.h>


METVariables::METVariables() :
  met_raw(0),
  phi_raw(0),
  met_original(0),
  phi_original(0),
  met(0),
  phi(0),
  met_METup(0),
  phi_METup(0),
  met_METdn(0),
  phi_METdn(0),
  met_JECup(0),
  phi_JECup(0),
  met_JECdn(0),
  phi_JECdn(0),
  met_JERup(0),
  phi_JERup(0),
  met_JERdn(0),
  phi_JERdn(0),
  met_PUup(0),
  phi_PUup(0),
  met_PUdn(0),
  phi_PUdn(0)
{}
METVariables::METVariables(METVariables const& other) :
  met_raw(other.met_raw),
  phi_raw(other.phi_raw),
  met_original(other.met_original),
  phi_original(other.phi_original),
  met(other.met),
  phi(other.phi),
  met_METup(other.met_METup),
  phi_METup(other.phi_METup),
  met_METdn(other.met_METdn),
  phi_METdn(other.phi_METdn),
  met_JECup(other.met_JECup),
  phi_JECup(other.phi_JECup),
  met_JECdn(other.met_JECdn),
  phi_JECdn(other.phi_JECdn),
  met_JERup(other.met_JERup),
  phi_JERup(other.phi_JERup),
  met_JERdn(other.met_JERdn),
  phi_JERdn(other.phi_JERdn),
  met_PUup(other.met_PUup),
  phi_PUup(other.phi_PUup),
  met_PUdn(other.met_PUdn),
  phi_PUdn(other.phi_PUdn)
{}
void METVariables::swap(METVariables& other){
  std::swap(met_raw, other.met_raw);
  std::swap(phi_raw, other.phi_raw);
  std::swap(met_original, other.met_original);
  std::swap(phi_original, other.phi_original);
  std::swap(met, other.met);
  std::swap(phi, other.phi);
  std::swap(met_METup, other.met_METup);
  std::swap(phi_METup, other.phi_METup);
  std::swap(met_METdn, other.met_METdn);
  std::swap(phi_METdn, other.phi_METdn);
  std::swap(met_JECup, other.met_JECup);
  std::swap(phi_JECup, other.phi_JECup);
  std::swap(met_JECdn, other.met_JECdn);
  std::swap(phi_JECdn, other.phi_JECdn);
  std::swap(met_JERup, other.met_JERup);
  std::swap(phi_JERup, other.phi_JERup);
  std::swap(met_JERdn, other.met_JERdn);
  std::swap(phi_JERdn, other.phi_JERdn);
  std::swap(met_PUup, other.met_PUup);
  std::swap(phi_PUup, other.phi_PUup);
  std::swap(met_PUdn, other.met_PUdn);
  std::swap(phi_PUdn, other.phi_PUdn);
}
METVariables& METVariables::operator=(const METVariables& other){
  METVariables tmp(other);
  swap(tmp);
  return *this;
}


METObject::METObject() :
  extras()
{}
METObject::METObject(const METObject& other) :
  extras(other.extras)
{}
void METObject::swap(METObject& other){
  extras.swap(other.extras);
}
METObject& METObject::operator=(const METObject& other){
  METObject tmp(other);
  swap(tmp);
  return *this;
}
METObject::~METObject(){}
