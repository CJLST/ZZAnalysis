#ifndef CutSet_h
#define CutSet_h

/** \class CutSet
 *
 *  Convert a ParameterSet containing cut strings into a map<name, StringCutObjectSelector> for 
 *  the specified object type.
 *  
 *
 *  $Date: 2012/05/02 22:17:08 $
 *  $Revision: 1.1 $
 *  \author N. Amapane - Torino
 */
#include <map>
#include <string>
#include <CommonTools/Utils/interface/StringCutObjectSelector.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>


template <typename T>
class CutSet : public std::map<std::string,  StringCutObjectSelector<T, true>* >{
public:

  /// Constructor: parse the ParameterSet and create the selectors.
  CutSet(const edm::ParameterSet& ps) {
    std::vector<std::string> cutNames = ps.getParameterNamesForType<std::string>();
    for( unsigned i=0; i<cutNames.size(); ++i) {    
      (*this)[cutNames[i]] = new StringCutObjectSelector<T, true>(ps.getParameter<std::string>(cutNames[i]));
    }
  }

  /// Destructor: delete the selectors.
  ~CutSet() {
    for (typename CutSet<T>::const_iterator cut = this->begin(); cut != this->end(); ++cut) {
      delete (cut->second);
    }
  }
};

#endif
