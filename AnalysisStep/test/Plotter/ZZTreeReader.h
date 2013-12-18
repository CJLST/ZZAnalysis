#include "../Macros/HZZ4lBase.h"

class ZZTreeReader: public HZZ4lBase {
  
public:
  ZZTreeReader(TChain *tree) : HZZ4lBase(tree){
  }

  ~ZZTreeReader(){};

  void Loop(bool useWeight=true);

};

#if defined(__MAKECINT__)  
#pragma link C++ class ZZTreeReader;
#endif
