#ifndef Cosmetics_h
#define Cosmetics_h

// C++
#include <iostream>

// ROOT
#include "TColor.h"


using namespace std;

class Cosmetics
{
   
public:
   Cosmetics ();
   ~Cosmetics();
   
   struct VBF
   {
      int fill_color = TColor::GetColor("#ff9b9b");
      int line_color = TColor::GetColor("#cc0000");      
   };
   
   struct VH
   {
      int fill_color = TColor::GetColor("#ff9b9b");
      int line_color = TColor::GetColor("#cc0000");      
   };
   
   struct ttH
   {
      int fill_color = TColor::GetColor("#ff9b9b");
      int line_color = TColor::GetColor("#cc0000");      
   };
   
   struct Higgs_other
   {
      int fill_color = TColor::GetColor("#ffdcdc");  
      int line_color = TColor::GetColor("#cc0000");  
   };
   
   struct Higgs_all
   {
      int fill_color = TColor::GetColor("#ff9b9b");  
      int line_color = TColor::GetColor("#cc0000");  
   };
   
   struct qqZZ
   {
      int fill_color = TColor::GetColor("#99ccff");  
      int line_color = TColor::GetColor("#000099");  
   };
   
   struct ggZZ
   {
      int fill_color = TColor::GetColor("#4b78ff");  
      int line_color = TColor::GetColor("#000099");  
   };
   
   struct ZX
   {
      int fill_color = TColor::GetColor("#669966");  
      int line_color = TColor::GetColor("#003300");  
   };
   
private:

};

#endif