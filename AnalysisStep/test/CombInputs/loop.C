#include "TROOT.h"
#include "TSystem.h"

#include <iostream>

void loop(){
	for(int i=0; i<3; i++)
		for(int j=0; j<6; j++)
			all(1,i,j,1);
}
void loop();
