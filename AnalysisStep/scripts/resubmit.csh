#!/bin/tcsh -f
foreach x (*Chunk*) 
 cd $x
 bsub -q 8nh < ./batchScript.sh
 cd -
end
