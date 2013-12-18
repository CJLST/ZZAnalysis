#!/bin/tcsh -f

foreach step (Gen Skim Trigger BestCand CandMZ2 CandpT MAllComb M70 M100: 100_1Jet 100_2Jet _MELA With1FSR With2FSR pseudoMELA graviMELA VBF)    
 set mmmm=`grep -a 4mu_Nevt $1 | grep $step | awk '{print $NF}' | sed s/"\.000"//g`
 set eeee=`grep -a 4e_Nevt $1 | grep $step | awk '{print $NF}' | sed s/"\.000"//g`
 set eemm=`grep -a 2e2mu_Nevt $1 | grep $step | awk '{print $NF}' | sed s/"\.000"//g`
 @ tot = $eeee + $mmmm + $eemm

 if ( $step != "Gen" && $step != "Skim" && $step != "Trigger" ) then
    echo $step "    \t" $tot/$eeee/$mmmm/$eemm
 else
    echo $step "    \t" $eeee/$mmmm/$eemm
 endif
end

