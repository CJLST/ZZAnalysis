#!/bin/tcsh -f

# Make the grid proxy available, if existing and valid
set proxy_valid=`voms-proxy-info --timeleft`
if ($proxy_valid > 10 ) then
   echo "GRID proxy found, validity: $proxy_valid s"
   if ($?X509_USER_PROXY) then 
    if ($X509_USER_PROXY != ~/x509up_u${uid}) cp $X509_USER_PROXY ~/x509up_u${uid}
   else if (-e /tmp/x509up_u${uid} ) then
     cp /tmp/x509up_u${uid} ~
   endif
else
   echo "Note: no valid GRID proxy found."
endif

foreach x (*Chunk*) 
 cd $x
 bsub -q 8nh < ./batchScript.sh
 cd -
end
