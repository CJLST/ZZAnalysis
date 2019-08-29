#!/bin/tcsh
# Search for missing root files
# parameter to be specified: mf = move failed jobs
#                            lf = link failed jobs (then you can run cleanup.csh; resubmit_Condor.csh in AAAFAIL)

#set echo

set opt=$1

if ( $#argv >= 2 ) then
  set gooddir=$2
else
  set gooddir=AAAOK
endif

if ( $#argv >= 3 ) then
  set baddir=$3
else
  set baddir=AAAFAIL
endif

foreach chunk ( *Chunk* )
 set fail="false"


 # Check that root file is existing and not empty
 set filename=${chunk}/ZZ4lAnalysis.root
 if ( ! -e $filename ) then
#   echo "Missing root file in " ${chunk}
   set fail="true"
 else if ( -z $filename ) then
#   echo "Empty file: " $filename
   set fail="true"
 endif

 # Check job exit status
 set exitStatus = 0
 if ( -es ${chunk}/exitStatus.txt ) then
   set exitStatus=`cat ${chunk}/exitStatus.txt`
   set fail="true"
 else if ( ! -e ${chunk}/exitStatus.txt ) then
   set fail="true"
 endif
 # Check for CPU time exceeded (this can also be reported explicitly with exit status 152, which is then catched above)

 set nonomatch
 set jobRep = ( ${chunk}/job_*.txt )
 if ( -e $jobRep[1] ) then
   if ( `grep -c -e "Exited with exit code 1" -e "CPU time limit exceeded" $jobRep[$#jobRep]` != 0 ) then
      set exitStatus=152
      set fail="true"
   endif
 endif

 set logFile = ( ${chunk}/log/*.log )
 if ( -e $logFile[1] ) then
   if ( `grep -c -e "The job attribute PeriodicRemove expression.*evaluated to TRUE" $logFile[$#logFile]` != 0 ) then
      set exitStatus=153
      set fail="true"
   else if ( `grep -c -e "Job was aborted by the user" $logFile[$#logFile]` != 0 ) then
      set exitStatus=154
      set fail="true"
   endif
 endif
 unset nonomatch

 # Archive succesful jobs, or report failure
 if ( $fail == "false" ) then
  mkdir -p $gooddir
  mv $chunk $gooddir/
 else
  set description=""
   if ( $exitStatus == 0 ) then
     echo $chunk ": still running (or unknown failure)"
   else
     if ( $exitStatus == 84 ) set description="(missing input file)"
     if ( $exitStatus == 134 ) set description="(Crashed)"
     if ( $exitStatus == 152 ) set description="(Exceeded CPU time)"
     if ( $exitStatus == 153 ) set description="(Condor crashed, see log file)"
     if ( $exitStatus == 154 ) set description="(You cancelled the job)"
    echo $chunk ": failed, exit status = " $exitStatus $description
   endif
   if ( $opt == "mf" && $exitStatus != 0 ) then
     mkdir -p $baddir
     mv $chunk $baddir/
   else if ( $opt == "lf" && $exitStatus != 0 ) then
     mkdir -p $baddir
     ln -s ../$chunk $baddir/
     ln -sf ../condor.sub $baddir
   endif
 endif
end
