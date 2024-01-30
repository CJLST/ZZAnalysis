#!/bin/tcsh

set nonomatch
set list = ( */*.root */*.corrupted */*.recovered */*.gz */*.txt */core* */jobid  */LSFJOB*/ */log/* */output/* */error/* */*.DAT */*.cc */br.sm?)

foreach f ( ${list} )
    if ( -e $f ) then
	rm -r $f
    endif
end

