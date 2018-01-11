program BigGamma
implicit none
character(len=40) :: arg
real(8) :: EHat, output

if (command_argument_count().ne.1) then
  print *, "Wrong number of arguments", command_argument_count()
  stop 1
endif

call GET_COMMAND_ARGUMENT(1,arg)
read(arg, *) EHat

call HTO_gridHt(EHat,output)
print *, output


end program BigGamma
