#/bin/tcsh -f

echo "NCU hours:"
echo "Sample                       min    avg    max  njobs"
echo "-----------------------------------------------------"

foreach x ( `ls -1 ./*Chunk*/log*.gz | sed s/"_Chunk.*"// | sort -u` )

 set sum=0
 set count=0
 set min=99999999999
 set max=0
 foreach CPU ( `grep "NCU hours" ${x}*/*.txt | awk '{print $9}' | sed s/"("// | sort -g` )
#   echo $CPU
   @ sum = $sum + $CPU
   @ count++
   if ( $CPU > $max ) set max=$CPU
   if ( $CPU < $min ) set min=$CPU
 end


    printf "%-25s %6.1f %6.1f %6.1f %6i\n" $x `dc -e "2 k $min 3600 / p"` `dc -e "2 k $sum $count / 3600 / p"` `dc -e "2 k $max 3600 / p"` $count
# echo $x `dc -e "2 k $min 3600 / p"` `dc -e "2 k $sum $count / 3600 / p"` `dc -e "2 k $max 3600 / p"` $count

end

