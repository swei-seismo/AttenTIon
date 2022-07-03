#!/bin/bash

for alpha in 0.27
do
for fcps in 1.0
do
      if [ $fcps == 1.0 ];then
            workdir="/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/tstar_fcp0.5-20MPa"$alpha"site"
      elif [ $fcps == 1.5 ];then
            workdir="/Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/tstar_fcs0.5-20MPa"$alpha"site"
      fi
      if ( ! [ -d $workdir ] );then
            mkdir $workdir
#      else
#            rm -r $workdir
#            mkdir $workdir
      fi
      for subdb in 4 1 11 2 22 3 33 44
#     for subdb in 3
      do
            if ([ $alpha == 0.27 ] && [ $fcps == 1.5 ] && [ $subdb == 4 ]);then
                  plotspec=0
            else
                  plotspec=0
            fi
#            echo $plotspec
            #rm $workdir"/eventfocal027_sub"$subdb".log"
            cd /Users/sowei/GoogleDriveMSU/Work/Lau/Qtomo/program
./main_tstar_MoPresspec.py<<EOF
$subdb $workdir $alpha $fcps $plotspec
EOF
      done

      cd $workdir
      cat result_sub*/*_pstar*.dat > p_raw.dat
      cat result_sub*/*_sstar*.dat > s_raw.dat
      awk '{if ($5 > 0) print $0}' p_raw.dat > p_tstarhigh.dat
      awk '{if ($5 > 0) print $0}' s_raw.dat > s_tstarhigh.dat
      awk '{if (($5 > 0) && ($8 < 30)) print $0}' p_raw.dat > p_tstar.dat
#      awk '{if (($5 > 0) && ($9 > 1 )) print $0}' s_raw.dat > s_tstar.dat
      awk '{if (($5 > 0) && ((($9/$10 >= 1) && ($10 > 1 ) && ($10 < 3)) || ($9/$10 < 1))) print $0}' s_raw.dat > s_tstar.dat
      awk '{if (($5 > 0) && ($10 > 1 ) && ($10 < 3)) print $0}' s_raw.dat > s_tstargood.dat
      for subdb in 1 11 2 22 3 33 4 44
      do
            cat result_sub$subdb/*_perr*.dat > "p_sub"$subdb"_err.dat"
            cat result_sub$subdb/*_serr*.dat > "s_sub"$subdb"_err.dat"
      done
done

done
