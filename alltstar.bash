#!/bin/bash

indir=/P/weisq/attentomo/tstarPS
outdir=$indir/alltstar
if (! [ -d $outdir ]);then
      mkdir $outdir
else
      rm $outdir/*
fi

cd $indir

for tdir in $(ls -d tstar*MPa*)
do
      para=$(echo $tdir | awk -F"_" '{print $2}')
      newpfl=$para"_p_tstar.dat"
      newsfl=$para"_s_tstar.dat"
      echo $newpfl
      cp $tdir/p_tstar.dat $outdir/$newpfl
      cp $tdir/s_tstar.dat $outdir/$newsfl
      newpfl=$para"_p_tstarhigh.dat"
      newsfl=$para"_s_tstarhigh.dat"
      cp $tdir/p_tstarhigh.dat $outdir/$newpfl
      cp $tdir/s_tstarhigh.dat $outdir/$newsfl
      newpfl=$para"_p_raw.dat"
      newsfl=$para"_s_raw.dat"
      cp $tdir/p_raw.dat $outdir/$newpfl
      cp $tdir/s_raw.dat $outdir/$newsfl
      newsfl=$para"_s_tstargood.dat"
      cp $tdir/s_tstargood.dat $outdir/$newsfl
      for subdb in 1 2 3 4 11 22 33 44
      do
            logfl=$(ls $tdir"/eventfocal"*"_sub"$subdb".log")
            newlogfl="eventfocal_"$para"_sub"$subdb".log"
            echo $newlogfl
            cp $logfl $outdir/$newlogfl
            newpfl=$para"_p_sub"$subdb"_err.dat"
            newsfl=$para"_s_sub"$subdb"_err.dat"
            cp $tdir"/p_sub"$subdb"_err.dat" $outdir/$newpfl
            cp $tdir"/s_sub"$subdb"_err.dat" $outdir/$newsfl
      done
done
