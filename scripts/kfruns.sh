#!/bin/bash
rm -rfv KF*
start_dir=$(pwd)
echo "Starting directory: "$start_dir
#kf=1.0
kf="0.1 0.3 1.0 3.0 10.0 30.0 100.0 300.0"

for k in $kf
do
   run_dir="KF"$k
   echo "Run directory: "$run_dir
   mkdir -v $run_dir
   sed {s/YYY/$k/} issams.yyy > $run_dir/issams.dat
   cd $run_dir
   nice issams  -noinit -iprint 1 -nr 10 #&> issams.out
   echo "iasams exit code = "$?" for KF = "$k
   cd $start_dir
done
echo "Done!"
