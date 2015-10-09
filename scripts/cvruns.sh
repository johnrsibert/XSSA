#!/bin/bash
rm -rfv CV*
start_dir=$(pwd)
echo "Starting directory: "$start_dir
#std="0.2"
std="0.05 0.1 0.15 0.201 0.3 0.4 0.5"

for s in $std
do
   run_dir="CV"$s
   echo "Run directory: "$run_dir
   mkdir -v $run_dir
   sed {s/XXX/$s/} issams.xxx > $run_dir/issams.dat
   cd $run_dir
   nice issams  -noinit -iprint 1 -nr 10 &> issams.out
   echo "iasams exit code = "$?" for CV = "$s
   cd $start_dir
done
echo "Done!"
