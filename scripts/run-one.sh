#!/bin/bash
if [ $# -lt 2 ]; then
   echo "Error: model-name and directory not specified on command line"
   echo "Usage: model path"
   echo " e.g.: run-one.sh issams r2/Q0"
   exit 1
fi
WD=$PWD
echo $WD
#
cd $HOME/xssa"/run-"$1/$2
echo "Entering "$PWD
#$WD/ADMB/admb-clean.sh $1  
$HOME/xssa/ADMB/admb-clean.sh $1  
$1 -noinit -shess -iprint 1 -nr 10 &> $1.out
#ls -lrt
cd $WD

