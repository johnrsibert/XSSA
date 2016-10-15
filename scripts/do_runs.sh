#!/bin/bash

do_issams()
{
   cd $d
   pwd
   issams -noinit -iprint 1 -shess &> issams.out
   cd ..
   pwd
}

DIRS="r2 r4 r0 r0-sdrprior"
echo $DIRS

for d in $DIRS
do
   do_issams
done

# cd r2
#  issams -noinit -iprint 1 -shess &> issams.out
#  cd ..
# cd r4
#  issams -noinit -iprint 1 -shess &> issams.out
#  cd ..
# cd r0
#  issams -noinit -iprint 1 -shess &> issams.out
#  cd ..
# cd r0-sdrprior
#  issams -noinit -iprint 1 -shess &> issams.out
#  cd ..
