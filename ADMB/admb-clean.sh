#!/bin/bash
if [ $# -lt 1 ]; then
   echo "Error: model name not specifid on command line"
   echo "Usage: admb-clean model"
   echo "       for example, admb-clean issams-msy"
   exit 1
fi

rm -fv $1.cor
rm -fv $1.eva
rm -fv $1.log
rm -fv $1.luu
rm -fv $1.out
rm -fv $1.p??
rm -fv $1.b??
rm -fv $1.r??
rm -fv $1.rhes
rm -fv $1.std
rm -fv hessian.diag
rm -fv admodel.cov
rm -fv fmin.log
rm -fv nf1b2list13
rm -fv nf1b2list12
rm -fv nf1b2list1
rm -fv f1b2list13
rm -fv f1b2list12
rm -fv f1b2list1
rm -fv cmpdiff.tmp
rm -fv gradfil1.tmp
rm -fv gradfil2.tmp
rm -fv varssave.tmp
rm -fv hessian.bin
rm -fv hesscheck
rm -fv admodel.hes
rm -fv admodel.dep
