#!/bin/bash
#prints normalized spectrum dN(E)/dE * E^2
inidir=$PWD

if [ $# -lt 1 ]
then
    echo "usage ./makeSpecE2.sh file"
    exit 1
fi

if [ -f $1 ]
then
	ofile=$1.spec
	cat $1 | awk '!/^($|[[:space:]]*#)/{ printf("%.1f %g\n",log($1)/log(10)+0.05,$2)}' | sort -g -k 1 | awk '{ if($1!=PREV&&NR>1) printf("%g %g %d\n",PREV,SUM,SUMN); if($1!=PREV||NR==1) {SUM=0; SUMN=0;}  SUM=SUM+$2; SUMN=SUMN+1; PREV=$1; } END {printf("%g %g %d\n",PREV,SUM,SUMN)}' | awk 'NF==3 { E=10^($1-0.05); printf("%g %g %d\n", E, $2*E/(10^0.05-10^(-0.05)),$3)}'>$ofile
else
	echo $1 not found
	exit 1
fi

