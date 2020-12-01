#!/bin/bash
    HHDIR=$PWD
    DATADIR=/$1

    echoerr() { echo "$@" 1>&2; }

    for dir in `ls $HHDIR$DATADIR`; do
        IFS='.' 
        read -a PIECES <<< "$dir"
        IFS=$'\n'
        OUT=''
        for (( i=2; i<${#PIECES[@]}; i++ )); do  
            OUT=$OUT.${PIECES[i]}
        done
        ami show dataset info ${OUT/./} >> $HHDIR$DATADIR/$dir/meta.txt
    done
