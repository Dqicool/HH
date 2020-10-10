#!/bin/bash
    HHDIR='/mnt/SSD/HH/'
    EXE='/build/selection'
    FLDIR='/fileLists/'
    DATA='/data/'
    OUT='/output/sel_out/'
    rm -rf $HHDIR/$OUT/*
    rm -rf $HHDIR/$DATA/$FLDIR/*
    mkdir -p $HHDIR/$OUT
    declare arr
    i=0
    j=0
    for dir in `ls $HHDIR$DATA`; do
        IFS='.' 
        read -a ADDR <<< "$dir"
        if [ ${ADDR[3]} = 'mc' ]; then
            IFS=$'\n'
            for entry in `ls $HHDIR$DATA$dir`; do
                echo $HHDIR$DATA$dir/$entry >> $HHDIR$DATA$FLDIR${ADDR[4]}.txt
            done
        else
            IFS=$'\n'
            for entry in `ls $HHDIR$DATA$dir`; do
                echo $HHDIR$DATA$dir/$entry >> $HHDIR$DATA$FLDIR${ADDR[3]}.txt
            done
        fi

    done
