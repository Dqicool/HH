#!/bin/bash
    HHDIR=$PWD
    DATADIR=/$1
    OUT=''
    echoerr() { echo "$@" 1>&2; }

    for dir in `ls $HHDIR$DATADIR`; do
        IFS='.' 
        read -a PIECES <<< "$dir"
        if [ ${PIECES[3]} == 'mc' ]; then
            IFS='_'
            read -a PIECES2 <<< "${PIECES[5]}"
            # echo ${PIECES2[0]}
            IFS=$'\n'
            if [ ${PIECES2[0]} = 'Sh221' ]; then
                OUT='Sherpa_221'
                if [ ${PIECES2[1]} = 'PDF30' ]; then
                    OUT=${OUT}_NNPDF30NNLO
                    for (( i=2; i<${#PIECES2[@]}; i++ )); do
                        if [ ${PIECES2[$i]} = 'Ztt' ]; then 
                            PIECES2[$i]=Ztautau
                        fi
                        if [[ "${PIECES2[$i]}" == *"MV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/MV/MAXHTPTV}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CV/CVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BV/BVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BF/BFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CF/CFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"ttb"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/ttb/ttbar_hdamp258p75}
                        fi
                        if [[ "${PIECES2[$i]}" == *"nonallh"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/nonallh/nonallhad}
                        fi
                        if [[ "${PIECES2[$i]}" == *"atop"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/atop/antitop}
                        fi
                        if [[ "${PIECES2[$i]}" == "st" ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/st/singletop}
                        fi
                        OUT=${OUT}_${PIECES2[$i]}
                    done
                else 
                    echoerr NOMATCH
                fi
            elif [ ${PIECES2[0]} = 'Sh222' ]; then
                OUT='Sherpa_222'
                if [ ${PIECES2[1]} = 'PDF30' ]; then
                    OUT=${OUT}_NNPDF30NNLO
                    for (( i=2; i<${#PIECES2[@]}; i++ )); do  
                        if [ ${PIECES2[$i]} = 'Ztt' ]; then 
                            PIECES2[$i]=Ztautau
                        fi
                        if [[ "${PIECES2[$i]}" == *"MV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/MV/MAXHTPTV}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CV/CVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BV/BVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BF/BFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CF/CFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"ttb"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/ttb/ttbar_hdamp258p75}
                        fi
                        if [[ "${PIECES2[$i]}" == *"nonallh"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/nonallh/nonallhad}
                        fi
                        if [[ "${PIECES2[$i]}" == *"atop"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/atop/antitop}
                        fi
                        if [[ "${PIECES2[$i]}" == "st" ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/st/singletop}
                        fi
                        OUT=${OUT}_${PIECES2[$i]}
                    done
                else 
                    echoerr NOMATCH
                fi
            elif [ ${PIECES2[0]} = 'PoPy8' ]; then
                OUT='PowhegPythia8EvtGen'
                if [ ${PIECES2[1]} = 'A14' ]; then
                    OUT=${OUT}_A14
                    for (( i=2; i<${#PIECES2[@]}; i++ )); do  
                        if [ ${PIECES2[$i]} = 'Ztt' ]; then 
                            PIECES2[$i]=Ztautau
                        fi
                        if [[ "${PIECES2[$i]}" == *"MV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/MV/MAXHTPTV}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CV/CVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BV/BVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BF/BFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CF/CFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"ttb"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/ttb/ttbar_hdamp258p75}
                        fi
                        if [[ "${PIECES2[$i]}" == *"nonallh"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/nonallh/nonallhad}
                        fi
                        if [[ "${PIECES2[$i]}" == *"atop"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/atop/antitop}
                        fi
                        if [[ "${PIECES2[$i]}" == "st" ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/st/singletop}
                        fi
                        OUT=${OUT}_${PIECES2[$i]}
                    done
                else
                    OUT=${OUT}_AZNLOCTEQ6L1
                    for (( i=1; i<${#PIECES2[@]}; i++ )); do  
                        if [ ${PIECES2[$i]} = 'Ztt' ]; then 
                            PIECES2[$i]=Ztautau
                        fi
                        if [[ "${PIECES2[$i]}" == *"MV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/MV/MAXHTPTV}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CV/CVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BV/BVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BF/BFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CF/CFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"ttb"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/ttb/ttbar_hdamp258p75}
                        fi
                        if [[ "${PIECES2[$i]}" == *"nonallh"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/nonallh/nonallhad}
                        fi
                        if [[ "${PIECES2[$i]}" == *"atop"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/atop/antitop}
                        fi
                        if [[ "${PIECES2[$i]}" == "st" ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/st/singletop}
                        fi
                        OUT=${OUT}_${PIECES2[$i]}
                    done
                fi
            elif [ ${PIECES2[0]} = 'PhPy8' ]; then
                OUT='PhPy8EG'
                if [ ${PIECES2[1]} = 'A14' ]; then
                    OUT=${OUT}_A14
                    for (( i=2; i<${#PIECES2[@]}; i++ )); do  
                        if [ ${PIECES2[$i]} = 'Ztt' ]; then 
                            PIECES2[$i]=Ztautau
                        fi
                        if [[ "${PIECES2[$i]}" == *"MV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/MV/MAXHTPTV}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CV/CVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BV"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BV/BVeto}
                        fi
                        if [[ "${PIECES2[$i]}" == *"BF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/BF/BFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"CF"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/CF/CFilter}
                        fi
                        if [[ "${PIECES2[$i]}" == *"ttb"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/ttb/ttbar_hdamp258p75}
                        fi
                        if [[ "${PIECES2[$i]}" == *"nonallh"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/nonallh/nonallhad}
                        fi
                        if [[ "${PIECES2[$i]}" == *"atop"* ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/atop/antitop}
                        fi
                        if [[ "${PIECES2[$i]}" == "st" ]]; then 
                            TMP=${PIECES2[$i]}
                            PIECES2[$i]=${TMP/st/singletop}
                        fi
                        OUT=${OUT}_${PIECES2[$i]}
                    done                    
                else 
                    echoerr NOMATCH
                fi
            else 
                echoerr NOMATCH
            fi
            # echo $OUT
            mv $HHDIR/$DATADIR/$dir $HHDIR/$DATADIR/qidong.v24.mc16_13TeV.${PIECES[4]}.$OUT.deriv.DAOD_STDM4.${PIECES[7]}
        else
            IFS=$'\n'
            mv $HHDIR/$DATADIR/$dir $HHDIR/$DATADIR/qidong.v24.${PIECES[3]}.${PIECES[4]}.${PIECES[5]}.PhysCont.DAOD_STDM4.${PIECES[7]} >> $HHDIR/$DATADIR/$dir/meta.txt
        fi
    done
