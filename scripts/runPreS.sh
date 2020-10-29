SECONDS=0
    HHDIR='/mnt/NVME/HH/'
    EXE='/build/presel'
    DATA='/output/01_convert_out/'
    OUT='/output/02_presel_out/'

    rm -rf $HHDIR/$OUT/*
    mkdir -p $HHDIR/$OUT
    make presel
    for entry in `ls $HHDIR$DATA$FLDIR`; do
        $HHDIR$EXE $HHDIR$DATA$FLDIR/$entry $HHDIR$OUT/$entry &
    done
    wait
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED