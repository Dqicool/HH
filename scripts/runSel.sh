SECONDS=0
    HHDIR='/mnt/NVME/HH/'
    EXE='/build/selection'
    DATA='/output/convert_out/'
    OUT='/output/sel_out/'

    rm -rf $HHDIR/$OUT/*
    mkdir -p $HHDIR/$OUT
    make selection
    for entry in `ls $HHDIR$DATA$FLDIR`; do
        $HHDIR$EXE $HHDIR$DATA$FLDIR/$entry $HHDIR$OUT/$entry &
    done
    wait
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED