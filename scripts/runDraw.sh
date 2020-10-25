SECONDS=0
    HHDIR='/mnt/NVME/HH/'
    EXE='/build/draw'
    DATA='/output/sel_out/'
    OUT='/output/plot_out/'

    rm -rf $HHDIR/$OUT/*
    mkdir -p $HHDIR/$OUT
    make draw
    for entry in `ls $HHDIR$DATA$FLDIR`; do
        $HHDIR$EXE $HHDIR$DATA$FLDIR/$entry $HHDIR$OUT/$entry &
    done
    wait
    make stack
    $HHDIR/build/stack
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED