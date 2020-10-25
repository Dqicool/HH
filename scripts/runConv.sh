SECONDS=0
    HHDIR='/mnt/NVME/HH/'
    EXE='/build/convert'
    FLDIR='/fileLists/'
    DATA='/data/'
    OUT='/output/convert_out/'

    rm -rf $HHDIR/$OUT/*
    mkdir -p $HHDIR/$OUT
    make convert
    for entry in `ls $HHDIR$DATA$FLDIR`; do
        $HHDIR$EXE $HHDIR$DATA$FLDIR/$entry $HHDIR$OUT/$entry.root &
    done
    wait
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED