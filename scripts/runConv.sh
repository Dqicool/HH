SECONDS=0
    HHDIR='/mnt/NVME/HH/'
    EXE='/build/convert'
    DATA='/data/'
    OUT='/output/convert_out/'

    # rm -rf $HHDIR/$OUT/*
    mkdir -p $HHDIR/$OUT
    make convert
    for entry in `ls $HHDIR$DATA`; do
        $HHDIR$EXE $HHDIR$DATA$entry $HHDIR$OUT/$entry.root &
    done
    wait
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED