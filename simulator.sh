#!/usr/bin/env bash
FUN=${1:-"run"}
shift

function run
{
    simulate $*
}

function set_anaconda
{
    if ! which conda >/dev/null 2>&1 ; then source ~/anaconda3/bin/activate ; fi
    if ! which conda >/dev/null 2>&1 ; then return 1 ; fi
}

function latest_file
{
    local ext=${1#.}
    [ x"$ext" != x ] && dotext=".$ext"
    ls -t *$dotext | head -1
}

function ask_latest_file
{
    filename=$(latest_file $*)
    if [ x"$filename" != x ] ; then
        ls -l $filename ; date ; echo
        echo "Is $filename the right file? [y/n]"
        read reply
        if [ $reply = Y ] || [ $reply = y ] || [ $reply = Yes ] || [ $reply = yes ] ; then
            LATEST_FILE=$filename
        fi
    fi
}


function simulate
{
    NAME=${1:-"input-1027p"}
    NAME=${NAME%%.csv}
    [ x"$2" != x ] && KVEL_ARG="--kvel $2"
    [ x"$3" != x ] && FUDGE_ARG="--fudge-pct $3"
    if ! set_anaconda ; then echo "can't find anaconda" >&2 ; return 1 ; fi
    python3 ocean_particle_simulator.py ${NAME}.csv --plot-all-lonlat $KVEL_ARG $FUDGE_ARG | tee results_${NAME}.log
}

function help
{
    echo "$0 function <...>"
}

function toZdrive
{
    local NAME=${1}
    [ x"$NAME" = x ] && ask_latest_file csv && NAME=$LATEST_FILE
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <name>" && exit 1
    NAME=${NAME%%".csv"}
    NAME=${NAME%%.mp4}

    local DST0='/Volumes/Anna/2021 Science Fair/results'
    if [ x"$2" != x ] ; then
        DST="$DST0/sample$2"
    else
        for i in {0..1000} ; do
            [ -d "$DST0/sample$i" ] && continue
            DST="$DST0/sample$i"
            break
        done
    fi
    [ x"DST" = x ] && echo "Cannot find dst" >&2 && exit 1

    echo "$DST"
    [ ! -d "$DST" ] && mkdir "$DST"
    cp ${NAME}* "$DST"
    for tile in {0..12} ; do
        if [ -f $tile/${NAME}.mp4 ] ; then
            cp $tile/${NAME}.mp4 "${DST}/${NAME}_t${tile}.mp4"
        fi
    done
}

function plot_tiles()
{
    NAME=${1}
    [ x"$NAME" = x ] && ask_latest_file csv && NAME=$LATEST_FILE
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <csvfile> " && exit 1
    NAME=${NAME%%.csv}
    # echo $NAME

    if ! set_anaconda ; then echo "cannot find anaconda" >&2 ; return 1 ; fi

    for i in {0..12} ; do
        [ ! -d $i ] && mkdir $i
        cp $NAME.csv $i/$NAME.csv
        python3 ocean_particle_simulator.py $i/$NAME.csv --only-plot --plot-1tile $i &
        # mv $i/$NAME.mp4 ${NAME}_t$i.mp4
    done

}

function plot_globe()
{
    local NAME=${1}
    [ x"$NAME" = x ] && ask_latest_file csv && NAME=$LATEST_FILE
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <name>" && exit 1
    NAME=${NAME%%".csv"}
    if ! set_anaconda ; then echo "cannot find anaconda" >&2 ; return 1 ; fi

    python3 ocean_particle_simulator.py $NAME.csv --only-plot --plot-all-lonlat
}

function plot_png_test()
{
    local NAME=${1}
    [ x"$NAME" = x ] && ask_latest_file csv && NAME=$LATEST_FILE
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <name>" && exit 1
    NAME=${NAME%%".csv"}
    if ! set_anaconda ; then echo "cannot find anaconda" >&2 ; return 1 ; fi

    python3 ocean_particle_simulator.py $NAME.csv --only-plot --plot-all-lonlat --png-ym 2005:10
}

function tarit
{
    NAME=${1%%.csv}
    NAME=${NAME%%.mp4}
    [ x"$NAME" = x ] && echo "$0 tarit <name>" && exit 1

    for tile in {0..12} ; do
        mv $tile/${NAME}.mp4 ${NAME}_t${tile}.mp4
    done

    tar cvfz ${NAME}.tgz ${NAME}.mp4 ${NAME}.csv ${NAME}_t*.mp4
}

$FUN $*
