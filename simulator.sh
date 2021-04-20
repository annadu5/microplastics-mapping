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

function simulate
{
    NAME=${1:-"input-1027p"}
    NAME=${NAME%%.csv}
    KVEL=${2:-"0.2"}
    FUDGE=${3:-"30"}
    if ! set_anaconda ; then echo "can't find anaconda" >&2 ; return 1 ; fi
    python3 ocean_particle_simulator.py ${NAME}.csv --plot-all-lonlat --kvel $KVEL --fudge-pct $FUDGE | tee results_${NAME}.log
}

function help
{
    echo "$0 function <...>"
}

function toZdrive
{
    NAME=${1%%.csv}
    NAME=${NAME%%.mp4}
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <name>" && exit 1
    DST='/Volumes/Anna/2021 Science Fair/results'
    for i in {0..100} ; do
        [ -d "$DST/sample$i" ] && continue
        DST="$DST/sample$i"
        echo "$DST"
        mkdir "$DST"
        cp ${NAME}* "$DST"
        break
    done
}

function plot_tiles()
{
    NAME=${1:-"results/results_input-1515p_f30kv0.25"}
    NAME=${1%%.csv}
    if ! set_anaconda ; then echo "can't find anaconda" >&2 ; return 1 ; fi
    for i in {0..12} ; do
        [ ! -d $i ] && mkdir $i
        cp $NAME.csv $i/$NAME.csv
        python3 ocean_particle_simulator.py $i/$NAME.csv --only-plot --plot-1tile $i &
        mv $i/$NAME.mp4 ${NAME}_t$i.mp4
    done

}

function plot_globe()
{
    local NAME=${1%%.csv}
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <name>" && exit 1
    if ! set_anaconda ; then echo "can't find anaconda" >&2 ; return 1 ; fi

    python3 ocean_particle_simulator.py $NAME.csv --only-plot --plot-all-lonlat
}

function plot_png_test()
{
    local NAME=${1%%.csv}
    [ x"$NAME" = x ] && echo "$0 ${FUNCNAME[0]} <name>" && exit 1
    if ! set_anaconda ; then echo "can't find anaconda" >&2 ; return 1 ; fi

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
