#!/bin/sh -e

expected=$1

arch=$(uname -m)
case $arch in
    arm*|aarch*)  actual=arm ;;
    x86*)         actual=intel ;;
    *)
	echo Unrecognized uname result $arch >&2
	exit 2
	;;
esac

if test $actual = $expected
then
    echo Running on $arch as expected
    exit 0
else
    echo Platform $arch is not the expected $expected >&2
    exit 1
fi
