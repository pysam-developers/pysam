#!/bin/sh -e

case $1 in
    ubuntu-*-arm)  expected=arm ;;
    macos-13)      expected=x86_64 ;;
    ubuntu-*)      expected=x86_64 ;;
    macos-*)       expected=arm ;;
    windows-*)     expected=x86_64 ;;
    *)
	echo Unknown platform $1 >&2
	exit 2
	;;
esac

arch=$(uname -m)
case $arch in
    arm*|aarch*)  actual=arm ;;
    x86*)         actual=x86_64 ;;
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
