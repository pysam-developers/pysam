#!/bin/sh -e

if test -x /usr/bin/dnf; then
    echo Installing prerequisites via dnf...
    dnf -y install epel-release
    dnf -y install zlib-devel bzip2-devel xz-devel curl-devel openssl-devel samtools bcftools htslib-tools

elif test -x /usr/bin/yum; then
    if yum -y install epel-release; then
        echo Installing prerequisites via yum...
        yum -y install zlib-devel bzip2-devel xz-devel curl-devel openssl-devel samtools bcftools htslib-tools
    else
        echo Installing non-test prerequisites via yum...
        yum -y install zlib-devel bzip2-devel xz-devel curl-devel openssl-devel
        emulate=yes
    fi

elif test -d /etc/dpkg; then
    echo Installing prerequisites via apt-get...
    apt-get update
    apt-get install -y --no-install-recommends --no-install-suggests libcurl4-openssl-dev libssl-dev zlib1g-dev libbz2-dev liblzma-dev samtools bcftools tabix

elif test -x /sbin/apk; then
    echo Installing non-test prerequisites via apk...
    apk update
    apk add zlib-dev bzip2-dev xz-dev curl-dev openssl-dev
    emulate=yes

elif test -x ${HOMEBREW_PREFIX-/usr/local}/bin/brew; then
    echo Installing prerequisites via brew...
    HOMEBREW_NO_AUTO_UPDATE=1 brew install -q samtools bcftools
    brew unlink xz || true

elif test -x /usr/sbin/pkg; then
    echo Installing prerequisites via pkg...
    pkg update
    pkg install -y bcftools gmake py311-cython py311-mypy py311-pytest samtools

elif test -x /usr/pkg/bin/pkgin; then
    echo Installing prerequisites via pkgin...
    pkgin update
    pkgin -y install bcftools gmake py314-cython py314-mypy py314-setuptools py314-test samtools

else
    echo No package manager detected
fi

if test -n "$emulate" && test $# -ge 2; then
    emulator=$1
    bindir=$2
    echo Creating symlinks to $emulator in $bindir...
    mkdir -p $bindir
    ln -s $emulator $bindir/samtools
    ln -s $emulator $bindir/bcftools
    ln -s $emulator $bindir/bgzip
    ln -s $emulator $bindir/tabix
fi
