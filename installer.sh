#!/bin/sh

CDIR="$(pwd)"
mkdir -p ~/bin

if test -d ~/bin/KhufuEnv; then rm -rf ~/bin/KhufuEnv; fi

echo "copying over files"
mkdir -p ~/bin/KhufuEnv/
cp -RP $CDIR/src/* ~/bin/KhufuEnv/

# perms
chmod 555 ~/bin/KhufuEnv/*.sh

## C++ compiling
g++ ~/bin/KhufuEnv/utilities/phred_score_stat.cpp -o ~/bin/KhufuEnv/utilities/phred_score_stat

echo "Successfully installed"
