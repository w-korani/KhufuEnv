#!/bin/sh

CDIR="$(pwd)"

if test -d ~/bin/KhufuEnv; then 
rm -rf ~/bin/KhufuEnv
echo "Successfully uninstalled; please reload ~/.bashrc"
else
echo "KhufuEnv is already uninstalled"
fi
