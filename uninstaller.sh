#!/bin/sh

CDIR="$(pwd)"

if test -d /etc/KhufuEnv; then 
rm -rf /etc/KhufuEnv
echo "Successfully uninstalled; please reload ~/.bashrc"
else
echo "KhufuEnv is already uninstalled"
fi
