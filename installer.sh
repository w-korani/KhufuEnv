#!/bin/sh

CDIR="$(pwd)"

if test -d /etc/KhufuEnv; then rm -rf /etc/KhufuEnv; fi

echo "copying over files"
mkdir -p /etc/KhufuEnv/
cp -RP $CDIR/src/* /etc/KhufuEnv/

# perms
chmod 555 /etc/KhufuEnv/*.sh

echo "Successfully installed"
