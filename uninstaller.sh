#!/bin/sh

CDIR="$(pwd)"

if test -d /etc/KhufuEnv; then rm -rf /etc/KhufuEnv; fi
sed -i "/^source \/etc\/KhufuEnv\/call.sh$/d"  ~/.bashrc

echo "Successfully uninstalled; please reload ~/.bashrc"
