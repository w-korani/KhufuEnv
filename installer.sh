#!/bin/bash

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

check_R_pkg() {
    pkg="$1"
    Rscript -e "if (nzchar(system.file(package='$pkg'))) quit(status=0) else quit(status=1)" \
        > /dev/null 2>&1
}

if command -v Rscript >/dev/null 2>&1 ; then 
    echo "Checking required R packages..."

    req_pkgs="data.table ggplot2 Rcpp"
    for i in $req_pkgs; do 
        if check_R_pkg "$i"; then
            echo "$i is already installed."
        else 
            echo "$i package is missing."
            echo "Installing $i..."
            Rscript -e "install.packages('$i')"
        fi
    done

else
    echo "Warning: Rscript not found. Please install R. Skipping R package checks."
fi
