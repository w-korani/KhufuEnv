khufu_dir="$HOME/bin/KhufuEnv"
source ""$khufu_dir"/01.KhufuEnvHapmap.sh"
source ""$khufu_dir"/02.KhufuEnvPanmap.sh"
source ""$khufu_dir"/03.KhufuEnvStats.sh"
source ""$khufu_dir"/04.KhufuEnvFasta.sh"
source ""$khufu_dir"/05.KhufuEnvFastq.sh"
source ""$khufu_dir"/06.KhufuEnvDataSet.sh"
PATH=$PATH:"$khufu_dir"
export PATH
KhufuEnvHelp(){ fun=$1; if [ -z $fun ]; then echo "";  echo "KhufuEnv.Ver1.0.0";  cat "$khufu_dir"/[0-9][0-9].KhufuEnv*.txt | awk 'BEGIN{X=1}{if($0 ~ "^[A-Za-z]+" ){printf "%03d: %s\n", X,$0; X+=1}else{print $0} }'  ; echo "" ; 
else if ! grep -q "###BEGIN_${fun}###" "$khufu_dir"/[0-9][0-9].KhufuEnv*.doc; then
    echo ""
    echo "Error: '$fun' not found. Tool does not exist."
    echo ""
    return
fi 
cat "$khufu_dir"/[0-9][0-9].KhufuEnv*.doc | awk -v fun=$fun 'BEGIN{x=0;print "\n#####################\nKhufuEnv.Ver1.0.0: "fun} {if(x==1 && $0!="##END_"fun"##" ){print $0}; if($0=="###BEGIN_"fun"###"){x=1}; if($0=="##END_"fun"##"){x=0} }END{print "#####################\n"}'  | tr '\t' ' '; fi; }
