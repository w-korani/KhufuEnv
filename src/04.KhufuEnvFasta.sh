# source /cluster/projects/khufu/korani_projects/KhufuEnv/04.KhufuEnvFasta.sh
############################################################################################################################################
###fasta processing
############################################################################################################################################
#Sec1#######################################################################################################################################
############################################################################################################################################
fasta2txt(){
	cat $1 | awk '{if($0~"^>") { if(NR>1) {print id"\t"seq}; id=$0; gsub(" ",";",id) ;seq=""} else {seq=seq$0}  }END {print id"\t"seq}' | sed "s:^>::g"
}
############################
txt2fasta(){
	width=$2
	if [[ -z $width ]] ; then width=120; fi
	cat $1  | awk -v width=$width '{print ">"$1; if(width==""){width=120} ;for(x=1;x<=length($2);x+=width) {print substr($2,x,width)}  }'
}
############################
fastaUnWrap (){
	cat $1 |awk '{if($0~"^>") { if(NR>1) {print id"\n"seq}; id=$0;seq=""} else {seq=seq$0}  }END {print id"\n"seq}'
}
############################
fastaWrap(){
	width=$2
	if [[ -z $width ]] ; then width=120; fi
	cat $1  | awk -v width=$width '{if($0~"^>") {print $0} else {if(width==""){width=120} ; for(x=1;x<=length($0);x+=width) {print substr($0,x,width)}  }  }'
}
############################
fastaSeqLen(){
   cat $1 |awk '{if($0~"^>") { if(NR>1) {print id"\t"length(seq)}; id=$0;seq=""} else {seq=seq$0}  }END {print id"\t"length(seq)}' | sed "s:^>::g"
}
############################
fastaLen(){
fastaSeqLen $1 | cut -f 2 | sumCol
}
############################
fasta2kmer(){
width=$2
move=$3
if [[ -z $width ]] ; then width=120; fi
if [[ -z $move ]]; then move=1; fi
cat $1 | fasta2fa | awk -v width=$width -v move=$move  '{if(NR%2==1){id=$0} ;if(NR%2==0) {X=1;for(x=1;x<=length($0);x+=move) {S=substr($0,x,width); if(length(S) == width) {printf "%s_%020d\n%s\n" , id, X, S ; X+=1}  } } }'
}
############################
fastaCombine(){
cat $1 |  fasta2txt | awk '{if($1==X) {S=S","$2} else { if(NR>1) print X"\t"S; S=$2} ; X=$1; chr=$1}END{print X"\t"S}' | sed "s:,::g"  | txt2fasta
}
############################
fastaAlleleCombine(){
 cat $1 | fasta2txt | awk '{if($1==X) {S=S","$2} else { if(NR>1) print X"\t"S; S=$2} ; X=$1; chr=$1}END{print X"\t"S}' |  txt2fasta
}
############################
fastatxtCombine(){
cat $1 | awk '{if($1==X) {S=S","$2} else { if(NR>1) print X"\t"S; S=$2} ; X=$1; chr=$1}END{print X"\t"S}' | sed "s:,::g" 
}
############################
fastatxtAlleleCombine(){
 cat $1  | awk '{if($1==X) {S=S","$2} else { if(NR>1) print X"\t"S; S=$2} ; X=$1; chr=$1}END{print X"\t"S}' 
}
############################
fastaSelID(){
tmpFile0401=$(mktemp "./KhufuEnviron.XXXXXXXXX")
lst=$2
cat $1 | fasta2txt > "$tmpFile0401"
merge $lst "$tmpFile0401"  | txt2fasta - 120
rm "$tmpFile0401"
trap "rm -f $tmpFile0401" EXIT
}
############################
fastaFilterLen(){
cat $1 | fasta2txt | awk -v len=$2 '{if(length($2) >= len) print $0}' | txt2fasta - 120
}
############################
fastaSplit(){
if [ -d "fastaCHRs" ]; then
   echo "ERROR: fastaCHRs already exists"
else
	mkdir fastaCHRs
	cat $1 | awk '{if($0~">"){id=$0; gsub(">","",id)}; print $0 > "fastaCHRs/"id".fa" }'
  echo "fastaCHRs was generated"
fi
}
###############################################################################################################################################


