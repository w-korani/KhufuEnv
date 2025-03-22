# source /cluster/projects/khufu/korani_projects/KhufuEnv/05.KhufuEnvFastq.sh
############################################################################################################################################
###fastq processing
############################################################################################################################################
#Sec1#######################################################################################################################################
############################################################################################################################################
fastq(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq
else
   cat $fq
fi
}
############################
fastq2txt(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | paste -d"\t" - - - - 
else
   cat $fq | paste -d"\t" - - - - 
fi
}
############################
txt2fastq(){
txt=$1
cat $txt | awk '{print $1"\n"$2"\n"$3"\n"$4}' 
}
############################ 
fastqLen(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | awk '{if(NR%4==2) {print length($0)} }'
else
   cat $fq | awk '{if(NR%4==2) {print length($0)} }'
fi
}
############################
fastqAllLen(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | awk '{if(NR%4==2) {print length($0)} }' | sumCol
else
   cat $fq | awk '{if(NR%4==2) {print length($0)} }' | sumCol
fi
}
############################
fastqFilterLen(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | awk -v len=$2 '{ if(NR%4==1){S=$0}; if(NR%4==2 &&  length($0) >= len) {S=S"\n"$0; s=1}; if(NR%4==0){S=S"\n+\n"$0; if(s==1){print S; s=0} }; }'
else
   cat $fq | awk -v len=$2 '{ if(NR%4==1){S=$0}; if(NR%4==2 &&  length($0) >= len) {S=S"\n"$0; s=1}; if(NR%4==0){S=S"\n+\n"$0; if(s==1){print S; s=0} }; }'
fi
}
############################
fastq2fasta(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | awk '{if(NR%4==1) {print ">"substr($0,2) } else if(NR%4==2) {print $0}  }'
else
   cat $fq | awk '{if(NR%4==1) {print ">"substr($0,2) } else if(NR%4==2) {print $0}  }'
fi
}
############################
fasta2fastq(){
   fasta2fa $1 | awk '{if(NR%2==0) {seq=$0; Q=$0; gsub(".","~",Q); print id;print seq; print "+"; print Q } else {id=$0;gsub(">","@",id)} }'
}
############################
fastq2kmer(){
fq=$1
width=$2
move=$3
if [[ -z $width ]] ; then width=20; fi
if [[ -z $move ]]; then move=1; fi
fastq2txt $fq| sed "1i\ " | awk -v width=$width -v move=$move 'FS="\t" {X=1;for(x=1;x<=length($2);x+=move) { if(length(substr($2,x,width)) == width)  { printf "%s_%020d\n%s\n+\n%s\n", $1,X,substr($2,x,width),substr($4,x,width)  }; X+=1  } }' 
}
############################
fastqSubSampling(){
fq1=$1
fq2=$2
len=$3
paste <(fastq $fq1) <(fastq $fq2) | awk -v len=$len -v fq1=$(echo $fq1 | sed "s:.*/::g" ) -v fq2=$(echo $fq2 | sed "s:.*/::g" ) '{if(X<=len) { print $1 |  "gzip > sub_"fq1".gz" ; print $2 |  "gzip > sub_"fq2".gz" } else if(X>len){exit;} ; if(NR%4==0) { X+=length($1)+length($2) } } '
echo "sub_"$fq1" & sub_"$fq2" was generated"
}
############################
fastqAverage(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | awk '{if(NR%4==2) {print length($0)} }' | averageCol
else
   cat $fq | awk '{if(NR%4==2) {print length($0)} }' | averageCol
fi
}
############################ 
fastqStd(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | awk '{if(NR%4==2) {print length($0)} }' | stdCol
else
   cat $fq | awk '{if(NR%4==2) {print length($0)} }' | stdCol
fi
}
############################
fastqStat(){
fq=$1
if (file $fq | grep -q compressed ) ; then
   zcat $fq | "$khufu_dir"/utilities/phred_score_stat 
else
   cat $fq | "$khufu_dir"/utilities/phred_score_stat 
fi
}
############################
phred33(){
phred=33
for i in $(seq 0 $((126-phred)) ); do 
   echo $i | awk -v phred=$phred '{ printf("%s\t%s\t%c\t%f\n",$0,$0+phred, $0+phred, 10^(-($0/10)) ); }' 
done | sed "1iQ\tint\tchar\tp"
}
############################
phred64(){
phred=64
for i in $(seq 0 $((126-phred)) ); do 
   echo $i | awk -v phred=$phred '{ printf("%s\t%s\t%c\t%f\n",$0,$0+phred, $0+phred, 10^(-($0/10)) ); }' 
done | sed "1iQ\tint\tchar\tp"
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

