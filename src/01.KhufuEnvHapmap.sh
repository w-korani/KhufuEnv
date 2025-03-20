# source /cluster/projects/khufu/korani_projects/KhufuEnv/01.KhufuEnvhapmap.sh
################################################################################################################################
###hapmap processing:
############################################################################################################################################
#Sec1#######################################################################################################################################
############################################################################################################################################
hapmapConcatenate(){
hap1=$1
hap2=$2
tmpDir0101=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0101"/hapmap1
cat $2 > "$tmpDir0101"/hapmap2
cat "$tmpDir0101"/hapmap1 | head -1 | cut -f 3-  | tr '\t' '\n' > "$tmpDir0101"/ls1
cat "$tmpDir0101"/hapmap2 | head -1 | cut -f 3-  | tr '\t' '\n' > "$tmpDir0101"/ls2
cat "$tmpDir0101"/ls1 "$tmpDir0101"/ls2 | sort | uniq -c | sed -E "s:^ +::g" | tr ' ' '\t'  | awk '{if($1==2) {print $2}}' | sed "1ichr\npos"  > "$tmpDir0101"/ls12
reorderDSfromList "$tmpDir0101"/hapmap1 "$tmpDir0101"/ls12 | hapmap2txt > "$tmpDir0101"/X1.txt
reorderDSfromList "$tmpDir0101"/hapmap2 "$tmpDir0101"/ls12 | hapmap2txt > "$tmpDir0101"/X2.txt
cat <(cat "$tmpDir0101"/X1.txt | sed 1d) <(cat "$tmpDir0101"/X2.txt | sed 1d )  | tr '_' '\t' | sort -k1,1 -k2,2n  > "$tmpDir0101"/X12.txt
cat "$tmpDir0101"/X12.txt  | awk '{ if($2 == x){ split(X,A,"\t");split($0,B,"\t"); str=$1"\t"$2; for (i=3;i<=NF;++i){ ab=A[i]","B[i]; str=str"\t"ab } ;X=str }else{print X; X=$0} ;x=$2}END{print X}' | sed 1d | awk 'function uniq (str){nA=split(str,B,",") ; asort(B, A); S=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s) {S=S","A[a]} ; s=A[a]};  if(S!="-"){gsub("-,","",S)} ; return S }{X=$1"\t"$2; for(i=3;i<=NF;++i){ X=X"\t"uniq($i)}; print X }' > "$tmpDir0101"/X12.txt2
cat <(cat "$tmpDir0101"/X1.txt | head -1 | tr '_' '\t' ) "$tmpDir0101"/X12.txt2 | uniq 
rm -rf $tmpDir0101
trap "rm -rf $tmpDir0101" EXIT
}
############################
hapmapMerge(){
#hapmap1=$1
#hapmap2=$2
tmpDir0102=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0102"/hapmap1
cat $2 > "$tmpDir0102"/hapmap2
###
i=1; for chr in $(cat $"$tmpDir0102"/hapmap1 | cut -f1 | uniq)
do
   Indx=$(printf "%05d\n" $i)
   mkdir "$tmpDir0102"/"$Indx"
   cat "$tmpDir0102"/hapmap1 | awk -v chr=$chr '{if($1==chr) print $0}' | cut -f 2- > "$tmpDir0102"/"$Indx"/X1.txt
   cat "$tmpDir0102"/"$Indx"/X1.txt | cut -f 1 > "$tmpDir0102"/"$Indx"/X2.txt
   cat "$tmpDir0102"/hapmap2 | awk -v chr=$chr '{if($1==chr) print $0}' | cut -f 2- > "$tmpDir0102"/"$Indx"/Y1.txt
   merge "$tmpDir0102"/"$Indx"/X2.txt "$tmpDir0102"/"$Indx"/Y1.txt  > "$tmpDir0102"/"$Indx"/XY1.txt
   paste "$tmpDir0102"/"$Indx"/X1.txt <(cat "$tmpDir0102"/"$Indx"/XY1.txt | cut -f 2- ) | sed "s:^:$chr\t:g" > "$tmpDir0102"/"$Indx"/XY2.txt
   i=$((i+1))
done
cat "$tmpDir0102"/*/XY2.txt > "$tmpDir0102"/XY2.txt
cat "$tmpDir0102"/XY2.txt | awk 'OFS="\t"{for (i=1;i<=NF;++i) {if($i=="NA") { gsub("NA","-",$i) } } }'1
rm -rf $tmpDir0102
trap "rm -rf $tmpDir0102" EXIT
}
############################
hapmapExtractSubtract()
{
hapmap=$1
bed=$2 # bed or hapmap
ExtractSubtract=$3 #1 extract, 0 subtract
tmpFile0101=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $hapmap > "$tmpFile0101"
cat "$tmpFile0101" | cut -f 1-2 | Rscript --vanilla -e 'args = commandArgs(trailingOnly=TRUE)
A = as.data.frame(read.table("stdin",header=FALSE))
B = as.data.frame(data.table::fread(args[1],header=FALSE,, select = c(1:2)))
A2 = as.data.frame(apply(A,1,function(x) paste0(x,collapse="_")))
colnames(A2) = "chr_pos"
if(B$V1[1]!="chr"){B = rbind.data.frame(c("chr","pos"),B) }
B2 = cbind.data.frame(chr_pos=apply(B,1,function(x) paste0(x,collapse="_")),bed="bed")
AB = plyr::join(A2,B2,by="chr_pos")
AB$bed[which(is.na(AB$bed))] <- 0
AB$bed[which(AB$bed == "bed")] <- 1
write.table(cbind.data.frame(AB$bed),"",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
'  $bed | paste - "$tmpFile0101" | awk -v ExtractSubtract=$ExtractSubtract '{if(NR==1 || $1==ExtractSubtract) print $0 }' | cut -f 2-
rm "$tmpFile0101"
trap "rm -f $tmpFile0101" EXIT
}
############################
hapmapExtract(){
   hapmapExtractSubtract $1 $2 1
}
############################
hapmapSubtract(){
   hapmapExtractSubtract $1 $2 0
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec2#######################################################################################################################################
############################################################################################################################################
hapmapSNP2SV(){
   cat $1 | awk '{if(NR == 1){print $0} else {S=$1"\t"$2; for(i=3;i<=NF;++i){ nA=split($i,A,""); s=A[1]; for(a=2;a<=nA;++a){s=s","A[a]}; S=S"\t"s   } ; print S} }'  
}

############################
hapmapGetSNP(){
cat $1 | awk 'function uniq (str){nA=split(str,B,",") ; asort(B, A); S=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s) {S=S","A[a]} ; s=A[a]};  if(S!="-"){gsub("-,","",S)} ; return S } {X="-";for(i=3;i<=NF;++i){X=X","$i}; X=uniq(X); x=0; nA=split(X,A,","); for(a=1;a<=nA;++a){if(length(A[a]) > 1){x=1}}; if(NR==1 || x == 0){print $0}  }'
}
############################
hapmapGetSV(){
cat $1 | awk 'function uniq (str){nA=split(str,B,",") ; asort(B, A); S=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s) {S=S","A[a]} ; s=A[a]};  if(S!="-"){gsub("-,","",S)} ; return S } {X="-";for(i=3;i<=NF;++i){X=X","$i}; X=uniq(X); x=0; nA=split(X,A,","); for(a=1;a<=nA;++a){if(length(A[a]) > 1){x=1}}; if(NR==1 || x == 1){print $0}  }'
}
############################
hapmapGetDiAlleles(){
cat $1 |awk 'function num_uniq_alleles (str){nA=split(str,B,",") ; asort(B, A) ; S=A[1]; s=A[1]; for(a=2;a <= nA ;++a){if(A[a] != s && A[a] != "-") {S=S","A[a]} ; s=A[a]}; gsub("^-,","",S);if(S=="-"){S=""}; return split(S,C,",") } {if(NR==1){print "2\t"$0} else {S="-"; for(i=3;i<=NF;++i){S=S","$i}; S=num_uniq_alleles(S); print S"\t"$0  }}' | awk '{if($1==2){print $0} }'  | cut -f 2-
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec3#######################################################################################################################################
############################################################################################################################################
hapmap2txt(){
cat $1 | awk '{print $1"_"$2"\t"$0}' | cut -f 1,4- 
}
############################
hapmap2stdhapmap(){
hapmap=$1
ref=$2
tmpDir0103=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 | sed "s:,::g" | hapmapGetSNP  > "$tmpDir0103"/hapmap
hapmapGetAlleles "$tmpDir0103"/hapmap | sed "1d" > "$tmpDir0103"/alleles
cat "$tmpDir0103"/hapmap | sed "1d" | cut -f 1-2 | awk '{print $1"\t"$2-1"\t"$2}' > "$tmpDir0103"/bed
bedtools getfasta -tab -fi $ref -bed "$tmpDir0103"/bed | cut -f 2 > "$tmpDir0103"/alleles_ref
paste "$tmpDir0103"/alleles_ref "$tmpDir0103"/alleles | awk '{X=$1;nA=split($2,A,","); for(a=1;a<=nA;++a){if($1 != A[a]){X=X"/"A[a]} }; print X }' | sed "1ialleles" > "$tmpDir0103"/alleles_refOrdered_unique
paste <(cat "$tmpDir0103"/hapmap | cut -f 3- | sed '1d' | tr '\t' ',' | awk 'function sortuniq (str){if(str==0 || str =="-"  ){return "-"}; nA=split(str,A,",");S=""; asort(A, B) ; for(a=1;a<=nA;++a){ if(B[a] !=0 && B[a] != s && B[a] != "-") {S=S","B[a]; s=B[a] }; s=B[a]};S=substr(S,2);if(S==""){S=0} ; return S } {print sortuniq($0) }' | tr ',' '/'  | paste - <(cat "$tmpDir0103"/hapmap | cut -f 3- | sed "s:[\t-]::g" | sed "1d" | perl -F -lane '{print sort@F}' | tr -s '[ACGT]' ) | awk '{if(length($2) ==2 ) {print substr($2,1,1)"/"substr($2,2,1)} else {print $1} }'|  sed '1iII') <(cat "$tmpDir0103"/hapmap | cut -f 1-2)        <(cat "$tmpDir0103"/hapmap | cut -f 3- | awk '{for( i=1;i<=NF;++i ){ s="";if($i=="-"){s="NN"} else if(length($i) >= 2){s=$i} else {s=$i$i} ; S=S","s   };S=substr(S,2); print S; S="" }'  | tr ',' '\t' )     | awk '{if(NR==1) {print "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t"$0} else {print $2"_"$3"\t"$1"\t"$2"\t"$3"\t+\tNA\tNA\tNA\tNA\tNA\tNA\t"$0} }' | cut -f 1-11,15- > "$tmpDir0103"/hapmap2
paste "$tmpDir0103"/alleles_refOrdered_unique "$tmpDir0103"/hapmap2 | awk '{print $2"\t"$0}' | cut -f 1-2,5-
rm -rf "$tmpDir0103"
trap "rm -rf $tmpDir0103" EXIT
}
############################
stdhapmap2hapmap(){
cat $1 | cut -f 3,4,12- | sed "s:NN:-:g" | tr '\t' ',' | sed "s:,AA:,A:g;s:,TT:,T:g;s:,CC:,C:g;s:,GG:,G:g" | tr ',' '\t' | sed "s:chrom\tpos\t:chr\tpos\t:g"
}
############################
hapmap2vcf(){
hapmap=$1
ref=$2
tmpDir0104=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $hapmap | sed "1d" | cut -f 3- > "$tmpDir0104"/txt2
cat "$tmpDir0104"/txt2|  tr '\t' ',' | sed "s:-,::g" | sed "s:,-$::g" > "$tmpDir0104"/seq
cat "$tmpDir0104"/seq | sortUniqueValuesDS > "$tmpDir0104"/alleles
cat $hapmap| sed "1d" | cut -f 1-2 | awk '{print $1"\t"$2-1"\t"$2}' > "$tmpDir0104"/bed
bedtools getfasta -tab -fi $ref -bed "$tmpDir0104"/bed | cut -f 2 | paste - "$tmpDir0104"/alleles |  awk '{nA=split($2,A,","); S=""; for(a=1;a<=nA;++a){if(A[a]!=$1){S=S","A[a]} }; S=substr(S,2); if(S==""){S="."} ; print $1"\t"S }'  > "$tmpDir0104"/refalt
#paste <(paste "$tmpDir0104"/refalt | tr '\t' ',' ) "$tmpDir0104"/seq | awk '{nA=split($1,A,",");S1=""; S2=$2;X=0; for(a=1;a<=nA;++a){nB=split(S2,B,","); S2="";x=0; for(b=1;b<=nB;++b){ if(B[b] == A[a] ){x+=1;X+=1}else{S2=S2","B[b]}  }; S1=S1","x  }; S1=substr(S1,2); print ".\t.\tDP4="S1":DP="X"\tGT:GP"  }' > "$tmpDir0104"/dep
paste <(paste "$tmpDir0104"/refalt | tr '\t' ',' ) "$tmpDir0104"/seq | awk '{nA=split($1,A,",");S1=""; S2=$2;X=0; for(a=1;a<=nA;++a){nB=split(S2,B,","); S2="";x=0; for(b=1;b<=nB;++b){ if(B[b] == A[a] ){x+=1;X+=1}else{S2=S2","B[b]}  }; S1=S1","x  }; S1=substr(S1,2); print ".\t.\tDP4="S1":DP="X"\tGT"  }' > "$tmpDir0104"/dep
cat "$tmpDir0104"/bed | awk '{print $1"\t"$3"\t."}' | paste - "$tmpDir0104"/refalt "$tmpDir0104"/dep | sed "1i#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" > "$tmpDir0104"/info
mkdir "$tmpDir0104"/Geno
num=$(cat $hapmap | head -1 | wc -w)
for i in $(seq 3 $num)
do
   #echo $i
   paste <(paste "$tmpDir0104"/refalt | tr '\t' ',' | sed "1ialleles") <(cat $hapmap | cut -f "$i" ) | awk '{if(NR==1){print $2}else if($2=="-"){ print "." } else {nA=split($1,A,",");nB=split($2,B,","); S=""; for(a=1;a<=nA;++a){ for(b=1;b<=nA;++b){ if(A[a] == B[b] ){S=S","a-1} }  }; S=substr(S,2) ; print S }  }' > "$tmpDir0104"/Geno/$(echo $i | awk '{printf "%020d\n", $0}' )
done

paste  "$tmpDir0104"/info "$tmpDir0104"/Geno/*
rm -rf "$tmpDir0104"
trap "rm -rf $tmpDir0104" EXIT
}
############################
vcf2hapmap(){
   cat $1 | grep -v "##" | awk 'function uniq (str){nA=split(str,B,",") ; asort(B, A) ; S0=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s && A[a] != "-") {S0=S0","A[a]} ; s=A[a]}; gsub("-,","",S0) ; return S0 } {if(NR==1){   nH=split($0,H,"\t"); S="chr\tpos"; for(i=10;i<=NF;++i){ S=S"\t"$i };inx=0   }else{ if(NR==2){nA=split($9,A,":");for(a=1;a<=nA;++a){if(A[a]=="GT"){inx=a}} }  ;  S=$1"\t"$2; alleles=$4","$5; split(alleles,Alleles,",") ; ; ; for(i=10;i<=NF;++i){ split($i,B,":"); allele=B[inx]; s="-"; nC=split(allele,C,"[|/,]"); for(c=1;c<=nC;++c){if(C[c]!="."){s=s","Alleles[C[c]+1]}}; S=S"\t"uniq(s) }  }; print S  }'
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec4#######################################################################################################################################
############################################################################################################################################
hapmapSplit(){
if [ -d "hapmapCHRs" ]; then
   echo "ERROR: hapmapCHRs is exist"
else
   cat $1 | awk '{if(NR==1){system("mkdir hapmapCHRs"); H=$0}else{if($1!=x){print H > "hapmapCHRs/"$1".hapmap"; print $0 > "hapmapCHRs/"$1".hapmap" }else{print $0 > "hapmapCHRs/"x".hapmap" } ;x=$1 } }'
   echo "hapmapCHRs was generated"
fi
}
############################
hapmapSampleIDsfromList(){
tmpFile0102=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $2 | sed "1ichr\npos" > "$tmpFile0102"
reorderDSfromList  $1 "$tmpFile0102"
rm "$tmpFile0102"
trap "rm -f $tmpFile0102" EXIT
}
############################
hapmapReAssignSampleIDsfromList(){
hapmap=$1
list=$2
tmpDir0105=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0105"/hapmap
head1=$(cat  "$tmpDir0105"/hapmap | head -1 | tr '\t' ',')
cat $list | awk -v h1=$head1 'BEGIN{nA=split(h1,A,",")}{ for(a=1;a<=nA;++a) { if(A[a]==$1) {print a"\t"$2}  }  }' | sed "1i1\tchr\n2\tpos" > "$tmpDir0105"/TXT
reord=$(cat "$tmpDir0105"/TXT | cut -f 1 | tr '\n' ',' | sed "s:,$:\n:g" )
reordHead=$(cat "$tmpDir0105"/TXT | cut -f 2 | tr '\n' '\t' | sed "s:\t$:\n:g" )
cat  "$tmpDir0105"/hapmap | awk -v reord=$reord 'BEGIN{nA=split(reord,A,",")} {for (a=1;a<=nA;++a) {S=S"\t"$A[a]}; S=substr(S,2); print S; S="" }' | sed "s:\t$::g" | sed '1d' | sed "1i$reordHead" 
rm -rf $tmpDir0105
trap "rm -rf $tmpDir0105" EXIT
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec5#######################################################################################################################################
############################################################################################################################################
hapmapCalculateMissingVariant () {
cat $1 | sed -E "1 s:[^\t]+:-:g" | cut -f 3- | awk '{if(NR==1){len=NF; x="miss" }else{x=gsub("-","",$0)/len }; print x }'
}
############################
hapmapCalculateMissingSample(){
cat $1 | awk '{if(NR==1){split($0,H,"\t"); str="" }else{for(i=3;i<=NF;++i){ if($i=="-"){S[i]+=1}  } } }END{ for(i=3;i<=NF;++i){ str=str","S[i]}  }END{str=substr(str,2); nA=split(str,A,","); for(a=1;a<=nA;++a){printf "%s\t%0.2f\n", H[a+2],A[a]/(NR-1) } } '
}
###########################
hapmapCalculateMAF(){
# cat $1 | awk  'function uniq (str){nA=split(str,B,",") ; asort(B, A) ; S0=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s && A[a] != "-") {S0=S0","A[a]} ; s=A[a]}; gsub("-,","",S0) ; return S0 } {if(NR==1) {print "MAF"} else { S=$3; for(i=4;i<=NF;++i){S=S","$i}; S2=S; gsub(",","",S2); gsub("-","",S2); len=length(S2) ; S3=uniq(S) ; S4="";z=1; nA=split(S,A,",")  ;for(a=1;a<=nA;++a){x=gsub(A[a],"",S2); S4=S4","x/len; if(x/len <= z){z=x/len} } ; S4=substr(S4,2) ;print z ; } }'
tmpFile0103=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0103"
cat "$tmpFile0103" | hapmapAlleleFreqStats | cut -f 1-2 | sed 1d| awk '{nA=split($1,A,",");split($2,B,","); x=A[1];y=B[1]; for(a=2;a<=nA;++a){ if(B[a] < y){x=A[a]; y=B[a]}  }; print x"\t"y }' | sed "1iallele\tmaf"
rm "$tmpFile0103"
trap "rm -f $tmpFile0103" EXIT
}
############################
hapmapFilterMissingVariant (){
   miss=$2
   cat $1 | awk -v miss=$miss '{ if(NR==1){print $0;N=NF-2}else{S=$0; x=gsub("-","",S); if(x/N <= miss){print $0}  } }'
}
############################
hapmapFilterMissingSample(){
Smiss=$2
# cat $1 | hapmap2txt  | transverseDS  | hapmapMiss - $Smiss | transverseDS | tr '_' '\t'
tmpDir0106=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0106"/hapmap
hapmapCalculateMissingSample "$tmpDir0106"/hapmap | awk -v Smiss=$Smiss '{if($2<=Smiss){print $1}}' | sed "1ichr\npos" > "$tmpDir0106"/sel.list
reorderDSfromList "$tmpDir0106"/hapmap "$tmpDir0106"/sel.list
rm -rf "$tmpDir0106"
trap "rm -rf $tmpDir0106" EXIT
}
############################
hapmapFilterMAF(){
maf=$2
# cat $1 | awk -v maf=$maf 'function uniq (str){nA=split(str,B,",") ; asort(B, A) ; S0=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s && A[a] != "-") {S0=S0","A[a]} ; s=A[a]}; gsub("-,","",S0) ; return S0 } {if(NR==1) {print $0} else { S=$3; for(i=4;i<=NF;++i){S=S","$i}; S2=S; gsub(",","",S2); gsub("-","",S2); len=length(S2) ; S3=uniq(S) ; S4="";z=1; nA=split(S,A,",")  ;for(a=1;a<=nA;++a){x=gsub(A[a],"",S2); S4=S4","x/len; if(x/len <= z){z=x/len} } ; S4=substr(S4,2) ;  if(z >= maf) {print $0} } }'
tmpFile0104=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0104"
cat "$tmpFile0104" | hapmapCalculateMAF | paste - "$tmpFile0104" | awk -v maf=$maf '{if(NR==1 || $2 >= maf){ print $0} }' | cut -f 3-
rm "$tmpFile0104"
trap "rm -f $tmpFile0104" EXIT
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



#Sec6#######################################################################################################################################
############################################################################################################################################
############################################################################################################################################
hapmapGetAlleles () {
  cat $1 | hapmapAlleleFreqStats | cut -f 1 | sed "s:-::g"
}
############################
hapmapAlleleFreqStats(){
cat $1 | awk 'function uniq (str){nA=split(str,B,",") ; asort(B, A) ; S0=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s && A[a] != "-") {S0=S0","A[a]} ; s=A[a]}; gsub("-,","",S0) ; return S0 } {if(NR==1) {print "alleles\tfreq\tcnt\tsum"} else { S=$3; for(i=4;i<=NF;++i){S=S","$i}; S2=uniq(S); if(S2=="-"){print S2"\t0\t0\t0"}else{sum=0; nA=split(S2,A,","); nB=split(S,B,",");S3=""; for(a=1;a<=nA;++a){s=0;for(b=1;b<=nB;++b){ if(A[a]==B[b]){s+=1}  }; sum+=s; S3=S3","s } ; S3=substr(S3,2); nS=split(S3,S4,","); S5=S4[1]; S6=S4[1]/sum; for(x=2;x<=nS;++x){S5=S5","S4[x]; S6=S6","S4[x]/sum} ; print S2,S6,S5,sum } } }' | tr ' ' '\t'
}
############################
hapmapAlleleTypeFreq(){
tmpFile0105=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0105"
hapmapGetAlleles "$tmpFile0105" | paste - "$tmpFile0105" | sed 1d | awk '{nA=split($1,A,","); if(nA==2){sum=0; homo1=0; homo2=0; het=0; for(i=4;i<=NF;++i){ if($i==A[1]){homo1+=1;sum+=1};if($i==A[2]){homo2+=1;sum+=1}; if($i==A[1]","A[2] || $i==A[2]","A[1]){het+=1;sum+=1}  }; if(sum>0) { printf "%s\t%s\t%s\t%0.2f,%0.2f,%0.2f\n",  $2,$3, $1, homo1/sum,het/sum,homo2/sum} } }' | sed "1ichr\tpos\tallele1,allele2\tA,B,H"
rm "$tmpFile0105"
trap "rm -f $tmpFile0105" EXIT
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec7#######################################################################################################################################
############################################################################################################################################
hapmapGetIDs(){
cat $1 | head -1 | cut -f 3- | tr '\t' '\n'
}
############################
hapmapGetFreqSampleAlleles(){
 cat $1 | awk -v id=$2 '{if(NR==1){inx=0; for(i=3;i<=NF;++i){if($i==id){inx=i}} }else{S=$1"\t"$2"\t"$inx; X=0; Y=0; for(i=3;i<=NF;++i){ if(i!=inx && $i != "-"){Y+=1;if($i==$inx){X+=1}} }; if(Y>0){freq=X/Y}else{freq=0}; print S,freq,Y }}' | sed "1ichr\tpos\tallele\tfreq\tdep"
# cat $1 | awk -v id=$2 '{if(NR==1){inx=0; for(i=3;i<=NF;++i){if($i==id){inx=i}} }else{S=$1"\t"$2"\t"$inx; str="" ; for(i=3;i<=NF;++i){ if(i!=inx){str=str""$i} } ; gsub("-","",str); gsub(",","",str); len=length(str);cnt=0; if(length($inx)==1 && len>0) { cnt=gsub($inx,"",str)  ; print S"\t"cnt/len"\t"len}}}' | sed "1ichr\tpos\tallele\tfreq\tdep"
}
############################
hapmapGetUniqueSampleAlleles(){
id=$2
cat $1  | awk -v id=$id '{if(NR==1){ for(i=3;i<=NF;++i){if($i==id){indx=i}} }else{S=$1"\t"$2; s=$indx; nA=split(s,A,","); for(i=3;i<=NF;++i) { if(i!=indx) {nB=split($i,B,","); for(b=1;b<=nB;++b){ for(a=1;a<=nA;++a){ if(A[a] == B[b] ){A[a] = ""}   }  }   } } ; s2=""; for(a=1;a<=nA;++a){if(A[a] != "" ) {s2=s2","A[a]}}; s2=substr(s2,2); if(s2==""){s2="-"} ; print S"\t"s2}}' | grep -v "-" | sed "1ichr\tpos\t$id"
}
############################
hapmapGetROHforSingleSample(){
cat $1 | grep -vw "-" | awk '{S=$1"\t"$2"\t"; if($0~","){S=S"0"}else{S=S"1"}; print S}' | awk '{if(NR==1){X=$2} else{ if($3!=homo){ print chr"\t"X1"\t"X"\t"X-X1+1"\t"cnt"\t"homo; X1=$2; cnt=1 }else{X0=$2;cnt+=1};homo=$3; X=$2; chr=$1} }END{ print chr"\t"X1"\t"X"\t"X-X1+1"\t"cnt"\t"homo }'  | sed 1d | awk '{if($6==1 && $5 >= 50) print $0}' |  cut -f 1-5 | grep -v "-"
}
############################
hapmapGetROH(){
echo "id chr pos1 pos2 len cnt" | tr ' ' '\t'
tmpFile0108=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0108"
head=($(cat "$tmpFile0108" | head -1 | tr '\t' ' ' ))
for i in $(seq 3 $(cat "$tmpFile0108" | head -1 | wc -w))
do
   cat "$tmpFile0108" | cut -f 1-2,$i  | grep -vw "-" | hapmapGetROHforSingleSample | sed "s:^:${head[$((i-1))]}\t:g"
done
rm "$tmpFile0108"
trap "rm -f $tmpFile0108" EXIT
}
############################
hapmapLowFreqAlleleMask(){
mask=$2
tmpDir0107=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0107"/hapmap
cat "$tmpDir0107"/hapmap | hapmapAlleleFreqStats | cut -f 1-2 | awk -v mask=$mask '{if(NR==1){print "-"}else{S="-"; nA=split($1,A,",");split($2,B,","); for(a=1;a<=nA;++a){if(B[a] < mask){S=S","A[a]} }; S=substr(S,3); if(S==""){S="-"} ; print S  } }' > "$tmpDir0107"/maskedAL
paste "$tmpDir0107"/maskedAL | paste - "$tmpDir0107"/hapmap |  awk 'function mask(alleles,str){ nA=split(alleles,A,","); nB=split(str,B,","); S=""; for(b=1;b<=nB;++b){s=B[b]; for(a=1;a<=nA;++a){ if(s==A[a]){s="-"}}; if(s!="-"){S=S","s}}; S=substr(S,2); if(S==""){S="-"} ; return S } {S=$2"\t"$3; for(i=4;i<=NF;++i){if(NR==1){S=S"\t"$i} else {S=S"\t"mask($1,$i)}}; print S }'
rm -rf $tmpDir0107
trap "rm -rf $tmpDir0107" EXIT
}
############################
hapmapGetHomoPolymorphic(){
tmpFile0106=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > $tmpFile0106
cat $tmpFile0106 | hapmapAlleleFreqStats | cut -f 1 | paste - <(cat $tmpFile0106 | cut -f 3- | tr '\t' ';' | sed 1d | awk '{nA=split($0,A,";"); homo=1; for(a=1;a<=nA;++a){ if(A[a] ~ ","){homo=0;break}  }; print homo }' | sed "1i1" ) $tmpFile0106 | awk '{if(NR==1 || ($1 ~ "," && $2 == 1) ) print $0 }' | cut -f 3-
rm "$tmpFile0106"
trap "rm -f $tmpFile0106" EXIT
}
############################
hapmapGetConsensusCall(){
tmpFile0107=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > $tmpFile0107
cat $tmpFile0107 | hapmapAlleleFreqStats | awk -v min=$2 -v max=$3 '{if(min==""){min=0.1}; if(max==""){max=0.9}; if(NR==1){print "ConsensusCall"} else {nA=split($1,A,",");split($2,B,","); s="-"; freq=0; for(a=1;a<=nA;++a){if(B[a] >= min && B[a] <= max && A[a] != "*"){s=s","A[a]; freq=B[a]} }; gsub("^-,","",s); print s }}' | paste <(cat $tmpFile0107 | cut -f 1-2) -

rm "$tmpFile0107"
trap "rm -f $tmpFile0107" EXIT
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
