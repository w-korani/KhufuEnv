# source /cluster/projects/khufu/korani_projects/KhufuEnv/06.KhufuEnvDataSet.sh
############################################################################################################################################
###data set processing
############################################################################################################################################
#Sec1#######################################################################################################################################
############################################################################################################################################
dimDS(){
cat $1 | awk 'END{print NR*NF": "NR"*"NF}'
}
############################
getIDsDS(){
cat $1 | head -1 | tr '\t' '\n' | grep -n '' | tr ':' '\t'
}
############################
revertDS(){
cat $1 |  awk '{for(i=1;i<=NF;++i){ Z[i][NR] = $i}}END{ for(x=NR;x>=1;--x) { for(y=1;y<=NF;++y){S=S"\t"Z[y][x]};S=substr(S,2);print S; S;S=""  } }'
}
############################
transverseDS(){
cat $1 | awk '{for(i=1;i<=NF;++i){ Z[NR][i] = $i}}END{ for(x=1;x<=NF;++x) { for(y=1;y<=NR;++y){S=S"\t"Z[y][x]};S=substr(S,2); print S;S=""  } }'  
}
############################
long2wideDS(){
cat $1 | Rscript -e 'args=commandArgs(TRUE); A=read.table("stdin",header=T); B= tidyr::spread(A,colnames(A)[2],colnames(A)[3]) ; write.table(B,stdout(),col.names=T,row.names=F,quote=F,sep="\t") '
}
############################
wide2longDS(){
cat $1 | awk '{if(NR==1) {head=$0; nA=split(head,A,"\t")} else { for(i=2;i<=NF;++i){ print $1"\t"A[i]"\t"$i}  } }' | sed "1ifactor1\tfactor2\tdata"
}
############################
sortValuesDS(){
   DS=$1
   cat $DS | awk 'function sort(str){ nA=split(str,A,","); asort(A, B); S=B[1] ; for(a=2;a<=nA;++a){ S=S","B[a] } return S } 
   {
      S=""; for(i=1;i<=NF;++i){S=S";"sort($i)}; S=substr(S,2); print S 
   }' | tr ';' '\t'
}
############################
sortReverseValuesDS(){
   cat $1 | awk 'function sort(str){ nA=split(str,A,","); asort(A, B); S=B[nA] ; for(a=(nA-1);a > 0 ; --a){ S=S","B[a] } return S } 
   {
      S=""; for(i=1;i<=NF;++i){S=S";"sort($i)}; S=substr(S,2); print S 
   }' | tr ';' '\t'
}
############################
sortUniqueValuesDS(){
   cat $1 | sortValuesDS | awk '
      function uniq (str){nA=split(str,A,",") ; S=A[1]; s=A[1]; for(a=2;a <= nA ;++a){ if(A[a] != s && A[a] != "-") {S=S","A[a]} ; s=A[a]}; ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"uniq($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}
############################
maxValuesDS(){
   cat $1 | sortValuesDS | awk '
      function max (str){nA=split(str,A,","); S=A[1] ; for(a=2;a<=nA;++a){ if(A[a] > S ) {S=A[a] } } ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"max($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}
############################
minValuesDS(){
   cat $1 | sortReverseValuesDS | awk '
      function min (str){nA=split(str,A,","); S=A[nA] ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"min($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}
############################
minSkipZeroValuesDS(){
   cat $1 | sortReverseValuesDS | awk '
      function min (str){nA=split(str,A,","); S=A[1] ; for(a=2;a<=nA;++a){ if(A[a] < S && A[a] != 0 ) {S=A[a] } } ; if(S==0){S="-"} ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"min($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}
############################
maxSkipZeroValuesDS(){
   cat $1 | sortValuesDS | awk '
      function max (str){nA=split(str,A,","); S=A[1] ; for(a=2;a<=nA;++a){ if(A[a] > S ) {S=A[a] } } ; if(S==0){S="-"} ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"max($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}
############################
maxLenDS(){
   cat $1 | sortValuesDS | awk '
      function max (str){nA=split(str,A,","); S=A[1] ; for(a=2;a<=nA;++a){ if(length(A[a]) > length(S) ) {S=A[a] } } ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"max($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}
############################
minLenDS(){
   cat $1 | sortValuesDS | awk '
      function min (str){nA=split(str,A,","); S=A[1] ; for(a=2;a<=nA;++a){ if(length(A[a]) < length(S) ) {S=A[a] } } ; return S }
      { 
         S=""; for(i=1;i<=NF;++i){S=S";"min($i)}; S=substr(S,2); print S 
      }' | tr ';' '\t'
}

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################


############################################################################################################################################
#Sec2#######################################################################################################################################
############################################################################################################################################
merge () {
Rscript -e 'args = commandArgs(trailingOnly=TRUE)
A = as.data.frame(data.table::fread(args[1],header=FALSE,sep="\t"))
B = as.data.frame(data.table::fread(args[2],header=FALSE,sep="\t"))  
AB = plyr::join(A,B,by="V1")
write.table(AB,"",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)' $1 $2 ;
}
############################
extract(){
tmpDir0601=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0601"/A
cat $2 > "$tmpDir0601"/B
Rscript -e 'args = commandArgs(trailingOnly=TRUE)
A = as.data.frame(data.table::fread(args[1],header=FALSE))
B = as.data.frame(data.table::fread(args[2],header=FALSE))
if(ncol(B)==1){B=cbind.data.frame(B,"B")}
AB = plyr::join(A,B,by="V1")
write.table(AB,"",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)' "$tmpDir0601"/A "$tmpDir0601"/B |  grep -v "NA$" | awk -v Bcol=$(cat "$tmpDir0601"/B | head -1 | wc -w) '{if(Bcol==1) {NF--} }'1 | tr ' ' '\t'
rm -rf $tmpDir0201
trap "rm -rf $tmpDir0601" EXIT
}
############################
subtract (){
tmpDir0601=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0601"/A
cat $2 > "$tmpDir0601"/B
Rscript -e 'args = commandArgs(trailingOnly=TRUE)
A = as.data.frame(data.table::fread(args[1],header=FALSE))
B = as.data.frame(data.table::fread(args[2],header=FALSE))
if(ncol(B)==1){B=cbind.data.frame(B,"B")}
AB = plyr::join(A,B,by="V1")
write.table(AB,"",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)' "$tmpDir0601"/A "$tmpDir0601"/B | 
grep NA | sed "s:$:\tNA\tNA:g" | sed "s:\tNA.*::g" ;
rm -rf $tmpDir0201
trap "rm -rf $tmpDir0601" EXIT
}
############################
reorderDSfromDS(){
tmpFile0601=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > $tmpFile0601
head1=$(cat $tmpFile0601 | head -1 | tr '\t' ',')
reord=$(cat $2 | head -1 | tr '\t' '\n' | awk -v h1=$head1 'BEGIN{nA=split(h1,A,",")}{ for(a=1;a<=nA;++a){ if(A[a]==$0){S=S","a} } }END{S=substr(S,2); print S}')
cat $tmpFile0601 | awk -v reord=$reord 'BEGIN{nA=split(reord,A,",")} {for (a=1;a<=nA;++a) {S=S"\t"$A[a]}; S=substr(S,2); print S; S="" }' 
rm "$tmpFile0601"
trap "rm -f $tmpFile0601" EXIT
}
############################
reorderDSfromList(){
tmpFile0602=$(mktemp "./KhufuEnviron.XXXXXXXXX")  
cat $1 > $tmpFile0602
head1=$(cat $tmpFile0602 | head -1 | tr '\t' ';')
reord=$(cat $2 | awk -v h1=$head1 'BEGIN{nA=split(h1,A,";")}{ for(a=1;a<=nA;++a){ if(A[a]==$0){S=S","a} } }END{S=substr(S,2); print S}')
cat $tmpFile0602 | awk -v reord=$reord 'BEGIN{nA=split(reord,A,",")} {for (a=1;a<=nA;++a) {S=S"\t"$A[a]}; S=substr(S,2); print S; S="" }' | sed "s:\t$::g"
rm "$tmpFile0602"
trap "rm -f $tmpFile0602" EXIT
}
############################
combineSimilarColumnDS(){
tmpDir0603=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0603"/ds.txt
cat "$tmpDir0603"/ds.txt | head -1 | tr '\t' '\n' | sort -V | uniq > "$tmpDir0603"/ds.list
reorderDSfromList "$tmpDir0603"/ds.txt "$tmpDir0603"/ds.list | awk '{
  if(NR==1){ nH=split($0,H,"\t")} ;nA=split($0,A,"\t"); S=""; skip=0 ; s="-";
  for(a=2;a<=nA;++a)
  {
     if(H[a]==H[a-1] )
     {
         if(skip == 1) { s=s","A[a] } else { s=A[a-1]","A[a] }; skip=1;
     }
     else
     {
        if(skip==1){skip=0;S=S";"s;s="-"}else{S=S";"A[a-1]}
     }
     if( a == NF ) { if(skip==1){S=S";"s}else{S=S";"A[a]} }
   }
   S=substr(S,2); print S
}' | tr ';' '\t'
rm -rf $tmpDir0603
trap "rm -rf $tmpDir0603" EXIT
}
############################
uniqueDS (){
cat $1 | awk 'function sortuniq (str){if(str =="-"  ){return "-"}; nA=split(str,A,"\t");S=""; asort(A, B) ; for(a=1;a<=nA;++a){ if( B[a] != s && B[a] != "-") {S=S","B[a]; s=B[a] }; s=B[a]};S=substr(S,2);if(S==""){S="-"} ; s=""; return S } {print sortuniq($0) }' 
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
