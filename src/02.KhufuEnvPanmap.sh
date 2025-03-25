# source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnvpanmap.sh
################################################################################################################################
###papmap processing:
############################################################################################################################################
#Sec1#######################################################################################################################################
############################################################################################################################################
panmapConcatenate(){
pan1=$1
pan2=$2
if [[ $(cat $pan1 | head -1 | cut -f 4) != $(cat $pan2 | head -1 | cut -f 4) ]]
then 
   echo "ERROR: the parents of the two panmaps (column 4) do not match"
else
   tmpDir0201=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
   panmap2hapmap $pan1 > "$tmpDir0201"/hap1
   panmap2hapmap $pan2 > "$tmpDir0201"/hap2
   id1=$(echo $pan1 | sed "s:.*/::g" | sed "s:.panmap::g" )
   id2=$(echo $pan2 | sed "s:.*/::g" | sed "s:.panmap::g" )
   hapmapConcatenate "$tmpDir0201"/hap1 "$tmpDir0201"/hap2 > "$tmpDir0201"/"$id1"."$id2".concatenated
   hapmap2panmap "$tmpDir0201"/"$id1"."$id2".concatenated $(cat $pan1 | head -1 | cut -f 4)
   rm -rf $tmpDir0201
   trap "rm -rf $tmpDir0201" EXIT
fi
}
############################
panmapMerge(){
pan1=$1
pan2=$2
if [[ $(cat $pan1 | head -1 | cut -f 4) != $(cat $pan2 | head -1 | cut -f 4) ]]
then 
   echo "ERROR: the parents of the two panmaps (column 4) do not match"
else
   tmpDir0202=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
   panmap2hapmap $pan1 > "$tmpDir0202"/hap1
   panmap2hapmap $pan2 | cut -f $(cat $pan1 | head -1 | cut -f 4 | tr ',' ' ' | wc -w | awk '{print "1-2,"$0+2+1"-"}') > "$tmpDir0202"/hap2
   id1=$(echo $pan1 | sed "s:.*/::g" | sed "s:.panmap::g" )
   id2=$(echo $pan2 | sed "s:.*/::g" | sed "s:.panmap::g" )
   hapmapMerge "$tmpDir0202"/hap1 "$tmpDir0202"/hap2 > "$tmpDir0202"/"$id1"."$id2".merged
   hapmap2panmap "$tmpDir0202"/"$id1"."$id2".merged $(cat $pan1 | head -1 | cut -f 4)
   rm -rf $tmpDir0202
   trap "rm -rf $tmpDir0202" EXIT
fi
}
############################
panmapExtract(){
pan=$1
bed=$2
id1=$(echo $pan | sed "s:.*/::g" | sed "s:.panmap::g" )
id2=$(echo $bed | sed "s:.*/::g" | sed "s:.bed::g" )
hapmapExtract $pan $bed > "$id1".extracted."$id2".panmap
echo ""$id1".extracted."$id2".panmap was generated"
panmapGetFasta "$id1".extracted."$id2".panmap $pan 
}
############################
panmapSubtract(){
pan=$1
bed=$2
id1=$(echo $pan | sed "s:.*/::g" | sed "s:.panmap::g" )
id2=$(echo $bed | sed "s:.*/::g" | sed "s:.bed::g" )
hapmapSubtract $pan $bed > "$id1".subtracted."$id2".panmap
echo ""$id1".subtracted."$id2".panmap was generated"
panmapGetFasta "$id1".subtracted."$id2".panmap $pan 
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec2#######################################################################################################################################
############################################################################################################################################
panmapGetFasta(){
panmap=$1
fa="$2".fa
tmpDir0203=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
cat $panmap | cut -f 1-2 | sed '1d' | tr '\t' '_' > "$tmpDir0203"/hap
fasta2txt $fa | tr '_' '\t' | awk '{print $1"_"$2"\t"$3"\t"$4}' > "$tmpDir0203"/fa
merge "$tmpDir0203"/hap "$tmpDir0203"/fa | awk '{print $1"_"$2"\t"$3}' | txt2fasta > "$panmap".fa
echo ""$panmap".fa was generated"
rm -rf $tmpDir0203
trap "rm -rf $tmpDir0203" EXIT
}
############################
panmapGetSNP(){
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap$::g" )
cat $1 | awk '{if(NR==1){print $0}else{nA=split($3,A,","); x=0; for(a=1;a<=nA;++a){if(A[a] > 1 ){x=1}}; if(x==0){print $0}   }}' > "$id".SNP.panmap
echo ""$id".SNP.panmap was generated"
panmapGetFasta "$id".SNP.panmap $1 
}
############################
panmapGetSV(){
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap$::g" )
cat $1 | awk '{if(NR==1){print $0}else{nA=split($3,A,","); x=0; for(a=1;a<=nA;++a){if(A[a] > 1 ){x=1}}; if(x==1){print $0}   }}' > "$id".SV.panmap
echo ""$id".SV.panmap was generated"
panmapGetFasta "$id".SV.panmap $1
}
############################
panmapGetDiAlleles(){
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap$::g" )
cat $1 | awk '{if(NR==1){print $0}else{nA=split($3,A,","); if(nA==2){print $0} }}' > "$id".DiAllelic.panmap
echo ""$id".DiAllelic.panmap was generated"
panmapGetFasta "$id".DiAllelic.panmap $1
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec3#######################################################################################################################################
############################################################################################################################################
panmap2txt(){
cat $1 | awk '{print $1"_"$2"\t"$0}' | cut -f 1,4- 
}
############################
panmap2hapmap(){
tmpFile0201=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0201"
cat "$1".fa | fasta2txt | tr '_' '\t' | awk '{print $1"_"$2"\t"$4}' | fastatxtAlleleCombine | cut -f 2 | sed "1ialleles" | paste - $tmpFile0201 | awk '{if(NR==1){S=$2"\t"$3; for(i=5;i<=NF;++i){S=S"\t"$i};gsub(",","\t",S); print S}else{S=$2"\t"$3; nA=split($1,A,","); nB=split($5,B,","); for(b=1;b<=nB;++b){s=A[B[b]]; if(s==""){s="-"} S=S"\t"s}; for(i=6;i<=NF;++i){ s="-"; nB=split($i,B,","); s="-"; for(b=1;b<=nB;++b){ if(A[B[b]]!=""){s=s","A[B[b]]} };gsub("^-,","",s) ; S=S"\t"s  } ; print S  } }'
rm "$tmpFile0201"
trap "rm -f $tmpFile0201" EXIT
}
############################
hapmap2panmap(){
parents=$2
tmpDir0204=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
id=$(echo $1 | sed "s:.*/::g" )
cat $1 > "$tmpDir0204"/hapmap
inds=$(echo $parents  | tr ',' '\n' | while read -r id;do  ind=$(cat "$tmpDir0204"/hapmap | head -1 | tr '\t' '\n' | grep -wn $id | sed "s:\:.*::g"); echo $ind  ;done | sort -n | tr '\n' ','| sed "s:,$:\n:g")
cat "$tmpDir0204"/hapmap | sed -E "s:[A-Z]+,[A-Z]+:-:g" | cut -f 1,2,"$inds" | hapmapGetAlleles |  sed "s:^$:-:g" | paste - "$tmpDir0204"/hapmap | awk '{if($1!="-") {print $0} }' > "$tmpDir0204"/hapmap2
cat "$tmpDir0204"/hapmap2 | awk '{if(NR==1){print $0}else {nA=split($1,A,","); S=$1"\t"$2"\t"$3; for(i=4;i<=NF;++i){ s="-"; nB=split($i,B,","); for(b=1;b<=nB;++b){ for(a=1;a<=nA;++a){if(B[b] == A[a]){s=s","a} };  }; gsub("^-,","",s) ;S=S"\t"s  }  ; print S } }' > "$tmpDir0204"/hapmap3
echo $parents | tr ',' '\n' > "$tmpDir0204"/list1
hapmapSampleIDsfromList "$tmpDir0204"/hapmap3 "$tmpDir0204"/list1 | cut -f 3- | tr '\t' ',' | sed "s:-:0:g" > "$tmpDir0204"/hapmap4B
inds1=$(echo $parents  | tr ',' '\n' | while read -r id;do  ind=$(cat "$tmpDir0204"/hapmap3 | head -1 | tr '\t' '\n' | grep -wn $id | sed "s:\:.*::g"); echo $ind  ;done | sort -n | tr '\n' ','| sed "s:,$:\n:g")
cat "$tmpDir0204"/hapmap3 | cut -f 1-3 | sed 1d| awk -v id=$id 'BEGIN{print "chr\tpos\tlen"}{nA=split($1,A,",");S=$2"\t"$3"\t"length(A[1]);print  ">"$2"_"$3"_1\n"A[1] > id".panmap.fa"; for(a=2;a<=nA;++a){ S=S","length(A[a]); print ">"$2"_"$3"_"a"\n"A[a] > id".panmap.fa" } ; print S} '  > "$tmpDir0204"/hapmap4A
paste "$tmpDir0204"/hapmap4A "$tmpDir0204"/hapmap4B <(cat "$tmpDir0204"/hapmap3 | cut -f 1-3,"$inds1" --complement ) > "$id".panmap
rm -rf $tmpDir0204
trap "rm -rf $tmpDir0204" EXIT
echo ">>> "$id".panmap & "$id".panmap.fa were generated"
}
############################
panmap2vcf(){
panmap=$1
tmpDir0205=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat "$panmap" > "$tmpDir0205"/panmap
cat "$tmpDir0205"/panmap | cut -f 4- | awk '{if(NR==1){print $0}else{ nA=split($1,A,","); S=A[1]-1; for(a=2;a<=nA;++a){S=S","A[a]-1}; for(i=2;i<=NF;++i){ nB=split($i,B,","); S=S"\t"B[1]-1; for(b=2;b<=nB;++b){S=S","B[b]-1}   } ; gsub("-1",".",S) ; print S}  }'  > "$tmpDir0205"/panmap2
cat "$tmpDir0205"/panmap2 | awk '{if(NR==1){S="-"; split($0,H,"\t"); for(i=1;i<=NF;++i); for(i=1;i<=NF;++i){ X=$i; gsub("[.][0-9]+$","",X); if(X!=Y){S=S"\t"X} ; Y=X};S=substr(S,3); print S }else{S=$1;X=H[1]; gsub("[.][0-9]+$","",X); for(i=2;i<=NF;++i){ Y=H[i]; gsub("[.][0-9]+$","",Y);  if(H[i]~"[.][0-9]+$" && X==Y){ S=S"|"$i} else {S=S"\t"$i}; X=H[i]; gsub("[.][0-9]+$","",X)  }; print S;} }' | awk '{X=$1;gsub(",","\t",X); for(i=2;i<=NF;++i){X=X"\t"$i} ; print X}' > "$tmpDir0205"/panmap2B
cat "$tmpDir0205"/panmap | cut -f 1-3 | paste <(cat "$panmap".fa | fasta2txt | tr '_' '\t' | awk '{print $1"_"$2"\t"$4}' | fastatxtAlleleCombine |cut -f 2 | sed "1ialleles") - | sed 1d | awk '{S=$2"\t"$3"\t.\t";nA=split($1,A,","); S=S"\t"A[1]; s="-"; for(a=2;a<=nA;++a){s=s","A[a]} ; s=substr(s,3); if(s==""){s="."} ; S=S"\t"s"\t.\t.\tMeth=Khufu\tGT" ; print S}' | sed "1i#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" > "$tmpDir0205"/panmap2A
paste "$tmpDir0205"/panmap2A "$tmpDir0205"/panmap2B | tr '\t' ';' | sed "s:;;:;:g" | tr ';' '\t'
rm -rf $tmpDir0205
trap "rm -rf $tmpDir0205" EXIT
}
############################
vcf2panmap(){
parents=$2
tmpDir0206=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
id=$(echo $1 | sed "s:.*/::g" )
cat $1 | grep -v "##" | awk '{if(NR>1){if($2 != x) {print X} };x=$2; X=$0  }END{ if($2 != $x) {print $0}}' | cut -f 1,2,4,5,9- > "$tmpDir0206"/X1
cat "$tmpDir0206"/X1 | head -2 | transverseDS | awk '{if($2~"[|]"){nA=split($2,A,"|"); for(a=1;a<=nA;++a){ print $1"."a}  }else{print $1} }' | transverseDS | cat - <(cat "$tmpDir0206"/X1 | sed 1d | tr '|' '\t' ) > "$tmpDir0206"/X2
cat "$tmpDir0206"/X2 | cut -f 1-4 | sed "1d" | awk '{print $1"\t"$2"\t"$3","$4}' | awk -v id=$id  '{nA=split($3,A,",");s="-"; for(a=1;a<=nA;++a){s=s","length(A[a]); print ">"$1"_"$2"_"a"\n"A[a] > id".panmap.fa" };gsub("^-,","",s); print $0"\t"s }' | cut -f 1-2,4 | sed "1ichr\tpos\tlen" > "$tmpDir0206"/X3A
cat "$tmpDir0206"/X2 | cut -f 5- | awk 'function GTinx(str){nA=split(str,A,":");X=0; for(a=1;a<=nA;++a){if(A[a]=="GT"){X=a;break}}; return X }{if(NR==1){print $0} else{ S=$1; for(i=2;i<=NF;++i){s="-"; nB=split($i,B,","); for(b=1;b<=nB;++b){if(B[b]=="."){s=s",-"}else{s=s","B[b]+1}; }; gsub("^-,","",s); S=S"\t"s   }; print S }}' | cut -f 2- > "$tmpDir0206"/X3B
inx=$(cat "$tmpDir0206"/X3B | head  | awk -v parents=$parents '{ if(NR==1){nA=split(parents,A,","); for(a=1;a<=nA;++a) { for(i=1;i<=NF;++i){X=$i;gsub("[.][0-9]+$","",X) ;  if(A[a]==X){ S=S","i} } }; S=substr(S,2);  print S  }}')
paste "$tmpDir0206"/X3A <(cat "$tmpDir0206"/X3B | cut -f $inx | tr '\t' ',') <(cat "$tmpDir0206"/X3B | cut -f $inx --complement ) > "$id".panmap
rm -rf $tmpDir0206
trap "rm -rf $tmpDir0206" EXIT
echo ">>> "$id".panmap & "$id".panmap.fa were generated"
}
############################
panmap2heatmap()
{
panmap=$1
cat $panmap |  awk 'function sortuniq (str){nA=split(str,A,","); if(nA==1){return str} ;S=""; asort(A, B) ; for(a=1;a<=nA;++a){ if(B[a] !="-" && B[a] != s) {S=S","B[a]; s=B[a] }; s=B[a]};S=substr(S,2);if(S==""){S=0} ; return S } { if(NR==1){H=$0; split(H,h,"\t"); ID=h[4]; nZ=split(ID,Z,","); for (z=1;z<=nZ;++z){z2=z2"\t"Z[z]} ; print H } else { 
nid=split(ID,id,","); 
for(i=1;i<=nid;++i) {parents=parents"\t"i};
nA = split($4,A,","); 
for(i=1;i<=NF;++i) {C[i]=$i}; 
for(i=5;i<=NF;++i) { nB=split(C[i],B,","); C[i] = "" ; for(b=1;b<=nB;++b) { for(a=1;a<=nA;++a) { if(B[b] == A[a] ) {C[i]=C[i]","a}  } }   };
for(i=5;i<=NF;++i) { if(C[i]=="") { C[i] = "0" } else {C[i] = substr(C[i],2) } ; if(C[i] ~ "," ) {D=D"\t"sortuniq(C[i])} else {D=D"\t"C[i]}  }; print $1"\t"$2"\t"$3"\t"$4D; D=""; parents=""
} }' | sed "s:\t\t:\t:g"
}
############################
Heatmap2panmap()
{
   heatmap=$1
   cat $heatmap | awk 'function sortuniq (str){nA=split(str,A,","); if(nA==1){return str} ;S=""; asort(A, B) ; for(a=1;a<=nA;++a){ if(B[a] !=0 && B[a] != s) {S=S","B[a]; s=B[a] }; s=B[a]};S=substr(S,2);if(S==""){S="-"} ; return S } 
{
   if(NR==1)
   {
      print $0
   }
   else
   {
      for(i=1;i<=NF;++i)
      {
         if(i==1){S=$i}
         else if(i==2 || i ==3 ) {S=S"\t"$i}
         else if(i==4) { S=S"\t"$i; nParents=split($4,Parents,",") }
         else
         {
            s=""
            nA=split($i,A,",")
            for(a=1;a<=nA;++a)
            {
               s=s","Parents[A[a]]
            }
            s=substr(s,2)
            if(s==""){s="-"}
            S=S"\t"sortuniq(s)
         }
      }
      print S
   }
}'
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec4#######################################################################################################################################
############################################################################################################################################
panmapSplit(){
if [ -d "panmapCHRs" ]; then
   echo "ERROR: panmapCHRs already exists"
else
   cat $1 | awk '{if(NR==1){system("mkdir panmapCHRs"); H=$0}else{if($1!=x){print H > "panmapCHRs/"$1".panmap"; print $0 > "panmapCHRs/"$1".panmap" }else{print $0 > "panmapCHRs/"x".panmap" } ;x=$1 } }'
   cat "$1".fa | fasta2txt | tr '_' '\t' | awk '{print ">"$1"_"$2"_"$3"\n"$4 > "panmapCHRs/"$1".panmap.fa" }'
   echo "panmapCHRs was generated"
fi
}
############################
panmapSampleIDsfromList(){
tmpFile0202=$(mktemp "./KhufuEnviron.XXXXXXXXX")
id1=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
id2=$(echo $2 | sed "s:.*/::g" | sed "s:.list::g" )
cat $1 | head -1 | cut -f 1-4 | tr '\t' '\n' | cat - <(cat $2) > "$tmpFile0202"
reorderDSfromList  $1 "$tmpFile0202" > "$id1".reordered.from"$id2".panmap
cat "$1".fa > "$id1".reordered.from"$id2".panmap.fa
echo ""$id1".reordered.from"$id2".panmap & "$id1".reordered.from"$id2".panmap.fa were generated"
rm "$tmpFile0202"
trap "rm -f $tmpFile0202" EXIT
}
############################
panmapReAssignSampleIDsfromList(){
hapmap=$1
list=$2
tmpDir0207=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
id1=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
id2=$(echo $2 | sed "s:.*/::g" | sed "s:.list::g" )
cat $1 | head -1 | cut -f 1-4 | tr '\t' '\n' | cat - <(cat $2 | cut -f 1 ) > "$tmpDir0207"/list1
cat $1 | head -1 | cut -f 1-4 | tr '\t' '\n' | cat - <(cat $2 | cut -f 2 ) | tr '\n' '\t' | sed "s:$:\n:g" > "$tmpDir0207"/list2
reorderDSfromList  $1 "$tmpDir0207"/list1 | sed 1d | sed "1i$(cat "$tmpDir0207"/list2 )" > "$id1".reassigned.from"$id2".panmap
cat "$1".fa > "$id1".reassigned.from"$id2".panmap.fa
echo ""$id1".reassigned.from"$id2".panmap & "$id1".reassigned.from"$id2".panmap.fa were generated"
rm -rf $tmpDir0207
trap "rm -rf $tmpDir0207" EXIT
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec5#######################################################################################################################################
############################################################################################################################################
panmapCalculateMissingVariant(){
cat $1 | sed -E "1 s:[^\t]+:-:g" | cut -f 5- | awk '{if(NR==1){len=NF; x="miss" }else{x=gsub("-","",$0)/len }; print x }'
}
############################
panmapCalculateMissingSample(){
cat $1 | awk '{if(NR==1){split($0,H,"\t"); str="" }else{for(i=5;i<=NF;++i){ if($i=="-"){S[i]+=1}  } } }END{ for(i=5;i<=NF;++i){ str=str","S[i]}  }END{str=substr(str,2); nA=split(str,A,","); for(a=1;a<=nA;++a){printf "%s\t%0.2f\n", H[a+4],A[a]/(NR-1) } } '
}
############################
panmapCalculateMAF(){
tmpFile0203=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0203"
panmap2hapmap $1 | cut -f $(cat "$tmpFile0203" | head -1 | cut -f 4 | tr ',' ' ' | wc -w | awk '{print "1-2,"$0+2+1"-" }') | hapmapCalculateMAF
rm "$tmpFile0203"
trap "rm -f $tmpFile0203" EXIT
}
############################
panmapFilterMissingVariant(){
miss=$2
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
tmpFile0204=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0204"
cat "$tmpFile0204" | panmapCalculateMissingVariant | paste - "$tmpFile0204" | awk -v miss=$miss '{if(NR==1 || $1 <= miss) print $0}' | cut -f 2- > "$id".miss"$miss".panmap
echo ""$id".miss"$miss".panmap was generated"
panmapGetFasta "$id".miss"$miss".panmap $1
rm "$tmpFile0204"
trap "rm -f $tmpFile0204" EXIT
}
############################
panmapFilterMissingSample(){
Smiss=$2
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
tmpDir0208=$(mktemp -d "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0208"/panmap
panmapCalculateMissingSample "$tmpDir0208"/panmap | awk -v Smiss=$Smiss '{if($2<=Smiss){print $1}}'  > "$tmpDir0208"/Smiss"$Smiss"
cat "$tmpDir0208"/panmap | head -1 | cut -f 1-4 | tr '\t' '\n' | cat - <(cat "$tmpDir0208"/Smiss"$Smiss") > "$tmpDir0208"/Smiss"$Smiss".2
reorderDSfromList  "$tmpDir0208"/panmap "$tmpDir0208"/Smiss"$Smiss".2  > "$id".Smiss"$Smiss".panmap
cat $1".fa" > "$id"Smiss"$Smiss".panmap.fa 
echo ""$id".Smiss"$Smiss".panmap & "$id".Smiss"$Smiss".panmap.fa were generated"
rm -rf "$tmpDir0208"
trap "rm -rf $tmpDir0208" EXIT
}
############################
panmapFilterMAF(){
maf=$2
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
tmpFile0204=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0204"
panmapCalculateMAF "$1" | paste - "$tmpFile0204" | awk -v maf=$maf '{if(NR==1 || $2 >= maf){ print $0} }' | cut -f 3- > "$id".MAF"$maf".panmap
echo ""$id".MAF"$maf".panmap was generated"
panmapGetFasta "$id".MAF"$maf".panmap $1
rm "$tmpFile0204"
trap "rm -f $tmpFile0204" EXIT
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec6#######################################################################################################################################
############################################################################################################################################
panmapGetAlleles(){
panmapAlleleFreqStats $1 | cut -f 1 
}
############################
panmapAlleleFreqStats(){
cat "$1".fa | fasta2txt | tr '_' '\t' | awk '{print $1"_"$2"\t"$4}' | fastatxtAlleleCombine | cut -f 2 | sed "1ialleles" | paste  - <(cat $1 | cut -f 5- | tr '\t' ',' ) | sed 1d |awk '{nA=split($1,A,",");nB=split($2,B,","); sum=0;S="-" ; for(a=1;a<=nA;++a){s=0;for(b=a;b<=nB;++b){ if(a==B[b]){s+=1; sum+=1};  };S=S","s }; S=substr(S,3); print $1"\t"S"\t"sum }' | awk '{nA=split($2,A,",");S="-"; for(a=1;a<=nA;++a){ if($3>0){s=A[a]/$3}else(s=0); S=S","s }; S=substr(S,3); print  $1"\t"S"\t"$2"\t"$3  }' | sed "1ialleles\tfreq\tcnt\tsum" 
}
############################
panmapAlleleTypeFreq(){
tmpFile0205=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > $tmpFile0205
panmapAlleleFreqStats $1 | paste - $tmpFile0205 | awk '{nA=split($1,A,","); if(nA==2){print $0} }' | cut -f 1,5,6,9- | awk '{sum=0; homo1=0;homo2=0;het=0; for(i=4;i<=NF;++i){ if($i==1){homo1+=1;sum+=1}; if($i==2){homo2+=1;sum+=1}; if($i~"[0-9]+,[0-9]+"){het+=1;sum+=1}  }; if(sum>0){s=homo1/sum","het/sum","homo2/sum}else{s="0,0,0"}; print $2"\t"$3"\t"$1"\t"s  }' |  sed "1ichr\tpos\tallele1,allele2\tA,B,H"
rm "$tmpFile0205"
trap "rm -f $tmpFile0205" EXIT
}
############################
panmapNumParents(){
cat $1 | head -1 | cut -f 4 | tr ',' '\t' | wc -w
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################



############################################################################################################################################
#Sec7#######################################################################################################################################
############################################################################################################################################
panmapGetIDs(){
cat $1 | head -1 | cut -f 5- | tr '\t' '\n'
}
############################
panmapGetFreqSampleAlleles(){
 cat $1 | cut -f 1-2,5- | awk -v id=$2 '{if(NR==1){inx=0; for(i=3;i<=NF;++i){if($i==id){inx=i}} }else{S=$1"\t"$2"\t"$inx; X=0; Y=0; for(i=3;i<=NF;++i){ if(i!=inx && $i != "-"){Y+=1;if($i==$inx){X+=1}} }; if(Y>0){freq=X/Y}else{freq=0}; print S,freq,Y }}' | sed "1ichr\tpos\tallele\tfreq\tdep"
}
############################
panmapGetUniqueSampleAlleles(){
tmpFile0206=$(mktemp "./KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpFile0206"
panmap2hapmap $1 | cut -f $(cat "$tmpFile0206"  | head -1 | cut -f 4 | tr ',' '\t' | wc -w | awk '{print "1-2,"$0+2+1"-"}') | hapmapGetUniqueSampleAlleles - $2
rm "$tmpFile0206"
trap "rm -f $tmpFile0206" EXIT
}
############################
panmapGetROH(){
cat $1 | cut -f 1-2,5- | hapmapGetROH - $2
}
############################
panmapLowFreqAlleleMask(){
mask=$2
panmap=$1
tmpDir0209=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
id=$(echo $panmap | sed "s:.*/::g" | sed "s:.panmap::g" )
cat $panmap > "$tmpDir0209"/panmap
panmap2hapmap $panmap > "$tmpDir0209"/hapmap
inx=$(cat "$tmpDir0209"/panmap  | head -1 | cut -f 4 | tr ',' '\t' | wc -w )
cat "$tmpDir0209"/hapmap | cut -f 3-$((inx+2)) --complement > "$tmpDir0209"/hapmap2A
cat "$tmpDir0209"/hapmap | cut -f 3-$((inx+2)) > "$tmpDir0209"/hapmap2B
hapmapLowFreqAlleleMask "$tmpDir0209"/hapmap2A $mask > "$tmpDir0209"/hapmap3
paste "$tmpDir0209"/hapmap3 "$tmpDir0209"/hapmap2B > "$tmpDir0209"/"$id".masked"$mask"
hapmap2panmap "$tmpDir0209"/"$id".masked"$mask" $(cat "$tmpDir0209"/panmap | cut -f 4 | head -1 )
rm -rf $tmpDir0209
trap "rm -rf $tmpDir0209" EXIT
}
############################
panmapGetHomoPolymorphic(){
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
tmpDir0210=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
cat $1 > "$tmpDir0210"/panmap
cat "$tmpDir0210"/panmap | cut -f 5- | awk '{if($0~","){print 0}else{print 1}}' | paste - "$tmpDir0210"/panmap | awk '{if($1==1) print $0}' | cut -f 2- > "$tmpDir0210"/panmap2
cat "$tmpDir0210"/panmap2 | awk '{if(NR==1){print $0}else {x=0; for(i=5;i<=NF;++i){if(X!= "" && $1!="-" && X!=$i ){x=1}; X=$i}; if(x==1){print $0} } }' > "$id".HomoPolymorphic.panmap
echo ""$id".HomoPolymorphic.panmap was generated"
panmapGetFasta "$id".HomoPolymorphic.panmap $1
rm -rf $tmpDir0210
trap "rm -rf $tmpDir0210" EXIT
}
############################
panmapSelParents(){
parents=$2
id=$(echo $1 | sed "s:.*/::g" | sed "s:.panmap::g" )
tmpDir0211=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
N=$(echo $2 | tr ',' '\t' | wc -w)
panmap2hapmap $panmap > "$tmpDir0211"/"$id".Sel_"$N"
hapmap2panmap "$tmpDir0211"/"$id".Sel_"$N" $2
rm -rf $tmpDir0211
trap "rm -rf $tmpDir0211" EXIT
}
############################
panmapGetConsensusCall(){
   panmap2hapmap $1 | hapmapGetConsensusCall
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
