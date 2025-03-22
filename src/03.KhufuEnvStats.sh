# source /cluster/projects/khufu/korani_projects/KhufuEnv/KhufuEnvStats.sh
################################################################################################################################
###general stats:
############################################################################################################################################
#Sec1#######################################################################################################################################
############################################################################################################################################
avg(){
indx=$2
if [ -z "$indx" ]; then indx=1; fi
cat $1 | Rscript -e 'args=commandArgs(T);A = data.table::data.table(read.table("stdin",header=TRUE)); indx=as.numeric(args[1]); cols = seq(1,ncol(A)); cols = cols[-indx]; fac=colnames(A)[indx]; B = A[, sapply(.SD, function(x) list(mean(x))), .SDcols = cols, by = fac]; C = subset(B,select=match(colnames(A),colnames(B))); write.table(C,stdout(),col.names=T,row.names=F,sep="\t", quote=F)' $indx
}
############################
std(){
indx=$2
if [ -z "$indx" ]; then indx=1; fi
cat $1 | Rscript -e 'args=commandArgs(T);A = data.table::data.table(read.table("stdin",header=TRUE)); indx=as.numeric(args[1]); cols = seq(1,ncol(A)); cols = cols[-indx]; fac=colnames(A)[indx]; B = A[, sapply(.SD, function(x) list(sd(x))), .SDcols = cols, by = fac]; C = subset(B,select=match(colnames(A),colnames(B))); write.table(C,stdout(),col.names=T,row.names=F,sep="\t", quote=F)' $indx
}
############################
sum(){
indx=$2
if [ -z "$indx" ]; then indx=1; fi
cat $1 | Rscript -e 'args=commandArgs(T);A = data.table::data.table(read.table("stdin",header=TRUE)); indx=as.numeric(args[1]); cols = seq(1,ncol(A)); cols = cols[-indx]; fac=colnames(A)[indx]; B = A[, sapply(.SD, function(x) list(sum(x))), .SDcols = cols, by = fac]; C = subset(B,select=match(colnames(A),colnames(B))); write.table(C,stdout(),col.names=T,row.names=F,sep="\t", quote=F)' $indx
}
############################
max(){
indx=$2
if [ -z "$indx" ]; then indx=1; fi
cat $1 | Rscript -e 'args=commandArgs(T);A = data.table::data.table(read.table("stdin",header=TRUE)); indx=as.numeric(args[1]); cols = seq(1,ncol(A)); cols = cols[-indx]; fac=colnames(A)[indx]; B = A[, sapply(.SD, function(x) list(max(x))), .SDcols = cols, by = fac]; C = subset(B,select=match(colnames(A),colnames(B))); write.table(C,stdout(),col.names=T,row.names=F,sep="\t", quote=F)' $indx
}
############################
min(){
indx=$2
##
if [ -z "$indx" ]; then indx=1; fi
cat $1 | Rscript -e 'args=commandArgs(T);A = data.table::data.table(read.table("stdin",header=TRUE)); indx=as.numeric(args[1]); cols = seq(1,ncol(A)); cols = cols[-indx]; fac=colnames(A)[indx]; B = A[, sapply(.SD, function(x) list(min(x))), .SDcols = cols, by = fac]; C = subset(B,select=match(colnames(A),colnames(B))); write.table(C,stdout(),col.names=T,row.names=F,sep="\t", quote=F)' $indx
}
############################
med(){
indx=$2
if [ -z "$indx" ]; then indx=1; fi
cat $1 | Rscript -e 'args=commandArgs(T);A = data.table::data.table(read.table("stdin",header=TRUE)); indx=as.numeric(args[1]); cols = seq(1,ncol(A)); cols = cols[-indx]; fac=colnames(A)[indx]; B = A[, sapply(.SD, function(x) list(median(x))), .SDcols = cols, by = fac]; C = subset(B,select=match(colnames(A),colnames(B))); write.table(C,stdout(),col.names=T,row.names=F,sep="\t", quote=F)' $indx
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
avgCol () {
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,2,mean, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' | tr '\n' '\t' | sed "s:$:\n:g" | sed "s:\t$::g"
}
############################
stdCol () {
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,2,sd, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' | tr '\n' '\t' | sed "s:$:\n:g" | sed "s:\t$::g"
}
############################
sumCol () {
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,2,sum, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' | tr '\n' '\t' | sed "s:$:\n:g"  | sed "s:\t$::g"
}
############################
maxCol(){
cat $1 | awk '{for(i=1;i<=NF;++i){ if(NR==1){Z[i] = $i} else if($i > Z[i]) {Z[i] = $i} } }END{ for(x=1;x<=NF;++x) {S=S"\t"Z[x]} ;S=substr(S,2); print(S); S=""  }' | sed "s:\t$::g"
}
############################
minCol(){
cat $1 | awk '{for(i=1;i<=NF;++i){ if(NR==1){Z[i] = $i} else if($i < Z[i]) {Z[i] = $i}  } }END{ for(x=1;x<=NF;++x) {S=S"\t"Z[x]} ;S=substr(S,2); print(S); S=""  }' | sed "s:\t$::g"
}
############################
medCol(){
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,2,median, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' | tr '\n' '\t' | sed "s:$:\n:g" | sed "s:\t$::g"
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
avgRow () {
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,1,mean, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' 
}
############################
stdRow () {
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,1,sd, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' 
}
############################
sumRow () {
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,1,sum, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' 
}
############################
maxRow(){
cat $1 | awk '{z=$1; for(i=2;i<=NF;++i){ if($i > z){z=$i} }; print z}'
}
############################
minRow(){
cat $1 | awk '{z=$1; for(i=2;i<=NF;++i){ if($i < z){z=$i} }; print z}'
}
############################
medRow(){
cat $1 | Rscript -e 'args = commandArgs(trailingOnly=TRUE); A=read.table("stdin",header=F); B=apply(A,1,median, na.rm = TRUE); write.table(B,stdout(),col.names=F,row.names=F,quote=F) ' 
}
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
