recluster.haplotypes<-function(fasta,pop=NULL){
number<-c(1:length(fasta))
len<-NULL
for (n in 1:length(fasta)){
	leng1<-length(which(fasta[[n]]=="-"))
	leng2<-length(which(fasta[[n]]=="N"))
	leng3<-length(which(fasta[[n]]=="n"))
	leng4<-length(which(fasta[[n]]=="f0"))
	len[n]<-length(fasta[[n]])-(leng1+leng2+leng3+leng4)
	}
matrdataor<-number[order(len,decreasing = T)]
fastaor<-fasta[order(len,decreasing = T)]
res<-NULL
matrdata<-matrix(0,length(fasta),6)
rownames(matrdata)<-names(fastaor)
matrdata[,1]<-matrdataor
sequ<-fastaor[1]
vector<-NULL
vectnames<-NULL
vecthaplotypes<-NULL
vectnames<-names(fastaor[1])
vecthaplotypes[1]<-1
pos<-1
for (n in 2:length(matrdataor)){
	check<-0
	for (giro in 1:(n-1)){
	fastacheck<-c(fastaor[giro],fastaor[n])
	distn<-dist.dna(fastacheck, model = "raw",pairwise.deletion = TRUE)
	if(distn==0){
	vecthaplotypes[n]<-vecthaplotypes[giro]
	vectnames[n]<-vectnames[giro]
	check<-1
	break
	}
	}
	if(check==0){
	pos<-pos+1
	sequ<-c(sequ, fastaor[n])
	vecthaplotypes[n]<-pos
	vectnames[n]<-names(fastaor[n])
	#print(n)
	}
}
res$seqs<-sequ
res$data_final<-cbind(matrdata[,1],vecthaplotypes,vectnames)
res$frequency<-aggregate(rep(1,length(vecthaplotypes))~vecthaplotypes,FUN="sum")[,2]
if(!is.null(pop)){
	popu<-unique(pop)
	haplo<-unique(vecthaplotypes)
	matrpop<-matrix(0,length(popu),length(haplo))
	for(spec in 1:length(vecthaplotypes)){
		#spec<-1
		matrpop[which(popu==pop[spec]),which(haplo==vecthaplotypes[spec])]<-matrpop[which(popu==pop[spec]),which(haplo==vecthaplotypes[spec])]+1
		rownames(matrpop)<-popu
		colnames(matrpop)<-haplo
	}
res$populations<-matrpop
}
return(res)
}