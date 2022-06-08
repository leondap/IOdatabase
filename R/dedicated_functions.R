
genetic.regionalization<-function(fasta,spec,longi,lati, minimum=3, minarea=2, nidw=15, gsttab=NULL, GST="all", alpha=0.1, power=1, powerW=0, method="mean", methodclust="ward.D"){
res<-NULL
species<-unique(spec)
specloc<-aggregate(rep(1,length(longi))~longi+lati+spec,FUN="sum")
specarea<-aggregate(rep(1,nrow(specloc))~specloc[,3],FUN="sum")
speciessel<-specarea[which(specarea[,2]>=minarea),1]
sel<-which(spec%in%speciessel)

longi<-longi[sel]
lati<-lati[sel]
fasta<-fasta[sel]
species<-speciessel
spec<-spec[sel]

loca<-paste(longi,lati)
areeun<-unique(loca)
areetot<-unique(cbind(longi,lati))
aree<-cbind(c(1:nrow(areetot)),areetot)
loc<-NULL
for(ar in 1:nrow(aree)){
#ar<-1
wh<-which(loca==areeun[ar])
loc[wh]<-ar
}


matrixGst<-array(NA,dim=c(nrow(aree),nrow(aree),length(species)))
areatabledissGst<-array(NA,dim=c(nrow(aree),nrow(aree)))
for (sp in 1:length(species)){
#sp<-2
	sel<-which(spec==species[sp])
	fastared<-fasta[sel]
	datamatsp<-loc[sel]
	#cbind(longi,lati)[sel,]
	dismatsp<-as.matrix(dist.dna(fastared, model = "raw",pairwise.deletion = TRUE))
	areamatsp<-unique(datamatsp)
	if(is.null(gsttab)){
	gstsp<-recluster.fst(as.dist(dismatsp),datamatsp,setzero=T,setnazero=T)$Gst
	}else{
	gstsp<-gsttab[which(names(gsttab)==species[sp])]
	}
	if(GST=="all"){
		wei<-1
	}
	if(GST=="gd"){
		wei<-(1-gstsp)
	}
	if(GST=="st"){
		wei<-(gstsp)
	}
	if (length(areamatsp)>1){
	for (n in 1: (length(areamatsp)-1)){
		#n<-1
		for (m in (n+1):length(areamatsp)){
		#m<-7			
			dista<-dismatsp[c(which(datamatsp==areamatsp[n]), which(datamatsp==areamatsp[m])),c(which(datamatsp==areamatsp[n]), which(datamatsp==areamatsp[m]))]
			popu<-c(rep(1,length(which(datamatsp==areamatsp[n]))),rep(2,length(which(datamatsp==areamatsp[m]))))
			value<-((recluster.fst(dista,popu,setzero=T,setnazero=T)$Ht)^power)*(wei)
			matrixGst[areamatsp[n], areamatsp[m],sp]<-value
			matrixGst[areamatsp[m], areamatsp[n],sp]<-value
			}
		}
	}
}

count<-areatabledissGst
for (r in 1: nrow(aree)){
	for (c in 1: nrow(aree)){
		count[r,c]<-length(which(!is.na(matrixGst[r,c,])))
		if(count[r,c]>1){
			if(method=="mean"){
			areatabledissGst[r,c]<-mean(matrixGst[r,c,],na.rm=T)
			areatabledissGst[c,r]<-mean(matrixGst[r,c,],na.rm=T)
			}
			if(method=="median"){
			areatabledissGst[r,c]<-median(matrixGst[r,c,],na.rm=T)
			areatabledissGst[c,r]<-median(matrixGst[r,c,],na.rm=T)
			}
		}
	}
}

matnew<-areatabledissGst

rich<-aggregate(rep(1,length(loc))~loc+spec,FUN="sum")
richness<-aggregate(rep(1,nrow(rich))~rich[,1],FUN="sum")
sitiguida<-which(richness[,2]>=minimum)
guida<-rep(0,nrow(richness))
guida[sitiguida]<-1
#siti<-aree[sitiguida[,1],]
areeg<-cbind(aree,richness[,2],guida)
distance<-rep(1,nrow(areeg))
aree2<-NULL
k<-1
for (qua in 1:nrow(matnew)){
#qua<-1
tabella<-cbind(areeg,matnew[,qua])
tabella<-tabella[-qua,]
tabella<-tabella[which(tabella[,5]==1),]
if(length(which(!is.na(tabella[,6])))>=nidw){
tabella<-tabella[complete.cases(tabella),]
res1 <- idw(tabella[,6], tabella[,2:3],areeg[,2:3],p=2)
distance<-cbind(distance,res1)
aree2[k]<-as.character(areeg[qua,1])
k<-k+1
}
}

aree2<-as.numeric(aree2)
distance<-distance[aree2,]
distance<-distance[,-1]
areefin<-sitiguida[which(sitiguida %in% aree2)]
#nrow(distance)
#rownames(distance)
sel<-which(rownames(distance)%in%areefin)
dist<-distance[sel,sel]
mat<-matrix(NA,nrow(dist),nrow(dist))
for(row in 1:nrow(mat)){
#row<-1
#col<-2
for(col in 1:ncol(mat)){
mat[row,col]<-(dist[row,col]+dist[col,row])/2
}}

dista<-as.dist(mat)
select<-as.numeric(rownames(dist))

coord<-areeg[select,]


distgeo<-earth.dist(coord[,2:3])


res$richness<-coord[,4]
clust<-hclustgeo2(dista,distgeo,alpha=alpha,w=coord[,4]^powerW,method=methodclust)
res$tree <- clust$solution
res$weightmatrix<-clust$weighteddiss
res$dist<-dista
res$distgeo<-distgeo
res$coord<-coord[,2:3]
res$originalcoord<-cbind(longi,lati)
res$originalsequences<-fasta
res$originalspec<-spec
res$minimum<-minimum
res$minarea<-minarea
res$gsttab<-gsttab
res$GST<-GST
res$alpha<-alpha
res$power<-power
res$powerW<-powerW
res$method<-method
return(res)
}


recluster.boot.regionalisation<-function (reg, boot=100,maxcl=8,nidw=13) {
minimum<-reg$minimum
minarea<-reg$minarea
gsttab<-reg$gsstab
alpha<-reg$alpha
GST=reg$GST
sequences<-reg$originalsequences
spec<-reg$originalspec
species<-unique(spec)
tree<-reg$tree
power<-reg$power
powerW<-reg$powerW
method<-reg$method
names<-as.data.frame(paste(reg$coord[,1],reg$coord[,2]))
tab<-names
restab<-NULL

for(c in 2:maxcl){
#c<-2
restab[[c-1]]<-matrix(NA, boot,c)
classification<-cutree(tree,c)
tab<-cbind(tab,classification)
}


for (i in 1 : boot){
xs<-species[sample(1:length(species),replace=T)]
select<-which(spec %in% xs)
sequences2<-sequences[select]
membership2<-spec [select]
lat2<-lat[select]
long2<-long[select]
treernd<-genetic.regionalization(sequences2,membership2,long2,lat2,minarea=minarea,minimum=minimum,nidw=nidw, GST=GST,gsttab=gsttab,alpha=alpha,powerW=powerW, power=power, method=method)
names<-paste(treernd$coord[,1],treernd$coord[,2])
match<-match(names,tab[,1])
nogood<-which(is.na(match))
if (length(nogood)>0){
match2<-match[-nogood]
names2<-names[-nogood]
}else{
match2<-match
names2<-names
}
tabrnd<-matrix(NA,nrow(tab),ncol(tab))
for(c in 2:maxcl){
#c<-2
classrnd<-cutree(treernd$tree,c)
if (length(nogood)>0){
tabrnd[match2,c]<-classrnd[-nogood]
}else{
tabrnd[match2,c]<-classrnd
}
}
tabsel<-tab[which(tab[,1] %in% names2),]
for(v in 2:maxcl){
#v<-7
tabrnd[match2,1]<-names2
somma<-aggregate(rep(1,nrow(tabsel))~tabsel[,v],FUN="sum")
vect<-NULL
for (gr in 1:nrow(somma)){
#gr<-1
quali<-tabsel[which(tabsel[,v]==somma[gr,1]),1]
tutti<-tabrnd[which(tabrnd[,1]%in%quali),v]
tot<-aggregate(rep(1,length(tutti))~tutti, FUN="sum")
vect[gr]<-max(tot[,2])/length(which(tabsel[,v]==somma[gr,1]))
}
restab[[v-1]][i,]<-vect
}
}
final<-NULL
for(giro in 1:(maxcl-1)){
final[[giro]]<-colMeans(restab[[giro]])
}
return(final)
}

hclustgeo2<-function (D0, D1 = NULL, alpha = 0, scale = TRUE, wt = NULL, method="ward.D") 
{
res<-NULL    
if (class(D0) != "dist") 
        stop("DO must be of class dist (use as.dist)", 
            call. = FALSE)
    if (!is.null(D1) && (class(D1) != "dist")) 
        stop("D1 must be of class dist (use as.dist)", 
            call. = FALSE)
    n <- as.integer(attr(D0, "Size"))
    if (is.null(n)) 
        stop("invalid dissimilarities", call. = FALSE)
    if (is.na(n) || n > 65536L) 
        stop("size cannot be NA nor exceed 65536", call. = FALSE)
    if (n < 2) 
        stop("must have n >= 2 objects to cluster", call. = FALSE)
    if (!is.null(D1) && length(D0) != length(D1)) 
        stop("the two dissimilarity structures must have the same size", 
            call. = FALSE)
    if ((max(alpha) > 1) || (max(alpha) < 0)) 
        stop("Values alpha must be in [0,1]", call. = FALSE)
    if ((scale == TRUE) && (!is.null(D1))) {
        D0 <- D0/max(D0)
        D1 <- D1/max(D1)
    }
    delta0 <- wardinit(D0, wt)
    if (!is.null(D1)) 
        delta1 <- wardinit(D1, wt)
    else delta1 <- 0
    delta <- (1 - alpha) * delta0 + alpha * delta1
    res$solution <- hclust(delta, method = method, members = wt)
	res$weighteddiss<-delta
    return(res)
}







slicesfun<-function(location, identification, sequences, coordinates, populations, range=c(0,90), GST=NULL, partition=c("ALL","GD","ST"), correct=F, width=3,asympt=FALSE, mini_sp=5, iter=10){
res<-NULL
interval<-width/2
start<-range[1]
slices<-round(((range[2]-range[1])/interval)+1,0)
matrixsp<-array(NA, dim=c(100000,14,iter))
dataset<-cbind(location,identification,coordinates,populations)
species<-unique(identification)
matricerichness<-array(0,dim=c(length(species),slices+1,iter))
rownames(matricerichness)<-species
matricehaplotypes<-array(0,dim=c(length(species),slices+1,iter))
rownames(matricehaplotypes)<-species
matricehaplotypeslow<-array(0,dim=c(length(species),slices+1,iter))
rownames(matricehaplotypeslow)<-species
matricedistances<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matricedistances)<-species
matriceGD<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matriceGD)<-species
matriceGDred<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matriceGD)<-species
matriceST<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matriceGD)<-species
matriceSTred<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matriceGD)<-species
matriceALL<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matriceGD)<-species
matriceALLred<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matriceGD)<-species
matricedistancesGD<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matricedistancesGD)<-species
matricedistancesST<-array(NA,dim=c(length(species),slices+1,iter))
rownames(matricedistancesST)<-species

for (it in 1:iter){
#it<-1
for (cut in 0:slices){
#cut<-1
lati<-start+interval*cut
stripe<-which(dataset[,1]>=lati-interval & dataset[,1]<lati+interval)
if(length(stripe)>0){
datalat<-dataset[stripe,]
speciessp<-aggregate(rep(1,nrow(datalat))~ datalat[,2], FUN="sum")
for(k in 1:nrow(speciessp)){
	#k<-1
spe<-which(rownames(matricerichness)==speciessp[k,1])
matricerichness[spe,cut+1,it]<-speciessp[k,2]
}
}
}
line<-1
for (spe in 1:length(species)){
#spe<-2
dataspeq<-dataset[which(dataset[,2]==species[spe]),]
seqspeq<-sequences[which(dataset[,2]==species[spe])]
stripes<-which(matricerichness[spe,,it]>=mini_sp)
if(length(stripes)>0){
minimum<-min(matricerichness[spe,,it][stripes])
max<-max(matricerichness[spe,,it])
for (strip in 1:length(stripes)){
#strip<-2
lati<-start+interval*(stripes[strip]-1)
stripe<-which(dataspeq[,1]>=lati-interval & dataspeq[,1]<lati+interval)
seqslat<-seqspeq[stripe]
datalat<-dataspeq[stripe,]

baricenter<-c(mean(datalat[,3]),mean(datalat[,4]))
coordin<-rbind(baricenter,datalat[,3:4])
qualiseqs<-sample(1:length(seqslat))[1:minimum]
baricenter<-c(mean(datalat[,3]),mean(datalat[,4]))
coordin<-rbind(baricenter,datalat[,3:4])
molt<-1

same<-as.vector(dist(datalat$populations))
same[which(same>0)]<-1

regi<-unique(datalat$populations)
regi<-regi[!is.na(regi)]

GDdis<-NULL
numero<-0
for(re in 1:length(regi)){
#re<-1
piglia<-datalat[which(datalat$populations==regi[re]),]
if(nrow(piglia)>1){
numero<-numero+nrow(piglia)
baricenterGD<-c(mean(piglia[,3]),mean(piglia[,4]))
coordinGD<-rbind(baricenterGD,piglia[,3:4])
distanzGD<-as.matrix(earth.dist(coordinGD, dist = TRUE))[1,]
distanzGD<-distanzGD[2:length(distanzGD)]
GDdis<-c(GDdis,distanzGD)
}
}


if(!is.null(GST)){
gst<-GST[which(names(GST)==species[spe])]
if(partition=="ALL"){
molt<-1
}
if(partition=="GD"){
molt<-1-gst
}
if(partition=="ST"){
molt<-gst
}
}

tutte<-dist.dna(seqslat, model = "raw",pairwise.deletion = TRUE)
tuttered<-(dist.dna(seqslat[qualiseqs], model = "raw",pairwise.deletion = TRUE))

matriceALL[spe,stripes[strip],it]<-(mean(tutte))*molt
matriceALLred[spe,stripes[strip],it]<-(mean(tuttered))*molt
if (correct==T){
matriceALL[spe,stripes[strip],it]<-sqrt((tutte)*(length(seqslat)/(length(seqslat)-1)))*molt
matriceALLred[spe,stripes[strip],it]<-(mean(tuttered)*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
}

intra<-which(same==0)
matriceGD[spe,stripes[strip],it]<-(mean(as.vector(tutte)[intra]))*molt
if (correct==T){
matriceGD[spe,stripes[strip],it]<-(mean(as.vector(tutte)[intra])*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
}

distanz<-as.matrix(earth.dist(coordin, dist = TRUE))[1,]
distanz<-distanz[2:length(distanz)]

matricedistancesGD[spe,stripes[strip],it]<-sqrt((sum(GDdis^2)/numero))

inter<-which(same==1)
matriceST[spe,stripes[strip],it]<-(mean(as.vector(tutte)[inter])*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
if (correct==T){
matriceST[spe,stripes[strip],it]<-(mean(as.vector(tutte)[inter])*(length(qualiseqs)/(length(qualiseqs)-1)))*molt
}

matricehaplotypeslow[spe,stripes[strip],it]<-length(recluster.haplotypes(seqslat[qualiseqs])$frequency)*molt
matricedistances[spe,stripes[strip],it]<-sqrt(mean(distanz^2))
matrixsp[line,1,it]<-species[spe]
matrixsp[line,2,it]<-matricedistances[spe,stripes[strip],it]
matrixsp[line,3,it]<-matriceALL[spe,stripes[strip],it]
matrixsp[line,4,it]<-matricehaplotypeslow[spe,stripes[strip],it]
matrixsp[line,11,it]<-lati
matrixsp[line,7,it]<-length(seqslat)
matrixsp[line,8,it]<-length(qualiseqs)
matrixsp[line,9,it]<-baricenter[1]
matrixsp[line,10,it]<-baricenter[2]
matrixsp[line,6,it]<-matriceALLred[spe,stripes[strip],it]
matrixsp[line,12,it]<-matricedistancesGD[spe,stripes[strip],it]
matrixsp[line,13,it]<-matriceGD[spe,stripes[strip],it]
matrixsp[line,14,it]<-matriceST[spe,stripes[strip],it]


if(asympt){
if(length(freq)>1){
accu<-iNEXT(freq, q=0, datatype="abundance", size=NULL, endpoint=max, knots=40, se=TRUE, conf=0.95, nboot=50)
matricehaplotypes[spe,stripes[strip],it]<-accu$AsyEst[1,2]*molt
matrixsp[line,5]<-matricehaplotypes[spe,stripes[strip]]
}
if(length(freq)==1){
matricehaplotypes[spe,stripes[strip],it]<-1*molt
matrixsp[line,5,it]<-matricehaplotypes[spe,stripes[strip],it]
}
}
line<-line+1
}
}
}
}

diversity<-matrix(NA,iter,ncol(matricehaplotypes))
speciesnum<-diversity
latitude<-diversity
distances<-diversity
reducedhap<-diversity
quanteslice<-diversity
reducedGD<-diversity
reducedhapsd<-diversity
reducedGDsd<-diversity
bar_long<-diversity
bar_lat<-diversity
centroid<-diversity
reducedALL<-diversity
reducedALLred<-diversity
reducedALLredsd<-diversity
reducedALLsd<-diversity
reducedST<-diversity
reducedSTsd<-diversity
distancesGD<-diversity




quanteslice<-NULL
for (slice in 1:ncol(matricehaplotypeslow)){
#slice<-1
	quali<-which(matricehaplotypeslow[,slice,1]>0)
	quanteslice[slice]<-length(quali)
	if(quanteslice[slice]>=mini_sp){
	for(it in 1:iter){
	diversity[it,slice]<-mean(matricehaplotypes[quali,slice,it],na.rm=T)
	speciesnum[it,slice]<-length(quali)
	latitude[it,slice]<-start+interval*(slice-1)
	distances[it,slice]<-mean(matricedistances[quali,slice,it],na.rm=T)
	reducedhap[it,slice]<-mean(matricehaplotypeslow[quali,slice,it],na.rm=T)
	reducedGD[it,slice]<-mean(matriceGD[quali,slice,it],na.rm=T)
	reducedhapsd[it,slice]<-sd(matricehaplotypeslow[quali,slice,it],na.rm=T)
	reducedGDsd[it,slice]<-sd(matriceGD[quali,slice,it],na.rm=T)
	reducedALL[it,slice]<-mean(matriceALL[quali,slice,it],na.rm=T)
	reducedALLsd[it,slice]<-sd(matriceALL[quali,slice,it],na.rm=T)
	reducedALLred[it,slice]<-mean(matriceALLred[quali,slice,it],na.rm=T)
	reducedALLredsd[it,slice]<-sd(matriceALLred[quali,slice,it],na.rm=T)
	reducedST[it,slice]<-mean(matriceST[quali,slice,it],na.rm=T)
	reducedSTsd[it,slice]<-sd(matriceST[quali,slice,it],na.rm=T)
	distancesGD[it,slice]<-mean(matricedistancesGD[quali,slice,it],na.rm=T)



}
}}

ncol(matrixsp)
colnames(matrixsp)<-c("species","Distances","IGV","haplotypes_red","haplotype_asy","IGV_red","Specimens","Specimens_red",
"Bar_long","Bar_lat","slice","Distances_GD","GD","ST")
rowgood<-which(complete.cases(matrixsp[,1:2,1]))
head(matrixsp)
dista<-as.numeric(matrixsp[rowgood,2,1])
haplotypes_red<-as.numeric(matrixsp[rowgood,4,1])
IGV_red<-as.numeric(matrixsp[rowgood,6,1])
Specimens_red<-as.numeric(matrixsp[rowgood,8,1])
Bar_long<-as.numeric(matrixsp[rowgood,9,1])
Bar_lat<-as.numeric(matrixsp[rowgood,10,1])

for(it in 2:iter){
dista<-cbind(dista,as.numeric(matrixsp[rowgood,2,it]))
haplotypes_red<-cbind(haplotypes_red,as.numeric(matrixsp[rowgood,4,it]))
IGV_red<-cbind(IGV_red,as.numeric(matrixsp[rowgood,6,it]))
Specimens_red<-cbind(Specimens_red,as.numeric(matrixsp[rowgood,8,it]))
Bar_long<-cbind(Bar_long,as.numeric(matrixsp[rowgood,9,it]))
Bar_lat<-cbind(Bar_lat,as.numeric(matrixsp[rowgood,10,it]))
}

res$separate<-as.data.frame(cbind(matrixsp[rowgood,1,1],rowMeans(dista),matrixsp[rowgood,3,1],rowMeans(haplotypes_red),matrixsp[rowgood,5,1],rowMeans(IGV_red),matrixsp[rowgood,7,1],matrixsp[rowgood,8,1],
rowMeans(Bar_long),rowMeans(Bar_lat),matrixsp[rowgood,11,1],matrixsp[rowgood,12,1],matrixsp[rowgood,13,1],matrixsp[rowgood,14,1]))
colnames(res$separate)<-c("species","distances","IGV","haplotypes_red","haplotype_asy","IGV_red","Specimens","Specimens_red",
"Bar_long","Bar_lat","slice","Distances_GD","GD","ST")
final<-as.data.frame(cbind(colMeans(diversity), colMeans(reducedALL), colMeans(reducedALLsd), colMeans(reducedALLred), colMeans(reducedALLredsd), 
colMeans(reducedhap), colMeans(reducedhapsd),colMeans(latitude),colMeans(speciesnum),colMeans(distances),colMeans(reducedGD),colMeans(reducedGDsd) ,
colMeans(reducedST),colMeans(reducedSTsd),colMeans(distancesGD)))
res$aggregate<-final[complete.cases(final[,1:2]),]
return(res)
}






recluster.landscape.dist<-function(mat,units=NULL,dist,transcorr,map_alt,map1,minsea=3){
	gendistances<-as.matrix(dist)
	if(!is.null(units)){
		tabuni<-gendistances
		for(r in 1:nrow(gendistances)){
			for(c in 1:nrow(gendistances)){
				if(units[c]==units[r]){
					tabuni[c,r]<-0
				}else{
					tabuni[c,r]<-1
				}
			}
		}
	gendistances<-gendistances*tabuni
	}
	
	mat[,1]<-mat[,1]+runif(nrow(mat), min = 0.001, max = 0.002)
	mat[,2]<-mat[,2]+runif(nrow(mat), min = 0.001, max = 0.002)
	image(map_alt)

	
#create the Delaunay triangles
try <- deldir(mat[,1],mat[,2])
maxdat<-max(max(try$delsgs[,5]),max(try$delsgs[,6]))
if(maxdat<nrow(mat))
{stop("Not all objects in triangulation consider collinearity")}
plot(try,wlines="triang",add=T)
#attribute to each edge its distance value
midpointx<-NULL
midpointy<-NULL
gendist<-NULL
lengthpath<-NULL
lengthsea<-NULL
seapointx<-NULL
seapointy<-NULL
altimin<-NULL
altimax<-NULL

for (i in 1:nrow(try$delsgs)){
#i<-20

gendist[i]<-gendistances[try$delsgs[i,5],try$delsgs[i,6]]

if(abs(try$delsgs[i,1]-try$delsgs[i,3])<0.043 & abs(try$delsgs[i,2]-try$delsgs[i,4])<0.043){
	try$delsgs[i,1]<-try$delsgs[i,1]+((0.043-abs(try$delsgs[i,1]-try$delsgs[i,3]))*sign(try$delsgs[i,1]-try$delsgs[i,3]))
	try$delsgs[i,2]<-try$delsgs[i,2]+((0.043-abs(try$delsgs[i,2]-try$delsgs[i,4]))*sign(try$delsgs[i,2]-try$delsgs[i,4]))

}

patch<-shortestPath(transcorr, c(try$delsgs[i,1],try$delsgs[i,2]),c(try$delsgs[i,3],try$delsgs[i,4]),output="SpatialLines")

lines(patch,col="red")
pa<-as(patch, "sf")
nume<-round(st_length(pa)*30.36,0)
punti<-st_sample(pa,nume,type = "regular")
coords<-(st_coordinates(punti)[,1:2])
coords<-rbind(c(10,10),coords)
value<-c(as.vector(extract(map1,coords)))
altitudes<-as.data.frame(as.vector(extract(map_alt,coords)))
val<-cbind(coords,value,altitudes)
altimin[i]<-min(val[,4],na.rm=T)
altimax[i]<-max(val[,4],na.rm=T)

lengthpath[i]<-nume


if(nume==1){
midpointx[i]<-(try$delsgs[i,1]+try$delsgs[i,3])/2
midpointy[i]<-(try$delsgs[i,2]+try$delsgs[i,4])/2
seapointx[i]<-(try$delsgs[i,1]+try$delsgs[i,3])/2
seapointy[i]<-(try$delsgs[i,2]+try$delsgs[i,4])/2
lengthsea[i]<-0
}

if(nume>1){
val<-val[-1,]
tabe<-with(rle(val[,3]), {
  ok <- values==1
  ends <- cumsum(lengths)[ok]
  starts <- ends - lengths[ok] + 1
  resi<-as.data.frame(cbind(starts, ends))
return(resi)
})

tabe<-cbind(tabe,(tabe[,2]-tabe[,1]+1))
tabe<-tabe[which(tabe[,3]>=minsea),]
line<-which.max(tabe[,3])[1]


midpointx[i]<-as.numeric(val[round((lengthpath[i]+0.0001)/2,0),1])
midpointy[i]<-as.numeric(val[round((lengthpath[i]+0.0001)/2,0),2])

lengthsea[i]<-sum(tabe[,3])
if(lengthsea[i]>0){
seapointx[i]<-as.numeric(val[round((tabe[line,1]+tabe[line,2])/2,0),1])
seapointy[i]<-as.numeric(val[round((tabe[line,1]+tabe[line,2])/2,0),2])


}


if(lengthsea[i]==0){
seapointx[i]<-midpointx[i]
seapointy[i]<-midpointy[i]
}
}
}
table<-cbind(try$delsgs,gendist,midpointx,midpointy,lengthpath,seapointx,seapointy,lengthsea,altimin,altimax)
return(table)
}



select.maximum.dist<-function(datasel,dist,size=0.045,max=3){
res<-NULL
lat<-(floor(datasel[,3]/size)*size)+(size/2)
long<-(floor(datasel[,4]/size)*size)+(size/2)
val<-paste(lat,long)

quant1<-aggregate(rep(1,nrow(datasel))~ val, FUN="sum")
	datasel2<-datasel[1,]
	for(loc in 1 : nrow(quant1)){
		#loc<-1
		which<-which(val==quant1[loc,1])
		if(quant1[loc,2]<(max+1)){
			datasel2<-rbind(datasel2,datasel[which,])
		}
		if(quant1[loc,2]>max){
			dataselt<-datasel[which,]
			sam<-sample(c(1:length(which)))
			dataselt<-dataselt[sam,]
			datasel2<-rbind(datasel2,dataselt[1:max,])
					}
	}


res$data<-datasel2[-1,]	
matchseqs<-match(res$data[,1],datasel[,1])
dista<-as.matrix(dist)[matchseqs,matchseqs]
res$dist<-as.dist(dista)
return(res)
}


genetic.regionalization.lambert<-function(fasta,spec,longi,lati, minimum=3, minarea=2, nidw=15, gsttab=NULL, GST="all", alpha=0.1, power=2, powerW=1, method="dist", methodclust="ward.D"){
res<-NULL
species<-unique(spec)
specloc<-aggregate(rep(1,length(longi))~longi+lati+spec,FUN="sum")
specarea<-aggregate(rep(1,nrow(specloc))~specloc[,3],FUN="sum")
speciessel<-specarea[which(specarea[,2]>=minarea),1]
sel<-which(spec%in%speciessel)

longi<-longi[sel]
lati<-lati[sel]
fasta<-fasta[sel]
species<-speciessel
spec<-spec[sel]

loca<-paste(longi,lati)
areeun<-unique(loca)
areetot<-unique(cbind(longi,lati))
aree<-cbind(c(1:nrow(areetot)),areetot)
loc<-NULL
for(ar in 1:nrow(aree)){
#ar<-1
wh<-which(loca==areeun[ar])
loc[wh]<-ar
}


matrixGst<-array(NA,dim=c(nrow(aree),nrow(aree),length(species)))
areatabledissGst<-array(NA,dim=c(nrow(aree),nrow(aree)))
for (sp in 1:length(species)){
#sp<-2
	sel<-which(spec==species[sp])
	fastared<-fasta[sel]
	datamatsp<-loc[sel]
	#cbind(longi,lati)[sel,]
	dismatsp<-as.matrix(dist.dna(fastared, model = "raw",pairwise.deletion = TRUE))
	areamatsp<-unique(datamatsp)
	
	if(GST=="all"){
		wei<-1
	}
	if(GST=="gd"){
	if(is.null(gsttab)){
	gstsp<-recluster.fst(as.dist(dismatsp),datamatsp,setzero=T,setnazero=T)$Gst
	}else{
	gstsp<-gsttab[which(names(gsttab)==species[sp])]
	}
		wei<-(1-gstsp)
	}
	if(GST=="st"){
	if(is.null(gsttab)){
	gstsp<-recluster.fst(as.dist(dismatsp),datamatsp,setzero=T,setnazero=T)$Gst
	}else{
	gstsp<-gsttab[which(names(gsttab)==species[sp])]
	}
		wei<-(gstsp)
	}
	if (length(areamatsp)>1){
	for (n in 1: (length(areamatsp)-1)){
		#n<-1
		for (m in (n+1):length(areamatsp)){
		#m<-7			
			dista<-dismatsp[c(which(datamatsp==areamatsp[n]), which(datamatsp==areamatsp[m])),c(which(datamatsp==areamatsp[n]), which(datamatsp==areamatsp[m]))]
			popu<-c(rep(1,length(which(datamatsp==areamatsp[n]))),rep(2,length(which(datamatsp==areamatsp[m]))))
			value<-((recluster.fst(dista,popu,setzero=T,setnazero=T)$Ht)^power)*(wei)
			matrixGst[areamatsp[n], areamatsp[m],sp]<-value
			matrixGst[areamatsp[m], areamatsp[n],sp]<-value
			}
		}
	}
}

count<-areatabledissGst
for (r in 1: nrow(aree)){
#r<-2
	for (c in 1: nrow(aree)){
#c<-2
		count[r,c]<-length(which(!is.na(matrixGst[r,c,])))
		if(count[r,c]>1){
			if(method=="dist"){
			areatabledissGst[r,c]<-(sum(matrixGst[r,c,],na.rm=T)^(1/power))/count[r,c]
			areatabledissGst[c,r]<-(sum(matrixGst[r,c,],na.rm=T)^(1/power))/count[r,c]
			}
			if(method=="mean"){
			areatabledissGst[r,c]<-mean(matrixGst[r,c,],na.rm=T)
			areatabledissGst[c,r]<-mean(matrixGst[r,c,],na.rm=T)
			}
			if(method=="median"){
			areatabledissGst[r,c]<-median(matrixGst[r,c,],na.rm=T)
			areatabledissGst[c,r]<-median(matrixGst[r,c,],na.rm=T)
			}
		}
	}
}

matnew<-areatabledissGst

rich<-aggregate(rep(1,length(loc))~loc+spec,FUN="sum")
richness<-aggregate(rep(1,nrow(rich))~rich[,1],FUN="sum")
sitiguida<-which(richness[,2]>=minimum)
guida<-rep(0,nrow(richness))
guida[sitiguida]<-1
#siti<-aree[sitiguida[,1],]
areeg<-cbind(aree,richness[,2],guida)
distance<-rep(1,nrow(areeg))
aree2<-NULL
k<-1
for (qua in 1:nrow(matnew)){
#qua<-1
tabella<-cbind(areeg,matnew[,qua])
tabella<-tabella[-qua,]
tabella<-tabella[which(tabella[,5]==1),]
if(length(which(!is.na(tabella[,6])))>=nidw){
tabella<-tabella[complete.cases(tabella),]
res1 <- idw(tabella[,6], tabella[,2:3],areeg[,2:3],p=2)
#
#colidw<-grey(distance[,2]/max(distance[,2]))
#plot(areeg[,2:3],col=colidw)
#hpts <-convhulln(tabella[,2:3])
#inhull<-which(!(inhulln(hpts ,areeg[,2:3])))
#head(res1)
#res1[inhull,]<-NA
distance<-cbind(distance,res1)

aree2[k]<-as.character(areeg[qua,1])


k<-k+1
}
}

aree2<-as.numeric(aree2)
distance<-distance[aree2,]
distance<-distance[,-1]
areefin<-sitiguida[which(sitiguida %in% aree2)]
#nrow(distance)
#rownames(distance)
sel<-which(rownames(distance)%in%areefin)
dist<-distance[sel,sel]
mat<-matrix(NA,nrow(dist),nrow(dist))
for(row in 1:nrow(mat)){
#row<-1
#col<-2
for(col in 1:ncol(mat)){
mat[row,col]<-(dist[row,col]+dist[col,row])/2
}}

dista<-as.dist(mat)
select<-as.numeric(rownames(dist))

coord<-areeg[select,]


distgeo<-dist(coord[,2:3])/1000


res$richness<-coord[,4]
clust<-hclustgeo2(dista,distgeo,alpha=alpha,w=coord[,4]^powerW,method=methodclust)
res$tree <- clust$solution
res$weightmatrix<-clust$weighteddiss
res$dist<-dista
res$distgeo<-distgeo
res$coord<-coord[,2:3]
res$originalcoord<-cbind(longi,lati)
res$originalsequences<-fasta
res$originalspec<-spec
res$minimum<-minimum
res$minarea<-minarea
res$gsttab<-gsttab
res$GST<-GST
res$alpha<-alpha
res$power<-power
res$powerW<-powerW
res$method<-method
return(res)
}


recluster.boot.regionalisation<-function (reg, boot=100,maxcl=8,nidw=13) {
minimum<-reg$minimum
minarea<-reg$minarea
gsttab<-reg$gsstab
alpha<-reg$alpha
GST=reg$GST
sequences<-reg$originalsequences
spec<-reg$originalspec
species<-unique(spec)
tree<-reg$tree
power<-reg$power
powerW<-reg$powerW
method<-reg$method
names<-as.data.frame(paste(reg$coord[,1],reg$coord[,2]))
tab<-names
restab<-NULL

for(c in 2:maxcl){
#c<-2
restab[[c-1]]<-matrix(NA, boot,c)
classification<-cutree(tree,c)
tab<-cbind(tab,classification)
}


for (i in 1 : boot){
xs<-species[sample(1:length(species),replace=T)]
select<-which(spec %in% xs)
sequences2<-sequences[select]
membership2<-spec [select]
lat2<-lat[select]
long2<-long[select]
treernd<-genetic.regionalization(sequences2,membership2,long2,lat2,minarea=minarea,minimum=minimum,nidw=nidw, GST=GST,gsttab=gsttab,alpha=alpha,powerW=powerW, power=power, method=method)
names<-paste(treernd$coord[,1],treernd$coord[,2])
match<-match(names,tab[,1])
nogood<-which(is.na(match))
if (length(nogood)>0){
match2<-match[-nogood]
names2<-names[-nogood]
}else{
match2<-match
names2<-names
}
tabrnd<-matrix(NA,nrow(tab),ncol(tab))
for(c in 2:maxcl){
#c<-2
classrnd<-cutree(treernd$tree,c)
if (length(nogood)>0){
tabrnd[match2,c]<-classrnd[-nogood]
}else{
tabrnd[match2,c]<-classrnd
}
}
tabsel<-tab[which(tab[,1] %in% names2),]
for(v in 2:maxcl){
#v<-7
tabrnd[match2,1]<-names2
somma<-aggregate(rep(1,nrow(tabsel))~tabsel[,v],FUN="sum")
vect<-NULL
for (gr in 1:nrow(somma)){
#gr<-1
quali<-tabsel[which(tabsel[,v]==somma[gr,1]),1]
tutti<-tabrnd[which(tabrnd[,1]%in%quali),v]
tot<-aggregate(rep(1,length(tutti))~tutti, FUN="sum")
vect[gr]<-max(tot[,2])/length(which(tabsel[,v]==somma[gr,1]))
}
restab[[v-1]][i,]<-vect
}
}
final<-NULL
for(giro in 1:(maxcl-1)){
final[[giro]]<-colMeans(restab[[giro]])
}
return(final)
}

hclustgeo2<-function (D0, D1 = NULL, alpha = 0, scale = TRUE, wt = NULL, method="ward.D") 
{
res<-NULL    
if (class(D0) != "dist") 
        stop("DO must be of class dist (use as.dist)", 
            call. = FALSE)
    if (!is.null(D1) && (class(D1) != "dist")) 
        stop("D1 must be of class dist (use as.dist)", 
            call. = FALSE)
    n <- as.integer(attr(D0, "Size"))
    if (is.null(n)) 
        stop("invalid dissimilarities", call. = FALSE)
    if (is.na(n) || n > 65536L) 
        stop("size cannot be NA nor exceed 65536", call. = FALSE)
    if (n < 2) 
        stop("must have n >= 2 objects to cluster", call. = FALSE)
    if (!is.null(D1) && length(D0) != length(D1)) 
        stop("the two dissimilarity structures must have the same size", 
            call. = FALSE)
    if ((max(alpha) > 1) || (max(alpha) < 0)) 
        stop("Values alpha must be in [0,1]", call. = FALSE)
    if ((scale == TRUE) && (!is.null(D1))) {
        D0 <- D0/max(D0)
        D1 <- D1/max(D1)
    }
    delta0 <- wardinit(D0, wt)
    if (!is.null(D1)) 
        delta1 <- wardinit(D1, wt)
    else delta1 <- 0
    delta <- (1 - alpha) * delta0 + alpha * delta1
    res$solution <- hclust(delta, method = method, members = wt)
	res$weighteddiss<-delta
    return(res)
}





