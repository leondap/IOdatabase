
define.areas<-function(coord,areas, square,areascoll=c("Africa", "Eurasia", "Britain", "Ireland")){
  res<-NULL
  val<-NULL
  islandlat<-aggregate(coord[,2]~ areas, FUN="mean")
  islandlong<-aggregate(coord[,1]~ areas, FUN="mean")
  coord<-cbind(coord,rep(0,nrow(coord)))
  coord<-cbind(coord,rep(0,nrow(coord)))
  for (i in 1:nrow(islandlong)){
    coord[which(areas==islandlong[i,1]),3]<-islandlong[i,2]
  }
  for (i in 1:nrow(islandlat)){
    coord[which(areas==islandlat[i,1]),4]<-islandlat[i,2]
  }
lat<-(floor(coord[,2]/square)*square)+(square/2)
long<-(floor(coord[,1]/square)*square)+(square/2)

main<-coord[which(areas%in%areascoll),]
mainar<-areas[which(areas%in%areascoll)]

if(nrow(main)>0){
  mainlat<-aggregate(main[,1]~ long[which(areas%in%areascoll)]+lat[which(areas%in%areascoll)]+mainar, FUN="mean")
  mainlat<-cbind(mainlat,aggregate(main[,2]~long[which(areas%in%areascoll)]+lat[which(areas%in%areascoll)]+mainar, FUN="mean")[,4])
  dat<-paste(mainlat[,1],mainlat[,2],mainlat[,3])
  sitid<-paste(long,lat,areas)

  for(k in 1 : length(dat)){
  qual<-which(sitid==dat[k])
  coord[qual,3]<-mainlat[k,1]
  coord[qual,4]<-mainlat[k,2]
  }
}

sites<-as.data.frame(unique(cbind(areas,coord[,3:4])))
sites<-cbind(sites,c(1:nrow(sites)))
siti<-paste(sites[,1],sites[,2],sites[,3])
sitid<-paste(areas,coord[,3],coord[,4])

for(ar in 1:nrow(sites)){
  wh<-which(sitid==siti[ar])
  val[wh]<-ar
}
res$val<-val
res$data<-cbind(coord,areas)
res$sites<-sites
return(res)
}
