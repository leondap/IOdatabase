

harmonise_taxonomy<-function(list, checkl=checklist[,5], conversiontab=conversion){
res<-NULL
absent<-which(!(list %in% checkl))
namesabs<-cbind(list[absent],absent)
correggibili<-which(namesabs[,1]%in%conversiontab[,1])
corretti<-match(namesabs[correggibili,1],conversiontab[,1])
nuovi<-conversiontab[corretti,2]
new<-namesabs
new[correggibili,1]<-nuovi
listnew<-list
listnew[as.numeric(new[,2])]<-new[,1]
absent2<-which(!(listnew %in% checkl))
res$listnew<-listnew
res$failed<-absent2
res$taxa_failed<-unique(listnew[absent2])
return(res)
}

tax<-harmonise_taxonomy(lis)