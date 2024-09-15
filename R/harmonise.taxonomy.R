harmonise_taxonomy<-function(list, checkl=checklist[,5], conversiontab=conversion){
res<-NULL
absent<-which(!(list %in% checkl))
#list[absent]
namesabs<-cbind(list[absent],absent)
correggibili<-which(namesabs[,1]%in%conversiontab[,1])
non_correggibili<-which(!(namesabs[,1]%in%conversiontab[,1]))
res$failed<-namesabs[non_correggibili,1]
corretti<-match(namesabs[correggibili,1],conversiontab[,1])
nuovi<-conversiontab[corretti,2]
new<-namesabs
new[correggibili,1]<-nuovi
listnew<-list
listnew[as.numeric(new[,2])]<-new[,1]
listnew[which(!(listnew%in%checkl))]<-NA
res$listnew<-listnew
res$taxa_failed<-unique(res$failed)
return(res)
}

