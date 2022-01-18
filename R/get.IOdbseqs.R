get.IOdbseqs<-function(families=NULL, subfamilies=NULL, genera=NULL, species=NULL,
                   countries=NULL, xlim=NULL, ylim=NULL){
  spetot<-checklist$Scientific_Name
  speciesusefam<-NULL
  speciesusesub<-NULL
  speciesusegen<-NULL
  speciesusesp<-NULL
  uselong<-NULL
  uselat<-NULL
  usecou<-NULL
  res<-NULL
  if(!is.null(families)){
    check<-which(families %in% checklist$Family)
    failed<-families[-check]
      if(length(failed)>0){
      print(paste("The following families were not found in the checklist:",failed))
    }
    speciesusefam<-spetot[which(checklist$Family %in% families)]
  }
  if(!is.null(subfamilies)){
    check<-which(subfamilies %in% checklist$Subfamily)
    failed<-subfamilies[-check]
    if(length(failed)>0){
      print(paste("The following subfamilies were not found in the checklist:",failed))
    }
    speciesusesub<-spetot[which(checklist$Subfamily %in% subfamilies)]
  }
  if(!is.null(genera)){
    check<-which(genera %in% checklist$Genus)
    failed<-genera[-check]
    if(length(failed)>0){
      print(paste("The following genera were not found in the checklist:",failed))
    }
    speciesusegen<-spetot[which(checklist$Genus %in% genera)]
  }
  if(!is.null(species)){
    check<-which(species %in% spetot)
    failed<-species[-check]
    if(length(failed)>0){
      print(paste("The following species were not found in the checklist:",failed))
    }
    speciesusesp<-spetot[which(checklist$Scientific_Name %in% species)]
  }
speciesuse<-unique(c(speciesusefam,speciesusesub,speciesusegen,speciesusesp))
select<-which(specimens$Scientific_Name %in% speciesuse)

if(!is.null(xlim)){
  uselong<-which(specimens$Longitude>=xlim[1] & specimens$Longitude<=xlim[2])
}
if(!is.null(ylim)){
  uselat<-which(specimens$Latitude>=ylim[1] & specimens$Latitude>=xlim[2])
}
if(!is.null(countries)){
  usecoun<-which(specimens$Country %in% countries)
}
useall<-unique(c(uselong,uselat,usecou))
if (!is.null(useall)){
  select<-select[which(select %in% useall)]
}

res$metadata<-specimens[select,]
res$fasta<-sequences[select]
return(res)
}



