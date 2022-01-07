get.IOdbseqs<-function(family=NULL, subfamily=NULL, genus=NULL, species=NULL,
                   country=NULL, xlim=NULL, ylim=NULL){
res<-NULL
check<-which(species %in% checklist$Taxa.name)
failed<-species[-check]
print(paste("The following species were not found in the checklist:",failed))
speciesuse<-species[check]
select<-which(specimens$SpeciesName %in% speciesuse)
res$species<-species[check]
res$metadata<-specimens[select,]
res$fasta<-sequences[select]
return(res)
}
