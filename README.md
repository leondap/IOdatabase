<img src="https://github.com/leondap/images/blob/main/io_database.jpg?raw=true" width="400" img align="right">

# IOdatabase
### Integrated and Open Butterly Database

The IOdatabase integrates several datasets of butterflies from Western Palaearctic interconnected with a harmonised taxonomy. IOdatabase includes some main datasets provided in recent paper like a revised version of the checklist and species features provided by Middleton Welling et al (2020) a series of more than 30000 COI sequences for 510 species (XXXX., 2021) some basic indexes of genetic diversity (number of haplotypes, GST, DST, haplotype diversity and nucleotide diversity)(XXXX et al 2022) and a revised version of the climatic variables introduced by Schweiger et al 2014 and updated by Platania et al (2020).

The IOdatabase project is currently contributed by four research units located in Florence University (Italy), Institut de Biologia Evolutiva (Spain), Institu Botànic de Barcelona (Spain) and University of Oulu (Finland).


To install IOdatabase you also need to install recluster. Use:
```
install.packages("remotes")
remotes::install_github("leondap/recluster")
remotes::install_github("leondap/iodatabase")
```

The function get.IOdbseqs extracts COI sequences and metadata for specimens. In particular the following query extracts all Satyrinae and Pieris napi and Pieris rapae while advertising that "Charles_Darwin" is not an available species in the checklist. Moreover, the xlim query limits the data to the specimens collected between 10 and 15 degrees of longitude.

```
library(rworldmap)
library(rworldxtra)
map<-getMap(resolution = "high")

library(iodatabase)
mydata<-get.IOdbseqs(subfamilies=c("Satyrinae"),species=c("Pieris_napi","Pieris_rapae","Charles_Darwin"),xlim=c(10,15))
```

The sequences can be obtained as a DNAbin object in mydata$fasta and the metadata as a table in mydata$metadata

Create a map for distribution of lineages. First select a species
```
mydata<-get.IOdbseqs(species="Lasiommata_megera")
metadata<-mydata$metadata
fasta<-mydata$fasta
```
Define the areas to be separated in pies defined as islands, and as size of continental cells in decimal degrees of latitude and longitude (square). The areas to be considered as continental should be indicated in the areascoll vector.
```
sitescodes<-define.areas(coord=metadata[,c(10,9)],areas=metadata[,7], square=2, areascoll=c("Africa", "Eurasia", "Britain", "Ireland"))
```
Compute genetic distances among specimens
```
sp.gendists <- dist.dna(fasta, model = "raw", pairwise.deletion = TRUE)
```
Perform a PCoA and project the configuration into RGB space as done by recluster.col function

```
sp.pcoa <- stats::cmdscale(sp.gendists, k=2)
colours<-recluster.col(sp.pcoa)
```
Plot the sequences as pie in the map, the mnore similar two sequences in colour the more genetically close they are
```
recluster.plot.pie(long=metadata[,10],lat=metadata[,9], mat=colours, loc = sitescodes$val,minsize=0.4)
plot(map,add=T)
```



# [UNDER CONSTRUCTION]
