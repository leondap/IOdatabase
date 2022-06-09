<img src="https://github.com/leondap/images/blob/main/io_database.jpg?raw=true" width="350" img align="right">

# IOdatabase
### Integrated and Open Butterly Database

The IOdatabase is intended as a project integrating several datasets of butterflies from Western Palaearctic harmonised with a shared taxonomy. IOdatabase includes some main datasets provided in recent paper like a revised version of the species checklist from West Palerarctic (Middleton Welling et al. 2020) and more than 30000 COI sequences for 533 species. Metadata (coordinates, country, BOLD and Genbank codes are also available). Using these data we also calculated some basic indexes of genetic diversity (number of haplotypes observed and predicted based on rarefaction curves, GST, DST, haplotype diversity and nucleotide diversity). 


<img src="https://github.com/leondap/images/blob/main/cover.jpg?raw=true" width="180" img align="left">
The DNA-barcodes are also used to generate The Atlas of mitochondrial genetic diversity for Western Palearctic butterflies that can be freely downloaded from this link:<br>
https://www.dropbox.com/s/6hznae7zk3swuqc/Atlas_v1.0.pdf?dl=0
<br>
<br>
<br>
The IOdatabase project is currently contributed by four research units located in Florence University (Italy), Institut de Biologia Evolutiva (Spain), Institu Botànic de Barcelona (Spain) and University of Oulu (Finland).





To install IOdatabase you also need to install recluster. Use:
```
install.packages("remotes")
remotes::install_github("leondap/iodatabase")
```

The function get.IOdbseqs extracts COI sequences and metadata for specimens. In particular the following query extracts all Satyrinae and Pieris napi and Pieris rapae while advertising that "Charles_Darwin" is not an available species in the checklist. Moreover, the xlim query limits the data to the specimens collected between 10 and 15 degrees of longitude.

```
library(rworldmap)
library(rworldxtra)
map<-getMap(resolution = "high")

library(iodatabase)
# load the data
load("~/iodatabase/data/IODBdata.RData")
#Extract the data
mydata<-get.IOdbseqs(subfamilies=c("Satyrinae"),species=c("Pieris_napi","Pieris_rapae","Charles_Darwin"),xlim=c(10,15))
```

The sequences can be obtained as a DNAbin object in mydata$fasta and the metadata as a table in mydata$metadata

### Create a map for distribution of lineages. 
First select a species and obtain the metadata table (metadata) and the sequences (fasta) in the same order.
```
mydata<-get.IOdbseqs(species="Lasiommata_megera")
metadata<-mydata$metadata
fasta<-mydata$fasta
```
Define the areas to be separated in pies defined as islands, and as size of continental cells in decimal degrees of latitude and longitude (square). The areas to be considered as continental should be indicated in the areascoll vector.
```
sitescodes<-define.areas(coord=metadata[,c(10,9)],areas=metadata[,7], square=3, areascoll=c("Africa", "Eurasia", "Britain", "Ireland"))
```
Compute genetic distances among specimens
```
sp.gendists <- dist.dna(fasta, model = "raw", pairwise.deletion = TRUE)
```
Perform a PCoA and project the configuration into RGB space as done by recluster.col function

```
sp.pcoa <- stats::cmdscale(sp.gendists, k=2)
colours<-recluster.col2(sp.pcoa)
```
Adjust colours for the map
```
colfunc <- colorRampPalette(c("#C1DABB","#E1D9BB","#e0c9a7","#cbb491","#a78a5d")) # colors elevation palette 
pal <- colfunc(20) # make elevation palette
```

Plot the sequences as pie in the map, the mnore similar two sequences in colour the more genetically close they are
```
plot(cbind(range(metadata[,10]),range(metadata[,9])),type="n",xlab="",ylab="")
plot(elev,col = pal,add=T,legend=F)
plot(sea,add = T, col = "azure2", border=F)
plot(map, add=T)
arrows(-13.18, 27.67, -8.67, 27.67, length = 0, col="black") #Draw the correct Morocco boundary
recluster.plot.pie(long=metadata[,10],lat=metadata[,9], mat=colours, loc = sitescodes$val,minsize=0.3,add=T)

```
![](https://github.com/leondap/images/blob/main/Lasiommata_megera_map.png?raw=true)
