<img src="https://github.com/leondap/images/blob/main/io_database.jpg?raw=true" width="350" img align="right">

# IOdatabase
### Integrated and Open Butterly Database

The IOdatabase is intended as a project integrating several datasets of butterflies from Western Palaearctic harmonised with a shared taxonomy. IOdatabase now includes a revised version of the species checklist from West Palerarctic (Wiemers et al. 2016 and Middleton Welling et al. 2020) and more than 30000 COI sequences for 532 species. COI sequences belong to three main sources: 1) De novo sequencing mostly for Maghreb and Macaronesia 2) Published DNA-barcode libraries and public BOLD projects compiled at regional level 3) Studies providing COI of single species or genera. In the latter case, we checked whether haplotypes, instead of specimens, were included in repositories (Paz‐Vinas et al. 2021). If the number of specimens sharing a given haplotype in a particular location was reported, we replicated the haplotype sequences to obtain data at the specimen level, otherwise the sequences were excluded. We verified species identifications by building neighbor-joining trees for each genus. When morphology could not be verified, we removed sequences not clustering within the species to which they were attributed if the mismatch did not involve one of the 74 species showing DNA-barcode sharing.
Metadata (coordinates, country, BOLD and Genbank codes are also available). Using these data we also calculated some basic indexes of genetic diversity (number of haplotypes observed and predicted based on rarefaction curves, GST, DST, haplotype diversity and nucleotide diversity). 

![](https://github.com/leondap/images/blob/main/Figure-1-new.png?raw=true)
The study area with the representation of the “taxonomic area” (blue perimeter) from where a complete species checklist has been assessed and updated. COI sequences of species also occurring outside the taxonomic area are also included. (a) Number of sequences and (b) number of species for each 100x100 km squared cells. 





<img src="https://github.com/leondap/images/blob/main/cover2.jpg?raw=true" width="180" img align="left">
The DNA-barcodes are also used to generate The Atlas of mitochondrial genetic diversity for Western Palearctic butterflies that can be freely downloaded from this link:<br>
https://drive.google.com/file/d/1RrIEQOQq1ch70iggk2ARF9sIX7Vu8ABd/view?usp=sharing
<br>
The Atlas contains the map of COI variation for all European species for which sequences are available together with their indexes of intraspecific genetic differentiation
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

#Extract the data
mydata<-get.IOdbseqs(subfamilies=c("Satyrinae"),species=c("Pieris_napi","Pieris_rapae","Charles_Darwin"),xlim=c(10,15))
```

The sequences can be obtained as a DNAbin object in mydata$fasta and the metadata as a table in mydata$metadata

### Inspect the checklist.
```
head(checklist)
```

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

Plot the sequences as pie in the map, the mnore similar two sequences in colour the more genetically close they are
```
plot(cbind(range(metadata[,10]),range(metadata[,9])),type="n",xlab="",ylab="")
plot(map, add=T)
arrows(-13.18, 27.67, -8.67, 27.67, length = 0, col="black") #Draw the correct Morocco boundary
recluster.plot.pie(long=metadata[,10],lat=metadata[,9], mat=colours, loc = sitescodes$val,minsize=0.3,add=T)

```
![](https://github.com/leondap/images/blob/main/megera_amp.png?raw=true)


It is also possible to select any combination of colours (i.e a colour blind friendly palette) by creating RGB proiection with costum colours at te four corners

```
corners<-palette.col("#000000","#FFC107","#1E88E5","#E6E6E6",size=20)
coloCB <- recluster.col.palette(sp.pcoa,palette=corners,st=T)       	
recluster.plot.col(coloCB, text=F, cex=1.5)
```
![](https://github.com/leondap/images/blob/main/pcoacb.png?raw=true)

And plot the map

```
plot(cbind(range(metadata[,10]),range(metadata[,9])),type="n",xlab="",ylab="")
plot(map, add=T)
arrows(-13.18, 27.67, -8.67, 27.67, length = 0, col="black") #Draw the correct Morocco boundary
recluster.plot.pie(long=metadata[,10],lat=metadata[,9], mat=coloCB, loc = sitescodes$val,minsize=0.3,add=T)
```
![](https://github.com/leondap/images/blob/main/megera%20colour%20blind.png?raw=true)

The IOdatabse also has coordinates in the Lambert Azimuthal equal area projection to obtain maps where distances among points are preserved and to group specimens into areas of the same size

```
crs.laea <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
```
For this it is only necessary to transform the map

```
library(sp)
newmap_lamb <- spTransform(map, crs.laea)
```
And repeat the analysis using the coordinates stored in columns #16 and #17
in this case peis are grouped to areas 200,000 metres large (in define.areas) 
pies for singletons appear 35,000 metres large on the map (in recluster.plot.pie)

```
sitescodes<-define.areas(coord=metadata[,c(16,17)],areas=metadata[,7], square=200000, areascoll=c("Africa", "Eurasia", "Britain", "Ireland"))
sp.gendists <- dist.dna(fasta, model = "raw", pairwise.deletion = TRUE)
sp.pcoa <- stats::cmdscale(sp.gendists, k=2)
corners<-palette.col("#000000","#FFC107","#1E88E5","#E6E6E6",size=20)
coloCB <- recluster.col.palette(sp.pcoa,palette=corners,st=T)       	
recluster.plot.col(coloCB, text=F, cex=1.5)
plot(cbind(range(metadata[,16]),range(metadata[,17])),type="n",xlab="",ylab="")
plot(newmap_lamb, add=T)
arrows(2019182,867633.3,2456683,748569.0,length=0)
recluster.plot.pie(long=metadata[,16],lat=metadata[,17], mat=coloCB, loc = sitescodes$val,minsize=35000,add=T)

```

![](https://github.com/leondap/images/blob/main/rhamni_map.png?raw=true)

The graph axes are now in metres

To add a different colour to the sea revert the newmap_lamb obejct and plot it

```
outline <- maps::map(newmap_lamb, plot=FALSE) # returns a list of x/y coords
xrange <- range(outline$x, na.rm=TRUE) # get bounding box
yrange <- range(outline$y, na.rm=TRUE)
xbox <- xrange + c(-2, 2)
ybox <- yrange + c(-2, 2)
plot(cbind(range(metadata[,16]),range(metadata[,17])),type="n",xlab="",ylab="")
plot(newmap_lamb, add=T)
polypath(c(outline$x, NA, c(xbox, rev(xbox))), c(outline$y, NA, rep(ybox, each=2)),col="azure2", rule="evenodd")
arrows(2019182,867633.3,2456683,748569.0,length=0)
recluster.plot.pie(long=metadata[,16],lat=metadata[,17], mat=coloCB, loc = sitescodes$val,minsize=35000,add=T)
```
![](https://github.com/leondap/images/blob/main/gonepteryx_sea.png?raw=true)



References

Middleton-Welling J, Dapporto L, García-Barros E, et al (2020) A new comprehensive trait database of European and Maghreb butterflies, Papilionoidea. Sci Data 7:. https://doi.org/10.1038/s41597-020-00697-7

Paz‐Vinas I, Jensen EL, Bertola LD, et al (2021) Macrogenetic studies must not ignore limitations of genetic markers and scale. Ecol Lett

Wiemers M, Balletto E, Dincă V, et al (2018) An updated checklist of the European butterflies (Lepidoptera, papilionoidea). Zookeys 2018:9–45. https://doi.org/10.3897/zookeys.811.28712



