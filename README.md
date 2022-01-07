<img src="https://github.com/leondap/images/blob/main/io_database.png?raw=true" width="300" img align="right">

# IOdatabase
### Integrated and Open Butterly Database

The IOdatabase integrates several datasets of butterflies from Western Palaearctic interconnected with a harmonised taxonomy. IOdatabase includes some main datasets provided in recent paper like a revised version of the checklist and species features provided by Middleton Welling et al (2020) a series of more than 30000 COI sequences for 510 species (XXXX., 2021) some basic indexes of genetic diversity (number of haplotypes, GST, DST, haplotype diversity and nucleotide diversity)(XXXX et al 2022) and a revised version of the climatic variables introduced by Schweiger et al 2014 and updated by Platania et al (2020).


mydata<-get.IOdbseqs(subfamilies=c("Satyrinae"),species=c("Pieris_napi","Pieris_rapae","Charles_Darwin"),xlim=c(10,15))

To install IOdatabase you also need to install recluster. Use:
```
install.packages("remotes")
remotes::install_github("leondap/recluster")
remotes::install_github("leondap/iodatabase")
```

The function get.IOdbseqs extracts COI sequences and metadata for specimens. In particular the following query extracts all Satyrinae and Pieris napi and Pieris rapae while advertising that "Charles_Darwin" is not an available species in the checklist. Moreover, the xlim query limits the data to the specimens collected between 10 and 15 degrees of longitude.

```
libray(recluster)
library(iodatabase)
mydata<-get.IOdbseqs(subfamilies=c("Satyrinae"),species=c("Pieris_napi","Pieris_rapae","Charles_Darwin"),xlim=c(10,15))
```

The sequences can be obtained as a DNAbin object in mydata$fasta and the metadata as a table in mydata$metadata

# [UNDER CONSTRUCTION]
