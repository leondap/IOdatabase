\name{get.IOdbseqs}
\alias{get.IOdbseqs}
\title{
Obtain data from the iodatabase dataset
}
\description{
This function is designed to filter the data (COI sequences and metadata) comprised in the iodatabase dataset.
}
\usage{
get.IOdbseqs(families=NULL, subfamilies=NULL, genera=NULL, species=NULL, countries=NULL, xlim=NULL, ylim=NULL)
}
\arguments{
  \item{families}{
a string or a vector indicating the butterfly families to select.
}
  \item{subfamilies}{
a string or a vector indicating the butterfly subfamilies to select.
}
  \item{genera}{
a string or a vector indicating the butterfly genera to select.
}
  \item{species}{
a string or a vector indicating the butterfly species to select.
}
  \item{countries}{
a string or a vector to select specimens belonging to one or to a series of countries.
}
  \item{xlim}{
a vector to select specimens based on a range of longitude (decimal degrees).
}
  \item{ylim}{
a vector to select specimens based on a range of latitude (decimal degrees).
}
}
\value{
\item{fasta}{
The COI sequences in fasta format.
	}
\item{metadata}{
The metadata of specimens in the same order as in the fasta file.
	}
}

\author{
Leonardo Dapporto
}
\keyword{selection}
