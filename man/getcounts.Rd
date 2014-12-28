\name{getcounts}
\alias{getcounts}
\title{Make a count matrix}
\usage{
getcounts(regions, binsize, marks, mapqs = 0, pairedends = FALSE,
  shifts = 75)
}
\arguments{
  \item{regions}{GRanges object containing the genomic
  regions of interest.  Each region will be divided into
  non-overlapping bins of equal size.  The reads falling in
  each bin will be the elements of the final count matrix.}

  \item{binsize}{The size of each bin in basepairs. Each
  region must have a width multiple of \code{binsize},
  otherwise an error will be thrown.  Use the function
  \code{refineRegions} to make sure that your GRanges
  object satisfies this constraint.  value for the models.
  See the item \code{models} in the return values to see
  how the}

  \item{marks}{Named list where each item is a path to an
  indexed bam file.  The chromosome names in the bam file
  must match with the sequence names in the GRanges
  object.}

  \item{mapqs}{Threshold on the 'mapq' field in the bam
  format. Reads with a mapping quality strictly lower than
  'mapq' will be discarded.}

  \item{pairedends}{Set this option to TRUE or FALSE to
  activate or deactivate the paired-end mode. Only read
  pairs where both ends are mapped will be considered and
  assigned to the bin where the midpoint of the read pair
  is located.  If this flag is set, the 'shift' option will
  be ignored.}

  \item{shifts}{Shift the reads in the 5' direction by a
  fixed number of base pairs. The read will be assigned to
  the bin where the shifted 5' end of the read is located.
  This option is ignored when the paired-end mode is
  active. (default: 75)}
}
\value{
A count matrix where the rows are the different marks
labelled with the names of the list \code{marks}, and the
columns are the different bins in the regions.
}
\description{
Make a count matrix from a GRanges object and a list of bam
files
}
\details{
The \code{mapqs, pairedends, shifts} parameters can be
provided once and they will be used for all bam files, or
as many times as there are bam files, in which case the
i-th path on the list will be matched to the i-th element
in the vector.
}

