#' Make branchpoint window regions
#'
#' Genrate branchpoint window regions corresponding to annotated exon(s) within a
#' queried gene, transcript or exon id
#' @param id identifier for the query gene/transcript/exon id
#' @param idType type of id to match in the exon annotation file (\code{"gene_id"},
#' \code{"transcript_id"}, or \code{"exon_id"})
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return Granges with formatted query
#' @export
#' @import GenomicRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' windowquery <- makeRegions("ENSG00000139618", "gene_id", exons)
#' windowquery <- makeRegions("ENST00000357654", "transcript_id", exons)
#' windowquery <- makeRegions("ENSE00003518965", "exon_id", exons)
#' @author Beth Signal

makeRegions <- function(id, idType, exons) {

  validTypes <- c("gene_id", "transcript_id","exon_id")

  #missing or invalid types
  noType <- missing(idType)

  if(!noType){
    noType <-  noType | !(idType %in% validTypes)
  }


  if (!noType) {
    x <- which(colnames(mcols(exons)) == idType)
    y <- grep(id, mcols(exons)[,x])
  }else{
    y <- NA
  }

  #go through possible columns if no matches found
  if (is.na(y) | noType) {
    idType <- validTypes[1]
    x <- which(colnames(mcols(exons)) == idType)
    y <- grep(id, mcols(exons)[,x])

    if (length(y) == 0) {
      idType <- validTypes[2]
      x <- which(colnames(mcols(exons)) == idType)
      y <- grep(id, mcols(exons)[,x])
    }

    if (length(y) == 0) {
      idType <- validTypes[3]
      x <- which(colnames(mcols(exons)) == idType)
      y <- y <- grep(id, mcols(exons)[,x])
    }

    if (length(y) == 0) {
      stop(paste0("cannot find ", id," in the exon annotation"))
    }

  }

  #use a subset of the exon annotation for faster processing
  if (idType != "gene_id") {
    gene_id <- exons$gene_id[y[1]]
    y2 <- which(!is.na(match(exons$gene_id,gene_id)))
  }else{
    y2 <- y
  }

  exons.subset <- exons[y]

  #by definition first exons shouldn' have branchpoints
  keep <- which(exons.subset$exon_number > 1)

  if (as.logical(strand(exons.subset)[1] == "+")) {
    windowStarts <- (start(ranges(exons.subset)) - 50)[keep]
  }else{
    windowStarts <- (end(ranges(exons.subset)) + 10)[keep]
  }

  window <- exons.subset[keep]

  # giving errors trying to set as start(ranges(window))
  window@ranges@start <- as.integer(windowStarts)
  width(ranges(window)) <- 41

  return(getQueryLoc(window,queryType = "region",exons = exons[y2]))

}
