#' Make branchpoint window regions
#'
#' Generate branchpoint window regions corresponding to annotated exon(s) within a
#' queried gene, transcript or exon id
#' @param id identifier(s) for the query gene/transcript/exon id
#' @param idType type of id to match in the exon annotation file (\code{"gene_id"},
#' \code{"transcript_id"}, or \code{"exon_id"})
#' @param exons GRanges containing exon co-ordinates.
#' @param forceClosestExon Force branchpointer to find the closest exon and not the
#' exon annotated as 5' to the query
#' @return Granges with formatted query
#' @export
#' @import GenomicRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v26.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' windowquery <- makeBranchpointWindowForExons("ENSG00000139618.14", "gene_id", exons)
#' windowquery <- makeBranchpointWindowForExons("ENST00000357654.7", "transcript_id", exons)
#' windowquery <- makeBranchpointWindowForExons("ENSE00003518965.1", "exon_id", exons)
#' @author Beth Signal

makeBranchpointWindowForExons <- function(id, idType, exons, forceClosestExon = FALSE) {
    
  validTypes <- c("gene_id", "transcript_id","exon_id")

  #missing or invalid types
  noType <- missing(idType)

  if(!noType){
    noType <-  noType | !(idType %in% validTypes)
  }


  if (!noType) {
    x <- which(colnames(mcols(exons)) == idType)
    y <- which(!is.na(match(mcols(exons)[,x], id)))
    #y <- grep(id[1], mcols(exons)[,x])
  }else{
    y <- NA
  }

  #go through possible columns if no matches found
  if (all(is.na(y)) | noType) {
    idType <- validTypes[1]
    x <- which(colnames(mcols(exons)) == idType)
    y <- which(!is.na(match(mcols(exons)[,x], id)))
    #y <- grep(id[1], mcols(exons)[,x])

    if (length(y) == 0) {
      idType <- validTypes[2]
      x <- which(colnames(mcols(exons)) == idType)
      y <- which(!is.na(match(mcols(exons)[,x], id)))
      #y <- grep(id[1], mcols(exons)[,x])
    }

    if (length(y) == 0) {
      idType <- validTypes[3]
      x <- which(colnames(mcols(exons)) == idType)
      y <- which(!is.na(match(mcols(exons)[,x], id)))
      #y <- grep(id[1], mcols(exons)[,x])
    }

    if (length(y) == 0) {
      stop(paste0("cannot find ", id[1]," in the exon annotation"))
    }

  }

  #use a subset of the exon annotation for faster processing
  
  gene_id <- unique(mcols(exons)$gene_id[y])
  subsetInd <- which(mcols(exons)$gene_id %in% gene_id)
  exons.subset <- exons[subsetInd]
  
  keep <- which(mcols(exons[y])$exon_number > 1)
  window <- exons[y][keep]
  
  # manually generate windows

  # make window starts & ends
  window.ends <- start(window) - 18
  window.starts <- window.ends - 26
    
  negStrandIndex <- which(as.logical(strand(window) == "-"))
    
  window.starts[negStrandIndex] <- 
    end(window)[negStrandIndex] + 18
  window.ends[negStrandIndex] <- 
    window.starts[negStrandIndex] + 26
    
  window.IRanges <- IRanges(start=window.starts, end=window.ends)
  # replace ranges
  ranges(window) <- window.IRanges
    
  # skips getQueryLoc() -- which can be time limiting step
  if(forceClosestExon == FALSE){
    
    # to_3prime will be 18 for all windows
    mcols(window)$to_3prime <- 18
    
    # get 5' exon
    # Note that this is not neccesarily the 'closest' exon
    # closest can be enforced with forceClosestExon = TRUE
    transcriptExonNum.5prime <- paste0(mcols(window)$transcript_id,"_",
                                       as.numeric(mcols(window)$exon_number) -1)
    transcriptExonNum <- paste0(mcols(exons.subset)$transcript_id,"_",
                                as.numeric(mcols(exons.subset)$exon_number))
    
    m <- match(transcriptExonNum.5prime, transcriptExonNum)
    to_5prime <- end(window) - end(exons.subset[m])
    to_5prime[negStrandIndex] <- 
      (start(exons.subset[m]) - start(window))[negStrandIndex]
    mcols(window)$to_5prime <- to_5prime
    mcols(window)$same_gene <- TRUE
    mcols(window)$exon_3prime <- mcols(window)$exon_id
    mcols(window)$exon_5prime <- mcols(exons.subset[m])$exon_id
    
  }else{
      window <- getQueryLoc(window,queryType = "region",exons = exons.subset)
  }
    
  mcols(window)$id <- mcols(window)$exon_id
  
  windowNames <- paste0(seqnames(window), "_", 
                        start(window),"_", 
                        end(window), "_", 
                        strand(window))
  rm <- which(duplicated(windowNames))
  if(length(rm) > 0){
    window <- window[-rm]
  }

  if(length(window) > 0){
    return(window)
  }

}

