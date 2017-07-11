#' Convert exon annotation GRanges to intron locations
#'
#' Converts exon annotation to intron locations overlapping the branchpoint region
#' for exculsion of non-branchpoint region SNPs
#' Returns a character vector of chromosome locations
#' @param exons GRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist Maximum distance from the 3' exon to create the branchpoint region.
#' @return GRanges containing intron window co-ordinates
#' @import GenomicRanges
#' @keywords internal
#' @author Beth Signal
exonsToIntrons <- function(exons, maxDist = 50){
  
  introns <- exons
  
  posInd <- which(as.logical(strand(exons) == "+"))
  start(introns[posInd]) <- start(introns[posInd]) - (maxDist + 1)
  end(introns[posInd]) <- start(introns[posInd]) + maxDist
  
  negInd <- which(as.logical(strand(exons) == "-"))
  end(introns[negInd]) <- end(introns[negInd]) + (maxDist + 1)
  start(introns[negInd]) <- end(introns[negInd]) - maxDist
  
  return(introns)
}

#' Get the closest 3' and 5' exons
#'
#' Finds the closest annotated exons from a genomic co-ordinate.
#' Returns the distance to the 3' exon, distance to the 5' exon,
#' ids of the 3' and 5' exon,
#' and if the exons are from the same parent gene
#' @param query GenomicRangesquery
#' @param exons GenomicRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param queryType type of query. "SNP" or "region"
#' @return GenomicRanges with distance to the closest 3' and 5' exons,
#' whether these exons are part of the same gene (i.e is the location intronic),
#' and the identifiers for the 3' and 5' exons.
#' @import GenomicRanges
#' @keywords internal
#' @author Beth Signal

getExonDists <- function(query, exons, queryType){
  
  # if using regions, region start will be used for negative strand queries
  # & region end for positive strand queries
  
  negInd <- which(as.logical(strand(query) == '-'))
  posInd <- which(as.logical(strand(query) == '+'))
  
  queryPoint <- query
  
  width(ranges(queryPoint))[negInd] <- 1
  end <- end(ranges(queryPoint))[posInd]
  start(ranges(queryPoint[posInd])) <- end
  width(ranges(queryPoint[posInd])) <- 1
  
  f <- GenomicRanges::follow(queryPoint,exons)
  gene5 <- exons$gene_id[f]
  exon5 <- exons$exon_id[f]
  
  keep <- which(!is.na(f))
  f <- f[keep]
  to5prime <- rep(NA, length(queryPoint))
  to5prime[keep] <- GenomicRanges::distance(queryPoint[keep],exons[f]) + 1
  
  
  p <- GenomicRanges::precede(queryPoint, exons)
  gene3 <- exons$gene_id[p]
  exon3 <- exons$exon_id[p]
  
  keep <- which(!is.na(p))
  p <- p[keep]
  to3prime <- rep(NA, length(queryPoint))
  to3prime[keep] <- GenomicRanges::distance(queryPoint[keep],exons[p]) + 1
  
  sameGene <- gene3 == gene5
  
  query$to_3prime <- to3prime
  query$to_5prime <- to5prime
  query$same_gene <- sameGene
  query$exon_3prime <- exon3
  query$exon_5prime <- exon5
  
  return(query)
}

#' Find the closest 3' and 5' exons to a branchpointer query
#'
#' Finds the closest annotated exons from genomic co-ordinates in a branchpointer query GRanges
#' @param query branchpointer query GenomicRanges
#' must have chromosome at position 2, genomic co-ordinate at position 3,
#' and strand at position 4.
#' @param queryType type of query file (\code{"SNP"} or \code{"region"})
#' @param exons data.frame containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist maximum distance a SNP can be from an annotated 3' exon.
#' @param filter remove SNP queries prior to finding finding nearest exons.
#' @return GenomicRanges with the query and its location relative to the 3' and 5' exons
#' @keywords internal
#' @import GenomicRanges
#' @author Beth Signal

getQueryLoc <- function(query, queryType, maxDist=50, filter=TRUE, exons){
  
  if(missing(queryType) | !(queryType %in% c("SNP", "region"))){
    
    stop("please specify queryType as \"region\" or \"SNP\"")
    
  }
  
  if(missing(exons)){
    
    stop("please specify exon annotation object")
    
  }
  
  if(queryType=="SNP"){
    
    if(filter){
      message("filtering for SNPs in branchpoint windows")
      introns <- exonsToIntrons(exons, maxDist)
      o <- GenomicRanges::findOverlaps(query, introns)
      keep <- unique(o@from)
      query <- query[keep]
      
    }
  }
  
  nearestExons <- getExonDists(query, exons, queryType)
  
  #remove any points not in same gene
  rm <- which(nearestExons$same_gene == FALSE |
                is.na(nearestExons$exon_3prime) | 
                is.na(nearestExons$exon_5prime))
  
  if(length(rm) > 0){
    nearestExons <- nearestExons[-rm]
  }
  
  
  if(queryType=="region"){
    #if introns regions are closer to the 3'SS than the branchpoint region
    # move the region to cover the window of the nearest exon
    near3 <- which(nearestExons$to_3prime < 18)
    if(length(near3) > 0){
      addDistance <- 18 - nearestExons$to_3prime[near3]
      negStrand <- which(as.logical(strand(nearestExons)[near3] == "-"))
      posStrand <- which(as.logical(strand(nearestExons)[near3] == "+"))
      
      if(length(negStrand) > 0){
        start(ranges(nearestExons))[near3][negStrand] <- 
          start(ranges(nearestExons))[near3][negStrand] + addDistance[negStrand]
      }
      if(length(posStrand) > 0){
        end(ranges(nearestExons))[near3][posStrand] <- 
          end(ranges(nearestExons))[near3][posStrand] - addDistance[posStrand]
      }
      
    }
    
    #if introns regions are further from the 3'SS than the branchpoint region
    notNear3 <- which(nearestExons$to_3prime > 44)
    if(length(notNear3) > 0){
      nearestExons <- nearestExons[-notNear3]
    }
    
    nearestExons <- getExonDists(nearestExons,exons,queryType)
    
    #adjust query location to only cover the 27nt window
    
    #negStrand -- move end
    negStrand <- which(as.logical(strand(nearestExons) == "-"))
    if(length(negStrand) > 0){
      end(ranges(nearestExons[negStrand])) <- 
        start(ranges(nearestExons[negStrand])) - 
        nearestExons$to_3prime[negStrand] + 44
    }
    #posStrand -- move start
    posStrand <- which(as.logical(strand(nearestExons) == "+"))
    if(length(posStrand) > 0){
      start(ranges(nearestExons[posStrand])) <- 
        end(ranges(nearestExons[posStrand])) +
        nearestExons$to_3prime[posStrand] - 44
    }
    
  }else if(queryType=="SNP"){
    
    #if a 5' or 3' exon can't be found remove from analysis
    rm <- which(nearestExons$to_3prime==-1 | nearestExons$to_5prime==-1)
    if(length(rm) > 0){
      nearestExons <- nearestExons[-rm,]
    }
    
    rm <- which(nearestExons$to_3prime > maxDist)
    if(length(rm) > 0){
      nearestExons <- nearestExons[-rm,]
    }
    
    #check both exons come from the same parent gene 
    #(i.e query is in intronic region)
    rm <- which(!(nearestExons$same_gene))
    if(length(rm) > 0){
      nearestExons <- nearestExons[-rm,]
    }
  }
  return(nearestExons)
}

#' Read a query file
#'
#' Reads and formats a manually generated query file, and finds realtive locations of the closest annotated exons
#' Converts unstranded SNPs to two entries for each strand.
#' Checks for duplicate names and replaces if found.
#' @param queryFile tab delimited file containing query information.
#' For intronic regions should be in the format: region id, chromosome name, region start, region end, strand.
#' For SNP variants should be in the format: SNP id, chromosome name, SNP position, strand, reference allele (A/T/C/G), alternative allele (A/T/C/G)
#' @param queryType type of query file (\code{"SNP"} or \code{"region"})
#' @param exons GRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist maximum distance a SNP can be from an annotated 3' exon.
#' @param filter remove SNP queries prior to finding finding nearest exons.
#' @return Formatted query GRanges
#' @export
#' @import GenomicRanges
#' @importFrom data.table fread
#' @importClassesFrom data.table data.table
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v26.annotation.small.gtf", package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' 
#' querySNPFile <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' querySNP <- readQueryFile(querySNPFile, queryType = "SNP", exons)
#'
#' queryIntronFile <- system.file("extdata","intron_example.txt", package = "branchpointer")
#' queryIntron <- readQueryFile(queryIntronFile,queryType = "region", exons)
#' @author Beth Signal

readQueryFile <- function(queryFile, queryType, exons, maxDist=50, filter=TRUE){

  if(missing(queryType)){
    queryTest <- data.table::fread(
      queryFile,header = TRUE,
      nrows=1)

    if(ncol(queryTest) == 5){
      queryType <- "region"
      message("first line of file has 5 columns")
      message("using query_type = \"region\"")
    }else if(ncol(queryTest) == 6){
      queryType <- "SNP"
      message("first line of file has 6 columns")
      message("using query_type = \"SNP\"")
    }else{
      stop("please specify query_type and provide correctly formatted file")
    }
  }

  if (!(queryType %in% c("SNP", "region"))) {
    message("please specify query_type as \"region\" or \"SNP\"")
  }else{
    if (queryType == "SNP") {
      query <- data.table::fread(
        queryFile,header = TRUE,
        colClasses = c(
          "character", "character", "numeric",
          "character", "character", "character"
        ), data.table=FALSE)
      colnames(query)[1:6] <-
        c("id", "chromosome","chrom_start","strand","ref_allele", "alt_allele")

      #check single nts
      query$ref_allele <- toupper(query$ref_allele)
      query$alt_allele <- toupper(query$alt_allele)

      #check for unstranded queries & replace with positive & negative
      unstranded <- which(query$strand != "+" | query$strand != "-")
      if (length(unstranded) > 0) {
        query.pos <- query[unstranded,]
        query.pos$id <- paste0(query.pos$id, "_pos")
        query.pos$strand <- "+"
        query.neg <- query[unstranded,]
        query.neg$id <- paste0(query.neg$id, "_neg")
        query.neg$strand <- "-"
        query <- rbind(query[-unstranded,], query.pos,query.neg)
        }
      
      queryGRanges <- GRanges(seqnames=S4Vectors::Rle(query$chromosome),
                              ranges=IRanges::IRanges(start=query$chrom_start, width=1),
                              strand=query$strand,
                              id=query$id,
                              ref_allele=query$ref_allele,
                              alt_allele=query$alt_allele)
    }else if (queryType == "region") {
      query <- data.table::fread(
        queryFile, header = TRUE,colClasses = c(
          "character","character","numeric","numeric","character"
        ), data.table=FALSE)
      colnames(query) <-
        c("id","chromosome","chrom_start","chrom_end","strand")
      
      queryGRanges <- GenomicRanges::GRanges(
          seqnames = S4Vectors::Rle(query$chromosome),
          ranges = IRanges::IRanges(start = query$chrom_start,
                                    end = query$chrom_end),
          strand = query$strand,
          id = query$id)
    }

    #check for duplicated query ids
    if (any(duplicated(queryGRanges$id))) {
      message(paste0(length(which(
        duplicated(queryGRanges$id)
      ))," query ids are not unique"))
      message("Check output for new names or rename")
      queryGRanges$id = make.names(queryGRanges$id, unique = TRUE)
    }
    
    #find 3'/5'exons
    if(length(queryGRanges) > 0){
      queryGRanges.loc <- getQueryLoc(queryGRanges, queryType, maxDist = maxDist, filter = filter,
                                    exons = exons)
      return(queryGRanges.loc)
    }
    
  }
}
