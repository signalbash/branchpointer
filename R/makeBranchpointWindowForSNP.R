#' Makes a branchpointer formatted GRanges object from refsnp ids
#'
#' Searches Biomart for refsnp ids, and pulls genomic location and sequence identity information
#' Reformats alleles so each query has only one alternative allele
#' @param refSNP Vector of refsnp ids
#' @param mart.snp biomaRt mart object specifying the BioMart database and dataset to be used
#' @param exons GRanges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @param maxDist maximum distance a SNP can be from an annotated 3' exon.
#' @param filter remove SNP queries prior to finding finding nearest exons?
#' @return formatted SNP query GRanges
#' @export
#' @import biomaRt
#' @import GenomicRanges
#' @importFrom stringr str_split
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v26.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' 
#' mart.snp <- biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="www.ensembl.org")
#' query <- makeBranchpointWindowForSNP("rs587776767", mart.snp, exons)
#' @author Beth Signal

makeBranchpointWindowForSNP <- function(refSNP, mart.snp,exons,maxDist=50, filter=TRUE){
  snpInfo <- biomaRt::getBM(attributes = c("refsnp_id",'refsnp_source', "chr_name",
                                   "chrom_start", "allele"),
                    filters = "snp_filter", values = refSNP, mart = mart.snp)

  #make sure each SNP has only 1 ref and 1 alternate allele
  multiAlleles <- which(nchar(snpInfo$allele) != 3)

  if (length(multiAlleles) > 0) {
    snpInfo.remade <- snpInfo[-multiAlleles,]
    snpInfo <- snpInfo[multiAlleles,]

    for (i in seq_along(snpInfo$refsnp_id)) {
      nts <- unlist(stringr::str_split(snpInfo$allele[i],"/"))
      ref <- nts[1]
      alt <- nts[-1]
      alleles <- paste0(ref,"/",alt)

      remade <- snpInfo[c(rep(i, length(alleles))),]
      remade$allele <- alleles
      snpInfo.remade <- rbind(snpInfo.remade, remade)
    }
    snpInfo <- snpInfo.remade
  }

  
  
  queryGRanges <- GRanges(seqnames=S4Vectors::Rle(paste0("chr",snpInfo$chr_name)),
                          ranges=IRanges::IRanges(start=snpInfo$chrom_start, width=1),
                          strand="*",
                          id=snpInfo$refsnp_id,
                          ref_allele=stringr::str_sub(snpInfo$allele,1,1),
                          alt_allele=stringr::str_sub(snpInfo$allele,3,3))

  #check for unstranded queries & replace with positive & negative
  unstranded <- which(as.logical(strand(queryGRanges) == "*"))
  if(length(unstranded)>0){
    queryGRanges.pos <- queryGRanges[unstranded]
    queryGRanges.pos$id <- 
      paste0(queryGRanges.pos$id, "_pos")
    strand(queryGRanges.pos) <- "+"
    
    queryGRanges.neg <- queryGRanges[unstranded]
    queryGRanges.neg$id <- 
      paste0(queryGRanges.neg$id, "_neg")
    strand(queryGRanges.neg) <- "-"
    
    queryGRanges <- do.call("c", list(queryGRanges[-unstranded],
                                      queryGRanges.pos,
                                      queryGRanges.neg))
    }

  #check for duplicated query ids
  if(any(duplicated(queryGRanges$id))){
    message(paste0(length(which(duplicated(queryGRanges$id)))," query ids are not unique"))
    message("Check output for new names or rename")
    queryGRanges$id <- make.names(queryGRanges$id, unique=TRUE)
  }
  
  #find 3'/5'exons
  if(length(queryGRanges) > 0){
    queryGRanges.loc <- getQueryLoc(queryGRanges, queryType="SNP", 
                                    maxDist = maxDist, filter = filter,
                                    exons = exons)
    return(queryGRanges.loc)
  }

}
