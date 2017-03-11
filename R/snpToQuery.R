#' Gets a branchpointer formatted query from refsnp ids
#'
#' Searches Biomart for refsnp ids, and pulls genomic location and sequence identity information
#' Reformats alleles so each query has only one alternative allele
#' @param refSNP Vector of refsnp ids
#' @param mart.snp biomaRt mart object specifying the BioMart database and dataset to be used
#' @return formatted SNP query data.frame
#' @export
#' @import biomaRt
#' @import GenomicRanges
#' @importFrom stringr str_split
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @examples
#' mart.snp <- biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host="www.ensembl.org")
#' query <- snpToQuery("rs17000647", mart.snp)
#' @author Beth Signal

snpToQuery <- function(refSNP, mart.snp){
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

  return(queryGRanges)

}
