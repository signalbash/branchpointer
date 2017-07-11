#' Convert GTF file to exon location file
#'
#' Converts a GTF annotation to exon locations
#' @param gtf file containing the gtf annotation.
#' @return exon annotation GRanges
#' @export
#' @import stringr
#' @import GenomicRanges
#' @importFrom utils write.table
#' @importFrom rtracklayer import
#' @examples
#' smallExons <- system.file("extdata","gencode.v26.annotation.small.gtf",
#' package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' @author Beth Signal

gtfToExons <- function(gtf){

  gtf.rtrack <- rtracklayer::import(gtf)

  #keep exon annotations only
  gtf.rtrack <- gtf.rtrack[gtf.rtrack$type=="exon"]

  #change "biotype" to "type" in GTFs from Ensembl
  if("gene_biotype" %in% colnames(mcols(gtf.rtrack))){
    colnames(mcols(gtf.rtrack)) <- gsub("biotype",
                                              "type",
                                        colnames(mcols(gtf.rtrack)))
  }

  if(!all(c("exon_id", "exon_number") %in% (colnames(mcols(gtf.rtrack))))){
    mcols(gtf.rtrack)$exon_number <- NA
    gtf.rtrack$exon_id <- gtf.rtrack$transcript_id

    #make exon names
    n <- 1
    while(any(is.na(gtf.rtrack$exon_number))){
      transcript_ids <- unique(gtf.rtrack$transcript_id)
      m <- match(transcript_ids,
                 gtf.rtrack$transcript_id[which(
                   is.na(gtf.rtrack$exon_number))])
      gtf.rtrack$exon_number[which(
        is.na(gtf.rtrack$exon_number))[m]] <- n
      n <- n+1
    }
  }

  exons <- gtf.rtrack
  mcols(exons) <- mcols(exons)[,c("gene_id","gene_type",
                                 "transcript_id","transcript_type",
                                 "exon_id","exon_number")]
  # close connections
  gc()

  return(exons)
}
