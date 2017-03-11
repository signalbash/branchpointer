#' Convert SNP branchpoint predictions across the branchpoint window to an intronic summary
#'
#' Takes predictions of branchpoint probabilities from a reference
#' and alternative SNP and summarises the effect(s) of the SNP.
#' @param predictions site-wide branchpoint proability predictions
#' produced from predictBranchpoints()
#' @param query query GRanges containing all SNP ids to be summarised
#' @param probabilityCutoff Value to be used as the cutoff for
#' discriminating branchpoint sites from non-branchpoint sites
#' (default = \code{0.5})
#' @param probabilityChange Minimum probability score change
#' required to call a branchpoint site as deleted or created by
#' a SNP (default = \code{0.2})
#' @return GRanges with summarised branchpoint changes
#' occuring within the intron due to a SNP.
#' @export
#' @import GenomicRanges
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",
#' package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' querySNP <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(querySNP,queryType = "SNP")
#' query <- getQueryLoc(query,queryType = "SNP",exons = exons, filter = FALSE)
#' predictions <- predictBranchpoints(query,queryType = "SNP",BSgenome = genome)
#' snpStats <- predictionsToStats(query,predictions)
#' @author Beth Signal
#'
predictionsToStats <- function(query,
                               predictions,
                               probabilityCutoff = 0.5,
                               probabilityChange = 0.2){

  snpIDs <- query$id
  
  mcols(query)$BP_num_REF = NA
  mcols(query)$BP_num_ALT = NA
  mcols(query)$deleted_n = NA
  mcols(query)$created_n = NA
  mcols(query)$dist_to_BP_REF = NA
  mcols(query)$dist_to_BP_ALT = NA
  mcols(query)$max_prob_REF = NA
  mcols(query)$max_prob_ALT = NA
  mcols(query)$max_U2_REF = NA
  mcols(query)$max_U2_ALT = NA

    for(z in seq(along = snpIDs)){
      #match snp id to predictions
      predictions.snp <- predictions[predictions$id == snpIDs[z]]

      #match snp id to query attributes
      queryIndex <- which(query$id == snpIDs[z])

      #find location of the SNP relative to the annotated 3'exon
       #elementMeta$to_3prime

      # find location of the SNP relative to the predicted BPs
      # in reference and alternative sequences
      
      branchpointIndex.ref <- which(predictions.snp$status == "REF" & 
                                      predictions.snp$branchpoint_prob >= probabilityCutoff)
      
      diffs.ref <- predictions.snp$to_3prime[branchpointIndex.ref] - 
        predictions.snp$to_3prime_point[branchpointIndex.ref] 
      query$BP_num_REF[queryIndex] <- length(diffs.ref)
      
      if(length(diffs.ref) > 0){
        query$dist_to_BP_REF[queryIndex] <- diffs.ref[which.min(abs(diffs.ref))]
      }

      branchpointIndex.alt <- which(predictions.snp$status == "ALT" & 
                                      predictions.snp$branchpoint_prob >= probabilityCutoff)
      
      diffs.alt <- predictions.snp$to_3prime[branchpointIndex.alt] - 
        predictions.snp$to_3prime_point[branchpointIndex.alt] 
      query$BP_num_ALT[queryIndex] <- length(diffs.alt)
      
      if(length(diffs.alt) > 0){
        query$dist_to_BP_ALT[queryIndex] <- diffs.alt[which.min(abs(diffs.alt))]
      }
          
      query$max_prob_REF[queryIndex] <- 
        max(predictions.snp$branchpoint_prob[predictions.snp$status == "REF"])
      query$max_prob_ALT[queryIndex] <- 
        max(predictions.snp$branchpoint_prob[predictions.snp$status == "ALT"])
      
      # maximum U2 binding energy of all branchpoint sites
      # above the probability cutoff
      if(length(diffs.ref) > 0){
        query$max_U2_REF[queryIndex] <- 
          max(predictions.snp$U2_binding_energy[branchpointIndex.ref])
      }
      
      if(length(diffs.alt) > 0){
        query$max_U2_ALT[queryIndex] <- 
          max(predictions.snp$U2_binding_energy[branchpointIndex.alt])
      }

      # all sites (ref & alt) called as branchpoints
      branchpointLocs <- unique(predictions.snp$to_3prime_point[c(branchpointIndex.ref, branchpointIndex.alt)])

      #probabiltiy scores for ref/alt at all BP sites
      prob.ref <- predictions.snp$branchpoint_prob[which(predictions.snp$status == "REF" & 
                                                           predictions.snp$to_3prime_point %in% branchpointLocs)]
      prob.alt <- predictions.snp$branchpoint_prob[which(predictions.snp$status == "ALT" & 
                                                           predictions.snp$to_3prime_point %in% branchpointLocs)]

      #change must be of sufficient magnitude to be called as created or deleted
      query$deleted_n[queryIndex] <- length(which((prob.ref - prob.alt) > probabilityChange &
                                                          (prob.ref < probabilityCutoff | prob.alt < probabilityCutoff)))
      query$created_n[queryIndex] <- length(which((prob.ref - prob.alt) < (probabilityChange * -1) &
                                                          (prob.ref < probabilityCutoff | prob.alt < probabilityCutoff)))
    }
  
  return(query)
}
