#' Plots transcript structures
#'
#' Plots transcript structures
#' @param exonID id of the exon to plot
#' @param exons Granges containing exon co-ordinates.
#' @param keepTranscripts which transcripts to plot (\code{"overlapping"} or \code{"withExon"})
#' \code{"overlapping"} will plot all transcripts overlapping the exon, 
#' whereas \code{"withExon"} will plot all transcripts containing the exon.
#' @return ggplot2 plot transcript structures
#' @import ggplot2
#' @export
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",
#' package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' plotStructure(exonID = "ENSE00001184784.4", exons)
#' @author Beth Signal
plotStructure <- function(exonID, exons, keepTranscripts="overlapping"){
  
  geneID <- exons$gene_id[match(exonID, exons$exon_id)]
  
  exonsForPlot <- exons[exons$gene_id == geneID]  
  
  exonsForPlot <- cbind(chromosome=as.character(seqnames(exonsForPlot)), 
                        as.data.frame(ranges(exonsForPlot)),
                        strand=as.character(strand(exonsForPlot)),
                        as.data.frame(mcols(exonsForPlot)))
  
  exonsForPlot$start <- as.numeric(exonsForPlot$start)
  exonsForPlot$end <- as.numeric(exonsForPlot$end)
  
  transcriptsForPlot <- exonsForPlot
  transcriptsForPlot$transcript_id <- transcriptsForPlot$gene_id
  exonsForPlot <- rbind(exonsForPlot, transcriptsForPlot)  
  
  exonsForPlot$length <- exonsForPlot$end - exonsForPlot$start
  
  
  if (as.logical(exonsForPlot$strand[1] == "+")) {
    BPexons = which(exonsForPlot$start ==
                      exonsForPlot$start[which(exonsForPlot$exon_id ==
                                                 exonID)][1])
  }else{
    BPexons <- which(exonsForPlot$end ==
                       exonsForPlot$end[which(exonsForPlot$exon_id ==
                                                exonID)][1])
  }
  
  exonsForPlot$SNP_exon <- 0
  exonsForPlot$SNP_exon[BPexons] <- 1
  exonsForPlot$SNP_exon[exonsForPlot$transcript_id == exonsForPlot$gene_id[1]] <- 2
  
  #make smaller data.frame for plotting connecting lines
  exonsForPlotLines <- exonsForPlot[!duplicated(exonsForPlot$transcript_id),]
  
  exonsForPlotLines <- stats::aggregate(start ~ transcript_id, exonsForPlot, max)
  exonsForPlotLines <- cbind(exonsForPlotLines, 
                             stats::aggregate(end ~ transcript_id, exonsForPlot, min)[,2],
                             stats::aggregate(end ~ transcript_id, exonsForPlot, max)[,2],
                             stats::aggregate(start ~ transcript_id, exonsForPlot, min)[,2])
  
  colnames(exonsForPlotLines) <- c("transcript_id", "maxs","mine","maxe","mins")
  
  #only plot exonsForPlot overlapping the query
  if(keepTranscripts == "overlapping"){
    keep <- which((exonsForPlotLines$maxe >= max(exonsForPlot$end[exonsForPlot$SNP_exon == 1])) &
                    (exonsForPlotLines$mins <= min(exonsForPlot$start[exonsForPlot$SNP_exon == 1])))
    exonsForPlotLines <- exonsForPlotLines[keep,]
    exonsForPlot <- exonsForPlot[exonsForPlot$transcript_id %in% exonsForPlotLines$transcript_id,]
  }else if(keepTranscripts == "withExon"){
    keep <- which(exonsForPlotLines$transcript_id %in% exonsForPlot$transcript_id[exonsForPlot$SNP_exon == 1])
    exonsForPlotLines <- exonsForPlotLines[keep,]
    exonsForPlot <- exonsForPlot[exonsForPlot$transcript_id %in% exonsForPlotLines$transcript_id,]
  }
  
  exonsForPlot$transcript_id_num <- as.numeric(as.factor(exonsForPlot$transcript_id)) * -1
  exonsForPlotLines$transcript_id_num <- as.numeric(as.factor(exonsForPlotLines$transcript_id)) * -1
  
  if (as.logical(exonsForPlot$strand[1] == "+")) {
    exonsForPlotLines$transcript_id <- paste0(exonsForPlotLines$transcript_id, " (+)")
  }else{
    exonsForPlotLines$transcript_id <- paste0(exonsForPlotLines$transcript_id, " (-)")
  }
  exonsForPlotLines$SNP_exon <- 0
  
  plot.structure <- ggplot() +
    geom_segment(data = exonsForPlotLines, aes_string(x = 'mins',xend = 'maxe', y = 'transcript_id_num',
                                                      yend = 'transcript_id_num')) +
    geom_rect(data=exonsForPlot, aes_string(xmin = 'start', xmax = 'start + length',
                                            ymin = 'transcript_id_num - 0.4',
                                            ymax = 'transcript_id_num + 0.4', fill = 'factor(SNP_exon)')) +
    scale_fill_manual(values = c("black","blue","grey60")) +
    scale_y_continuous(breaks = seq(-1, min(exonsForPlotLines$transcript_id_num),-1),
                       labels = exonsForPlotLines$transcript_id[
                         match(seq(-1, min(exonsForPlotLines$transcript_id_num),-1),
                               exonsForPlotLines$transcript_id_num)]) +
    theme(panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none")
  
  return(plot.structure)
}

#' Plots branchpointer predictions
#'
#' Plots branchpointer predictions
#' @param queryName query id used to identify the SNP or region
#' @param predictions Granges object generated by predictBranchpoints()
#' @param probabilityCutoff probability score cutoff value for displaying U2 binding energy
#' @param plotStructure plot structures for gene and 3' exon containing and skipping isoforms
#' @param plotMutated plot alternative sequence predicitons alongside reference sequence predictions
#' @param exons Granges containing exon co-ordinates.
#' Should be produced by gtfToExons()
#' @return ggplot2 plot with branchpoint features in the specified intronic region
#' @export
#' @import ggplot2
#' @importFrom cowplot ggdraw
#' @importFrom cowplot draw_plot
#' @examples
#' smallExons <- system.file("extdata","gencode.v24.annotation.small.gtf",package = "branchpointer")
#' exons <- gtfToExons(smallExons)
#' genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'
#' querySNP <- system.file("extdata","SNP_example.txt", package = "branchpointer")
#' query <- readQueryFile(querySNP,queryType = "SNP")
#' query <- getQueryLoc(query,queryType = "SNP",exons = exons, filter = FALSE)
#' 
#' predictions <- predictBranchpoints(query,queryType = "SNP", BSgenome = genome)
#' plotBranchpointWindow(query$id[1], predictions,
#' plotMutated = TRUE, exons = exons)
#' @author Beth Signal

plotBranchpointWindow <- function(queryName,
                                 predictions,
                                 probabilityCutoff = 0.5,
                                 plotMutated = FALSE,
                                 plotStructure = TRUE,
                                 exons) {
  #find the query
  ind <- which(!is.na(match(predictions$id, queryName)))
  predictions.snp <- predictions[which(predictions$id == queryName)]

  predictions.snp$to_3prime_point <- predictions.snp$to_3prime_point * -1
  
  predictions.snp.ref <- as.data.frame(mcols(predictions.snp)[predictions.snp$status == "REF",])
  colnames(predictions.snp.ref) <- gsub("to_3prime_point","distance",colnames(predictions.snp.ref))
  colnames(predictions.snp.ref) <- gsub("seq_pos0","nucleotide",colnames(predictions.snp.ref))
  predictions.snp.ref$nucleotide <- as.character(predictions.snp.ref$nucleotide)
  
  #make plots for reference sequence
  
  #probability score by sequence position
  plot.prob.ref <- ggplot(predictions.snp.ref, aes_string(x = 'distance', y = 'branchpoint_prob', fill = 'nucleotide',
                                           col = 'nucleotide',alpha = 'U2_binding_energy')) +
    geom_bar(stat = "identity") +
    scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                       breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                          "0.60","0.70","0.80","0.90","1.00")) +
    scale_fill_manual(values = mercer_nt_cols, drop = TRUE,
                      limits = c("A","C","G","T")) +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  #U2 binding energy by sequence position
  plot.U2.ref <- ggplot(predictions.snp.ref[predictions.snp.ref$branchpoint_prob >= probabilityCutoff,],
              aes_string(x = 'distance', y = 'U2_binding_energy')) +
    geom_bar(stat = "identity",width = 1) +
    scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                       breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
    theme(legend.position = "none") +
    scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                       breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

  #Sequence identity
  plot.seq.ref <- ggplot(predictions.snp.ref, aes_string(x = 'distance', y = 1, col = 'nucleotide',label = 'nucleotide')) +
    geom_text(size = 4, family = "Courier") +
    scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                       limits = c("A","C","G","T")) +
    theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
          axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
          panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
          axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(-45,-17))

  if (plotMutated == TRUE) {
    
    predictions.snp.alt <- as.data.frame(mcols(predictions.snp)[predictions.snp$status == "ALT",])
    colnames(predictions.snp.alt) <- gsub("to_3prime_point","distance",colnames(predictions.snp.alt))
    colnames(predictions.snp.alt) <- gsub("seq_pos0","nucleotide",colnames(predictions.snp.alt))
    predictions.snp.alt$nucleotide <- as.character(predictions.snp.alt$nucleotide)
    
    #probability score by sequence position
    plot.prob.alt <- ggplot(predictions.snp.alt, aes_string(x = 'distance', y = 'branchpoint_prob', fill = 'nucleotide',
                                             col = 'nucleotide',alpha = 'U2_binding_energy')) +
      geom_bar(stat = "identity") +
      scale_y_continuous(name = "branchpointer probability score",limits = c(0,1),
                         breaks = seq(0,1,0.10), labels = c("0.00","0.10","0.20","0.30","0.40","0.50",
                                                            "0.60","0.70","0.80","0.90","1.00")) +
      scale_fill_manual(values = mercer_nt_cols, drop = TRUE,
                        limits = c("A","C","G","T")) +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = c("A","C","G","T")) +
      theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17)) + ggtitle("Alternative")

    plot.prob.ref <- plot.prob.ref + ggtitle("Reference")

    #U2 binding energy by sequence position
    plot.U2.alt <- ggplot(predictions.snp.alt[predictions.snp.alt$branchpoint_prob >= probabilityCutoff,],
                aes_string(x = 'distance', y = 'U2_binding_energy')) +
      geom_bar(stat = "identity",width = 1) +
      scale_x_continuous(limits = c(-45,-17), labels = seq(45,17,-5),
                         breaks = seq(-45,-17, 5), name ="Distance to 3' exon (nt)") +
      theme(legend.position = "none") +
      scale_y_continuous(name = "U2 binding energy",limits = c(0,10),
                         breaks = seq(0,9,3), labels = c("0.00","3.00","6.00","9.00"))

    #Sequence identity
    plot.seq.alt <- ggplot(predictions.snp.alt, aes_string(x = 'distance', y = 1, col = 'nucleotide',label = 'nucleotide')) +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values = mercer_nt_cols,drop = TRUE,
                         limits = c("A","C","G","T")) +
      theme(legend.position = "none",axis.text.y = element_text(colour = "white"),
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-45,-17))

    ######Get reference vs. alt sequence for 0-50 window#####
    
    seqSplit.ref <- unlist(stringr::str_split(predictions.snp.ref$seq[1], ""))
    seqSplit.alt <- unlist(stringr::str_split(predictions.snp.alt$seq[1], ""))
    
    wholeSeq <- data.frame(nt = c(seqSplit.ref, seqSplit.alt),
                           pos = rep(seq(-50,-1,1),2),
                           set = rep(c("REF","ALT"), each =50))
    mutPosition <- wholeSeq$pos[which(wholeSeq$nt[wholeSeq$set == "REF"] !=
                                    wholeSeq$nt[wholeSeq$set == "ALT"])]

    plot.seq.comparison <- ggplot(wholeSeq, aes_string(x = 'pos', y = 'set', col = 'nt',label = 'nt')) +
      geom_rect(xmin = mutPosition - 0.5,xmax = mutPosition + 0.5, ymin = 0.5, ymax = 0.5 + 2,
                col = NA, fill = "grey90") +
      geom_text(size = 4, family = "Courier") +
      scale_color_manual(values = mercer_nt_cols) +
      theme(legend.position = "none",
            axis.title.y = element_text(colour ="white"), axis.ticks.y = element_line(color = "white"),
            panel.background = element_rect(fill = "white"), axis.text.x = element_blank(),
            axis.title.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_continuous(limits = c(-50,-1)) +
      scale_y_discrete(labels = c("Alternative", "Reference")) +
      geom_rect(xmin = -44.5,xmax = -17.5, ymin = 0.5, ymax = 0.5 + 2, fill = NA, col ="black")

  }

  #plot gene/transcript structure
  if (plotStructure == TRUE) {
    exonID <- predictions.snp$exon_3prime[1]
    plot.structure <- plotStructure(exonID, exons)
  }
    
  theme_set(theme_gray())
  if (plotStructure == TRUE & plotMutated == TRUE) {
    ggdraw() +
      draw_plot(plot.structure,0,0.775,1,0.225) + 
      draw_plot(plot.seq.comparison,0,0.675,1,0.1) +
      draw_plot(plot.prob.ref,0,.275,0.5,.40) + draw_plot(plot.seq.ref,0,.2,0.5,.075) + draw_plot(plot.U2.ref,0,0,0.5,.2) +
      draw_plot(plot.prob.alt,0.5,.275,0.5,.40) + draw_plot(plot.seq.alt,0.5,.2,0.5,.075) + draw_plot(plot.U2.alt,0.5,0,0.5,.2)

  }else if (plotStructure == TRUE & plotMutated == FALSE) {
    ggdraw() +
      draw_plot(plot.structure,0,0.775,1,0.225) +
      draw_plot(plot.prob.ref,0,.325,1,.45) + draw_plot(plot.seq.ref,0,.25,1,.075) + draw_plot(plot.U2.ref,0,0,1,.25)
  }else if (plotStructure == FALSE & plotMutated == FALSE) {
    ggdraw() +
      draw_plot(plot.prob.ref,0,.325,1,.45) + draw_plot(plot.seq.ref,0,.25,1,.075) + draw_plot(plot.U2.ref,0,0,1,.25)
  }else if (plotStructure == FALSE & plotMutated == TRUE) {
    ggdraw() +
      draw_plot(plot.seq.comparison,0,0.9,1,0.1) +
      draw_plot(plot.prob.ref,0,.4,0.5,.5) + draw_plot(plot.seq.ref,0,.3,0.5,.1) + draw_plot(plot.U2.ref,0,0,0.5,.3) +
      draw_plot(plot.prob.alt,0.5,.4,0.5,.5) + draw_plot(plot.seq.alt,0.5,.3,0.5,.1) + draw_plot(plot.U2.alt,0.5,0,0.5,.3)
  }

}
