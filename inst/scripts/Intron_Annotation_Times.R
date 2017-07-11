library(branchpointer)

gtf <- rtracklayer::import("gencode.v26.annotation.gtf")

pc_level1_transcripts <- gtf$transcript_id[which(gtf$exon_number == 2 &
                                                     gtf$type=="exon" &
                                                     seqnames(gtf) == "chr22" &
                                                     gtf$transcript_type=="protein_coding" &
                                                     gtf$transcript_support_level==1)]

linc_transcripts <- gtf$transcript_id[which(gtf$exon_number == 2 &
                                                gtf$type=="exon" &
                                                seqnames(gtf) == "chr22" &
                                                gtf$gene_type=="lincRNA")]


t1 <- system.time({g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38})
t2 <- system.time({exons <- gtfToExons("gencode.v26.annotation.gtf")})

t3 <- system.time({query_chr22.pc <- makeBranchpointWindowForExons(id=pc_level1_transcripts, idType = "transcript_id", exons)})
t4 <- system.time({predictions_chr22 <- predictBranchpoints(query_chr22.pc, queryType = "region",BSgenome = g)})


t5 <- system.time({query_chr22.linc <- makeBranchpointWindowForExons(id=linc_transcripts, idType = "transcript_id", exons)})
t6 <- system.time({predictions_chr22.linc <- predictBranchpoints(query_chr22.linc, queryType = "region",BSgenome = g)})
