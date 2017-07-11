# Test times for 'predictBranchpoints()' using parallelisation

library(branchpointer)
library(reshape2)
library(ggplot2)

exons <- gtfToExons("gencode.v26.annotation.gtf")
g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

gene_ids <- unique(exons@elementMetadata$gene_id[exons@elementMetadata$exon_number == 2])

times <- list()
times_par2 <- list()
times_par4 <- list()
times_par8 <- list()

gene_set_sizes <- rep(c(1,2,5,10,20,50,100,200,500), each = 5)

query <- makeBranchpointWindowForExons(gene_ids, "region", exons)

for(i in 1:length(gene_set_sizes)){
  
  index <- sample(1:length(gene_ids), gene_set_sizes[i])
  
  times_par2[[i]] <- system.time({
    preds <- predictBranchpoints(query[index], queryType = "region",
                                 BSgenome = g, useParallel = TRUE, cores=2)
  })
  times_par4[[i]] <- system.time({
    preds <- predictBranchpoints(query[index], queryType = "region",
                                 BSgenome = g, useParallel = TRUE, cores=4)
  })
  times_par8[[i]] <- system.time({
      preds <- predictBranchpoints(query[index], queryType = "region",
                                   BSgenome = g, useParallel = TRUE, cores=8)
  })
  times[[i]] <- system.time({
    preds <- predictBranchpoints(query[index], queryType = "region", BSgenome = g)
  })
  message(gene_set_sizes[i])
}

predictBranchpointsTimes <- data.frame(introns = rep(gene_set_sizes, 4), 
                                       times = c(unlist(lapply(times, "[[", 3)),
                                                 unlist(lapply(times_par2, "[[", 3)),
                                                 unlist(lapply(times_par4, "[[", 3)),
                                                 unlist(lapply(times_par8, "[[", 3))),
                                       useParallel = rep(factor(c("FALSE","TRUE","TRUE","TRUE"), levels=c("TRUE","FALSE")), each=length(gene_set_sizes)),
                                       parallel = rep(c(1,2,4,8), each=length(gene_set_sizes)))

Vignette_Figure_predictTime <- ggplot(predictBranchpointsTimes, aes(x=introns, y=times, col=factor(parallel), shape=useParallel)) + 
  geom_point() + 
  scale_x_continuous(name="Query regions processed") +
  scale_y_continuous(name="Time elapsed (seconds)") +
  scale_color_brewer(palette="Set1", name="Cores") +
  geom_smooth(method="lm")

pdf(file="vignettes/predictBranchpoints_time.pdf", width=4,height=3)
Vignette_Figure_predictTime
dev.off()

lm(times ~ introns, data=predictBranchpointsTimes[predictBranchpointsTimes$parallel == 1,])

