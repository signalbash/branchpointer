# Test times for 'readQueryFile()' using filtering

library(branchpointer)
library(data.table)
library(ggplot2)

file_name <- "snp_subset.txt"

gtex_vars <- fread("data/GTEX/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt")
gtex_vars <- as.data.frame(gtex_vars)

gtex_vars$Chr <- paste0("chr", gtex_vars$Chr)

norm <- nchar(gtex_vars$Ref_b37)
alt <- nchar(gtex_vars$Alt)

keep <- which(norm == 1 & alt == 1)
gtex_vars <- gtex_vars[keep,]

snp_info_query <- data.frame(id=gtex_vars[,3], chromosome=gtex_vars[,1], 
                             chrom_start=gtex_vars[,2], strand=2,
                             allele_ref=gtex_vars[,4],
                             allele_alt=gtex_vars[,5])

write.table(snp_info_query, file=file_name, row.names = F)

exons <- gtfToExons("gencode.v26.annotation.gtf")
query_SNP <- readQueryFile(file_name, queryType = "SNP", exons=exons)

gtex_snps_NOTwindow
snp_info_query
gtex_snps_window <- snp_info_query[snp_info_query$id %in% query]

#### times with 1% and 0.1% SNP windows

snp_reading_time.filter.1 <- list()
snp_reading_time.Nofilter.1 <- list()
snp_reading_time.filter.01 <- list()
snp_reading_time.Nofilter.01 <- list()

query_set_sizes <- rep(c(1,2,5,10,20,50,100,200,500, 1000,2000, 5000)*1000, each = 5)

for(i in seq_along(query_set_sizes)){
    
    message(query_set_sizes[i])
    
    # 1% SNP in window
    number_windowSNP <- query_set_sizes[i] * 0.01
    if(number_windowSNP < nrow(gtex_snps_window)){
        
        snp_file <- rbind(gtex_snps_NOTwindow[sample(1:(dim(gtex_snps_NOTwindow)[1]), query_set_sizes[i] - number_windowSNP), ],
                          gtex_snps_window[sample(1:(dim(gtex_snps_window)[1]), number_windowSNP),])
        write.table(snp_file, file=file_name, row.names = F)
        
        snp_reading_time.filter.1[[i]] <- system.time({
            query <- readQueryFile(file_name, queryType = "SNP",exons = exons, filter=TRUE)
        })
        message(round(snp_reading_time.filter.1[[i]][3], 2))
        snp_reading_time.Nofilter.1[[i]] <- system.time({
            query <- readQueryFile(file_name, queryType = "SNP",exons = exons, filter=FALSE)
        })
        message(round(snp_reading_time.Nofilter.1[[i]][3], 2))
    }else{
        snp_reading_time.filter.1[[i]] <- rep(NA,3)
        snp_reading_time.Nofilter.1[[i]] <- rep(NA,3)
    }
    
    # 0.1% SNP in window
    number_windowSNP <- query_set_sizes[i] * 0.001
    snp_file <- rbind(gtex_snps_NOTwindow[sample(1:(dim(gtex_snps_NOTwindow)[1]), query_set_sizes[i] - number_windowSNP), ],
                      gtex_snps_window[sample(1:(dim(gtex_snps_window)[1]), number_windowSNP),])
    write.table(snp_file, file=file_name, row.names = F)
    
    snp_reading_time.filter.01[[i]] <- system.time({
        query <- readQueryFile(file_name, queryType = "SNP",exons = exons, filter=TRUE)
    })
    message(round(snp_reading_time.filter.01[[i]][3], 2))
    snp_reading_time.Nofilter.01[[i]] <- system.time({
        query <- readQueryFile(file_name, queryType = "SNP",exons = exons, filter=FALSE)
    })
    message(round(snp_reading_time.Nofilter.01[[i]][3], 2))
    
}

getTime <- function(list, index = 3){
    return(unlist(lapply(list,"[[", index)))
}

snp_reading_times <- data.frame(query_size=rep(query_set_sizes, 4),
                                snps_percent=rep(c(1,0.1), each=length(query_set_sizes)*2),
                                filter=rep(rep(factor(c("TRUE","FALSE"), levels=c("TRUE", "FALSE")), each = length(query_set_sizes)),4),
                                time_elapsed=c(getTime(snp_reading_time.filter.1),getTime(snp_reading_time.Nofilter.1),
                                               getTime(snp_reading_time.filter.01),getTime(snp_reading_time.Nofilter.01)))
snp_reading_times <- snp_reading_times[!is.na(snp_reading_times$time_elapsed),]
snp_reading_times$snps <- snp_reading_times$query_size * snp_reading_times$snps_percent/100

snp_reading_times$set <- paste0(snp_reading_times$filter, "_", snp_reading_times$snps_percent)

Vignette_Figure_readingTime <- 
  ggplot(snp_reading_times[snp_reading_times$snps_percent==0.1,], 
       aes(x=query_size, y=time_elapsed, col=filter)) + 
    geom_point() + 
    scale_color_brewer(palette = "Set1", name="filter") +
    geom_smooth(method="loess") + 
    scale_y_log10(name="Time elapsed (seconds)") + scale_x_log10(name="Number of queries read") 

pdf(file="vignettes/readQueryFile_time.pdf", width=4,height=3)
Vignette_Figure_readingTime
dev.off()
